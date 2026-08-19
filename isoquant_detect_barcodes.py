#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2023-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Standalone barcode detection tool for single-cell and spatial transcriptomics.

This is a CLI wrapper around isoquant_lib.barcode_calling.detect_barcodes.
"""

import argparse
import logging
import os
import sys
from traceback import print_exc

from isoquant_lib.modes import (IsoQuantMode, BarcodeCorrectionMethod, LARGE_WHITELIST_SIZE,
                                AUTO_BARCODES, DEPRECATED_MODE_ALIASES, SPLIT_MOLECULES_CHOICES,
                                SPLIT_MOLECULES_TRUE, SPLIT_MOLECULES_FALSE, SPLIT_MOLECULES_AUTO)
from isoquant_lib.utils.error_codes import IsoQuantExitCode
from isoquant_lib.barcode_calling.detect_barcodes import (
    process_single_thread,
    process_in_parallel,
    get_barcode_length,
    BARCODE_CALLING_MODES,
)
from isoquant_lib.barcode_calling.detect_barcodes import detect_cell_barcode_list

logger = logging.getLogger('IsoQuant')


def set_logger(logger_instance, args):
    logger_instance.setLevel(logging.INFO)
    if args.debug:
        logger_instance.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    if args.debug:
        ch.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def parse_args(sys_argv):
    def add_hidden_option(*args, **kwargs):  # show command only with --full-help
        kwargs['help'] = argparse.SUPPRESS
        parser.add_argument(*args, **kwargs)

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", required=True)
    parser.add_argument("--barcodes", "-b", nargs='+', type=str, help="barcode whitelist(s)", required=False)
    mode_names = [x.name for x in BARCODE_CALLING_MODES.keys()]
    parser.add_argument("--mode", type=str, help="mode to be used",
                        choices=mode_names + list(DEPRECATED_MODE_ALIASES),
                        metavar="{%s}" % ",".join(mode_names))
    parser.add_argument("--split_molecules", type=str, choices=SPLIT_MOLECULES_CHOICES, default=None,
                        help="split reads containing several cDNA molecules ligated end to end: "
                             "%s splits wherever the protocol allows it, %s also fails for protocols "
                             "that cannot [%s]" % (SPLIT_MOLECULES_AUTO, SPLIT_MOLECULES_TRUE,
                                                   SPLIT_MOLECULES_AUTO))
    parser.add_argument("--molecule", type=str, help="MDF files with molecule description (for custom_sc mode only)")

    parser.add_argument("--input", "-i", nargs='+', type=str, help="input reads in [gzipped] FASTA, FASTQ, BAM",
                        required=True)
    parser.add_argument("--threads", "-t", type=int, help="threads to use [16]", default=16)
    parser.add_argument("--tmp_dir", type=str, help="folder for temporary files")
    parser.add_argument("--n_cells", type=str,
                        help='expected number of cell-associated barcodes, or "%s" to estimate it from the '
                             'barcode count distribution. When set, the whitelist is treated as a pool to '
                             'select cell barcodes from; when omitted, the whitelist is taken to be the '
                             'cell barcodes themselves' % AUTO_BARCODES)
    add_hidden_option("--n_cells_interval", type=int, default=25)
    add_hidden_option("--barcode_correction", type=str, choices=[e.name for e in BarcodeCorrectionMethod],
                      default=BarcodeCorrectionMethod.auto.name)
    add_hidden_option('--debug', action='store_true', default=False, help='debug log output.')

    args = parser.parse_args(sys_argv)
    args.out_fasta = None
    args.output_tsv = None
    return args


def resolve_deprecated_mode(args):
    """Translate a superseded mode name into a chemistry plus a --split_molecules default."""
    alias = DEPRECATED_MODE_ALIASES.get(args.mode)
    if not alias:
        return
    mode_name, splits = alias
    logger.warning("Mode %s is deprecated, use `--mode %s --split_molecules %s` instead" %
                   (args.mode, mode_name, SPLIT_MOLECULES_TRUE if splits else SPLIT_MOLECULES_FALSE))
    args.mode = mode_name
    if args.split_molecules is None:
        args.split_molecules = SPLIT_MOLECULES_TRUE if splits else SPLIT_MOLECULES_FALSE


def resolve_split_molecules(args):
    """Turn --split_molecules into a bool. See isoquant.py for the rules."""
    requested = args.split_molecules or SPLIT_MOLECULES_AUTO
    supported = args.mode.supports_molecule_splitting()
    if requested == SPLIT_MOLECULES_FALSE:
        args.split_molecules = False
        return
    if requested == SPLIT_MOLECULES_TRUE and not supported:
        logger.critical("Mode %s cannot split reads into separate molecules" % args.mode.name)
        sys.exit(IsoQuantExitCode.INCOMPATIBLE_OPTIONS)
    args.split_molecules = supported


def count_whitelist_barcodes(whitelist_files):
    import gzip
    total = 0
    for file_name in whitelist_files:
        handle = gzip.open(file_name, "rt") if file_name.endswith(("gz", "gzip")) else open(file_name)
        with handle:
            total += sum(1 for _ in handle)
    return total


def resolve_cell_barcode_detection(args):
    """Decide whether cell barcodes are supplied or detected. See isoquant.py for the rules."""
    args.detect_cell_barcodes = False
    whitelist_is_auto = args.barcodes == [AUTO_BARCODES]

    if args.n_cells is not None and args.n_cells != AUTO_BARCODES:
        try:
            args.n_cells = int(args.n_cells)
        except ValueError:
            logger.critical('--n_cells must be a positive integer or "%s"' % AUTO_BARCODES)
            sys.exit(IsoQuantExitCode.INVALID_PARAMETER)
        if args.n_cells <= 0:
            logger.critical("--n_cells must be positive")
            sys.exit(IsoQuantExitCode.INVALID_PARAMETER)

    requested = BarcodeCorrectionMethod[args.barcode_correction]
    if requested == BarcodeCorrectionMethod.whitelist:
        if whitelist_is_auto:
            logger.critical('--barcode_correction whitelist cannot be used with --barcodes %s' % AUTO_BARCODES)
            sys.exit(IsoQuantExitCode.INCOMPATIBLE_OPTIONS)
        return

    if whitelist_is_auto and args.n_cells is None:
        args.n_cells = AUTO_BARCODES

    if not (whitelist_is_auto or args.n_cells is not None or requested == BarcodeCorrectionMethod.graph):
        if args.barcodes and count_whitelist_barcodes(args.barcodes) > LARGE_WHITELIST_SIZE:
            logger.warning("Barcode whitelist is large and is treated as the list of cell barcodes; "
                           "matching every read against it effectively requires an exact match. "
                           'Set --n_cells (or --n_cells %s) to select cell barcodes instead.' % AUTO_BARCODES)
        return

    if not args.mode.supports_cell_barcode_detection():
        logger.critical("Detecting cell barcodes from the data is not supported for mode %s" % args.mode.name)
        sys.exit(IsoQuantExitCode.INCOMPATIBLE_OPTIONS)
    if requested == BarcodeCorrectionMethod.graph and args.n_cells is None:
        args.n_cells = AUTO_BARCODES
    args.detect_cell_barcodes = True


def check_args(args):
    """Set up output file lists based on input files."""
    # args.input is always a list (nargs='+')
    num_files = len(args.input)

    if args.output_tsv is None:
        if num_files == 1:
            args.output_tsv = [args.output + ".barcoded_reads.tsv"]
        else:
            args.output_tsv = [args.output + "_%d.barcoded_reads.tsv" % i for i in range(num_files)]

    if args.out_fasta is None and args.split_molecules:
        if num_files == 1:
            args.out_fasta = [args.output + ".split_reads.fasta"]
        else:
            args.out_fasta = [args.output + "_%d.split_reads.fasta" % i for i in range(num_files)]


def run_barcode_calling(args):
    if args.threads == 1 or args.mode.enforces_single_thread():
        process_single_thread(args)
    else:
        process_in_parallel(args)


def main(sys_argv):
    args = parse_args(sys_argv)
    set_logger(logger, args)
    # resolved here rather than in parse_args so the warnings reach the configured logger
    resolve_deprecated_mode(args)
    args.mode = IsoQuantMode[args.mode]
    resolve_split_molecules(args)
    check_args(args)
    resolve_cell_barcode_detection(args)

    out_dir = os.path.dirname(args.output)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    if args.detect_cell_barcodes:
        # first pass: extract barcodes only to count them. Nothing is written, and no FASTA
        # is produced -- the second pass does the extraction whose results are kept.
        final_output_tsv, final_out_fasta = args.output_tsv, args.out_fasta
        args.output_tsv, args.out_fasta = None, None
        args.whitelist_matching = False
        cell_barcodes = detect_cell_barcode_list(args, args.output + ".cell_barcodes.tsv",
                                                 get_barcode_length(args.mode),
                                                 args.n_cells, args.n_cells_interval,
                                                 args.output + ".cell_barcodes.stats")
        # second pass: ordinary barcode calling against the detected cell barcodes
        args.barcodes = [cell_barcodes]
        args.output_tsv, args.out_fasta = final_output_tsv, final_out_fasta

    args.whitelist_matching = True
    run_barcode_calling(args)


def main_entry():
    """Entry point for console_scripts (pip install)."""
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except Exception:
        print_exc()
        sys.exit(IsoQuantExitCode.UNCAUGHT_EXCEPTION)


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except Exception:
        print_exc()
        sys.exit(IsoQuantExitCode.UNCAUGHT_EXCEPTION)
