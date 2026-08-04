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

from isoquant_lib.modes import IsoQuantMode, BarcodeCorrectionMethod, LARGE_WHITELIST_SIZE
from isoquant_lib.utils.error_codes import IsoQuantExitCode
from isoquant_lib.barcode_calling.detect_barcodes import (
    process_single_thread,
    process_in_parallel,
    get_barcode_length,
    BARCODE_CALLING_MODES,
)
from isoquant_lib.barcode_calling.correct_barcodes import correct_barcodes

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
    parser.add_argument("--mode", type=str, help="mode to be used", choices=[x.name for x in BARCODE_CALLING_MODES.keys()])
    parser.add_argument("--molecule", type=str, help="MDF files with molecule description (for custom_sc mode only)")

    parser.add_argument("--input", "-i", nargs='+', type=str, help="input reads in [gzipped] FASTA, FASTQ, BAM",
                        required=True)
    parser.add_argument("--threads", "-t", type=int, help="threads to use [16]", default=16)
    parser.add_argument("--tmp_dir", type=str, help="folder for temporary files")
    parser.add_argument("--barcode_correction", type=str, choices=[e.name for e in BarcodeCorrectionMethod],
                        default=BarcodeCorrectionMethod.auto.name,
                        help="how to correct extracted barcodes: match each read against the whitelist, "
                             "or select actually sequenced barcodes and correct via an edit-distance "
                             "graph [%s: graph for large whitelists]" % BarcodeCorrectionMethod.auto.name)
    parser.add_argument("--n_cells", type=int,
                        help="expected number of cell-associated barcodes, used by "
                             "--barcode_correction graph [estimated from barcode counts]")
    add_hidden_option("--n_cells_interval", type=int, default=25)
    add_hidden_option("--barcode_graph_threshold", type=int, default=1)
    add_hidden_option("--barcode_graph_rounds", type=int, default=2)
    add_hidden_option("--barcode_graph_impl", type=str, choices=["centers", "full"], default="centers")
    add_hidden_option('--debug', action='store_true', default=False, help='debug log output.')

    args = parser.parse_args(sys_argv)
    args.mode = IsoQuantMode[args.mode]
    args.out_fasta = None
    args.output_tsv = None
    return args


def count_whitelist_barcodes(whitelist_files):
    import gzip
    total = 0
    for file_name in whitelist_files:
        handle = gzip.open(file_name, "rt") if file_name.endswith(("gz", "gzip")) else open(file_name)
        with handle:
            total += sum(1 for _ in handle)
    return total


def resolve_barcode_correction(args):
    """Turn --barcode_correction into args.whitelist_matching, validating the request."""
    args.whitelist_matching = True
    requested = BarcodeCorrectionMethod[args.barcode_correction]
    if requested == BarcodeCorrectionMethod.whitelist:
        return

    if not args.mode.supports_graph_correction():
        if requested == BarcodeCorrectionMethod.graph:
            logger.critical("Graph-based barcode correction is not supported for mode %s" % args.mode.name)
            sys.exit(IsoQuantExitCode.INCOMPATIBLE_OPTIONS)
        return
    if not args.barcodes:
        if requested == BarcodeCorrectionMethod.graph:
            logger.critical("Graph-based barcode correction requires --barcodes")
            sys.exit(IsoQuantExitCode.BARCODE_WHITELIST_MISSING)
        return

    if requested == BarcodeCorrectionMethod.auto:
        barcode_count = count_whitelist_barcodes(args.barcodes)
        if barcode_count <= LARGE_WHITELIST_SIZE:
            logger.info("Whitelist contains %d barcodes, using per-read whitelist matching" % barcode_count)
            return
        logger.info("Whitelist contains %d barcodes, per-read matching would degenerate into exact "
                    "matching; using graph-based barcode correction instead" % barcode_count)

    args.whitelist_matching = False


def check_args(args):
    """Set up output file lists based on input files."""
    # args.input is always a list (nargs='+')
    num_files = len(args.input)

    if args.output_tsv is None:
        if num_files == 1:
            args.output_tsv = [args.output + ".barcoded_reads.tsv"]
        else:
            args.output_tsv = [args.output + "_%d.barcoded_reads.tsv" % i for i in range(num_files)]

    if args.out_fasta is None and args.mode.produces_new_fasta():
        if num_files == 1:
            args.out_fasta = [args.output + ".split_reads.fasta"]
        else:
            args.out_fasta = [args.output + "_%d.split_reads.fasta" % i for i in range(num_files)]


def main(sys_argv):
    args = parse_args(sys_argv)
    set_logger(logger, args)
    check_args(args)
    resolve_barcode_correction(args)

    out_dir = os.path.dirname(args.output)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # with graph correction the detector writes uncorrected barcodes to an intermediate
    # table, and the correction stage produces the final barcoded read tables
    final_output_tsv = args.output_tsv
    if not args.whitelist_matching:
        args.output_tsv = [name.replace(".barcoded_reads.tsv", ".raw_barcodes.tsv")
                           for name in final_output_tsv]

    if args.threads == 1 or args.mode.enforces_single_thread():
        process_single_thread(args)
    else:
        process_in_parallel(args)

    if not args.whitelist_matching:
        correct_barcodes(args.output_tsv, final_output_tsv, args.barcodes,
                         barcode_length=get_barcode_length(args.mode),
                         n_cells=args.n_cells,
                         n_cells_interval=args.n_cells_interval,
                         threshold=args.barcode_graph_threshold,
                         rounds=args.barcode_graph_rounds,
                         threads=args.threads,
                         stats_file=args.output + ".barcode_correction.stats",
                         implementation=args.barcode_graph_impl)


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
