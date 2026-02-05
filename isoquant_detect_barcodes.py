#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2023-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Standalone barcode detection tool for single-cell and spatial transcriptomics.

This is a CLI wrapper around src.barcode_calling.detect_barcodes.
"""

import argparse
import logging
import os
import sys
from traceback import print_exc

from src.modes import IsoQuantMode
from src.error_codes import IsoQuantExitCode
from src.barcode_calling.detect_barcodes import (
    process_single_thread,
    process_in_parallel,
    get_umi_length,
    BARCODE_CALLING_MODES,
)

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

    parser.add_argument("--input", "-i", nargs='+', type=str, help="input reads in [gzipped] FASTA, FASTQ, BAM, SAM",
                        required=True)
    parser.add_argument("--threads", "-t", type=int, help="threads to use (16)", default=16)
    parser.add_argument("--tmp_dir", type=str, help="folder for temporary files")
    add_hidden_option('--debug', action='store_true', default=False, help='Debug log output.')

    args = parser.parse_args(sys_argv)
    args.mode = IsoQuantMode[args.mode]
    args.out_fasta = None
    args.output_tsv = None
    return args


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

    out_dir = os.path.dirname(args.output)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    if args.threads == 1 or args.mode.enforces_single_thread():
        process_single_thread(args)
    else:
        process_in_parallel(args)


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(IsoQuantExitCode.UNCAUGHT_EXCEPTION)
