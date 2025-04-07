#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import argparse
from traceback import print_exc
import logging

from umi_filtering import UMIFilter, filter_bam, load_barcodes, create_transcript_info_dict


logger = logging.getLogger('IsoQuant')


def set_logger(logger_instance):
    logger_instance.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", required=True)
    parser.add_argument("--barcodes", "-b", nargs='+', type=str, help="read - barcode - UMI table", required=True)
    parser.add_argument("--read_assignments", "-r", nargs='+', type=str, help="IsoQuant read assignments", required=True)
    parser.add_argument("--genedb", "-g", help="gene database in gffutils DB format (optional)", type=str)
    parser.add_argument("--bam", type=str, help="original BAM file, provide only if you need a BAM file"
                                                " with UMI-filtered alignments")

    parser.add_argument("--min_distance", type=int, help="minimal edit distance for UMIs to be considered distinct;"
                                                         "length difference is added to this values by default;"
                                                         "0 for equal UMIs, -1 for keeping only a single gene-barcode "
                                                         "read. By default will be process with -1, 2, 3", default=4)
    parser.add_argument("--untrusted_umis", action="store_true", help="allow untrusted UMIs to be used", default=False)
    parser.add_argument("--only_spliced", action="store_true", help="keep only spliced reads", default=False)
    parser.add_argument("--only_unique", action="store_true", help="keep only non-ambiguous reads", default=False)
    parser.add_argument("--disregard_length_diff", action="store_true", help="do not account for length difference "
                                                                             "when comparing UMIs", default=False)
    parser.add_argument("--downsample_fraction", type=float, help="downsample reads with specified fraction (0, 1)")
    parser.add_argument("--old_barcode_format", action="store_true",
                        help="use previous barcode formate (barcode column = 6)", default=False)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    set_logger(logger)

    if args.old_barcode_format:
        barcode_umi_dict = load_barcodes(args.barcodes, args.untrusted_umis, 6, 7, 8, 9)
    else:
        barcode_umi_dict = load_barcodes(args.barcodes, args.untrusted_umis)

    if args.genedb:
        transcript_info_dict = create_transcript_info_dict(args.genedb)
    else:
        transcript_info_dict = {}

    for d in sorted({4, args.min_distance}):
        logger.info("== Filtering by UMIs with edit distance %d ==" % d)
        output_prefix = args.output + (".ALL" if d < 0 else "ED%d" % d)
        logger.info("Results will be saved to %s" % output_prefix)
        umi_filter = UMIFilter(barcode_umi_dict, d, args.disregard_length_diff,
                               args.only_unique, args.only_spliced)
        umi_filter.process(args.read_assignments, output_prefix, transcript_info_dict, args.downsample_fraction)
        if args.bam:
            filter_bam(args.bam, output_prefix  + ".UMI_filtered.reads.bam", umi_filter.selected_reads)
        logger.info("== Done filtering by UMIs with edit distance %d ==" % d)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
