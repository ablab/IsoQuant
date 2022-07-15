#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import subprocess
import sys
import argparse
from collections import defaultdict
from traceback import print_exc
from collections import namedtuple

import pysam
import gffutils
from enum import Enum, unique
from common import *


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output folder", default="gtf_stats")
    parser.add_argument("--genedb", "-d", type=str, help="annotation for assessment")
    parser.add_argument("--gtf", "-g", type=str, help="gtf to assess")
    parser.add_argument("--tool", type=str, choices=['isoquant', 'talon', 'sqanti', 'flair', 'bambu', 'stringtie'],
                        help="tool used for generating GTF")

    args = parser.parse_args()
    if not args.genedb or not args.gtf or not args.tool:
        parser.print_usage()
        exit(-1)
    return args


def main():
    args = parse_args()
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    out_full_path = os.path.join(args.output, args.tool + ".full.gtf")
    out_known_path = os.path.join(args.output, args.tool + ".known.gtf")
    out_novel_path= os.path.join(args.output, args.tool + ".novel.gtf")
    print("Seprating known and novel transcripts")
    separator = SEPARATE_FUNCTORS[args.tool](args.gtf)
    split_gtf(args.gtf, separator, out_full_path, out_known_path, out_novel_path)
    print("Running gffcompare for entire GTF")
    expressed_gtf = args.genedb
    run_gff_compare(expressed_gtf, out_full_path, os.path.join(args.output, args.tool + ".full.stats"), additional_option="-R")
    print("Running gffcompare for known transcripts")
    run_gff_compare(expressed_gtf, out_known_path, os.path.join(args.output, args.tool + ".known.stats"), additional_option="-R")
    print("Running gffcompare for novel transcripts")
    run_gff_compare(expressed_gtf, out_novel_path, os.path.join(args.output, args.tool + ".novel.stats"), additional_option="-R")


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)

