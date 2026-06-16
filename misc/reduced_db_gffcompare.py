#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2020-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
from traceback import print_exc

from common import *


# Terminal-end tolerances (bp) for transcript-level matching. None = default
# end-agnostic match; the integers use the gffcompare fork's --terminal-delta.
TERMINAL_DELTAS = [None, 50, 10]


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output folder", default="gtf_stats")
    parser.add_argument("--genedb", "-d", type=str, help="prefix to reduced gene db")
    parser.add_argument("--gtf", "-g", type=str, help="gtf to assess")
    parser.add_argument("--tool", type=str, choices=SEPARATE_FUNCTORS.keys(),
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

    # (split, reference subset, query split GTF)
    splits = [
        ("full",  args.genedb + ".expressed.gtf",      out_full_path),
        ("known", args.genedb + ".expressed_kept.gtf", out_known_path),
        ("novel", args.genedb + ".excluded.gtf",       out_novel_path),
    ]
    # Score each split at several terminal-end tolerances. None = default
    # end-agnostic transcript match (-> "<tool>.<split>.stats"); the integer
    # deltas use the gffcompare fork's --terminal-delta so the transcript-level
    # metric becomes end-sensitive (-> "<tool>.<split>.td<delta>.stats").
    # See .claude/GFFCOMPARE.md. Requires the gffcompare fork for the deltas;
    # the default run works with stock gffcompare too.
    for split, reference_gtf, compared_gtf in splits:
        for delta in TERMINAL_DELTAS:
            suffix = "" if delta is None else (".td%d" % delta)
            option = "" if delta is None else ("--terminal-delta=%d" % delta)
            out_stats = os.path.join(args.output, "%s.%s%s.stats" % (args.tool, split, suffix))
            print("Running gffcompare for %s transcripts (terminal-delta=%s)" % (split, delta))
            run_gff_compare(reference_gtf, compared_gtf, out_stats, additional_option=option)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)

