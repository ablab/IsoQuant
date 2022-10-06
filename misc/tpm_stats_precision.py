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
from enum import Enum, unique
from common import *


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genedb", "-g", type=str, required=True, help="prefix to reference gene db")
    parser.add_argument("--tracking", "-t", type=str, required=True, help="gffcompare tracking wrt to the given annotation")
    parser.add_argument("--counts", "-c", type=str, help="TSV with isoform counts: id, count, tpm")
    parser.add_argument("--tool", type=str, required=True, choices=['stringtie', 'isoquant', 'talon', 'flair', 'bambu'],
                        help="tool used for generating annotation")

    args = parser.parse_args()
    return args


def default_split_func(tid):
    return tid


def flair_split(tid):
    new_tid = "_".join(tid.split("_")[:-1])
    return new_tid


def convert_dict_to_tpm(count_dict):
    scale_factor = 1000000.0 / sum(count_dict.values())
    for k in count_dict.keys():
        count_dict[k] *= scale_factor
    return  count_dict


def load_counts(counts_file, tid_column=1, count_column=3, convert_to_tpm=True, split_func=default_split_func):
    count_dict = {}
    for l in open(counts_file):
        if l.startswith("#") or l.startswith('TXNAME') or l.startswith('feature_id') or l.startswith('gene_ID') or l.startswith('ids'):
            continue
        t = l.strip().split()
        if len(t) < count_column:
            continue
        tid = split_func(t[tid_column-1])
        count_dict[tid] = float(t[count_column - 1])
    print("Loaded %d counts" % len(count_dict))
    if convert_to_tpm:
        return convert_dict_to_tpm(count_dict)
    return count_dict


def load_transcripts(gtf_file, load_tpm=False):
    transcripts = set()
    count_dict = {}
    for l in open(gtf_file):
        if l.startswith("#"):
            continue
        t = l.strip().split()
        if len(t) < 9:
            continue
        if t[2] not in {'transcript', 'mRNA'}:
            continue

        transcript = None
        tpm = 0.0
        for i, v in enumerate(t):
            if v == 'transcript_id':
                transcript = t[i+1][1:-2]
            elif load_tpm and v == 'TPM':
                tpm = float(t[i+1][1:-2])
        if transcript:
            transcripts.add(transcript)
            count_dict[transcript] = tpm

    print("Loaded %d reference transcripts" % len(transcripts))
    return transcripts, count_dict


def count_stats(tracking_file, count_dict, transcript_set, bins):
    transcript_hist = defaultdict(int)
    processed_transcripts = set()
    for l in open(tracking_file):
        t = l.strip().split()
        if t[2] == '-' or t[3] != '=' or t[4] == '-':
            continue
        transcript_id = t[4].split('|')[1]
        if transcript_id not in transcript_set:
            continue
        if transcript_id not in count_dict:
            sys.stderr.write("Counts for transcript %s is not found\n" % transcript_id)
            continue
        if transcript_id in processed_transcripts:
            continue

        tpm = count_dict[transcript_id]
        for b in reversed(bins):
            if tpm >= b:
                transcript_hist[b] += 1
                processed_transcripts.add(transcript_id)
                break
    return transcript_hist


def count_total(count_dict, transcript_set, bins):
    transcript_hist = defaultdict(int)
    for transcript_id in transcript_set:
        if transcript_id not in count_dict:
            sys.stderr.write("Counts for reference transcript %s is not found\n" % transcript_id)
            continue

        tpm = count_dict[transcript_id]
        for b in reversed(bins):
            if tpm >= b:
                transcript_hist[b] += 1
                break
    return transcript_hist


def main():
    args = parse_args()
    bins = [0, 5, 10, 50, 100]
    split_funct = flair_split if args.tool == 'flair' else default_split_func

    if args.tool == 'stringtie':
        transcript_set, count_dict = load_transcripts(args.genedb, load_tpm=True)
    else:
        transcript_set, count_dict = load_transcripts(args.genedb, load_tpm=False)
        normalize = args.tool in {'bambu', 'talon'}
        # 1 base column index
        if  args.tool in {'isoquant', 'flair'}:
            count_column = 2
            tid_column = 1
        elif args.tool in {'talon'}:
            count_column = 12
            tid_column = 4
        else:
            count_column = 3
            tid_column = 1
        split_funct = flair_split if args.tool == 'flair' else default_split_func

        count_dict = load_counts(args.counts, tid_column, count_column, convert_to_tpm=normalize, split_func=split_funct)

    total_hist = count_total(count_dict, transcript_set, bins)
    restored_hist = count_stats(args.tracking, count_dict, transcript_set, bins)

    print(args.tracking)
    print("\t".join(map(str, bins)))
    print("\t".join(map(str, [restored_hist[b] for b in bins])) + "\t%d" % sum(restored_hist.values()))
    print("\t".join(map(str, [total_hist[b] for b in bins])) + "\t%d" % sum(total_hist.values()))
    print("\t".join(["%.2f" % (restored_hist[b] / total_hist[b]) for b in bins]) + "\t%.2f" % (sum(restored_hist.values()) / sum(total_hist.values())))
    print(" ")


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)

