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
    parser.add_argument("--counts", "-c", type=str, required=True, help="TSV with isoform counts: id, count, tpm")

    args = parser.parse_args()
    return args


def load_counts(counts_file):
    count_dict = {}
    for l in open(counts_file):
        if l.startswith("#"):
            continue
        t = l.strip().split()
        if len(t) < 3:
            continue
        count_dict[t[0]] = float(t[2])
    print("Loaded %d counts" % len(count_dict))
    return count_dict


def load_transcripts(gtf_file):
    transcripts = set()
    for l in open(gtf_file):
        if l.startswith("#"):
            continue
        t = l.strip().split()
        if len(t) < 9:
            continue
        if t[2] not in {'transcript', 'mRNA'}:
            continue

        for i, v in enumerate(t):
            if v == 'transcript_id':
                transcripts.add(t[i+1][1:-2])
                break
    print("Loaded %d reference transcripts" % len(transcripts))
    return transcripts


def count_stats(tracking_file, count_dict, transcript_set, bins):
    transcript_hist = defaultdict(int)
    processed_transcripts = set()
    for l in open(tracking_file):
        t = l.strip().split()
        if t[2] == '-' or t[3] != '=':
            continue
        transcript_id = t[2].split('|')[1]
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
    count_dict = load_counts(args.counts)
    transcript_set = load_transcripts(args.genedb)
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

