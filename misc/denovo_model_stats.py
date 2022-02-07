############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import os
import shutil
import random
import argparse
from traceback import print_exc
from common import *

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gtf_list", "-g", help="list of GTFs to process", type=str, required=True)
    parser.add_argument("--output", "-o", help="output prefix", type=str, required=True)
    args = parser.parse_args()

    return args


def process_tracking(gffcompare_output_prefix, gtf_num):
    reliable_transcripts = 0
    almost_reliable = [0] * gtf_num # supported by gtf_num - 2 tools
    unique_transcripts = [0] * gtf_num
    missed_transcripts = [0] * gtf_num
    total_transcripts = [0] * gtf_num

    for l in open(gffcompare_output_prefix + ".tracking", "r"):
        vals = l.strip().split()[4:]
        assert len(vals) == gtf_num
        equality_vector = [0 if vals[i] == '-' else 1 for i in range(gtf_num)]
        for i, e in enumerate(equality_vector):
            if e == 1:
                total_transcripts[i] += 1

        matches_transcripts = equality_vector.count(1)
        assert matches_transcripts > 0

        if matches_transcripts == gtf_num:
            reliable_transcripts += 1
        elif matches_transcripts == gtf_num - 1:
            for i, e in enumerate(equality_vector):
                if e == 1:
                    almost_reliable[i] += 1
            index = equality_vector.index(0)
            missed_transcripts[index] += 1
        elif matches_transcripts == 1:
            index = equality_vector.index(1)
            unique_transcripts[index] += 1

    return reliable_transcripts, almost_reliable, total_transcripts, unique_transcripts, missed_transcripts


def print_stats(name, reliable_transcripts, almost_reliable, total_transcripts, unique_transcripts, missed_transcripts, out=sys.stdout):
    out.write(name + "\n")
    total_columns = len(total_transcripts)
    out.write("\t".join([str(reliable_transcripts)] * total_columns) + "\n")
    out.write("\t".join([str(almost_reliable[i]) for i in range(total_columns)]) + "\n")
    out.write("\t".join([str(total_transcripts[i] - unique_transcripts[i] - almost_reliable[i] - reliable_transcripts)
                         for i in range(total_columns)]) + "\n")
    out.write("\t".join([str(unique_transcripts[i]) for i in range(total_columns)]) + "\n")
    out.write("\t".join([str(-missed_transcripts[i]) for i in range(total_columns)]) + "\n")
    out.flush()


def split_all(args):
    full_gtfs = []
    known_gtfs = []
    novel_gtfs = []
    for l in open(args.gtf_list):
        v = l.strip().split()

        if len(v) < 2:
            continue
        gtf_path = v[0]
        tool = v[1]
        if len(v) > 2:
            gtf_id = v[2]
        else:
            gtf_id = tool
        #print(gtf_path, gtf_id)

        out_full_path = os.path.join(args.output, gtf_id + ".full.gtf")
        out_known_path = os.path.join(args.output, gtf_id + ".known.gtf")
        out_novel_path = os.path.join(args.output, gtf_id + ".novel.gtf")
        if not os.path.exists(out_full_path):
            print("Seprating known and novel transcripts for " + gtf_id)
            separator = SEPARATE_FUNCTORS[tool](gtf_path)
            split_gtf(gtf_path, separator, out_full_path, out_known_path, out_novel_path)
        full_gtfs.append(out_full_path)
        known_gtfs.append(out_known_path)
        novel_gtfs.append(out_novel_path)

    return {"full" : full_gtfs, "known" : known_gtfs, "novel" : novel_gtfs}


def main():
    args = parse_args()
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    gtfs = split_all(args)
    #print(gtfs)
    out_stats = open(os.path.join(args.output, "all_stats.tsv"), "w")
    for name, gtf_list in gtfs.items():
        gff_compare_out = os.path.join(args.output, name)
        if not os.path.exists(gff_compare_out + ".tracking"):
            print("Running gffcompare for " + name)
            run_gff_compare_noref(gtf_list, gff_compare_out)
        reliable_transcripts, almost_reliable, total_transcripts, unique_transcripts, missed_transcripts = \
            process_tracking(gff_compare_out, len(gtf_list))
        print_stats(name, reliable_transcripts, almost_reliable, total_transcripts, unique_transcripts, missed_transcripts, out=out_stats)
    out_stats.close()


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
