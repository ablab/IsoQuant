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
    parser.add_argument("--tracking", "-t", help="gffcompare tracking file", type=str, required=True)
    parser.add_argument("--output", "-o", help="output prefix", type=str, required=True)
    args = parser.parse_args()

    return args


def filter_gtf(ingtf_path, isoform_ids, out_full_path):
    print("Filtering %d transcripts from %s to %s" % (len(isoform_ids), ingtf_path, out_full_path))
    #print(list(isoform_ids)[:10])
    out_full = open(out_full_path, "w")
    count = 0
    for l in open(ingtf_path):
        if l.startswith("#"):
            continue

        tpos = l.find(' transcript_id')
        if tpos == -1:
            continue
        idpos = tpos + len(' transcript_id') + 2
        endpos = l.find(";", idpos)
        
        if endpos == -1:
            print("Warning, unable to find ;")
        tid = l[idpos:endpos-1]
        count += 1
        if tid not in isoform_ids:
            continue
        out_full.write(l)
    out_full.close()


def process_tracking(gtfs, tracking, outdir):
    gtf_num = len(gtfs)
    unique_transcripts = [set() for i in range(gtf_num)]
    missed_transcripts = [set() for i in range(gtf_num)]

    for l in open(tracking, "r"):
        vals = l.strip().split()[4:]
        assert len(vals) == gtf_num
        equality_vector = [0 if vals[i] == '-' else 1 for i in range(gtf_num)]
        matches_transcripts = equality_vector.count(1)
        assert matches_transcripts > 0

        if matches_transcripts == 1:
            index = equality_vector.index(1)
            unique_transcripts[index].add(vals[index].split("|")[1])
        elif matches_transcripts == gtf_num - 1:
            index = equality_vector.index(0)
            if index != gtf_num - 1:
                missed_transcripts[index].add(vals[-1].split("|")[1])
            else:
                missed_transcripts[index].add(vals[-2].split("|")[1])

    for i in range(gtf_num):
        filter_gtf(gtfs[i][0], unique_transcripts[i], os.path.join(outdir, gtfs[i][1] + ".unique.gtf"))
        if i != gtf_num - 1:
            filter_gtf(gtfs[-1][0], missed_transcripts[i], os.path.join(outdir, gtfs[i][1] + ".missed.gtf"))
        else:
            filter_gtf(gtfs[-2][0], missed_transcripts[i], os.path.join(outdir, gtfs[i][1] + ".missed.gtf"))


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


def main():
    args = parse_args()
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    gtfs = []
    for l in open(args.gtf_list):
        v = l.strip().split()

        if len(v) < 2:
            continue
        gtfs.append((v[0], v[1]))

    process_tracking(gtfs, args.tracking, args.output)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
