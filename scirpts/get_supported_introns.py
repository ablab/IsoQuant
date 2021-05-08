import sys
import gffutils
import argparse
from traceback import print_exc
from collections import defaultdict
import re

match_types = (["extra_intron_novel", "fake_terminal_exon_3", "fake_terminal_exon_5", "extra_intron_5", "extra_intron_3", "alternative_structure_novel"])


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--assignments1", help="read assignment tsv", type=str)
    parser.add_argument("--assignments2", help="read assignment tsv", type=str) ## SQANTI2 output
    parser.add_argument("--cutoff", "-c", help="min number of reads when single tsv specified", type=int, default=3)
    parser.add_argument("--gene_db", "-g", help="gene database", type=str)
    parser.add_argument("--output", "-o", help="gene database", type=str, default="supported_introns.tsv")

    args = parser.parse_args()

    if not args.gene_db or not args.assignments1:
        parser.print_usage()
        exit(-1)
    return args


def get_isoform_chromosomes(gene_db):
    print("Loading gene db from " + gene_db)
    gffutils_db = gffutils.FeatureDB(gene_db, keep_order=True)
    isoform_map = {}
    for g in gffutils_db.features_of_type('transcript', order_by=('seqid')):
        isoform_map[g.id] = g.seqid
    return isoform_map


def load_tsv(read_assignments, isoform_map):
    print("Loading assignments from " + read_assignments)
    support_map = defaultdict(set)
    for l in open(read_assignments):
        t = l.strip().split()
        if t[5] != "inconsistent":
            continue

        read_id = t[0]
        isoform_id = t[3]
        chrid = isoform_map[isoform_id]

        for mt in match_types:
            event_info = t[6]
            for occ in re.finditer(mt, event_info):
                event_pos = occ.start()
                next_event = event_info.find(",", event_pos)
                if next_event == -1:
                    intron_str = event_info[event_pos + len(mt) + 1:]
                else:
                    intron_str = event_info[event_pos + len(mt) + 1:next_event]
                coords = intron_str.split("-")
                intron = (chrid, int(coords[0]), int(coords[1]))
                support_map[intron].add(read_id)
    return support_map


def main():
    args = parse_args()
    isoform_map = get_isoform_chromosomes(args.gene_db)
    support_map = load_tsv(args.assignments1, isoform_map)

    outf = open(args.output, 'w')
    if args.assignments2:
        support_map2 = load_tsv(args.assignments2, isoform_map)
        print("Saving shared introns to " + args.output)
        for intron in support_map.keys():
            if intron in support_map2:
                outf.write("%s\t%d\t%d\n" % (intron[0], intron[1], intron[2]))
    else:
        print("Saving supported introns to " + args.output)
        for intron in support_map.keys():
            if len(support_map[intron]) >= args.cutoff:
                outf.write("%s\t%d\t%d\n" % (intron[0], intron[1], intron[2]))
    outf.close()
    print("Done")


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)