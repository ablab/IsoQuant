import numpy
import sys
from collections import defaultdict
import gffutils
import re


def print_stats(read_intron_dict, max_len = 5):
    # sum up everythin above max_len
    for intron_count in sorted(read_intron_dict.keys()):
        if intron_count <= max_len:
            continue
        for k in read_intron_dict[intron_count].keys():
            read_intron_dict[max_len][k] += read_intron_dict[intron_count][k]
            read_intron_dict[intron_count][k] = 0

    print("introns\ttotal\t" + "\t".join(match_types))
    for intron_count in sorted(read_intron_dict.keys()):
        if intron_count > max_len:
            continue
        print("%d\t%s" % (
        intron_count, "\t".join([str(read_intron_dict[intron_count][x]) for x in ["total_reads"] + match_types])))

    print("introns\t" + "\t".join(match_types))
    for intron_count in sorted(read_intron_dict.keys()):
        if intron_count > max_len:
            continue
        total_reads = read_intron_dict[intron_count]["total_reads"]
        fractions = [(100 * read_intron_dict[intron_count][x] / total_reads) for x in match_types]
        print("%d\t%s" % (intron_count, "\t".join(["%.2f" % x for x in fractions])))


def get_isoform_chromosomes(gene_db):
    gffutils_db = gffutils.FeatureDB(gene_db, keep_order=True)
    isoform_map = {}
    for g in gffutils_db.features_of_type('transcript', order_by=('seqid')):
        isoform_map[g.id] = g.seqid
    return isoform_map


match_types = (["extra_intron_novel", "fake_terminal_exon_3", "fake_terminal_exon_5", "extra_intron_5", "extra_intron_3", "alternative_structure_novel"])

isoform_map = get_isoform_chromosomes(sys.argv[3])

supported_introns = set()
for l in open(sys.argv[2]):
    t = l.strip().split()
    supported_introns.add((t[0], int(t[1]), int(t[2])))


read_intron_dict = {}
confirmed_reads_introns = {}
for l in open(sys.argv[1]):
    t = l.strip().split()
    if t[2] != "inconsistent":
        continue
    read_id = t[0]
    isoform_id = t[3]
    chrid = isoform_map[isoform_id]
    intron_count = t[7].count(",")
    if intron_count not in read_intron_dict:
        read_intron_dict[intron_count] = defaultdict(int)
    if intron_count not in confirmed_reads_introns:
        confirmed_reads_introns[intron_count] = defaultdict(int)

    for mt in match_types:
        event_info = t[3]
        event_found = 0
        event_supported = 0
        for occ in re.finditer(mt, event_info):
            event_pos = occ.start()
            event_found = 1
            next_event = event_info.find(",", event_pos)
            if next_event == -1:
                intron_str = event_info[event_pos+len(mt)+1:]
            else:
                intron_str = event_info[event_pos + len(mt) + 1:next_event]
            coords = intron_str.split("-")
            intron = (chrid, int(coords[0]), int(coords[1]))
            if intron in supported_introns:
                event_supported = 1
        read_intron_dict[intron_count][mt] += event_found
        confirmed_reads_introns[intron_count][mt] += event_supported

    read_intron_dict[intron_count]["total_reads"] += 1
    confirmed_reads_introns[intron_count]["total_reads"] += 1

print("Inconsistency event distribution")
print_stats(read_intron_dict)
print("Supported inconsistency event distribution")
print_stats(confirmed_reads_introns)




