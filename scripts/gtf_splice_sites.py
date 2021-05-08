import sys
import gffutils
from Bio import SeqIO
from collections import defaultdict


def junctions_from_blocks(sorted_blocks):
    junctions = []
    if len(sorted_blocks) >= 2:
        for i in range(0, len(sorted_blocks) - 1):
            if sorted_blocks[i][1] + 1 < sorted_blocks[i + 1][0]:
                junctions.append((sorted_blocks[i][1] + 1, sorted_blocks[i + 1][0] - 1))
    return junctions


def get_all_introns_sites(gene_db):
    print("Loading gene db from " + gene_db)
    gffutils_db = gffutils.FeatureDB(gene_db, keep_order=True)
    transcript_introns = defaultdict(list)
    for t in gffutils_db.features_of_type('transcript', order_by=('seqid', 'start')):
        exon_list = []
        for e in gffutils_db.children(t, order_by='start'):
            if e.featuretype == 'exon':
                exon_list.append((e.start, e.end))
        introns = junctions_from_blocks(exon_list)
        for intron in introns:
            transcript_introns[(t.seqid, t.strand)].append(intron)
    return transcript_introns


def check_canonical(intron, ref_region, strand=None):
    intron_left_pos = intron[0]
    intron_right_pos = intron[1]
    left_site = ref_region[intron_left_pos:intron_left_pos + 2]
    right_site = ref_region[intron_right_pos - 1:intron_right_pos + 1]
    if (left_site == "GT" and right_site == "AG") and strand != '-':
        return '+'
    elif (left_site == "CT" and right_site == "AC") and strand != '+':
        return '-'
    else:
        return None


def analyse_intron_sites(intron, ref_region, strand):
    seq_size = 10
    intron_left_pos = intron[0]
    intron_right_pos = intron[1]

    if strand not in ['+', '-']:
        return None, None, None, None

    left_upper = ref_region[intron_left_pos - seq_size:intron_left_pos]
    left_lower = ref_region[intron_left_pos + 2:intron_left_pos + seq_size + 2]
    right_upper = ref_region[intron_right_pos - seq_size - 1:intron_right_pos - 1]
    right_lower = ref_region[intron_right_pos + 1:intron_right_pos + seq_size + 1]

    # upstream and downstream here are relative to the genome
    if strand == "+":
        donor_upstream = left_upper.rfind("GT")
        donor_downstream = left_lower.find("GT")
        acc_upstream = right_upper.rfind("AG")
        acc_downstream = right_lower.find("AG")
    else:
        acc_upstream = left_upper.rfind("CT")
        acc_downstream = left_lower.find("CT")
        donor_upstream = right_upper.rfind("AC")
        donor_downstream = right_lower.find("AC")

    donor_upstream = seq_size - donor_upstream if donor_upstream != -1 else 0
    donor_downstream = 2 + donor_downstream if donor_downstream != -1 else 0
    acc_upstream = seq_size - acc_upstream if acc_upstream != -1 else 0
    acc_downstream = 2 + acc_downstream if acc_downstream != -1 else 0

    if strand == '+':
        return -donor_upstream, donor_downstream, -acc_upstream, acc_downstream
    else:
        return -donor_downstream, donor_upstream, -acc_downstream, acc_upstream


def count_splice_site_stats(ref_records, intron_map, allow_repeats=True):
    acceptor_set = set()
    donor_set = set()
    donor_counts = defaultdict(int)
    acc_counts = defaultdict(int)
    for chr_info in sorted(intron_map.keys()):
        chr_id = chr_info[0]
        annotation_strand = chr_info[1]
        chr_seq = "N" + str(ref_records[chr_id].seq)
        for intron in intron_map[chr_info]:
            strand = check_canonical(intron, chr_seq, annotation_strand)
            if not strand:
                continue
            donor_up, donor_down, acceptor_up, acceptor_down = analyse_intron_sites(intron, chr_seq, strand)
            donor_pos = intron[0] if strand == "+" else intron[1]
            acc_pos = intron[0] if strand == "-" else intron[1]
            if allow_repeats or (chr_id, strand, donor_pos) not in donor_set:
                if donor_up != 0 and donor_down != 0:
                    donor_counts[donor_up] += 1
                    donor_counts[donor_down] += 1
                elif donor_up != 0:
                    donor_counts[donor_up] += 1
                elif donor_down != 0:
                    donor_counts[donor_down] += 1
                else:
                    donor_counts[0] += 1
                donor_set.add((chr_id, strand, donor_pos))
            if allow_repeats or (chr_id, strand, acc_pos) not in acceptor_set:
                if acceptor_up != 0 and acceptor_down != 0:
                    acc_counts[acceptor_up] += 1
                    acc_counts[acceptor_down] += 1
                elif acceptor_up != 0:
                    acc_counts[acceptor_up] += 1
                if acceptor_down != 0:
                    acc_counts[acceptor_down] += 1
                else:
                    acc_counts[0] += 1
                acceptor_set.add((chr_id, strand, acc_pos))
    return donor_counts, acc_counts


def print_stats(donor_counts, acc_counts):
    shifts = [x for x in range(-10, 11)]
    print("shift\ttotal\t%s" % "\t".join(str(x) for x in shifts))
    print("d\t%d\t%s" % (sum(donor_counts.values()), "\t".join(str(donor_counts[x]) for x in shifts)))
    print("a\t%d\t%s" % (sum(acc_counts.values()), "\t".join(str(acc_counts[x]) for x in shifts)))


sys.stderr.write("Loading gene db from %s\n" % sys.argv[1])
transcript_introns = get_all_introns_sites(sys.argv[1])
sys.stderr.write("Loading genome from %s\n" % sys.argv[2])
reference_record_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[2], "fasta"))
sys.stderr.write("Counting splice sites\n")
donor_counts, acc_counts = count_splice_site_stats(reference_record_dict, transcript_introns)
print("= All transcripts =")
print_stats(donor_counts, acc_counts)

unique_introns = {}
for k in transcript_introns.keys():
    unique_introns[k] = sorted(set(transcript_introns[k]))
donor_counts, acc_counts = count_splice_site_stats(reference_record_dict, unique_introns)
print("= All introns =")
print_stats(donor_counts, acc_counts)

donor_counts, acc_counts = count_splice_site_stats(reference_record_dict, unique_introns, allow_repeats=False)
print("= All splice sites =")
print_stats(donor_counts, acc_counts)

