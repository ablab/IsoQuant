import sys
from collections import defaultdict
import pysam
from Bio import Seq
from Bio import SeqIO

bam_file_pb = sys.argv[1]
bam_file_ont = sys.argv[2]
tsv_file = sys.argv[3]


def get_sequence_to_check(alignment):
    if alignment.seq is None or alignment.cigartuples is None:
        return -1

    cigar_tuples = alignment.cigartuples
    clipped_len = -1
    sequence_to_check = ''
    if len(cigar_tuples) > 1 and cigar_tuples[0][0] == 5 and cigar_tuples[1][0] == 4:
        # hard clipped
        clipped_len = cigar_tuples[1][1]
    elif cigar_tuples[0][0] == 4:
        # soft clipped
        clipped_len = cigar_tuples[0][1]
    if (clipped_len != -1):
        sequence_to_check = alignment.seq[
                            :clipped_len]  # str(Seq.Seq(alignment.seq[:clipped_len]).reverse_complement()).upper()

    clipped_len = -1

    if len(cigar_tuples) > 1 and cigar_tuples[-1][0] == 5 and cigar_tuples[-2][0] == 4:
        # hard clipped
        clipped_len = cigar_tuples[-2][1]
    elif cigar_tuples[-1][0] == 4:
        # soft clipped
        clipped_len = cigar_tuples[-1][1]

    sequence_to_check_end = ''
    if (clipped_len != -1):
        sequence_to_check_end = str((alignment.seq[-clipped_len:]).upper())

    return sequence_to_check, sequence_to_check_end


def get_direction(alignment):
    if alignment.seq is None or alignment.cigartuples is None:
        return -1

    sequence_to_check_start, sequence_to_check_end = get_sequence_to_check(alignment)

    sequence_to_check_start = str(Seq.Seq(sequence_to_check_start).reverse_complement()).upper()

    if check_polya(sequence_to_check_end) != -1:
        return 1
    elif check_polya(sequence_to_check_start) != -1:
        return -1

    return 0


def check_polya(sequence):
    len_slice = 12
    n = len(sequence)
    if n < len_slice:
        return -1

    pos = 12
    for i in range(n - len_slice):
        slice_ = sequence[i: i + len_slice]
        k = slice_.count('A')

        if (k >= 0.75 * len(slice_)):
            return i + slice_.find("AA")

    return -1


def get_tss_polya(bam_file):
    read_tss = defaultdict(int)
    read_polya = defaultdict(int)
    polya_dist = []
    tss_dist = []
    wrong_reads = set()
    for r in pysam.AlignmentFile(bam_file, "rb").fetch():
        if r.is_supplementary or r.is_secondary or r.reference_id == -1 or not r.cigartuples:
            continue
        strand = get_direction(r)
        if not strand: continue
        intron_count = sum([1 if x[0] == 3 else 0 for x in r.cigartuples])
        #if intron_count:
        #    wrong_reads.add(r.query_name)
        #    continue

        min_dist_up = 1000000000
        min_tss_up = None

        min_dist_down = 1000000000
        min_tss_down = None
        strand = "+" if strand > 0 else "-"
        tss_read_coord = r.reference_start if strand == "+" else r.reference_end
        polya_read_coord = r.reference_start if strand == "-" else r.reference_end
        for i in range(tss_read_coord - max_dist, tss_read_coord + max_dist + 1):
            tss_key = r.reference_name + "_" + str(i) + "_" + strand
            if tss_key in tss_dict:
                if strand == "+" and i <= tss_read_coord and abs(tss_read_coord - i) < min_dist_up:
                    min_dist_up = tss_read_coord - i
                    min_tss_up = tss_dict[tss_key]
                elif strand == "+" and i >= tss_read_coord and abs(tss_read_coord - i) < abs(min_dist_down):
                    min_dist_down = tss_read_coord - i
                    min_tss_down = tss_dict[tss_key]
                elif strand == "-" and i >= tss_read_coord and abs(tss_read_coord - i) < min_dist_up:
                    min_dist_up = tss_read_coord - i
                    min_tss_up = tss_dict[tss_key]
                elif strand == "-" and i <= tss_read_coord and abs(tss_read_coord - i) < abs(min_dist_down):
                    min_dist_down = tss_read_coord - i
                    min_tss_down = tss_dict[tss_key]
        if abs(min_dist_up) < abs(min_dist_down):
            if min_dist_up < min_dist_down:
                read_tss[r.query_name] = min_tss_up
                tss_dist.append(min_dist_up)
            else:
                read_tss[r.query_name] = min_tss_down
                tss_dist.append(min_dist_down)
        min_dist_up = 1000000000
        min_polya_up = None
        min_dist_down = 1000000000
        min_polya_down = None
        for i in range(polya_read_coord - max_dist, polya_read_coord + max_dist + 1):
            polya_key = r.reference_name + "_" + str(i) + "_" + strand
            if polya_key in polya_dict:
                if strand == "+" and i <= polya_read_coord and abs(polya_read_coord - i) < abs(min_dist_up):
                    min_dist_up = polya_read_coord - i
                    min_polya_up = polya_dict[polya_key]
                elif strand == "+" and i >= polya_read_coord and abs(polya_read_coord - i) < abs(min_dist_down):
                    min_dist_down = polya_read_coord - i
                    min_polya_down = polya_dict[polya_key]
                elif strand == "-" and i >= polya_read_coord and abs(polya_read_coord - i) < abs(min_dist_up):
                    min_dist_up = polya_read_coord - i
                    min_polya_up = polya_dict[polya_key]
                elif strand == "-" and i <= polya_read_coord and abs(polya_read_coord - i) < abs(min_dist_down):
                    min_dist_down = polya_read_coord - i
                    min_polya_down = polya_dict[polya_key]
        if not min_polya_up and not min_polya_down: continue
        if abs(min_dist_up) < abs(min_dist_down):
            read_polya[r.query_name] = min_polya_up
            polya_dist.append(min_dist_up)
        else:
            read_polya[r.query_name] = min_polya_down
            polya_dist.append(min_dist_down)
    return read_tss, read_polya,wrong_reads

cage_file = "data/mm10_fair+new_CAGE_peaks_phase1and2.bed"
polya_file = "data/atlas.clusters_chr.mm10.2-0.bed"

cage_cutoff = 1
polya_cutoff = 0.0001

tss_dict = {}
with open(cage_file) as f:
    for line in f:
        fs = line.split()
        if int(fs[4]) < cage_cutoff: continue
        cage_pos = int(fs[1]) + 1
        cage_peak = "%s_%d_%s_%s" % (fs[0], cage_pos, fs[2], fs[5])
        for i in range(cage_pos, int(fs[2]) + 1):
            tss_key = "%s_%d_%s" % (fs[0], i, fs[5])
            tss_dict[tss_key] = cage_peak

polya_dict = {}
with open(polya_file) as f:
    for line in f:
        fs = line.split()
        if float(fs[4]) < polya_cutoff: continue
        polya_pos = int(fs[1]) + 1
        polya_peak = "%s_%d_%s_%s" % (fs[0], polya_pos, fs[2], fs[5])
        for i in range(polya_pos, int(fs[2]) + 1):
            polya_key = "%s_%d_%s" % (fs[0], i, fs[5])
            polya_dict[polya_key] = polya_peak

max_dist = 100
read_tss_pb, read_polya_pb,wrong_reads_pb = get_tss_polya(bam_file_pb)
read_tss_ont, read_polya_ont,wrong_reads_ont = get_tss_polya(bam_file_ont)

tss_agree, tss_pb, tss_ont, tss_disagree = 0, 0, 0, 0
polya_agree, polya_pb, polya_ont, polya_disagree = 0, 0, 0, 0
tss_assigned_pb,tss_assigned_ont = 0,0
polya_assigned_pb,polya_assigned_ont = 0,0
total_tss_pairs = 0
total_polya_pairs = 0
with open(tsv_file) as f:
    for line in f:
        _, _, pb_read, ont_read = line.split()
        if pb_read in wrong_reads_pb or ont_read in wrong_reads_ont: continue
        if pb_read in read_tss_pb or ont_read in read_tss_ont:
            total_tss_pairs +=1
            if read_tss_pb[pb_read]: tss_assigned_pb+=1
            if read_tss_ont[ont_read]: tss_assigned_ont+=1
            if not read_tss_ont[ont_read]:
                tss_pb += 1
            elif not read_tss_pb[pb_read]:
                tss_ont += 1
            elif read_tss_ont[ont_read] and read_tss_pb[pb_read] == read_tss_ont[ont_read]:
                tss_agree += 1
            elif read_tss_ont[ont_read] and read_tss_pb[pb_read] != read_tss_ont[ont_read]:
                tss_disagree += 1
        if pb_read in read_polya_pb or ont_read in read_polya_ont:
            total_polya_pairs+=1
            if read_polya_pb[pb_read]: polya_assigned_pb+=1
            if read_polya_ont[ont_read]: polya_assigned_ont+=1
            if not read_polya_ont[ont_read]:
                polya_pb += 1
            elif not read_polya_pb[pb_read]:
                polya_ont += 1
            elif read_polya_ont[ont_read] and read_polya_pb[pb_read] == read_polya_ont[ont_read]:
                polya_agree += 1
            elif read_polya_ont[ont_read] and read_polya_pb[pb_read] != read_polya_ont[ont_read]:
                polya_disagree += 1
#print(total_tss_pairs, tss_assigned_pb, tss_assigned_ont)
#print(total_polya_pairs, polya_assigned_pb, polya_assigned_ont)
print("TSS. Agree: ", tss_agree, ", PB only: ", tss_pb,  ", ONT only: ", tss_ont, ", disagree: ", tss_disagree)
print("PolyA. Agree: ", polya_agree,  ", PB only: ",polya_pb,  ", ONT only: ", polya_ont,  ", disagree: ", polya_disagree)
