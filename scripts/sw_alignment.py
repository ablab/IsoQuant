#!/usr/bin/end python
import sys
import numpy as np
from Bio import SeqIO, Seq


from ssw import (
    SSW,
    force_align,
    format_force_align
)

def rev_comp(seq):
    c = dict(zip('ATCGNatcgn', 'TAGCNtagcn'))
    return ''.join(c.get(nucleotide, '') for nucleotide in reversed(seq))

def check_polya(sequence):
    len_slice = 12
    n = len(sequence)
    if n < len_slice:
        return -1
    for i in range(n - len_slice):
        slice_ = sequence[i: i + len_slice]
        k = slice_.count('A')

        if (k >= 0.75 * len(slice_)):
            return i + slice_.find("AA")

    return -1


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


def is_reverse(sequence_to_check_start, sequence_to_check_end):
    sequence_to_check_start = str(Seq.Seq(sequence_to_check_start).reverse_complement()).upper()

    if check_polya(sequence_to_check_start) != -1:
        return 1
    elif check_polya(sequence_to_check_end) != -1:
        return 0

    return -1

bins = 20
all_qualities = [[] for i in range(bins*3)]

fastq_file_pb = sys.argv[1]
fastq_file_ont = sys.argv[2]
tsv_file = sys.argv[3]

pb_seqs = SeqIO.parse(fastq_file_pb, format="fastq")
ont_seqs = SeqIO.parse(fastq_file_ont, format="fastq")

reverse_reads = set()
pb_qual = {seq.id: (str(seq.seq), seq.letter_annotations["phred_quality"])  for seq in pb_seqs}
ont_qual = {seq.id: (str(seq.seq), seq.letter_annotations["phred_quality"]) for seq in ont_seqs}
all_qualities_ont = [[] for i in range(bins*3)]

a = SSW()
with open(tsv_file) as f:
    for line in f:
        _, _, pb_read, ont_read = line.split()
        if pb_read not in pb_qual or ont_read not in ont_qual: continue
        pb_s, ont_s = pb_qual[pb_read][0], ont_qual[ont_read][0]
        pb_q, ont_q = pb_qual[pb_read][1], ont_qual[ont_read][1]
        is_rev1 = is_reverse(pb_s[:100], pb_s[-100:])
        is_rev2 = is_reverse(ont_s[:100], ont_s[-100:])
        if is_rev1 == -1 or is_rev2 == -1: continue
        if len(pb_s) <= len(ont_s) and pb_s:
            if is_rev1 == is_rev2:
                a.setRead(pb_s)
            else:
                a.setRead(rev_comp(pb_s))
            a.setReference(ont_s)
            alignment = a.align()
            pb_start, pb_end = alignment.read_start, alignment.read_end
            ont_start, ont_end = alignment.reference_start, alignment.reference_end
        elif pb_s and ont_s:
            if is_rev1 == is_rev2:
                a.setRead(ont_s)
            else:
                a.setRead(rev_comp(ont_s))
            a.setReference(pb_s)
            alignment = a.align()
            pb_start,pb_end = alignment.reference_start, alignment.reference_end
            ont_start,ont_end = alignment.read_start, alignment.read_end
        prefixq,alnq,suffixq=pb_q[:pb_start], \
                             pb_q[pb_start:pb_end],\
                             pb_q[pb_end:]

        if is_rev1 == 1:
            prefixq,alnq,suffixq=suffixq[::-1],alnq[::-1],prefixq[::-1]
        step = len(prefixq) * 1.0 / bins
        for i in range(bins):
            all_qualities[i].append(np.mean(prefixq[int(i * step):int((i + 1) * step)]))
        step = max(1, len(alnq) * 1.0 / bins)
        for i in range(bins):
            all_qualities[bins + i].append(np.mean(alnq[int(i * step):int((i + 1) * step)]))
        step = max(1, len(suffixq) * 1.0 / bins)
        for i in range(bins):
            all_qualities[2 * bins + i].append(np.mean(suffixq[int(i * step):int((i + 1) * step)]))

        prefixq,alnq,suffixq=ont_q[:ont_start], \
                             ont_q[ont_start:ont_end],\
                             ont_q[ont_end:]

        if is_rev2 == 1:
            prefixq,alnq,suffixq=suffixq[::-1],alnq[::-1],prefixq[::-1]
        step = len(prefixq)*1.0/bins
        for i in range(bins):
            all_qualities_ont[i].append(np.mean(prefixq[int(i*step):int((i+1)*step)]))
        step = max(1,len(alnq)*1.0/bins)
        for i in range(bins):
            all_qualities_ont[bins+i].append(np.mean(alnq[int(i*step):int((i+1)*step)]))
        step = max(1,len(suffixq)*1.0/bins)
        for i in range(bins):
            all_qualities_ont[2*bins+i].append(np.mean(suffixq[int(i*step):int((i+1)*step)]))

data = dict()
data['x'] = list(range(bins*3)) + list(range(bins*3))
data['Quality'] = [np.nanmean(p) for p in all_qualities] + [np.nanmean(p) for p in all_qualities_ont]
data['Quality_std'] = [np.nanstd(p) for p in all_qualities] + [np.nanstd(p) for p in all_qualities_ont]
data['Platform'] = ['PacBio'] * len(all_qualities) + ['ONT'] * len(all_qualities_ont)

with open("metaread.txt", "w") as f:
    f.write(" ".join([str(s) for s in data['x']]))
    f.write("\n")
    f.write(" ".join([str(s) for s in data['Quality']]))
    f.write("\n")
    f.write(" ".join([str(s) for s in data['Quality_std']]))
    f.write("\n")
    f.write(" ".join([str(s) for s in data['Platform']]))
    f.write("\n")
