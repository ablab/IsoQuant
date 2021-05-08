#!/usr/bin/end python
import sys
import pysam
from Bio import SeqIO

samfile_pb = sys.argv[1]
samfile_ont = sys.argv[2]
tsv_file = sys.argv[3]

infile = pysam.AlignmentFile(samfile_pb, "r", check_sq=False)
len_dict = dict()
for aln in infile:
    if aln.is_secondary:continue
    if not aln.reference_name: continue
    read_name = aln.query_name
    len_dict[read_name] = aln.query_alignment_length

infile= pysam.AlignmentFile(samfile_ont, "r", check_sq=False)
for aln in infile:
    if aln.is_secondary:continue
    if not aln.reference_name: continue
    read_name = aln.query_name
    len_dict[read_name] = aln.query_alignment_length

pb_lens, ont_lens = [], []
with open(tsv_file) as f:
    for line in f:
        _, _, pb_read, ont_read = line.split()
        if pb_read not in len_dict or ont_read not in len_dict: continue
        pb_lens.append(len_dict[pb_read])
        ont_lens.append(len_dict[ont_read])

with open("aligned_len.tsv","w") as outf:
    outf.write(" ".join([str(s) for s in pb_lens]))
    outf.write("\n")
    outf.write(" ".join([str(s) for s in ont_lens]))
