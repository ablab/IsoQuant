import sys
from Bio import SeqIO
window_size = 16
def rev_comp(seq):
    c = dict(zip('ATCGNatcgn', 'TAGCNtagcn'))
    return ''.join(c.get(nucleotide, '') for nucleotide in reversed(seq))

barcodes = dict()

pb_file = sys.argv[1]
ont_file = sys.argv[2]
tsv_file = sys.argv[3]

with open(tsv_file) as f:
    for line in f:
        bc, umi, pb_read, ont_read = line.split()
        if ',' in line: continue
        barcodes[pb_read] = (bc,umi)
        barcodes[ont_read] = (bc,umi)
        umi_size=len(bc+umi)

pb_seqs = SeqIO.parse(pb_file, format="fastq")
trim_seqs = []
pb_lens = []
tail_lens = dict()
for seq in pb_seqs:
    s, qual = str(seq.seq), seq.letter_annotations["phred_quality"]
    if seq.id not in barcodes: continue
    bc, umi = barcodes[seq.id]
    umi_pos = s.find(bc+umi)
    is_reverse =False
    if umi_pos == -1:
        s=rev_comp(s)
        qual = qual[::-1]
        umi_pos = s.find(bc+umi)
        is_reverse = True
    if umi_pos == -1:
        continue
    polya_pos = umi_pos + umi_size
    mismatches = 0
    i = 0
    polya_startpos = polya_pos
    a_count = s[polya_pos:polya_pos+window_size].count('T')
    while True:
        if a_count < window_size*0.8:
            polya_lastpos = s[:polya_pos+window_size].rindex('TT') + 2
            break
        first_base_a = s[polya_pos] == 'T'
        new_base_a = polya_pos + window_size < len(s) and s[polya_pos + window_size] == 'T'
        if first_base_a and not new_base_a:
            a_count -= 1
        elif not first_base_a and new_base_a:
            a_count += 1
        polya_pos += 1

    s = s[polya_lastpos:]
    qual = qual[polya_lastpos:]
    trim_seq = seq[polya_lastpos:]
    if is_reverse:
        trim_seq = seq[:(len(seq.seq)-polya_lastpos)]
    trim_seqs.append(trim_seq)
    tail_lens[seq.id] = polya_lastpos-polya_startpos

ont_seqs = SeqIO.parse(ont_file, format="fastq")
trim_seqs = []
ont_lens = []
for seq in ont_seqs:
    s, qual = str(seq.seq), seq.letter_annotations["phred_quality"]
    if seq.id not in barcodes: continue
    bc, umi = barcodes[seq.id]
    umi_pos = s.find(bc+umi)
    is_reverse = False
    if umi_pos == -1:
        s=rev_comp(s)
        qual = qual[::-1]
        umi_pos = s.find(bc+umi)
        is_reverse = True
    if umi_pos == -1:
        continue
    polya_pos = umi_pos + umi_size
    mismatches = 0
    polya_startpos = polya_pos
    i = 0
    a_count = s[polya_pos:polya_pos+window_size].count('T')
    while True:
        if a_count < window_size*0.8:
            polya_lastpos = s[:polya_pos+window_size].rindex('TT') + 2
            break
        first_base_a = s[polya_pos] == 'T'
        new_base_a = polya_pos + window_size < len(s) and s[polya_pos + window_size] == 'T'
        if first_base_a and not new_base_a:
            a_count -= 1
        elif not first_base_a and new_base_a:
            a_count += 1
        polya_pos += 1

    s = s[polya_lastpos:]
    qual = qual[polya_lastpos:]
    tail = seq[polya_lastpos:]
    trim_seq = seq[polya_lastpos:]
    if is_reverse:
        trim_seq = seq[:(len(seq.seq)-polya_lastpos)]
    trim_seqs.append(trim_seq)
    tail_lens[seq.id] = polya_lastpos-polya_startpos

with open(tsv_file) as f:
    for line in f:
        bc, umi, pb_read, ont_read = line.split()
        if ',' in line: continue
        if pb_read not in tail_lens or ont_read not in tail_lens:
            continue
        pb_lens.append(tail_lens[pb_read])
        ont_lens.append(tail_lens[ont_read])
print(" ".join(str(s) for s in pb_lens))
print(" ".join(str(s) for s in ont_lens))
