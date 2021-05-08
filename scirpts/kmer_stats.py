import argparse

from Bio import SeqIO
from itertools import groupby
import gffutils
import sys
from Bio.Seq import Seq

is_hp = 1
k = 14

def parse_args(args=None, namespace=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--reference", "-r", help="reference genome in FASTA format", type=str)
    parser.add_argument("--genedb", "-g", help="gene database in gffutils DB format or GTF/GFF format", type=str)

    parser.add_argument('--assign_pb', type=str, help='IsoQuant read assignments file for PacBio')
    parser.add_argument('--assign_ont', type=str, help='IsoQuant read assignments file for ONT')

    parser.add_argument('--fastq_pb', type=str, help='FASTQ file for PacBio')
    parser.add_argument('--fastq_ont', type=str, help='FASTQ file for ONT')
    args = parser.parse_args(args)
    return args

def rev_comp(seq):
    c = dict(zip('ATCGNatcgn', 'TAGCNtagcn'))
    return ''.join(c.get(nucleotide, '') for nucleotide in reversed(seq))
def compress_hp(seq):
    if not is_hp: return seq
    return ''.join(x[0] for x in groupby(list(seq)))

args = parse_args(sys.argv[1:])
db = gffutils.FeatureDB(args.genedb, keep_order = True)
transcripts = SeqIO.parse(args.reference, format="fasta")
chr_seqs = {seq.id: str(seq.seq) for seq in transcripts}

pb_seqs = SeqIO.parse(args.fastq_pb, format="fastq")
ont_seqs = SeqIO.parse(args.fastq_ont, format="fastq")

pb_seq_dict = {seq.id: str(seq.seq) for seq in pb_seqs}
ont_seq_dict = {seq.id: str(seq.seq) for seq in ont_seqs}

def canon_kmer(kmer):
    kmer = ''.join(x[0] for x in groupby(list(kmer)))
    return min(kmer, rev_comp(kmer))

exon_lens = [21, 50, 100, 200, 20000]
exon_titles = ['<' + str(exon_lens[0])] + ['%d-%d' % (exon_lens[i]+1, exon_lens[i+1]) for i in range(len(exon_lens)-2)] + ['>' + str(exon_lens[-2]+1)]

ont_filt_data = []
ont_data = []

with open(args.assign_ont) as f:
    f.readline()
    for line in f:
        read_name, isoform_id, assignment_type, assignment_events = line.split()[:4]
        aligned_blocks = line.split()[5]
        if isoform_id == 'None' or isoform_id == 'noninformative': continue
        if 'unique' not in assignment_type: continue
        read_st, read_en = int(aligned_blocks.split('-')[0]), int(aligned_blocks.split('-')[-1])
        read_seq = ont_seq_dict[read_name]
        read_seq = compress_hp(read_seq)
        read_kmers = set([canon_kmer(read_seq[i:i+k]) for i in range(len(read_seq)-k+1)])
        t=db[isoform_id]
        for e in db.children(t, featuretype='exon', order_by='start'):
            st, en = e.start, e.end
            if en-st+1 < k: continue
            if read_en < en or read_st > st: continue
            exon_seq = chr_seqs[e.chrom][st:en+1]
            exon_seq = compress_hp(exon_seq)
            exon_len = len(exon_seq)
            exon_kmers = set([canon_kmer(exon_seq[i:i+k]) for i in range(len(exon_seq)-k+1)])
            if not exon_kmers: continue
            intsect = len(exon_kmers.intersection(read_kmers))*1.0/len(exon_kmers)
            ont_data.append((exon_len, intsect))

pb_data = []
with open(args.assign_pb) as f:
    f.readline()
    for line in f:
        read_name, isoform_id, assignment_type, assignment_events = line.split()[:4]
        aligned_blocks = line.split()[5]
        if isoform_id == 'None' or isoform_id == 'noninformative': continue
        if 'unique' not in assignment_type: continue
        read_st, read_en = int(aligned_blocks.split('-')[0]), int(aligned_blocks.split('-')[-1])
        read_seq = pb_seq_dict[read_name]
        read_seq = compress_hp(read_seq)
        read_kmers = set([canon_kmer(read_seq[i:i+k]) for i in range(len(read_seq)-k+1)])
        t=db[isoform_id]
        for e in db.children(t, featuretype='exon', order_by='start'):
            st, en = e.start, e.end
            if en-st+1 < k: continue
            if read_en < en or read_st > st: continue
            exon_seq = chr_seqs[e.chrom][st:en+1]
            exon_len = len(exon_seq)
            exon_seq = compress_hp(exon_seq)
            exon_kmers = set([canon_kmer(exon_seq[i:i+k]) for i in range(len(exon_seq)-k+1)])
            if not exon_kmers: continue
            intsect = len(exon_kmers.intersection(read_kmers))*1.0/len(exon_kmers)
            pb_data.append((exon_len, intsect))

pb_filt_data = []
pb_exons = []
for exon_len, kmer_pct in pb_data:
    for i, l in enumerate(exon_lens):
        if exon_len <= l:
            pb_filt_data.append(kmer_pct*100)
            pb_exons.append(exon_titles[i])
            break

ont_exons = []
for exon_len, kmer_pct in ont_data:
    for i, l in enumerate(exon_lens):
        if exon_len <= l: 
            ont_filt_data.append(kmer_pct*100)
            ont_exons.append(exon_titles[i])
            break

with open("kmer_stats_compress%d_k%d.txt" % (is_hp, k), "w") as f:
    f.write(" ".join([str(s) for s in (pb_filt_data+ont_filt_data)]))
    f.write("\n")
    f.write(" ".join(['PacBio'] * len(pb_filt_data) + ['ONT'] * len(ont_filt_data)))
    f.write("\n")
    f.write(" ".join([str(s) for s in (pb_exons+ont_exons)]))

