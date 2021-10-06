############################################################################
# Copyright (c) 2021 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import argparse
from Bio import SeqIO
import numpy
from traceback import print_exc
import gffutils


def read_probs(args):
    values = []
    probs = []
    for l in open(args.prob):
        t = l.strip().split()
        values.append(int(t[0]))
        probs.append(float(t[1]))
    return values, probs


def simulate_10x(args, isoforms):
    values, probs = read_probs(args)
    outf = os.path.join(args.output, 'tmp.fa')
    outfq_prefix = os.path.join(args.output, 'tmp')
    isoform_counts = {}

    for record in SeqIO.parse(args.fasta, 'fasta'):
        if len(record.seq) < 200:
            continue
        isoform_id = record.description.strip().replace('>', '')
        if len(isoforms) > 0 and isoform_id not in isoforms:
            continue

        if isoform_id not in isoform_counts:
            isoform_counts[isoform_id] = 0
        isoform_counts[isoform_id] += 1

        if os.path.exists(outf):
            os.remove(outf)
        SeqIO.write([record], outf, 'fasta')
        nreads = numpy.random.choice(values, p = probs)

        if os.path.exists(outfq_prefix + '_R1.fastq'):
            os.remove(outfq_prefix + '_R1.fastq')
        if os.path.exists(outfq_prefix + '_R2.fastq'):
            os.remove(outfq_prefix + '_R2.fastq')

        # install insilicoseq first e.g. `pip install InSilicoSeq`
        os.system('$(which iss) generate -p 6  --genomes ' + outf + ' --n_reads ' + str(nreads) + ' --model HiSeq  -a uniform -o ' + outfq_prefix + ' 2> /dev/null')

        isoform_name = isoform_id + '_' + str(isoform_counts[isoform_id])
        left_reads = []
        for read in SeqIO.parse(outfq_prefix + '_R1.fastq', 'fastq'):
            read.id = read.id.split('/')[0] + ':1___' + isoform_id + '_' + str(isoform_counts[isoform_id])
            read.description = ''
            left_reads.append(read)
        left_fq = os.path.join(args.output, isoform_name) + '_R1.fastq'
        SeqIO.write(left_reads, left_fq, 'fastq')

        right_reads = []
        for read in SeqIO.parse(outfq_prefix + '_R2.fastq', 'fastq'):
            read.id = read.id.split('/')[0] + ':2___' + isoform_id + '_' + str(isoform_counts[isoform_id])
            read.description = ''
            right_reads.append(read)
        right_fq = os.path.join(args.output, isoform_name) + '_R2.fastq'
        SeqIO.write(right_reads, right_fq, 'fastq')


def get_isoform_set(args):
    isoforms = set()
    if args.genedb is None or args.gene_id is None:
        return isoforms

    db = gffutils.FeatureDB(args.genedb)
    gene = db[args.gene_id]

    for t in db.children(gene, featuretype='transcript'):
        isoforms.add(t.id)
    return isoforms


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", "-f", help="initial file with sequences", type=str)
    parser.add_argument("--prob", help="probability values", type=str)
    parser.add_argument("--output", "-o", help="output folder", type=str)
    parser.add_argument("--genedb", "-g", help="gene database in gffutils db format", type=str)
    parser.add_argument("--gene_id", help="gene id to take isoforms from", type=str)


    args = parser.parse_args()

    if args.fasta is None or args.prob is None:
        parser.print_help()
        exit(-1)
    if args.output is None:
        print("No output folder specified, writing to current directory")
        args.output = './'

    if not os.path.exists(args.output):
        os.makedirs(args.output)
    return args


def main():
    args = parse_args()
    isoforms = get_isoform_set(args)
    simulate_10x(args, isoforms)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
