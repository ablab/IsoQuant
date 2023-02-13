#!/usr/bin/env python
#
# ############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import logging
import argparse
from traceback import print_exc

import gffutils
from Bio import SeqIO
from collections import defaultdict

from .common import junctions_from_blocks
from .gtf2db import convert_gtf_to_db


logger = logging.getLogger('IsoQuant')


def array_to_coutns(arr):
    counts = defaultdict(int)
    for e in arr:
        counts[e] += 1
    return counts


def dump_dict_to_tsv(table, outf):
    for k in sorted(table.keys()):
        outf.write(str(k) + '\t' + str(table[k]) + '\n')


class AnnotationStats:
    canonical_splic_sites = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
    def __init__(self, gene_db, reference_record_dict):
        self.gene_db = gene_db
        self.reference_record_dict = reference_record_dict
        self.exon_lengths = []
        self.intron_length = []
        self.exon_lengths_no_dup = []
        self.transcript_lengths = []
        self.exons_per_transcript = []
        self.exons_per_gene = []
        self.transcripts_per_gene = []
        self.noncanonical_splice_sites_per_isoform = defaultdict(int)
        self.splice_site_dict = defaultdict(int)
        self.transcripts_with_non_canonical = defaultdict(list)

    def add_exon(self, exon):
        elen = exon[1] - exon[0] + 1
        self.exon_lengths.append(elen)
        return elen

    # return whether splice sites are canonical
    def add_intron(self, intron, strand, reference_region, gene_start):
        self.intron_length.append(intron[1] - intron[0] + 1)
        # splice site counts
        splice_site = None
        if reference_region:
            donor_site = reference_region[intron[0] - gene_start:intron[0]-gene_start+2]
            acceptor_site = reference_region[intron[1] - gene_start - 1:intron[1] - gene_start + 1]
            if strand == '+':
                splice_site = (str(donor_site), str(acceptor_site))
                self.splice_site_dict[splice_site] += 1
            else:
                donor_site = str(donor_site.reverse_complement())
                acceptor_site = str(acceptor_site.reverse_complement())
                splice_site = (acceptor_site, donor_site)
                self.splice_site_dict[splice_site] += 1
        return splice_site in self.canonical_splic_sites

    def add_transcript(self, exons, strand, reference_region, gene_start, t_id):
        tlen = 0
        for e in exons:
            tlen += self.add_exon(e)
        introns = junctions_from_blocks(exons)
        non_canonical_splice_sites = 0
        for i in introns:
            is_canonical = self.add_intron(i, strand, reference_region, gene_start)
            if not is_canonical:
                non_canonical_splice_sites += 1
                self.transcripts_with_non_canonical[t_id].append(i)

        self.transcript_lengths.append(tlen)
        self.exons_per_transcript.append(len(exons))
        self.noncanonical_splice_sites_per_isoform[non_canonical_splice_sites] += 1

    def add_gene(self, exon_set, transcript_count):
        self.exons_per_gene.append(len(exon_set))
        for e in exon_set:
            self.exon_lengths_no_dup.append(e[1] - e[0] + 1)
        self.transcripts_per_gene.append(transcript_count)


    def count_gene_stats(self, gene_data, current_chr_record):
        exon_set = set()
        transcript_count = 0

        reference_region = None
        gene_start = gene_data.start
        if self.reference_record_dict:
            reference_region = current_chr_record[gene_start - 1:gene_data.end].seq

        for t in self.gene_db.children(gene_data, featuretype=('transcript', 'mRNA'), order_by='start'):
            exon_list = []
            transcript_count += 1
            for e in self.gene_db.children(t, order_by='start'):
                if e.featuretype == 'exon':
                    exon = (e.start, e.end)
                    exon_set.add(exon)
                    exon_list.append(exon)
            self.add_transcript(exon_list, t.strand, reference_region, gene_start, t.id)
        self.add_gene(exon_set, transcript_count)

    def count_stats(self):
        current_chr = ""
        current_chr_record = None
        for g in self.gene_db.features_of_type('gene', order_by=('seqid', 'start')):
            if current_chr != g.seqid:
                current_chr = g.seqid
                logger.info("Processing " + current_chr)
                current_chr_record = None if not self.reference_record_dict else self.reference_record_dict[current_chr]

            gene_name = g.id
            self.count_gene_stats(self.gene_db[gene_name], current_chr_record)

    def print_to_file(self, output):
        outf = open(output, "w")
        outf.write("Exon length distribution (same exons in different isoforms counted multiple times)\n")
        dump_dict_to_tsv(array_to_coutns(self.exon_lengths), outf)
        outf.write("Exon length distribution (each unique exon is counted once)\n")
        dump_dict_to_tsv(array_to_coutns(self.exon_lengths_no_dup), outf)
        outf.write("Intron length distribution\n")
        dump_dict_to_tsv(array_to_coutns(self.intron_length), outf)
        outf.write("Exons per transcripts\n")
        dump_dict_to_tsv(array_to_coutns(self.exons_per_transcript), outf)
        outf.write("Exons per gene\n")
        dump_dict_to_tsv(array_to_coutns(self.exons_per_gene), outf)
        outf.write("Transcripts per gene\n")
        dump_dict_to_tsv(array_to_coutns(self.transcripts_per_gene), outf)
        if len(self.splice_site_dict) > 0:
            outf.write("Splice site per gene\n")
            dump_dict_to_tsv(self.splice_site_dict, outf)
        if len(self.noncanonical_splice_sites_per_isoform) > 0:
            outf.write("Non-canonical splice sites per isoform\n")
            dump_dict_to_tsv(self.noncanonical_splice_sites_per_isoform, outf)
            outf.write("Transcripts with non-canonical splice sites:\n")
            for t in sorted(self.transcripts_with_non_canonical.keys()):
                outf.write(t + ":\t")
                outf.write("\t".join(list(map(str, self.transcripts_with_non_canonical[t]))) + "\n")
        outf.close()

def parse_args(args=None, namespace=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--output", "-o", help="output file prefix [default=gtf_stats]",
                        type=str, default="gtf_stats.tsv")

    # REFERENCE
    parser.add_argument("--genedb", "-g", help="gene database in gffutils DB format or GTF/GFF format", type=str)
    parser.add_argument('--complete_genedb', action='store_true', default=False,
                        help="use this flag if gene annotation contains transcript and gene metafeatures, "
                             "e.g. with official annotations, such as GENCODE; "
                             "speeds up gene database conversion")
    parser.add_argument("--reference", "-r", help="reference genome in FASTA format,"
                                                  "should be provided to compute splice site stats", type=str)
    parser.add_argument('--clean_start', action='store_true', default=False,
                        help='Do not use previously generated genee db')

    args = parser.parse_args(args, namespace)

    config_dir = os.path.join(os.environ['HOME'], '.config', 'IsoQuant')
    os.makedirs(config_dir, exist_ok=True)
    args.db_config_path = os.path.join(config_dir, 'db_config.json')

    if not args.genedb:
        parser.print_usage()
        exit(-1)
    return args


def run_pipeline(args):
    logger.info(" === Counting gene annotation statistics === ")

    if not args.genedb.endswith('db'):
        args.genedb = convert_gtf_to_db(args, output_is_dir=False)

    logger.info("Loading gene database from " + args.genedb)
    gffutils_db = gffutils.FeatureDB(args.genedb, keep_order=True)
    reference_record_dict = None
    if args.reference:
        logger.info("Loading reference genome from " + args.reference)
        reference_record_dict = SeqIO.index(args.reference, "fasta")

    stats = AnnotationStats(gffutils_db, reference_record_dict)
    logger.info("Counting stats, may take a while...")
    stats.count_stats()
    logger.info("Writing stats to " + args.output + ".tsv")
    stats.print_to_file(args.output + ".tsv")

    logger.info(" === Counting done === ")


def main(args):
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    args = parse_args(args)
    run_pipeline(args)


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)



