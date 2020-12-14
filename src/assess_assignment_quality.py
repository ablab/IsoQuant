############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import re
import sys
import argparse
from collections import defaultdict
from traceback import print_exc

import pysam
from Bio import SeqIO
import gffutils
from enum import Enum
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

from src.common import overlaps

id_pattern = re.compile("[A-Z]+\.?(\d+\.\d+)")


class ReadType(Enum):
    CORRECT = 1
    INCORRECT_SAME_GENE = 2
    INCORRECT_OTHER_GENE = 3
    NOT_ASSIGNED = 4


class MappingData:
    def __init__(self, args):
        self.seq_set = set()
        self.seqid_to_isoform = dict()
        self.parse_fasta(args.fasta, args.real_data)

        self.mapped_seqs = {}
        self.secondary_mapped_seqs = {}
        if not args.real_data:
            self.parse_bam(args.mapping)

    def parse_fasta(self, fasta, is_real_data):
        for record in SeqIO.parse(fasta, "fasta"):
            if is_real_data:
                #TODO: check SQANTI seq names
                seq_id = record.id
            else:
                seq_id = re.search(id_pattern, record.description).group(1)
                isoform_id = record.description.replace('_', ' ').split()[-1]
                self.seqid_to_isoform[seq_id] = isoform_id
            self.seq_set.add(seq_id)

    def parse_bam(self, bam_file):
        alignment_file_in = pysam.AlignmentFile(bam_file, "rb")
        for alignment in alignment_file_in:
            if alignment.reference_id == -1:
                continue
            blocks = sorted(alignment.get_blocks())
            if len(blocks) == 0:
                continue
            seq_id = id_pattern.search(alignment.query_name).group(1)
            if alignment.is_secondary or alignment.is_supplementary:
                self.secondary_mapped_seqs[seq_id] = (blocks[0][0], blocks[-1][1])
            else:
                self.mapped_seqs[seq_id] = (blocks[0][0], blocks[-1][1])


class AssignmentData:
    def __init__(self, tsv_file, is_real_data):
        self.assigned_isoforms = defaultdict(str)
        self.parse_tsv(tsv_file, is_real_data)

    def parse_tsv(self, tsv_file, is_real_data):
        with open(tsv_file) as f:
            for i,l in enumerate(f):
                if i == 0:
                    continue
                tokens = l.strip().split()
                seq_id = tokens[0] if is_real_data else id_pattern.search(tokens[0]).group(1)
                if len(tokens) > 10:  ## SQANTI2
                    self.assigned_isoforms[seq_id] = tokens[7]
                else:
                    self.assigned_isoforms[seq_id] = tokens[1]


class StatCounter:
    def __init__(self, prefix, mapping_data, assignment_data):
        self.prefix = prefix
        self.mapping_data = mapping_data
        self.assignment_data = assignment_data
        self.correct_seqs = set()
        self.mismapped_seqs = set()
        self.seq_assignments = defaultdict()

    def count_mapping_stats(self, db):
        for seq_id in self.mapping_data.seq_set:
            if seq_id in self.mapping_data.mapped_seqs:
                gene_name = db.isoform_to_gene_map[self.mapping_data.seqid_to_isoform[seq_id]]
                if overlaps(self.mapping_data.mapped_seqs[seq_id], db.get_gene_coords(gene_name)):
                    self.correct_seqs.add(seq_id)
                else:
                    self.mismapped_seqs.add(seq_id)

    def count_assignment_stats(self, db, use_mismapped=False):
        self.seq_assignments = defaultdict()
        c_fname = self.prefix + ".corr_isoforms.txt" if use_mismapped else self.prefix + ".full_corr_isoforms.txt"
        w_fname = self.prefix + ".wrong_isoforms.txt" if use_mismapped else self.prefix + ".full_wrong_isoforms.txt"
        u_fname = self.prefix + ".unassigned_isoforms.txt" if use_mismapped else self.prefix + ".full_unassigned_isoforms.txt"
        c_f = open(c_fname, "w")
        w_f = open(w_fname, "w")
        u_f = open(u_fname, "w")
        for seq_id in self.mapping_data.seq_set:
            if not use_mismapped and seq_id in self.mismapped_seqs:
                continue
            if seq_id in self.assignment_data.assigned_isoforms.keys():
                real_isoform_id = self.mapping_data.seqid_to_isoform[seq_id]
                assigned_isoform_id = self.assignment_data.assigned_isoforms[seq_id]
                if assigned_isoform_id == real_isoform_id:
                    self.seq_assignments[seq_id] = ReadType.CORRECT
                    c_f.write(seq_id + "\n")
                elif assigned_isoform_id == 'novel':
                    self.seq_assignments[seq_id] = ReadType.NOT_ASSIGNED
                    u_f.write(seq_id + "\n")
                else:
                    gene_name = db.isoform_to_gene_map[real_isoform_id]
                    if assigned_isoform_id in db.gene_to_isoforms_map[gene_name]:
                        self.seq_assignments[seq_id] = ReadType.INCORRECT_SAME_GENE
                        w_f.write(seq_id + "\n")
                    else:
                        self.seq_assignments[seq_id] = ReadType.INCORRECT_OTHER_GENE
                        w_f.write(seq_id + "\n")
            else:
                self.seq_assignments[seq_id] = ReadType.NOT_ASSIGNED
                u_f.write(seq_id + "\n")

    def get_read_counts(self, read_type):
        return sum(val == read_type for val in self.seq_assignments.values())

    def calc_precision(self, tp, fp):
        return tp * 100.0 / (tp + fp)

    def calc_recall(self, tp, fn):
        return tp * 100.0 / (tp + fn)

    def print_stats(self, tp, fp, fn, stream, name=""):
        stream.write("correct\t%d\nincorrect\t%d\nunmapped/unassigned\t%d\n" % (tp, fp, fn))
        precision, recall = self.calc_precision(tp, fp), self.calc_recall(tp, fn)
        stream.write("%sprecision\t%d\n%srecall\t%d\n" % (name, precision, name, recall))


class DbHandler:
    def __init__(self, gene_db):
        self.db = gffutils.FeatureDB(gene_db, keep_order=True)
        self.isoform_to_gene_map = {}
        self.gene_to_isoforms_map = {}

        self.parse_db()

    def parse_db(self):
        for g in self.db.features_of_type('gene'):
            gene_name = g.id
            self.gene_to_isoforms_map[gene_name] = set()
            gene_db = self.db[gene_name]
            for t in self.db.children(gene_db, featuretype='transcript'):
                self.isoform_to_gene_map[t.id] = gene_name
                self.gene_to_isoforms_map[gene_name].add(t.id)

    def get_gene_by_isoform(self, t_id):
        return self.isoform_to_gene_map[t_id]

    def get_isoforms_by_gene(self, gene_name):
        return self.gene_to_isoforms_map[gene_name]

    def get_gene_coords(self, gene_name):
        return self.db[gene_name].start, self.db[gene_name].end


def compare_stats(data_a, data_b):
    stats_a, label_a = data_a
    stats_b, label_b = data_b
    a, b, ab = 0, 0, 0
    for seq_id in stats_a.mapping_data.seq_set:
        correct_a = stats_a.seq_assignments[seq_id] == ReadType.CORRECT
        correct_b = stats_b.seq_assignments[seq_id] == ReadType.CORRECT
        if correct_a and correct_b:
            ab += 1
        elif correct_a:
            a += 1
        elif correct_b:
            b += 1

    venn2(subsets=(a, b, ab), set_labels=(label_a, label_b))
    plt.savefig("venn.png")
    print("%d correct reads by both methods, %d in %s only and %d in %s only" % (ab, a, label_a, b, label_b))
    print("Venn diagram saved to venn.png")
    return


def compare_real_results(data_a, data_b):
    stats_a, label_a = data_a
    stats_b, label_b = data_b
    a, b, ab, diff_ab, not_ab = 0, 0, 0, 0, 0
    for seq_id in stats_a.mapping_data.seq_set:
        isoform_a, isoform_b = stats_a.assignment_data.assigned_isoforms[seq_id], stats_b.assignment_data.assigned_isoforms[seq_id]
        if isoform_a and isoform_b:
            if stats_a.mapping_data.assigned_isoforms[seq_id] == isoform_b:
                ab += 1
            else:
                diff_ab += 1
        elif isoform_a:
            a += 1
        elif isoform_b:
            b += 1
        else:
            not_ab += 1
    print('Common assignments: %d, different assignments: %d, '
          'only %s assignments: %d, only %s assignments: %d, both not assigned: %d' %
          (ab, diff_ab, label_a, a, label_b, b, not_ab))

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="Output file name")
    parser.add_argument("--fasta", "-f", type=str, help="initial file with sequences")
    parser.add_argument("--mapping", "-m", type=str, help="mapped sequences (SAM or BAM format)") ## SQANTI2 output
    parser.add_argument("--tsv", "-t", type=str, nargs='+', help="assigned isoforms, max number of files to compare: 2")
    parser.add_argument("--gene_db", "-g", type=str, help="gene database")
    parser.add_argument("--real_data", default=False, action="store_true",
                        help="real data is used (correct isoforms for each read are not provided)")

    args = parser.parse_args()
    if len(args.tsv) > 2:
        print("ERROR! Maximum number of files to compare: 2")
        exit(-1)
    return args


def main():
    args = parse_args()
    db = DbHandler(args.gene_db)
    mapping_data = MappingData(args)
    all_stats = []
    output_file = sys.stdout if not args.output else open(args.output, "w")

    for tsv_file in args.tsv:
        #TODO: diff bams
        print()
        print("Calculating stats for %s" % tsv_file)
        assignment_data = AssignmentData(tsv_file, args.real_data)
        prefix = os.path.splitext(os.path.basename(tsv_file))[0]
        stat_counter = StatCounter(prefix, mapping_data, assignment_data)
        if args.real_data:
            all_stats.append((stat_counter, prefix))
            continue
        stat_counter.count_mapping_stats(db)

        total_reads = len(stat_counter.mapping_data.seq_set)
        correctly_mapped = len(stat_counter.correct_seqs)
        mismapped_reads = len(stat_counter.mismapped_seqs)
        unmapped_reads = total_reads - (correctly_mapped + mismapped_reads)
        output_file.write("# MAPPING STATS \n")
        stat_counter.print_stats(correctly_mapped, mismapped_reads, unmapped_reads, output_file, name="mapping_")

        # use only correctly mapped reads
        stat_counter.count_assignment_stats(db)
        correct_assignments = stat_counter.get_read_counts(ReadType.CORRECT)
        incorrect_assignments = stat_counter.get_read_counts(ReadType.INCORRECT_OTHER_GENE) + stat_counter.get_read_counts(ReadType.INCORRECT_SAME_GENE)
        unassigned_reads = stat_counter.get_read_counts(ReadType.NOT_ASSIGNED)
        print("# ASSIGNMENT STATS")
        stat_counter.print_stats(correct_assignments, incorrect_assignments, unassigned_reads, output_file, name="assignment_")

        # use all reads
        stat_counter.count_assignment_stats(db, use_mismapped=True)
        correct_assignments = stat_counter.get_read_counts(ReadType.CORRECT)
        incorrect_assignments = stat_counter.get_read_counts(ReadType.INCORRECT_OTHER_GENE) + stat_counter.get_read_counts(ReadType.INCORRECT_SAME_GENE)
        unassigned_reads = stat_counter.get_read_counts(ReadType.NOT_ASSIGNED)
        print("# ASSIGNMENT STATS (INCLUDING MISMAPPED READS) ")
        stat_counter.print_stats(correct_assignments, incorrect_assignments, unassigned_reads, output_file, name="overall_")
        all_stats.append((stat_counter, prefix))

    if len(all_stats) == 2:
        if args.real_data:
            compare_real_results(all_stats[0], all_stats[1])
        else:
            compare_stats(all_stats[0], all_stats[1])

    if args.output:
        output_file.close()


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)

