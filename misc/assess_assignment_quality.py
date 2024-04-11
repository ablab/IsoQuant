#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################
import gzip
import os
import re
import sys
import argparse
from collections import defaultdict
from traceback import print_exc
from collections import namedtuple

import pysam
from Bio import SeqIO
import gffutils
from enum import Enum
import logging


logger = logging.getLogger('IsoQuantQA')

id_pattern = re.compile("[A-Z]+\.?(\d+\.\d+)")


def overlaps(range1, range2):
    return not (range1[1] < range2[0] or range1[0] > range2[1])


def set_logger(logger_instance):
    logger_instance.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


class ReadType(Enum):
    CORRECT = 1
    INCORRECT_SAME_GENE = 2
    INCORRECT_OTHER_GENE = 3
    NOT_ASSIGNED = 4


class MappingData:
    def __init__(self, args):
        self.seq_set = set()
        self.seqid_to_isoform = dict()
        self.args = args
        self.parse_fasta(args.fasta, args.real_data)
        self.mapped_seqs = {}
        self.secondary_mapped_seqs = defaultdict(list)
        if not args.real_data:
            self.parse_bam(args.mapping)

    def parse_fasta(self, fasta, is_real_data=False):
        logger.info("Loading reads from %s" % fasta)
        basename_plus_inner_ext, outer_ext = os.path.splitext(fasta.lower())
        if outer_ext not in ['.zip', '.gz', '.gzip', '.bz2', '.bzip2']:
            basename_plus_inner_ext, outer_ext = fasta, ''  # not a supported archive
        basename, fasta_ext = os.path.splitext(basename_plus_inner_ext)
        if fasta_ext in {'.fastq', '.fq'}:
            data_type = 'fastq'
        elif fasta_ext in {'.fasta', '.fa', '.fna'}:
            data_type = 'fasta'
        else:
            logger.error("Unsupported extension: %s" % fasta_ext)
            return

        for record in SeqIO.parse(fasta, data_type):
            if is_real_data:
                #TODO: check SQANTI seq names
                seq_id = record.id
            else:
                seq_id = record.id
                tokens = record.id.split('_')
                if len(tokens) <= self.args.isoform_id_index:
                    logger.warning("Malformed read id %s" % seq_id)
                    continue
                isoform_id = tokens[self.args.isoform_id_index]
                if isoform_id.startswith("Pac"):
                    self.args.isoform_id_index = -1
                    isoform_id = tokens[self.args.isoform_id_index]
                if isoform_id.startswith("Y"):
                    isoform_id = tokens[-3]
                if isoform_id.startswith("E"):
                    self.seqid_to_isoform[seq_id] = isoform_id
                else:
                    logger.warning("Unexpected isoform id %s" % isoform_id)
            self.seq_set.add(seq_id)
        logger.info("Total %d sequences loaded" % len(self.seq_set))

    def parse_bam(self, bam_file):
        logger.info("Loading alignments from %s" % bam_file)
        if bam_file.endswith("sam"):
            alignment_file_in = pysam.AlignmentFile(bam_file, "r")
        else:
            alignment_file_in = pysam.AlignmentFile(bam_file, "rb")
        for alignment in alignment_file_in:
            if alignment.reference_id == -1:
                continue
            read_coords = (alignment.reference_start, alignment.reference_end)
            if not read_coords[1]:
                continue
            seq_id = alignment.query_name
            if alignment.is_secondary or alignment.is_supplementary:
                self.secondary_mapped_seqs[seq_id].append(read_coords)
            else:
                self.mapped_seqs[seq_id] = read_coords
        logger.info("Total alignments loaded: %d" % len(self.mapped_seqs))


class AssignmentData:
    # read id column, isoform id column, assignment type column
    ParseParams = namedtuple("ParseParams", ("read_id_column", "isoform_id_column", "assignment_type_column"))
    PRESETS = {'isoquant': ParseParams(0, 3, 5), 'talon': ParseParams(0, 12, 16), 'sqanti': ParseParams(0, 7, 5)}
    UNIQUE_ASSIGNMENTS_TYPES = {"unique", "unique_minor_difference", "Known", "ISM", "full-splice_match", "incomplete-splice_match"}
    AMB_ASSIGNMENTS_TYPES = {"unique", "unique_minor_difference", "ambiguous"}
    ALL_ASSIGNMENTS_TYPES = {"unique", "unique_minor_difference", "ambiguous", "inconsistent"}

    def __init__(self, tsv_file, preset='isoquant', assignment_types = UNIQUE_ASSIGNMENTS_TYPES):
        self.parse_params = self.PRESETS[preset]
        self.assigned_isoforms = defaultdict(str)
        self.assignment_types = assignment_types
        self.parse_tsv(tsv_file)

    def parse_tsv(self, tsv_file):
        logger.info("Reading assignments from %s" % tsv_file)
        if tsv_file.endswith(".gz") or tsv_file.endswith(".gzip"):
            handler = gzip.open(tsv_file, "rt")
        else:
            handler = open(tsv_file, "rt")
        for l in handler:
            if l.startswith('#'):
                continue
            tokens = l.strip().split()
            seq_id = tokens[self.parse_params.read_id_column] # if is_real_data else id_pattern.search(tokens[0]).group(1)
            assignment_type = tokens[self.parse_params.assignment_type_column]
            if assignment_type in self.assignment_types:
                isoform_id = tokens[self.parse_params.isoform_id_column]
                if not isoform_id[0] == 'E':
                    logger.warning("Assigned to weird isoform " + isoform_id)
                self.assigned_isoforms[seq_id] = isoform_id
        handler.close()
        logger.info("Total assignments loaded: %d" % len(self.assigned_isoforms))


class StatCounter:
    def __init__(self, prefix, mapping_data, assignment_data):
        self.prefix = prefix
        self.mapping_data = mapping_data
        self.assignment_data = assignment_data
        self.correct_seqs = set()
        self.mismapped_seqs = set()
        self.secondary_mapped = set()
        self.seq_assignments = defaultdict()

    def count_mapping_stats(self, db, read_list=None):
        reads_to_iterate = read_list if read_list is not None else self.mapping_data.mapped_seqs
        for seq_id in reads_to_iterate:
            isoform_id = self.mapping_data.seqid_to_isoform[seq_id]
            if isoform_id not in db.isoform_to_gene_map:
                self.mismapped_seqs.add(seq_id)
                continue

            gene_name = db.isoform_to_gene_map[isoform_id]
            if overlaps(self.mapping_data.mapped_seqs[seq_id], db.get_gene_coords(gene_name)):
                self.correct_seqs.add(seq_id)
            elif any(overlaps(sec, db.get_gene_coords(gene_name)) for sec in self.mapping_data.secondary_mapped_seqs[seq_id]):
                self.secondary_mapped.add(seq_id)
            else:
                self.mismapped_seqs.add(seq_id)

    def count_assignment_stats(self, db, use_mismapped=False, output_dir=None):
        self.seq_assignments = defaultdict()
        c_f = None
        w_f = None
        u_f = None
        if output_dir:
            c_fname = self.prefix + ".corr_isoforms.txt" if use_mismapped else self.prefix + ".full_corr_isoforms.txt"
            w_fname = self.prefix + ".wrong_isoforms.txt" if use_mismapped else self.prefix + ".full_wrong_isoforms.txt"
            u_fname = self.prefix + ".unassigned_isoforms.txt" if use_mismapped else self.prefix + ".full_unassigned_isoforms.txt"
            c_f = open(os.path.join(output_dir, c_fname), "w")
            w_f = open(os.path.join(output_dir, w_fname), "w")
            u_f = open(os.path.join(output_dir, u_fname), "w")
        for seq_id in self.mapping_data.seq_set:
            if not use_mismapped and seq_id not in self.correct_seqs:
                continue
            if seq_id in self.assignment_data.assigned_isoforms.keys():
                real_isoform_id = self.mapping_data.seqid_to_isoform[seq_id]
                assigned_isoform_id = self.assignment_data.assigned_isoforms[seq_id]
                if assigned_isoform_id == real_isoform_id:
                    self.seq_assignments[seq_id] = ReadType.CORRECT
                    if c_f:
                        c_f.write(seq_id + "\n")
                elif assigned_isoform_id == 'novel':
                    self.seq_assignments[seq_id] = ReadType.NOT_ASSIGNED
                    if u_f:
                        u_f.write(seq_id + "\n")
                else:
                    if real_isoform_id not in db.isoform_to_gene_map or \
                            assigned_isoform_id not in db.gene_to_isoforms_map[db.isoform_to_gene_map[real_isoform_id]]:
                        self.seq_assignments[seq_id] = ReadType.INCORRECT_OTHER_GENE
                        if w_f:
                            w_f.write(seq_id + "\n")
                    else:
                        self.seq_assignments[seq_id] = ReadType.INCORRECT_SAME_GENE
                        if w_f:
                            w_f.write(seq_id + "\n")

            else:
                self.seq_assignments[seq_id] = ReadType.NOT_ASSIGNED
                if w_f:
                    u_f.write(seq_id + "\n")

    def get_read_counts(self, read_type):
        return sum(val == read_type for val in self.seq_assignments.values())

    def calc_precision(self, tp, fp):
        if tp + fp == 0.0:
            return 0.0
        return tp * 100.0 / (tp + fp)

    def calc_recall(self, tp, fn):
        if tp + fn == 0.0:
            return 0.0
        return tp * 100.0 / (tp + fn)

    def print_stats(self, tp, fp, fn, stream, name=""):
        stream.write("correct\t%d\nincorrect\t%d\nunmapped/unassigned\t%d\n" % (tp, fp, fn))
        precision, recall = self.calc_precision(tp, fp), self.calc_recall(tp, fn)
        stream.write("%sprecision\t%2.4f\n%srecall\t%2.4f\n" % (name, precision, name, recall))


class DbHandler:
    def __init__(self, gene_db):
        self.db = gffutils.FeatureDB(gene_db, keep_order=True)
        self.isoform_to_gene_map = {}
        self.gene_to_isoforms_map = {}
        self.gene_coords = {}

        self.parse_db()

    def parse_db(self):
        logger.info("Loading gene database")
        for g in self.db.features_of_type('gene'):
            gene_name = g.id
            self.gene_coords[gene_name] = (self.db[gene_name].start, self.db[gene_name].end)
            self.gene_to_isoforms_map[gene_name] = set()
            gene_db = self.db[gene_name]
            for t in self.db.children(gene_db, featuretype=('transcript', 'mRNA')):
                self.isoform_to_gene_map[t.id] = gene_name
                self.gene_to_isoforms_map[gene_name].add(t.id)
        logger.info("Gene database loaded: %d genes, %d transcripts" % (len(self.gene_to_isoforms_map), len(self.isoform_to_gene_map)))

    def get_gene_by_isoform(self, t_id):
        return self.isoform_to_gene_map[t_id]

    def get_isoforms_by_gene(self, gene_name):
        return self.gene_to_isoforms_map[gene_name]

    def get_gene_coords(self, gene_name):
        return self.gene_coords[gene_name]


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

    # venn2(subsets=(a, b, ab), set_labels=(label_a, label_b))
    # plt.savefig("venn.png")
    logger.info("%d correct reads by both methods, %d in %s only and %d in %s only" % (ab, a, label_a, b, label_b))
    logger.info("Venn diagram saved to venn.png")
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
    logger.info('Common assignments: %d, different assignments: %d, '
                'only %s assignments: %d, only %s assignments: %d, both not assigned: %d' %
                (ab, diff_ab, label_a, a, label_b, b, not_ab))

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="Output file name")
    parser.add_argument("--detailed_output", type=str, help="Output folder for detailed information on every read")
    parser.add_argument("--fasta", "-f", type=str, help="initial file with sequences")
    parser.add_argument("--mapping", "-m", type=str, help="mapped sequences (SAM or BAM format)") ## SQANTI2 output
    parser.add_argument("--tsv", "-t", type=str, nargs='+', help="assigned isoforms, max number of files to compare: 2")
    parser.add_argument("--gene_db", "-g", type=str, help="gene database")
    parser.add_argument("--tool", type=str, choices=['isoquant', 'talon', 'sqanti'], default='isoquant',
                        help="tool used for generating TSV (isoquant, talon, sqanti)")
    parser.add_argument("--isoform_id_index", type=int, default=0,
                        help="read id is split by underscore, indicate where isoform id is stored")

    parser.add_argument("--real_data", default=False, action="store_true",
                        help="real data is used (correct isoforms for each read are not provided)")
    parser.add_argument('--additional_stats', action='store_true', default=False,
                        help="count additional stats including ambiguous/inconsistent reads")

    args = parser.parse_args()
    if len(args.tsv) > 2:
        logger.error("ERROR! Maximum number of files to compare: 2")
        exit(-1)
    return args


def main():
    set_logger(logger)
    args = parse_args()
    db = DbHandler(args.gene_db)
    mapping_data = MappingData(args)
    all_stats = []
    output_file = sys.stdout if not args.output else open(args.output, "w")
    detailed_output_dir = None
    if args.detailed_output:
        if not os.path.exists(args.detailed_output):
            os.makedirs(args.detailed_output)
        detailed_output_dir = args.detailed_output

    for tsv_file in args.tsv:
        #TODO: diff bams
        logger.info("Calculating stats for %s" % tsv_file)
        assignment_data = AssignmentData(tsv_file, args.tool)
        prefix = os.path.splitext(os.path.basename(tsv_file))[0]
        stat_counter = StatCounter(prefix, mapping_data, assignment_data)
        if args.real_data:
            all_stats.append((stat_counter, prefix))
            continue

        logger.info("   Counting mapping stats...")
        stat_counter.count_mapping_stats(db)
        #stat_counter.count_mapping_stats(db, assignment_data.assigned_isoforms.keys())
        total_reads = len(stat_counter.mapping_data.seq_set)
        correctly_mapped = len(stat_counter.correct_seqs)
        secondary_mapped = len(stat_counter.secondary_mapped)
        mismapped_reads = len(stat_counter.mismapped_seqs)
        unmapped_reads = total_reads - (correctly_mapped + mismapped_reads)
        output_file.write("# MAPPING STATS\n")
        logger.info("Secondary mapped %d" % secondary_mapped)
        stat_counter.print_stats(correctly_mapped, mismapped_reads + secondary_mapped, unmapped_reads, output_file, name="mapping_")
        logger.info("   Done")

        # use only correctly mapped reads
        logger.info("   Counting pure assignment stats...")
        stat_counter.count_assignment_stats(db, output_dir=detailed_output_dir)
        correct_assignments = stat_counter.get_read_counts(ReadType.CORRECT)
        incorrect_assignments = stat_counter.get_read_counts(ReadType.INCORRECT_OTHER_GENE) + stat_counter.get_read_counts(ReadType.INCORRECT_SAME_GENE)
        unassigned_reads = stat_counter.get_read_counts(ReadType.NOT_ASSIGNED)
        output_file.write("# ASSIGNMENT STATS\n")
        stat_counter.print_stats(correct_assignments, incorrect_assignments, unassigned_reads, output_file, name="assignment_")
        logger.info("   Done")

        # use all reads
        logger.info("   Counting overall assignment stats...")
        stat_counter.count_assignment_stats(db, use_mismapped=True, output_dir=detailed_output_dir)
        correct_assignments = stat_counter.get_read_counts(ReadType.CORRECT)
        incorrect_assignments = stat_counter.get_read_counts(ReadType.INCORRECT_OTHER_GENE) + stat_counter.get_read_counts(ReadType.INCORRECT_SAME_GENE)
        unassigned_reads = stat_counter.get_read_counts(ReadType.NOT_ASSIGNED)
        output_file.write("# ASSIGNMENT STATS (INCLUDING MISMAPPED READS)\n")
        stat_counter.print_stats(correct_assignments, incorrect_assignments, unassigned_reads, output_file, name="overall_")
        all_stats.append((stat_counter, prefix))
        logger.info("   Done")

    if len(all_stats) == 2:
        if args.real_data:
            compare_real_results(all_stats[0], all_stats[1])
        else:
            compare_stats(all_stats[0], all_stats[1])

    if args.output:
        logger.info("Quality report saved to " + args.output)
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

