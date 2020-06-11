############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import gffutils
import pysam
from collections import defaultdict

from src.input_data_storage import *
from src.alignment_processor import *
from src.assignment_io import *
from src.isoform_assignment import *
from src.gene_info import *
from src.long_read_counter import *
from src.multimap_resolver import *
from src.read_groups import *

logger = logging.getLogger('IsoQuant')

class GeneClusterConstructor:
    def __init__(self, gene_db):
        self.gene_db = gene_db
        self.gene_sets = None

    def get_gene_sets(self):
        if self.gene_sets is None:
            self.gene_sets = self.fill_gene_sets()
        return self.gene_sets

    def fill_gene_sets(self):
        gene_sets = []
        current_gene_db_list = []
        for g in self.gene_db.features_of_type('gene', order_by=('seqid', 'start')):
            gene_name = g.id
            gene_db = self.gene_db[gene_name]
            if len(current_gene_db_list) == 0 or any(genes_overlap(cg, gene_db) for cg in current_gene_db_list):
                current_gene_db_list.append(gene_db)
            else:
                gene_sets.append(current_gene_db_list)
                current_gene_db_list = [gene_db]

        gene_sets.append(current_gene_db_list)
        return gene_sets


class OverlappingExonsGeneClusterConstructor(GeneClusterConstructor):
    def get_gene_sets(self):
        if self.gene_sets is None:
            self.gene_sets = self.cluster_on_shared_exons()
        return self.gene_sets

    def cluster_on_shared_exons(self):
        gene_clusters = self.fill_gene_sets()
        new_gene_sets = []

        for cluster in gene_clusters:
            new_gene_sets += self.split_cluster(cluster)
        return new_gene_sets

    def split_cluster(self, gene_cluster):
        gene_sets = []
        gene_exon_sets = []

        for gene_db in gene_cluster:
            gene_exons = set()
            for e in self.gene_db.children(gene_db):
                if e.featuretype == 'exon':
                    gene_exons.add((e.start, e.end))

            overlapping_sets = []
            for i in range(len(gene_exon_sets)):
                if any(e in gene_exon_sets[i] for e in gene_exons):
                    overlapping_sets.append(i)

            if len(overlapping_sets) == 0:
                #non-overlapping gene
                gene_sets.append([gene_db])
                gene_exon_sets.append(gene_exons)
            elif len(overlapping_sets) == 1:
                #overlaps with 1 gene
                index = overlapping_sets[0]
                gene_sets[index].append(gene_db)
                gene_exon_sets[index].update(gene_exons)
            else:
                #merge all overlapping genes
                new_gene_set = [gene_db]
                new_exons_set = gene_exons
                for index in overlapping_sets:
                    new_gene_set += gene_sets[index]
                    new_exons_set.update(gene_exon_sets[index])
                for index in overlapping_sets[::-1]:
                    del gene_sets[index]
                    del gene_exon_sets[index]
                gene_sets.append(new_gene_set)
                gene_exon_sets.append(new_exons_set)

        return gene_sets

# Class for processing all samples against gene database
class DatasetProcessor:
    def __init__(self, args):
        self.args = args
        logger.info("Loading gene database from " + self.args.genedb)
        self.gffutils_db = gffutils.FeatureDB(self.args.genedb, keep_order=True)
        self.gene_cluster_constructor = GeneClusterConstructor(self.gffutils_db)
        self.gene_clusters = self.gene_cluster_constructor.get_gene_sets()
        self.read_grouper = create_read_grouper(args.read_group)

        #self.correct_assignment_checker = PrintOnlyFunctor([ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference])
        #self.novel_assignment_checker = PrintOnlyFunctor(ReadAssignmentType.contradictory)
        #self.rest_assignment_checker = PrintOnlyFunctor([ReadAssignmentType.empty, ReadAssignmentType.ambiguous])

    def process_all_samples(self, input_data):
        logger.info("Processing " + proper_plural_form("sample", len(input_data.samples)))
        for sample in input_data.samples:
            self.process_sample(sample)
        logger.info("Processed " + proper_plural_form("sample", len(input_data.samples)))

    # Run though all genes in db and count stats according to alignments given in bamfile_name
    def process_sample(self, sample):
        logger.info("Processing sample " + sample.label)
        logger.info("Sample has " + proper_plural_form("BAM file", len(sample.file_list)) + ": " + ", ".join(map(lambda x: x[0], sample.file_list)))

        self.tmp_dir = os.path.join(sample.out_dir, "tmp")
        if not os.path.isdir(self.tmp_dir):
            os.makedirs(self.tmp_dir)

        current_chromosome = ""
        counter = 0
        self.processed_reads = set()
        self.multimapped_reads = set()
        self.reads_assignments = []
        for g in self.gene_clusters:
            chr_id = g[0].seqid
            if chr_id != current_chromosome:
                logger.info("Processing chromosome " + chr_id)
                current_chromosome = chr_id

            gene_info = GeneInfo(g, self.gffutils_db, self.args.delta)
            bam_files = list(map(lambda x: x[0], sample.file_list))
            alignment_processor = LongReadAlignmentProcessor(gene_info, bam_files, self.args, self.read_grouper)
            assignment_storage = alignment_processor.process()
            self.dump_reads(assignment_storage, counter)
            counter += 1

        logger.info("Combining output")
        self.aggregate_reads(sample)
        os.rmdir(self.tmp_dir)
        logger.info("Processed sample " + sample.label)

    def dump_reads(self, read_storage, gene_counter):
        if self.args.memory_efficient:
            # TODO: dump to file
            out_tmp_tsv = os.path.join(self.tmp_dir, str(gene_counter) + ".processed_reads.tsv")
            assert (False)
        else:
            for read_assignment in read_storage:
                if read_assignment.assignment_type is not None:
                    if read_assignment.read_id in self.processed_reads:
                        self.multimapped_reads.add(read_assignment.read_id)
                    else:
                        self.processed_reads.add(read_assignment.read_id)
            self.reads_assignments.append(read_storage)

    def create_aggregators(self, sample):
        out_assigned_tsv = os.path.join(sample.out_dir, self.args.prefix + sample.label + ".read_assignments.tsv")
        self.basic_printer = BasicTSVAssignmentPrinter(out_assigned_tsv, self.args)
        out_alt_tsv = os.path.join(sample.out_dir, self.args.prefix + sample.label + ".SQANTI-like.tsv")
        self.sqanti_printer = SqantiTSVPrinter(out_alt_tsv, self.args)
        self.global_printer = ReadAssignmentCompositePrinter([self.basic_printer, self.sqanti_printer])

        out_gene_counts_tsv = os.path.join(sample.out_dir, self.args.prefix + sample.label + ".gene_counts.tsv")
        self.gene_counter = create_gene_counter(out_gene_counts_tsv)
        out_transcript_counts_tsv = os.path.join(sample.out_dir, self.args.prefix + sample.label + ".transcript_counts.tsv")
        self.transcript_counter = create_transcript_counter(out_transcript_counts_tsv)
        # TODO make optional ?
        out_exon_counts_tsv = os.path.join(sample.out_dir, self.args.prefix + sample.label + ".exon_counts.tsv")
        self.exon_counter = ExonCounter(out_exon_counts_tsv)
        out_intron_counts_tsv = os.path.join(sample.out_dir,
                                                 self.args.prefix + sample.label + ".intron_counts.tsv")
        self.intron_counter = IntronCounter(out_intron_counts_tsv)
        self.global_counter = CompositeCounter([self.gene_counter, self.transcript_counter, self.exon_counter, self.intron_counter])

    def pass_to_aggregators(self, read_assignment):
        self.global_printer.add_read_info(read_assignment)
        self.global_counter.add_read_info(read_assignment)

    def finalize_aggregators(self, sample):
        self.global_counter.dump()
        logger.info("Finished processing sample " + sample.label)
        logger.info("Gene counts are stored in " + self.gene_counter.output_file_name)
        logger.info("Transcript counts are stored in " + self.transcript_counter.output_file_name)

    def aggregate_reads(self, sample):
        self.create_aggregators(sample)

        multimap_reads_assignments = defaultdict(list)
        if self.args.memory_efficient:
            assert (False)
        else:
            for storage in self.reads_assignments:
                for read_assignment in storage:
                    read_id = read_assignment.read_id
                    if read_id in self.multimapped_reads:
                        multimap_reads_assignments[read_id].append(read_assignment)
                    else:
                        self.pass_to_aggregators(read_assignment)

        #  TODO: resolve multimappers
        multimap_resover = MultimapResolver(self.args.matching_stategy)
        for read_id in multimap_reads_assignments.keys():
            read_assignment = multimap_resover.resolve(multimap_reads_assignments[read_id])
            self.pass_to_aggregators(read_assignment)

        self.finalize_aggregators(sample)
