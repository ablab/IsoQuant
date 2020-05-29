############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import gffutils
import pysam

from src.input_data_storage import *
from src.alignment_processor import *
from src.assignment_io import *
from src.isoform_assignment import *
from src.gene_info import *
from src.long_read_counter import FeatureCounter

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
            for e in self.db.children(gene_db):
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

        self.correct_assignment_checker = PrintOnlyFunctor([ReadAssignmentType.unique, ReadAssignmentType.unique_minor_difference])
        self.novel_assignment_checker = PrintOnlyFunctor(ReadAssignmentType.contradictory)
        self.rest_assignment_checker = PrintOnlyFunctor([ReadAssignmentType.empty, ReadAssignmentType.ambiguous])

    def process_all_samples(self, input_data):
        logger.info("Processing " + proper_plural_form("sample", len(input_data.samples)))
        for sample in input_data.samples:
            self.process_sample(sample)
        logger.info("Processed " + proper_plural_form("sample", len(input_data.samples)))

    # Run though all genes in db and count stats according to alignments given in bamfile_name
    def process_sample(self, sample):
        logger.info("Processing sample " + sample.label)
        logger.info("Sample has " + proper_plural_form("BAM file", len(sample.file_list)) + ": " + ", ".join(map(lambda x: x[0], sample.file_list)))

        out_assigned_tsv = os.path.join(sample.out_dir, self.args.prefix + sample.label + ".assigned_reads.tsv")
        correct_printer = BasicTSVAssignmentPrinter(out_assigned_tsv, self.args,
                                                    assignment_checker=self.correct_assignment_checker)
        out_unmatched_tsv = os.path.join(sample.out_dir, self.args.prefix + sample.label + ".unmatched_reads.tsv")
        unmatched_printer = BasicTSVAssignmentPrinter(out_unmatched_tsv, self.args,
                                                      assignment_checker=self.rest_assignment_checker)
        out_alt_tsv = os.path.join(sample.out_dir, self.args.prefix + sample.label + ".altered_reads.tsv")
        alt_printer = BasicTSVAssignmentPrinter(out_alt_tsv, self.args,
                                                assignment_checker=self.novel_assignment_checker)
        global_printer = ReadAssignmentCompositePrinter([correct_printer, unmatched_printer, alt_printer])

        out_counts_tsv = os.path.join(sample.out_dir, self.args.prefix + sample.label + ".counts.tsv")
        feature_counter = FeatureCounter(self.gffutils_db, out_counts_tsv)

        current_chromosome = ""
        for g in self.gene_clusters:
            chr_id = g[0].seqid
            if chr_id != current_chromosome:
                logger.info("Processing chromosome " + chr_id)
                current_chromosome = chr_id

            gene_info = GeneInfo(g, self.gffutils_db)
            bam_files = list(map(lambda x: x[0], sample.file_list))
            alignment_processor = LongReadAlignmentProcessor(gene_info, bam_files, self.args, global_printer, feature_counter)
            alignment_processor.process()

        feature_counter.dump()
        logger.info("Finished processing sample " + sample.label)
        logger.info("Assigned reads are stored in " + out_assigned_tsv)
        logger.info("Unmatched reads are stored in " + out_unmatched_tsv)
        logger.info("Reads with alternative structure are stored in " + out_alt_tsv)
        logger.info("Read counts are stored in " + out_counts_tsv)
