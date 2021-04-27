############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import _pickle as pickle
import gc
import gzip
from multiprocessing import Pool
from collections import namedtuple

from src.input_data_storage import *
from src.alignment_processor import *
from src.assignment_io import *
from src.long_read_counter import *
from src.multimap_resolver import *
from src.read_groups import *
from src.transcript_model_constructor import *
from src.stats import *


logger = logging.getLogger('IsoQuant')


class GeneClusterConstructor:
    MAX_GENE_CLUSTER = 50
    MAX_GENE_LEN = 100000

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

            if len(current_gene_db_list) > 0 and \
                    (all(not genes_overlap(cg, gene_db) for cg in current_gene_db_list) or
                     (len(current_gene_db_list) > self.MAX_GENE_CLUSTER and
                      all(not genes_contain(cg, gene_db) for cg in current_gene_db_list))):
                gene_sets.append(current_gene_db_list)
                current_gene_db_list = []

            if g.end - g.start > self.MAX_GENE_LEN:
                gene_sets.append([gene_db])
            else:
                current_gene_db_list.append(gene_db)

        if current_gene_db_list:
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
                # non-overlapping gene
                gene_sets.append([gene_db])
                gene_exon_sets.append(gene_exons)
            elif len(overlapping_sets) == 1:
                # overlaps with 1 gene
                index = overlapping_sets[0]
                gene_sets[index].append(gene_db)
                gene_exon_sets[index].update(gene_exons)
            else:
                # merge all overlapping genes
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


def assign_reads_in_parallel(sample, chr_id, cluster, args, read_grouper, current_chr_record):
    tmp_printer = TmpFileAssignmentPrinter("{}_{}".format(sample.out_raw_file, chr_id), args)
    processed_reads = []
    bam_files = list(map(lambda x: x[0], sample.file_list))
    bam_files = [pysam.AlignmentFile(bam, "rb") for bam in bam_files]
    logger.info("Processing chromosome " + chr_id)
    gffutils_db = gffutils.FeatureDB(args.genedb, keep_order=True)
    for g in cluster:
        if len(g) > 100:
            logger.debug("Potential slowdown in %s due to large gene cluster of size %d" % (chr_id, len(g)))
        gene_info = GeneInfo(g, gffutils_db, args.delta)
        alignment_processor = LongReadAlignmentProcessor(gene_info, bam_files, args,
                                                         current_chr_record, read_grouper)
        assignment_storage = alignment_processor.process()
        gene_info.db = None
        tmp_printer.add_gene_info(gene_info)
        for read_assignment in assignment_storage:
            tmp_printer.add_read_info(read_assignment)
            processed_reads.append(ShortReadAssignment(read_assignment, gene_info))
    logger.info("Finished processing chromosome " + chr_id)
    return processed_reads


# Class for processing all samples against gene database
class DatasetProcessor:
    def __init__(self, args):
        self.args = args
        logger.info("Loading gene database from " + self.args.genedb)
        self.gffutils_db = gffutils.FeatureDB(self.args.genedb, keep_order=True)
        if self.args.needs_reference:
            logger.info("Loading reference genome from " + self.args.reference)
            _, outer_ext = os.path.splitext(self.args.reference)
            if outer_ext.lower() in ['.gz', '.gzip']:
                with gzip.open(self.args.reference, "rt") as handle:
                    self.reference_record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
            else:
                self.reference_record_dict = SeqIO.to_dict(SeqIO.parse(self.args.reference, "fasta"))
        else:
            self.reference_record_dict = None
        self.gene_cluster_constructor = GeneClusterConstructor(self.gffutils_db)
        self.gene_clusters = self.gene_cluster_constructor.get_gene_sets()
        self.read_grouper = create_read_grouper(args)
        self.io_support = IOSupport(self.args)
        self.multimapped_reads = defaultdict(list)

    def process_all_samples(self, input_data):
        logger.info("Processing " + proper_plural_form("sample", len(input_data.samples)))
        for sample in input_data.samples:
            self.process_sample(sample)
        logger.info("Processed " + proper_plural_form("sample", len(input_data.samples)))

    # Run though all genes in db and count stats according to alignments given in bamfile_name
    def process_sample(self, sample):
        logger.info("Processing sample " + sample.label)
        logger.info("Sample has " + proper_plural_form("BAM file", len(sample.file_list)) + ": " + ", ".join(map(lambda x: x[0], sample.file_list)))
        self.multimapped_reads = defaultdict(list)

        self.tmp_dir = os.path.join(sample.out_dir, "tmp")
        if not os.path.isdir(self.tmp_dir):
            os.makedirs(self.tmp_dir)

        if self.args.read_assignments:
            saves_file = sample.file_list[0][0]
            logger.info('Using read assignments from {}*'.format(saves_file))
            total_alignments, polya_found = self.load_reads(saves_file)
        else:
            self.assign_reads(sample)
            total_alignments, polya_found = self.load_reads(sample.out_raw_file)
            logger.info('Read assignments files saved to {}*.\nYou can reuse it later with option --read_assignments'.
                        format(sample.out_raw_file))

        intial_polya_required = self.args.require_polyA
        polya_fraction = polya_found / total_alignments if total_alignments > 0 else 0.0
        logger.info("Total alignments processed: %d, polyA tail detected in %d (%.1f%%)" %
                    (total_alignments, polya_found, polya_fraction * 100.0))
        if polya_fraction < 0.1 and self.args.require_polyA:
            logger.warning("PolyA tail found in less than 10% of the reads, "
                           "polyA tails will not be required for transcript model construction. "
                           "Set --polya_trimmed to avoid this warning in future.")
            self.args.require_polyA = False

        logger.info("Combining output for sample " + sample.label)
        self.aggregate_reads(sample)
        self.construct_models(sample)
        os.rmdir(self.tmp_dir)
        self.args.require_polyA = intial_polya_required
        logger.info("Processed sample " + sample.label)

    def assign_reads(self, sample):
        logger.info('Assigning reads to isoforms')
        chrom_clusters = self.get_chromosome_gene_clusters()
        pool = Pool(self.args.threads)
        processed_reads = pool.starmap(assign_reads_in_parallel, [(sample, chr_id, c, self.args, self.read_grouper,
                                                                   (self.reference_record_dict[
                                                                        chr_id] if self.reference_record_dict else None))
                                                                  for (chr_id, c) in chrom_clusters])
        pool.close()
        pool.join()

        logger.info("Resolving multimappers")
        all_processed_reads = set()
        self.multimapped_reads = defaultdict(list)
        for storage in processed_reads:
            for read_assignment in storage:
                if read_assignment.read_id in all_processed_reads:
                    self.multimapped_reads[read_assignment.read_id].append(read_assignment)
                else:
                    all_processed_reads.add(read_assignment.read_id)

        multimap_resolver = MultimapResolver(self.args.multimap_strategy)
        multimap_pickler = pickle.Pickler(open(sample.out_raw_file + "_multimappers", "wb"),  -1)
        multimap_pickler.fast = True
        for assignment_list in self.multimapped_reads.values():
            multimap_resolver.resolve(assignment_list)
            multimap_pickler.dump(assignment_list)
        logger.info('Finishing read assignment')

    def load_reads(self, dump_filename):
        self.reads_assignments = []

        read_storage = []
        total_assignments = 0
        polya_found = 0
        gene_info = None
        gc.disable()

        if not self.multimapped_reads:
            multimap_unpickler = pickle.Unpickler(open(dump_filename + "_multimappers", "rb"), fix_imports=False)
            while True:
                try:
                    obj = multimap_unpickler.load()
                    if isinstance(obj, list):
                        assignment_list = obj
                        read_id = assignment_list[0].read_id
                        self.multimapped_reads[read_id] = assignment_list
                    else:
                        raise ValueError("Multimap assignment file {} is corrupted!".format(dump_filename))
                except EOFError:
                    break

        for chr_id in self.get_chromosome_list():
            chr_dump_file = dump_filename + "_" + chr_id
            logger.info("Loading read assignments from " + chr_dump_file)
            if not os.path.exists(chr_dump_file): continue
            with open(chr_dump_file, "rb") as f:
                p = pickle.Unpickler(f, fix_imports=False)
                while True:
                    try:
                        obj = p.load()
                        if isinstance(obj, ReadAssignment):
                            read_assignment = obj
                            assert gene_info is not None
                            read_assignment.gene_info = gene_info
                            if read_assignment.read_id in self.multimapped_reads:
                                resolved_assignment = None
                                for a in self.multimapped_reads[read_assignment.read_id]:
                                    if a.start == read_assignment.start() and a.end == read_assignment.end() and \
                                            a.gene_id == gene_info.gene_db_list[0].id and a.chr_id == read_assignment.chr_id:
                                        if resolved_assignment is not None:
                                            logger.warning("Duplicate read: %s %s %s" % (read_assignment.read_id, a.gene_id, a.chr_id))
                                        resolved_assignment = a

                                if not resolved_assignment or resolved_assignment.assignment_type == ReadAssignmentType.noninformative:
                                    continue
                                read_assignment.assignment_type = resolved_assignment.assignment_type
                                read_assignment.multimapper = resolved_assignment.multimapper
                            read_storage.append(read_assignment)
                            total_assignments += 1
                            if read_assignment.polyA_found:
                                polya_found += 1
                        elif isinstance(obj, GeneInfo):
                            if gene_info and read_storage:
                                self.reads_assignments.append((gene_info, read_storage))
                                read_storage = []
                            gene_info = obj
                            gene_info.db = self.gffutils_db
                        else:
                            raise ValueError("Read assignment file {} is corrupted!".format(chr_dump_file))
                    except EOFError:
                        break
        gc.enable()
        if gene_info and read_storage:
            self.reads_assignments.append((gene_info, read_storage))
        return total_assignments, polya_found

    def get_chromosome_list(self):
        chromosomes = []
        current_chromosome = ""
        for g in self.gene_clusters:
            chr_id = g[0].seqid
            if chr_id != current_chromosome:
                chromosomes.append(chr_id)
                current_chromosome = chr_id
        return chromosomes

    def get_chromosome_gene_clusters(self):
        chrom_clusters = []
        cur_cluster = []
        current_chromosome = ""
        for g in self.gene_clusters:
            chr_id = g[0].seqid
            if chr_id != current_chromosome:
                if cur_cluster:
                    chrom_clusters.append((current_chromosome, cur_cluster))
                    cur_cluster = []
                current_chromosome = chr_id
            cur_cluster.append(g)
        if cur_cluster:
            chrom_clusters.append((current_chromosome, cur_cluster))

        # chromosomes with large clusters take more time and should go first
        chrom_clusters = sorted(chrom_clusters, key=lambda x: sum(len(cluster) ** 2 for cluster in x[1]), reverse=True)
        return chrom_clusters

    def create_aggregators(self, sample):
        self.read_stat_counter = EnumStats()
        self.basic_printer = BasicTSVAssignmentPrinter(sample.out_assigned_tsv, self.args, self.io_support)
        self.bed_printer = BEDPrinter(sample.out_mapped_bed, self.args)
        printer_list = [self.basic_printer, self.bed_printer]
        if self.args.sqanti_output:
            self.sqanti_printer = SqantiTSVPrinter(sample.out_alt_tsv, self.args, self.io_support)
            printer_list.append(self.sqanti_printer)
        self.global_printer = ReadAssignmentCompositePrinter(printer_list)

        self.gene_counter = create_gene_counter(sample.out_gene_counts_tsv, ignore_read_groups=True)
        self.transcript_counter = create_transcript_counter(sample.out_transcript_counts_tsv, ignore_read_groups=True)
        if self.args.count_exons:
            self.exon_counter = ExonCounter(sample.out_exon_counts_tsv, ignore_read_groups=True)
            self.intron_counter = IntronCounter(sample.out_intron_counts_tsv, ignore_read_groups=True)
            self.global_counter = CompositeCounter([self.gene_counter, self.transcript_counter,
                                                    self.exon_counter, self.intron_counter])
        else:
            self.global_counter = CompositeCounter([self.gene_counter, self.transcript_counter])

        if self.args.read_group:
            self.gene_grouped_counter = create_gene_counter(sample.out_gene_grouped_counts_tsv)
            self.transcript_grouped_counter = create_transcript_counter(sample.out_transcript_grouped_counts_tsv)
            if self.args.count_exons:
                self.exon_grouped_counter = ExonCounter(sample.out_exon_grouped_counts_tsv)
                self.intron_grouped_counter = IntronCounter(sample.out_intron_grouped_counts_tsv)
                self.global_counter.add_counters([self.gene_grouped_counter, self.transcript_grouped_counter,
                                                  self.exon_grouped_counter, self.intron_grouped_counter])
            else:
                self.global_counter.add_counters([self.gene_grouped_counter, self.transcript_grouped_counter])

    def pass_to_aggregators(self, read_assignment):
        if read_assignment is None:
            return
        self.read_stat_counter.add(read_assignment.assignment_type)
        self.global_printer.add_read_info(read_assignment)
        self.global_counter.add_read_info(read_assignment)

    def finalize_aggregators(self, sample):
        self.global_counter.dump()
        logger.info("Gene counts are stored in " + self.gene_counter.output_file_name)
        logger.info("Transcript counts are stored in " + self.transcript_counter.output_file_name)
        logger.info("Read assignments are stored in " + self.basic_printer.output_file_name)
        self.read_stat_counter.print_start("Read assignment statistics")

    def aggregate_reads(self, sample):
        self.create_aggregators(sample)

        if self.args.memory_efficient:
            assert False
        else:
            for storage in self.reads_assignments:
                gene_read_assignments = storage[1]
                for read_assignment in gene_read_assignments:
                    self.pass_to_aggregators(read_assignment)

        self.finalize_aggregators(sample)

    def construct_models(self, sample):
        logger.info("Constructing transcript models for " + sample.label)
        gff_printer = GFFPrinter(sample.out_dir, sample.label, self.io_support)
        transcript_stat_counter = EnumStats()

        for gene_info, assignment_storage in self.reads_assignments:
            transcript_generator = TranscriptModelConstructor(gene_info, assignment_storage, self.args)
            transcript_generator.process()
            gff_printer.dump(transcript_generator)
            for t in transcript_generator.transcript_model_storage:
                transcript_stat_counter.add(t.transcript_type)

        logger.info("Transcript model file " + gff_printer.model_fname)
        transcript_stat_counter.print_start("Transcript model statistics")
