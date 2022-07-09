############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import _pickle as pickle
import gc
import glob
import gzip
import os
from multiprocessing import Pool

from src.file_utils import *
from src.input_data_storage import *
from src.alignment_processor import *
from src.long_read_counter import *
from src.multimap_resolver import *
from src.read_groups import *
from src.transcript_printer import *
from src.stats import *
from src.graph_based_model_construction import *


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
    bam_file_pairs = [(pysam.AlignmentFile(bam, "rb"), bam) for bam in bam_files]
    logger.info("Processing chromosome " + chr_id)
    gffutils_db = gffutils.FeatureDB(args.genedb, keep_order=True)
    for g in cluster:
        if len(g) > 100:
            logger.debug("Potential slowdown in %s due to large gene cluster of size %d" % (chr_id, len(g)))
        gene_info = GeneInfo(g, gffutils_db, args.delta)
        alignment_processor = LongReadAlignmentProcessor(gene_info, bam_file_pairs, args,
                                                         current_chr_record, read_grouper)
        assignment_storage = alignment_processor.process()
        gene_info.db = None
        tmp_printer.add_gene_info(gene_info)
        for read_assignment in assignment_storage:
            tmp_printer.add_read_info(read_assignment)
            processed_reads.append(BasicReadAssignment(read_assignment, gene_info))
    logger.info("Finished processing chromosome " + chr_id)
    return processed_reads, read_grouper.read_groups


def collect_reads_in_parallel(sample, chr_id, args, read_grouper, current_chr_record):
    tmp_printer = TmpFileAssignmentPrinter("{}_{}".format(sample.out_raw_file, chr_id), args)
    processed_reads = []
    bam_files = list(map(lambda x: x[0], sample.file_list))
    bam_file_pairs = [(pysam.AlignmentFile(bam, "rb"), bam) for bam in bam_files]
    gffutils_db = gffutils.FeatureDB(args.genedb, keep_order=True) if args.genedb else None

    logger.info("Processing chromosome " + chr_id)
    alignment_collector = IntergenicAlignmentCollector(chr_id, bam_file_pairs, args, gffutils_db, current_chr_record, read_grouper)

    for gene_info, assignment_storage in alignment_collector.process():
        tmp_printer.add_gene_info(gene_info)
        for read_assignment in assignment_storage:
            tmp_printer.add_read_info(read_assignment)
            processed_reads.append(BasicReadAssignment(read_assignment, gene_info))
    logger.info("Finished processing chromosome " + chr_id)
    return processed_reads, read_grouper.read_groups


class ReadAssignmentLoader:
    def __init__(self, save_file_name, gffutils_db, multimappers_dict):
        logger.info("Loading read assignments from " + save_file_name)
        assert os.path.exists(save_file_name)
        self.save_file_name = save_file_name
        self.unpickler = pickle.Unpickler(open(save_file_name, "rb"), fix_imports=False)
        self.current_gene_info_obj = None
        self.is_updated = False
        self.gffutils_db = gffutils_db
        self.multimapped_reads = multimappers_dict


def load_assigned_reads(save_file_name, gffutils_db, multimapped_chr_dict):
    gc.disable()
    logger.info("Loading read assignments from " + save_file_name)
    assert os.path.exists(save_file_name)
    save_file_name = save_file_name
    unpickler = pickle.Unpickler(open(save_file_name, "rb"), fix_imports=False)
    read_storage = []
    current_gene_info = None

    while True:
        try:
            obj = unpickler.load()
            if isinstance(obj, ReadAssignment):
                read_assignment = obj
                assert current_gene_info is not None
                read_assignment.gene_info = current_gene_info
                if read_assignment.read_id in multimapped_chr_dict:
                    resolved_assignment = None
                    for a in multimapped_chr_dict[read_assignment.read_id]:
                        if a.start == read_assignment.start() and a.end == read_assignment.end() and \
                                a.gene_id == current_gene_info.gene_db_list[0].id and \
                                a.matches == (0 if not read_assignment.isoform_matches else len(read_assignment.isoform_matches)) and \
                                a.chr_id == read_assignment.chr_id:
                            if resolved_assignment is not None:
                                logger.debug("Duplicate read: %s %s %s" % (read_assignment.read_id, a.gene_id, a.chr_id))
                            resolved_assignment = a

                    if not resolved_assignment:
                        logger.warning("Incomplete information on read %s" % read_assignment.read_id)
                        continue
                    elif resolved_assignment.assignment_type == ReadAssignmentType.suspended:
                        continue
                    else:
                        read_assignment.assignment_type = resolved_assignment.assignment_type
                        read_assignment.multimapper = resolved_assignment.multimapper
                read_storage.append(read_assignment)
            elif isinstance(obj, GeneInfo):
                if current_gene_info and read_storage:
                    yield current_gene_info, read_storage
                read_storage = []
                current_gene_info = obj
                current_gene_info.db = gffutils_db
            else:
                raise ValueError("Read assignment file {} is corrupted!".format(save_file_name))
        except EOFError:
            break
    gc.enable()
    if current_gene_info and read_storage:
        yield current_gene_info, read_storage


def construct_models_in_parallel(sample, chr_id, dump_filename, args, multimapped_reads, read_grouper, current_chr_record,
                                 io_support):
    processor = DatasetProcessor(args, parallel=True)
    processor.read_grouper = read_grouper
    processor.create_aggregators(sample)
    if args.genedb:
        gffutils_db = gffutils.FeatureDB(args.genedb, keep_order=True)
    else:
        gffutils_db = None

    # debug info only
    expressed_db = None
    if args.expressed_genedb:
        logger.info("Loading expreseed genes from %s" % args.expressed_genedb)
        expressed_db = gffutils.FeatureDB(args.expressed_genedb, keep_order=True)

    logger.info("Processing chromosome " + chr_id)
    transcripts = []
    tmp_gff_printer = GFFPrinter(sample.out_dir, sample.label, io_support)
    chr_dump_file = dump_filename + "_" + chr_id
    for gene_info, assignment_storage in load_assigned_reads(chr_dump_file, gffutils_db, multimapped_reads):
        for read_assignment in assignment_storage:
            if read_assignment is None:
                continue
            processor.read_stat_counter.add(read_assignment.assignment_type)
            processor.global_printer.add_read_info(read_assignment)
            processor.global_counter.add_read_info(read_assignment)

        if not args.no_model_construction:
            expressed_gene_info = None
            if expressed_db:
                gene_list = []
                if gffutils_db:
                    for g in gene_info.gene_db_list:
                        try:
                            gene_list.append(expressed_db[g.id])
                        except gffutils.exceptions.FeatureNotFoundError:
                            pass
                else:
                    gene_list = list(expressed_db.region(region=(gene_info.chr_id, gene_info.start, gene_info.end),
                                                         completely_within=False, featuretype='gene'))
                if gene_list:
                    expressed_gene_info = GeneInfo(gene_list, expressed_db, delta=args.delta)


            model_constructor = GraphBasedModelConstructor(gene_info, current_chr_record, args, processor.transcript_model_global_counter,
                                                           expressed_gene_info=expressed_gene_info)
            model_constructor.process(assignment_storage)
            tmp_gff_printer.dump(model_constructor)
            for t in model_constructor.transcript_model_storage:
                transcripts.append(t.transcript_type)

    logger.info("Finished processing chromosome " + chr_id)
    processor.global_counter.dump()
    processor.transcript_model_global_counter.dump()
    return processor.read_stat_counter, transcripts


# Class for processing all samples against gene database
class DatasetProcessor:
    def __init__(self, args, parallel=False):
        self.args = args
        self.args.gunzipped_reference = None
        self.common_header = "# Command line: " + args._cmd_line + "\n# IsoQuant version: " + args._version + "\n"
        self.read_grouper = create_read_grouper(args)
        self.io_support = IOSupport(self.args)
        if parallel: return
        if args.genedb:
            logger.info("Loading gene database from " + self.args.genedb)
            self.gffutils_db = gffutils.FeatureDB(self.args.genedb, keep_order=True)
            self.gene_cluster_constructor = GeneClusterConstructor(self.gffutils_db)
            self.gene_clusters = self.gene_cluster_constructor.get_gene_sets()
        else:
            self.gffutils_db = None
            self.gene_clusters = None

        if self.args.needs_reference:
            logger.info("Loading reference genome from " + self.args.reference)
            ref_name, outer_ext = os.path.splitext(os.path.basename(self.args.reference))
            low_ext = outer_ext.lower()
            if low_ext in ['.gz', '.gzip', '.bgz']:
                with gzip.open(self.args.reference, "rt") as handle:
                    if self.args.low_memory:
                        if low_ext == '.bgz':
                            self.reference_record_dict = SeqIO.index(self.args.reference, "fasta")
                        else:
                            self.args.gunzipped_reference = os.path.join(args.output, ref_name)
                            with open(self.args.gunzipped_reference, "w") as outf:
                                shutil.copyfileobj(handle, outf)
                            logger.info("Loading uncompressed reference from " + self.args.gunzipped_reference)
                            self.reference_record_dict = SeqIO.index(self.args.gunzipped_reference, "fasta")
                    else:
                        self.reference_record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
            else:
                if self.args.low_memory:
                    self.reference_record_dict = SeqIO.index(self.args.reference, "fasta")
                else:
                    self.reference_record_dict = SeqIO.to_dict(SeqIO.parse(self.args.reference, "fasta"))
        else:
            self.reference_record_dict = None

        self.multimapped_reads = defaultdict(list)
        # chr_id -> read_id -> list of assignments
        self.multimapped_info_dict = defaultdict(lambda : defaultdict(list))

    def __del__(self):
        if not self.args.keep_tmp and self.args.gunzipped_reference:
            if os.path.exists(self.args.gunzipped_reference):
                os.remove(self.args.gunzipped_reference)

    def process_all_samples(self, input_data):
        logger.info("Processing " + proper_plural_form("sample", len(input_data.samples)))
        for sample in input_data.samples:
            self.process_sample(sample)
        logger.info("Processed " + proper_plural_form("sample", len(input_data.samples)))

    # Run though all genes in db and count stats according to alignments given in bamfile_name
    def process_sample(self, sample):
        logger.info("Processing sample " + sample.label)
        logger.info("Sample has " + proper_plural_form("BAM file", len(sample.file_list)) + ": " + ", ".join(map(lambda x: x[0], sample.file_list)))
        self.args.use_technical_replicas = self.args.read_group == "file_name" and len(sample.file_list) > 1
        self.multimapped_reads = defaultdict(list)

        if self.args.read_assignments:
            saves_file = self.args.read_assignments[0]
            logger.info('Using read assignments from {}*'.format(saves_file))
        else:
            self.collect_reads(sample)
            saves_file = sample.out_raw_file
            logger.info('Read assignments files saved to {}*. '.
                        format(sample.out_raw_file))
            if not self.args.keep_tmp:
                logger.info("To keep these intermediate files for debug purposes use --keep_tmp flag")

        total_alignments, polya_found, all_read_groups = self.load_read_info(saves_file)

        polya_fraction = polya_found / total_alignments if total_alignments > 0 else 0.0
        logger.info("Total alignments processed: %d, polyA tail detected in %d (%.1f%%)" %
                    (total_alignments, polya_found, polya_fraction * 100.0))
        self.args.needs_polya_for_construction = polya_fraction >= 0.7

        self.read_grouper.read_groups = all_read_groups
        self.process_assigned_reads(sample, saves_file)
        if not self.args.read_assignments and not self.args.keep_tmp:
            for f in glob.glob(saves_file + "_*"):
                os.remove(f)
        logger.info("Processed sample " + sample.label)

    def collect_reads(self, sample):
        logger.info('Collecting read alignments')
        pool = Pool(self.args.threads)
        chr_ids = sorted(self.reference_record_dict.keys(), key=lambda x: len(self.reference_record_dict[x]), reverse=True)
        #chr_ids = ['chr1']
        results = pool.starmap(collect_reads_in_parallel, [(sample, chr_id, self.args, self.read_grouper,
                                                            (self.reference_record_dict[chr_id] if self.reference_record_dict else None))
                                                           for chr_id in chr_ids], chunksize=1)
        pool.close()
        pool.join()
        processed_reads, read_groups = [x[0] for x in results], [x[1] for x in results]

        logger.info("Resolving multimappers")
        self.multimapped_reads = defaultdict(list)
        for storage in processed_reads:
            for read_assignment in storage:
                self.multimapped_reads[read_assignment.read_id].append(read_assignment)

        multimap_resolver = MultimapResolver(self.args.multimap_strategy)
        multimap_pickler = pickle.Pickler(open(sample.out_raw_file + "_multimappers", "wb"),  -1)
        multimap_pickler.fast = True
        total_assignments = 0
        polya_assignments = 0

        for read_id in list(self.multimapped_reads.keys()):
            assignment_list = self.multimapped_reads[read_id]
            if len(assignment_list) == 1:
                total_assignments += 1
                polya_assignments += 1 if assignment_list[0].polyA_found else 0
                del self.multimapped_reads[read_id]
                continue
            multimap_resolver.resolve(assignment_list)
            multimap_pickler.dump(assignment_list)
            for a in assignment_list:
                if a.assignment_type != ReadAssignmentType.suspended:
                    total_assignments += 1
                    if a.polyA_found:
                        polya_assignments += 1

        info_pickler = pickle.Pickler(open(sample.out_raw_file + "_info", "wb"),  -1)
        info_pickler.dump(total_assignments)
        info_pickler.dump(polya_assignments)
        all_read_groups = set([r for sublist in read_groups for r in sublist])
        info_pickler.dump(all_read_groups)
        if total_assignments == 0:
            logger.warning("No reads were assigned to isoforms, check your input files")
        else:
            logger.info('Finishing read assignment, total assignments %d, polyA percentage %.1f' %
                        (total_assignments, 100 * polya_assignments / total_assignments))

    def assign_reads(self, sample):
        logger.info('Assigning reads to isoforms')
        chrom_clusters = self.get_chromosome_gene_clusters()
        pool = Pool(self.args.threads)
        results = pool.starmap(assign_reads_in_parallel, [(sample, chr_id, c, self.args, self.read_grouper,
                                                                   (self.reference_record_dict[
                                                                        chr_id] if self.reference_record_dict else None))
                                                                  for (chr_id, c) in chrom_clusters], chunksize=1)
        pool.close()
        pool.join()
        processed_reads, read_groups = [x[0] for x in results], [x[1] for x in results]

        logger.info("Resolving multimappers")
        self.multimapped_reads = defaultdict(list)
        for storage in processed_reads:
            for read_assignment in storage:
                self.multimapped_reads[read_assignment.read_id].append(read_assignment)

        multimap_resolver = MultimapResolver(self.args.multimap_strategy)
        multimap_pickler = pickle.Pickler(open(sample.out_raw_file + "_multimappers", "wb"),  -1)
        multimap_pickler.fast = True
        total_assignments = 0
        polya_assignments = 0
        for read_id in list(self.multimapped_reads.keys()):
            assignment_list = self.multimapped_reads[read_id]
            if len(assignment_list) == 1:
                total_assignments += 1
                polya_assignments += 1 if assignment_list[0].polyA_found else 0
                del self.multimapped_reads[read_id]
                continue
            multimap_resolver.resolve(assignment_list)
            multimap_pickler.dump(assignment_list)
            for a in assignment_list:
                if a.assignment_type != ReadAssignmentType.suspended:
                    total_assignments += 1
                    if a.polyA_found:
                        polya_assignments += 1

        info_pickler = pickle.Pickler(open(sample.out_raw_file + "_info", "wb"),  -1)
        info_pickler.dump(total_assignments)
        info_pickler.dump(polya_assignments)
        all_read_groups = set([r for sublist in read_groups for r in sublist])
        info_pickler.dump(all_read_groups)
        if total_assignments == 0:
            logger.warning("No reads were assigned to isoforms, check your input files")
        else:
            logger.info('Finishing read assignment, total assignments %d, polyA percentage %.1f' %
                        (total_assignments, 100 * polya_assignments / total_assignments))

    def process_assigned_reads(self, sample, dump_filename):
        if self.gffutils_db:
            chrom_clusters = self.get_chromosome_gene_clusters()
            chr_ids = [chr_id for chr_id, c in chrom_clusters]
        else:
            chr_ids = sorted(self.reference_record_dict.keys(), key=lambda x: len(self.reference_record_dict[x]), reverse=True)
            #chr_ids = ['chr1'] #

        pool = Pool(self.args.threads)
        logger.info("Processing assigned reads " + sample.label)

        results = pool.starmap(construct_models_in_parallel, [(SampleData(sample.file_list, sample.label + "_" + chr_id,
                                                                          sample.out_dir), chr_id, dump_filename,
                                                               self.args, self.multimapped_info_dict[chr_id], self.read_grouper,
                                                               (self.reference_record_dict[chr_id] if self.reference_record_dict else None),
                                                               self.io_support) for chr_id in chr_ids], chunksize=1)
        pool.close()
        pool.join()

        read_stat_counters, transcript_model_storages = [x[0] for x in results], [x[1] for x in results]

        self.create_aggregators(sample)
        for read_stat_counter in read_stat_counters:
            for k, v in read_stat_counter.stats_dict.items():
                self.read_stat_counter.stats_dict[k] += v

        gff_printer = GFFPrinter(sample.out_dir, sample.label, self.io_support,
                                 header=self.common_header)
        transcript_stat_counter = EnumStats()
        self.merge_assignments(sample, chr_ids)
        if not self.args.no_model_construction:
            for storage in transcript_model_storages:
                for t in storage:
                    transcript_stat_counter.add(t)
            self.merge_transcript_models(sample.label, chr_ids, gff_printer)

        self.finalize_aggregators(sample)
        logger.info("Transcript model file " + gff_printer.model_fname)
        transcript_stat_counter.print_start("Transcript model statistics")

    def load_read_info(self, dump_filename):
        gc.disable()
        self.multimapped_info_dict = defaultdict(lambda : defaultdict(list))
        if self.multimapped_reads:
            for read_id in self.multimapped_reads:
                for a in self.multimapped_reads[read_id]:
                    self.multimapped_info_dict[a.chr_id][a.read_id].append(a)
            self.multimapped_reads = None
        else:
            multimap_unpickler = pickle.Unpickler(open(dump_filename + "_multimappers", "rb"), fix_imports=False)
            while True:
                try:
                    obj = multimap_unpickler.load()
                    if isinstance(obj, list):
                        assignment_list = obj
                        for a in assignment_list:
                            self.multimapped_info_dict[a.chr_id][a.read_id].append(a)
                    else:
                        raise ValueError("Multimap assignment file {} is corrupted!".format(dump_filename))
                except EOFError:
                    break

        info_unpickler = pickle.Unpickler(open(dump_filename + "_info", "rb"), fix_imports=False)
        total_assignments = info_unpickler.load()
        polya_assignments = info_unpickler.load()
        all_read_groups = info_unpickler.load()

        gc.enable()
        return total_assignments, polya_assignments, all_read_groups

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
        self.basic_printer = BasicTSVAssignmentPrinter(sample.out_assigned_tsv, self.args, self.io_support,
                                                       additional_header=self.common_header)
        self.corrected_bed_printer = BEDPrinter(sample.out_corrected_bed, self.args, print_corrected=True)
        printer_list = [self.basic_printer, self.corrected_bed_printer]
        if self.args.sqanti_output:
            self.sqanti_printer = SqantiTSVPrinter(sample.out_alt_tsv, self.args, self.io_support)
            printer_list.append(self.sqanti_printer)
        self.global_printer = ReadAssignmentCompositePrinter(printer_list)

        self.gene_counter = create_gene_counter(sample.out_gene_counts_tsv, self.args.gene_quantification,
                                                ignore_read_groups=True, output_zeroes=False)
        self.transcript_counter = create_transcript_counter(sample.out_transcript_counts_tsv,
                                                            self.args.transcript_quantification,
                                                            ignore_read_groups=True, output_zeroes=False)
        self.transcript_model_counter = create_transcript_counter(sample.out_transcript_model_counts_tsv,
                                                                  self.args.transcript_quantification,
                                                                  ignore_read_groups=True, output_zeroes=False)

        self.transcript_model_global_counter = CompositeCounter([self.transcript_model_counter])
        if self.args.count_exons:
            self.exon_counter = ExonCounter(sample.out_exon_counts_tsv, ignore_read_groups=True)
            self.intron_counter = IntronCounter(sample.out_intron_counts_tsv, ignore_read_groups=True)
            self.global_counter = CompositeCounter([self.gene_counter, self.transcript_counter,
                                                    self.exon_counter, self.intron_counter])
        else:
            self.global_counter = CompositeCounter([self.gene_counter, self.transcript_counter])

        if self.args.read_group:
            self.gene_grouped_counter = create_gene_counter(sample.out_gene_grouped_counts_tsv,
                                                            self.args.gene_quantification, self.read_grouper)
            self.transcript_grouped_counter = create_transcript_counter(sample.out_transcript_grouped_counts_tsv,
                                                                        self.args.transcript_quantification,
                                                                        self.read_grouper)
            self.transcript_model_grouped_counter = create_transcript_counter(sample.out_transcript_model_grouped_counts_tsv,
                                                                              self.args.transcript_quantification,
                                                                              self.read_grouper, output_zeroes=False)

            self.transcript_model_global_counter.add_counters([self.transcript_model_grouped_counter])
            if self.args.count_exons:
                self.exon_grouped_counter = ExonCounter(sample.out_exon_grouped_counts_tsv)
                self.intron_grouped_counter = IntronCounter(sample.out_intron_grouped_counts_tsv)
                self.global_counter.add_counters([self.gene_grouped_counter, self.transcript_grouped_counter,
                                                  self.exon_grouped_counter, self.intron_grouped_counter])
            else:
                self.global_counter.add_counters([self.gene_grouped_counter, self.transcript_grouped_counter])

    def merge_assignments(self, sample, chr_ids):
        merge_files([rreplace(sample.out_assigned_tsv, sample.label, sample.label + "_" + chr_id) for chr_id in chr_ids],
                    sample.out_assigned_tsv, copy_header=False)
        merge_files([rreplace(sample.out_corrected_bed, sample.label, sample.label + "_" + chr_id) for chr_id in chr_ids],
                    sample.out_corrected_bed, copy_header=False)
        for p in self.global_counter.counters:
            merge_files([rreplace(p.output_counts_file_name, sample.label, sample.label + "_" + chr_id) for chr_id in chr_ids],
                        p.output_counts_file_name,
                        stats_file_names=[rreplace(p.output_stats_file_name, sample.label, sample.label + "_" + chr_id) for chr_id in chr_ids]
                        if p.output_stats_file_name else None,
                        ignore_read_groups=p.ignore_read_groups)
            p.convert_counts_to_tpm()

    def merge_transcript_models(self, label, chr_ids, gff_printer):
        merge_files([rreplace(gff_printer.model_fname, label, label + "_" + chr_id) for chr_id in chr_ids],
                    gff_printer.model_fname, copy_header=False)
        merge_files([rreplace(gff_printer.r2t_fname, label, label + "_" + chr_id) for chr_id in chr_ids],
                    gff_printer.r2t_fname, copy_header=False)
        for p in self.transcript_model_global_counter.counters:
            merge_files([rreplace(p.output_counts_file_name, label, label + "_" + chr_id) for chr_id in chr_ids],
                        p.output_counts_file_name,
                        stats_file_names=[rreplace(p.output_stats_file_name, label, label + "_" + chr_id) for chr_id in chr_ids]
                        if p.output_stats_file_name else None,
                        ignore_read_groups=p.ignore_read_groups)
            p.convert_counts_to_tpm()

    def finalize_aggregators(self, sample):
        logger.info("Gene counts are stored in " + self.gene_counter.output_counts_file_name)
        logger.info("Transcript counts are stored in " + self.transcript_counter.output_counts_file_name)
        logger.info("Read assignments are stored in " + self.basic_printer.output_file_name)
        self.read_stat_counter.print_start("Read assignment statistics")

