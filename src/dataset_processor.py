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
import time
from multiprocessing import Pool
import itertools

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


def collect_reads_in_parallel(sample, chr_id, args, read_grouper, current_chr_record):
    lock_file = "{}_{}_collected".format(sample.out_raw_file, chr_id)
    save_file = "{}_{}".format(sample.out_raw_file, chr_id)
    group_file = "{}_{}_groups".format(sample.out_raw_file, chr_id)
    processed_reads = []

    if os.path.exists(lock_file) and args.resume:
        logger.info("Detected processed reads for " + chr_id)
        if os.path.exists(group_file) and os.path.exists(save_file):
            read_grouper.read_groups.clear()
            for g in open(group_file):
                read_grouper.read_groups.add(g.strip())
            for gene_info, assignment_storage in load_assigned_reads(save_file, None, None):
                for a in assignment_storage:
                    processed_reads.append(BasicReadAssignment(a, gene_info))
            logger.info("Loaded data for " + chr_id)
            return processed_reads, read_grouper.read_groups
        else:
            logger.warning("Something is wrong with save files for %s, will process from scratch " % chr_id)

    tmp_printer = TmpFileAssignmentPrinter(save_file, args)
    bam_files = list(map(lambda x: x[0], sample.file_list))
    bam_file_pairs = [(pysam.AlignmentFile(bam, "rb"), bam) for bam in bam_files]
    gffutils_db = gffutils.FeatureDB(args.genedb, keep_order=True) if args.genedb else None

    logger.info("Processing chromosome " + chr_id)
    alignment_collector = IntergenicAlignmentCollector(chr_id, bam_file_pairs, args, gffutils_db, current_chr_record, read_grouper)

    if args.low_memory:
        alignment_iterator = alignment_collector.process_slow()
    else:
        alignment_iterator = alignment_collector.process()
    for gene_info, assignment_storage in alignment_iterator:
        tmp_printer.add_gene_info(gene_info)
        for read_assignment in assignment_storage:
            tmp_printer.add_read_info(read_assignment)
            processed_reads.append(BasicReadAssignment(read_assignment, gene_info))
    with open(group_file, "w") as group_dump:
        for g in read_grouper.read_groups:
            group_dump.write("%s\n" % g)

    logger.info("Finished processing chromosome " + chr_id)
    open(lock_file, "w").close()

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
                if multimapped_chr_dict is not None and read_assignment.read_id in multimapped_chr_dict:
                    resolved_assignment = None
                    for a in multimapped_chr_dict[read_assignment.read_id]:
                        if a.start == read_assignment.start() and a.end == read_assignment.end() and \
                                a.chr_id == read_assignment.chr_id and \
                                (a.assignment_type in [ReadAssignmentType.intergenic, ReadAssignmentType.noninformative] or
                                 ((not current_gene_info.gene_db_list or a.gene_id == current_gene_info.gene_db_list[0].id) and
                                  a.matches == (0 if not read_assignment.isoform_matches else len(read_assignment.isoform_matches)))):
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
                if current_gene_info:
                    yield current_gene_info, read_storage
                read_storage = []
                current_gene_info = obj
                current_gene_info.db = gffutils_db
            else:
                raise ValueError("Read assignment file {} is corrupted!".format(save_file_name))
        except EOFError:
            break
    gc.enable()
    if current_gene_info:
        yield current_gene_info, read_storage


def construct_models_in_parallel(sample, chr_id, dump_filename, args, multimapped_reads, read_grouper, current_chr_record,
                                 io_support):
    chr_dump_file = dump_filename + "_" + chr_id
    lock_file = "{}_processed".format(chr_dump_file)
    read_stat_file = "{}_read_stat".format(chr_dump_file)
    transcript_stat_file = "{}_transcript_stat".format(chr_dump_file)

    if os.path.exists(lock_file) and args.resume:
        logger.info("Processed assignments from chromosome " + chr_id + " detected")
        return EnumStats(read_stat_file), EnumStats(transcript_stat_file)

    processor = DatasetProcessor(args, parallel=True)
    processor.read_grouper = read_grouper
    processor.create_aggregators(sample)
    if args.genedb:
        gffutils_db = gffutils.FeatureDB(args.genedb, keep_order=True)
    else:
        gffutils_db = None

    logger.info("Processing chromosome " + chr_id)
    transcript_stat_counter = EnumStats()
    tmp_gff_printer = GFFPrinter(sample.out_dir, sample.label, io_support)
    tmp_extended_gff_printer = None
    if gffutils_db:
        tmp_extended_gff_printer = GFFPrinter(sample.out_dir, sample.label, io_support,
                                             gtf_suffix=".extended_annotation.gtf", output_r2t=False)

    for gene_info, assignment_storage in load_assigned_reads(chr_dump_file, gffutils_db, multimapped_reads):
        logger.debug("Processing %d reads" % len(assignment_storage))
        for read_assignment in assignment_storage:
            if read_assignment is None:
                continue
            processor.read_stat_counter.add(read_assignment.assignment_type)
            processor.global_printer.add_read_info(read_assignment)
            processor.global_counter.add_read_info(read_assignment)

        if not args.no_model_construction:
            model_constructor = GraphBasedModelConstructor(gene_info, current_chr_record, args, processor.transcript_model_global_counter)
            model_constructor.process(assignment_storage)
            tmp_gff_printer.dump(model_constructor)
            if tmp_extended_gff_printer:
                tmp_extended_gff_printer.dump(model_constructor, model_constructor.extended_annotation_storage)
            for t in model_constructor.transcript_model_storage:
                transcript_stat_counter.add(t.transcript_type)

    processor.global_counter.dump()
    processor.transcript_model_global_counter.dump()
    processor.read_stat_counter.dump(read_stat_file)
    transcript_stat_counter.dump(transcript_stat_file)
    logger.info("Finished processing chromosome " + chr_id)
    open(lock_file, "w").close()
    # TODO: switch to returning transcript_stat_counter
    return processor.read_stat_counter, transcript_stat_counter


# Class for processing all samples against gene database
class DatasetProcessor:
    def __init__(self, args, parallel=False):
        self.args = args
        self.args.gunzipped_reference = None
        self.common_header = "# Command line: " + args._cmd_line + "\n# IsoQuant version: " + args._version + "\n"
#        self.read_grouper = create_read_grouper(args)
        self.io_support = IOSupport(self.args)
        if parallel: return
        self.read_grouper = create_read_grouper(args)
        if args.genedb:
            logger.info("Loading gene database from " + self.args.genedb)
            self.gffutils_db = gffutils.FeatureDB(self.args.genedb, keep_order=True)
        else:
            self.gffutils_db = None

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

    # Run through all genes in db and count stats according to alignments given in bamfile_name
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
        info_file = sample.out_raw_file + "_info"
        lock_file = sample.out_raw_file + "_lock"
        # TODO: load save only if command line was the same or continue option is enabled
        if os.path.exists(lock_file) and self.args.resume:
            logger.info("Collected reads detected, will not process")
            return

        chr_ids = sorted(self.reference_record_dict.keys(), key=lambda x: len(self.reference_record_dict[x]),
                         reverse=True)
        if self.args.threads == 1:
            results = itertools.starmap(collect_reads_in_parallel, [(sample, chr_id, self.args, self.read_grouper,
                                                                     (self.reference_record_dict[
                                                                          chr_id] if self.reference_record_dict else None))
                                                                    for chr_id in chr_ids])
        else:
            pool = Pool(self.args.threads)
            results = pool.starmap(collect_reads_in_parallel, [(sample, chr_id, self.args, self.read_grouper,
                                                                (self.reference_record_dict[chr_id] if self.reference_record_dict else None))
                                                               for chr_id in chr_ids], chunksize=1)
            pool.close()
            pool.join()
        processed_reads, read_groups = zip(*results)

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

        info_pickler = pickle.Pickler(open(info_file, "wb"),  -1)
        info_pickler.dump(total_assignments)
        info_pickler.dump(polya_assignments)
        all_read_groups = set([r for sublist in read_groups for r in sublist])
        info_pickler.dump(all_read_groups)
        open(lock_file, "w").close()
        if total_assignments == 0:
            logger.warning("No reads were assigned to isoforms, check your input files")
        else:
            logger.info('Finishing read assignment, total assignments %d, polyA percentage %.1f' %
                        (total_assignments, 100 * polya_assignments / total_assignments))

    def process_assigned_reads(self, sample, dump_filename):
        chr_ids = sorted(self.reference_record_dict.keys(), key=lambda x: len(self.reference_record_dict[x]), reverse=True)
        logger.info("Processing assigned reads " + sample.label)
        if self.args.threads == 1:
            results = itertools.starmap(construct_models_in_parallel, [(SampleData(sample.file_list, sample.label + "_" + chr_id,
                                                                                   sample.out_dir), chr_id, dump_filename,
                                                                        self.args, self.multimapped_info_dict[chr_id], self.read_grouper,
                                                                        (self.reference_record_dict[chr_id] if self.reference_record_dict else None),
                                                                        self.io_support) for chr_id in chr_ids])
        else:
            pool = Pool(self.args.threads)
            results = pool.starmap(construct_models_in_parallel, [(SampleData(sample.file_list, sample.label + "_" + chr_id,
                                                                              sample.out_dir), chr_id, dump_filename,
                                                                   self.args, self.multimapped_info_dict[chr_id], self.read_grouper,
                                                                   (self.reference_record_dict[chr_id] if self.reference_record_dict else None),
                                                                   self.io_support) for chr_id in chr_ids], chunksize=1)
            pool.close()
            pool.join()

        read_stat_counters, transcript_stat_counters = zip(*results)
        self.create_aggregators(sample)
        for read_stat_counter in read_stat_counters:
            for k, v in read_stat_counter.stats_dict.items():
                self.read_stat_counter.stats_dict[k] += v

        gff_printer = GFFPrinter(sample.out_dir, sample.label, self.io_support,
                                 header=self.common_header)
        extended_gff_printer = None
        if self.args.genedb:
            extended_gff_printer = GFFPrinter(sample.out_dir, sample.label, self.io_support,
                                              gtf_suffix=".extended_annotation.gtf", output_r2t=False,
                                              header=self.common_header)

        transcript_stat_counter = EnumStats()
        self.merge_assignments(sample, chr_ids)
        if self.args.sqanti_output:
            merge_files(
                [rreplace(sample.out_alt_tsv, sample.label, sample.label + "_" + chr_id) for chr_id in chr_ids],
                sample.out_alt_tsv, copy_header=False)

        if not self.args.no_model_construction:
            for tsc in transcript_stat_counters:
                for k, v in tsc.stats_dict.items():
                    transcript_stat_counter.stats_dict[k] += v

            self.merge_transcript_models(sample.label, chr_ids, gff_printer)
            logger.info("Transcript model file " + gff_printer.model_fname)
            if extended_gff_printer:
                merge_files(
                    [rreplace(extended_gff_printer.model_fname, sample.label, sample.label + "_" + chr_id) for chr_id in
                     chr_ids],
                    extended_gff_printer.model_fname, copy_header=False)
                logger.info("Extended annotation is saved to " + extended_gff_printer.model_fname)
            transcript_stat_counter.print_start("Transcript model statistics")

        self.finalize_aggregators(sample)

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

