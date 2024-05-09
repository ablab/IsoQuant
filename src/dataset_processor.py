############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import glob
import gzip
import itertools
import logging
import os
import shutil
from enum import Enum, unique
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import gffutils
import pysam
from pyfaidx import Fasta, Faidx, UnsupportedCompressionFormat

from .common import proper_plural_form
from .serialization import *
from .isoform_assignment import BasicReadAssignment, ReadAssignmentType
from .stats import EnumStats
from .file_utils import merge_files, merge_counts
from .input_data_storage import SampleData
from .alignment_processor import AlignmentCollector, AlignmentType
from .long_read_counter import (
    ExonCounter,
    IntronCounter,
    CompositeCounter,
    create_gene_counter,
    create_transcript_counter,
)
from .multimap_resolver import MultimapResolver
from .read_groups import (
    create_read_grouper,
    prepare_read_groups
)
from .assignment_io import (
    IOSupport,
    BEDPrinter,
    ReadAssignmentCompositePrinter,
    SqantiTSVPrinter,
    BasicTSVAssignmentPrinter,
    TmpFileAssignmentPrinter,
    TmpFileAssignmentLoader,
)
from .transcript_printer import GFFPrinter, VoidTranscriptPrinter, create_extended_storage
from .graph_based_model_construction import GraphBasedModelConstructor
from .gene_info import TranscriptModelType, get_all_chromosome_genes, get_all_chromosome_transcripts

logger = logging.getLogger('IsoQuant')


def reads_collected_lock_file_name(sample_out_raw, chr_id):
    return "{}_{}_collected".format(sample_out_raw, chr_id)


def reads_processed_lock_file_name(dump_filename, chr_id):
    chr_dump_file = dump_filename + "_" + chr_id
    return "{}_processed".format(chr_dump_file)


def read_group_lock_filename(sample):
    return sample.read_group_file + "_lock"


def clean_locks(chr_ids, base_name, fname_function):
    for chr_id in chr_ids:
        fname = fname_function(base_name, chr_id)
        if os.path.exists(fname):
            os.remove(fname)


@unique
class PolyAUsageStrategies(Enum):
    auto = 1
    never = 2
    always = 3


def set_polya_requirement_strategy(flag, polya_requirement_strategy):
    if polya_requirement_strategy == PolyAUsageStrategies.auto:
        return flag
    elif polya_requirement_strategy == PolyAUsageStrategies.never:
        return False
    else:
        return True


def collect_reads_in_parallel(sample, chr_id, args):
    current_chr_record = Fasta(args.reference, indexname=args.fai_file_name)[chr_id]
    if args.high_memory:
        current_chr_record = str(current_chr_record)
    read_grouper = create_read_grouper(args, sample, chr_id)
    lock_file = reads_collected_lock_file_name(sample.out_raw_file, chr_id)
    save_file = "{}_{}".format(sample.out_raw_file, chr_id)
    group_file = "{}_{}_groups".format(sample.out_raw_file, chr_id)
    bamstat_file = "{}_{}_bamstat".format(sample.out_raw_file, chr_id)
    processed_reads = []

    if os.path.exists(lock_file) and args.resume:
        logger.info("Detected processed reads for " + chr_id)
        if os.path.exists(group_file) and os.path.exists(save_file):
            read_grouper.read_groups.clear()
            for g in open(group_file):
                read_grouper.read_groups.add(g.strip())
            alignment_stat_counter = EnumStats(bamstat_file)
            loader = ReadAssignmentLoader(save_file, None, None, None)
            while loader.has_next():
                gene_info, assignment_storage = loader.get_next()
                for a in assignment_storage:
                    processed_reads.append(BasicReadAssignment(a))
            logger.info("Loaded data for " + chr_id)
            return processed_reads, read_grouper.read_groups, alignment_stat_counter
        else:
            logger.warning("Something is wrong with save files for %s, will process from scratch " % chr_id)

    tmp_printer = TmpFileAssignmentPrinter(save_file, args)
    bam_files = list(map(lambda x: x[0], sample.file_list))
    bam_file_pairs = [(pysam.AlignmentFile(bam, "rb", require_index=True), bam) for bam in bam_files]
    gffutils_db = gffutils.FeatureDB(args.genedb) if args.genedb else None
    illumina_bam = sample.illumina_bam

    logger.info("Processing chromosome " + chr_id)
    alignment_collector = \
        AlignmentCollector(chr_id, bam_file_pairs, args, illumina_bam, gffutils_db, current_chr_record, read_grouper)

    for gene_info, assignment_storage in alignment_collector.process():
        tmp_printer.add_gene_info(gene_info)
        for read_assignment in assignment_storage:
            tmp_printer.add_read_info(read_assignment)
            processed_reads.append(BasicReadAssignment(read_assignment))
    with open(group_file, "w") as group_dump:
        for g in read_grouper.read_groups:
            group_dump.write("%s\n" % g)
    alignment_collector.alignment_stat_counter.dump(bamstat_file)

    logger.info("Finished processing chromosome " + chr_id)
    open(lock_file, "w").close()
    for bam in bam_file_pairs:
        bam[0].close()

    return processed_reads, read_grouper.read_groups, alignment_collector.alignment_stat_counter


class ReadAssignmentLoader:
    def __init__(self, save_file_name, gffutils_db, chr_record, multimapped_chr_dict):
        logger.info("Loading read assignments from " + save_file_name)
        assert os.path.exists(save_file_name)
        self.save_file_name = save_file_name
        self.unpickler = TmpFileAssignmentLoader(save_file_name, gffutils_db, chr_record)
        self.multimapped_chr_dict = multimapped_chr_dict

    def has_next(self):
        return self.unpickler.has_next()

    def get_next(self):
        if not self.unpickler.has_next():
            return None, None

        assert self.unpickler.is_gene_info()
        gene_info = self.unpickler.get_object()
        assignment_storage = []
        while self.unpickler.is_read_assignment():
            read_assignment = self.unpickler.get_object()
            if self.multimapped_chr_dict is not None and read_assignment.read_id in self.multimapped_chr_dict:
                resolved_assignment = None
                for a in self.multimapped_chr_dict[read_assignment.read_id]:
                    if a.assignment_id == read_assignment.assignment_id and a.chr_id == read_assignment.chr_id:
                        if resolved_assignment is not None:
                            logger.info("Duplicate read: %s %s %s" % (read_assignment.read_id, a.gene_id, a.chr_id))
                        resolved_assignment = a

                if not resolved_assignment:
                    logger.warning("Incomplete information on read %s" % read_assignment.read_id)
                    continue
                elif resolved_assignment.assignment_type == ReadAssignmentType.suspended:
                    continue
                else:
                    read_assignment.assignment_type = resolved_assignment.assignment_type
                    read_assignment.gene_assignment_type = resolved_assignment.gene_assignment_type
                    read_assignment.multimapper = resolved_assignment.multimapper
            assignment_storage.append(read_assignment)

        return gene_info, assignment_storage


def construct_models_in_parallel(sample, chr_id, dump_filename, args, read_groups):
    logger.info("Processing chromosome " + chr_id)
    construct_models = not args.no_model_construction
    current_chr_record = Fasta(args.reference, indexname=args.fai_file_name)[chr_id]
    multimapped_reads = defaultdict(list)
    multimap_loader = open(dump_filename + "_multimappers_" + chr_id, "rb")
    list_size = read_int(multimap_loader)
    while list_size != TERMINATION_INT:
        for i in range(list_size):
            a = BasicReadAssignment.deserialize(multimap_loader)
            if a.chr_id == chr_id:
                multimapped_reads[a.read_id].append(a)
        list_size = read_int(multimap_loader)

    chr_dump_file = dump_filename + "_" + chr_id
    lock_file = reads_processed_lock_file_name(dump_filename, chr_id)
    read_stat_file = "{}_read_stat".format(chr_dump_file)
    transcript_stat_file = "{}_transcript_stat".format(chr_dump_file)

    if os.path.exists(lock_file) and args.resume:
        logger.info("Processed assignments from chromosome " + chr_id + " detected")
        read_stat = EnumStats(read_stat_file)
        transcript_stat = EnumStats(transcript_stat_file) if construct_models else EnumStats()
        return read_stat, transcript_stat


    if args.genedb:
        gffutils_db = gffutils.FeatureDB(args.genedb)
    else:
        gffutils_db = None
    aggregator = ReadAssignmentAggregator(args, sample, read_groups, gffutils_db, chr_id)

    transcript_stat_counter = EnumStats()
    io_support = IOSupport(args)
    tmp_gff_printer = GFFPrinter(sample.out_dir, sample.prefix, check_canonical=args.check_canonical) \
        if construct_models else VoidTranscriptPrinter()
    tmp_extended_gff_printer = GFFPrinter(sample.out_dir, sample.prefix, gtf_suffix=".extended_annotation.gtf",
                                          output_r2t=False, check_canonical=args.check_canonical) \
        if (construct_models and args.genedb) else VoidTranscriptPrinter()

    sqanti_t2t_printer = SqantiTSVPrinter(sample.out_t2t_tsv, args, IOSupport(args)) \
        if args.sqanti_output else VoidTranscriptPrinter()
    novel_model_storage = []

    loader = ReadAssignmentLoader(chr_dump_file, gffutils_db, current_chr_record, multimapped_reads)
    while loader.has_next():
        gene_info, assignment_storage = loader.get_next()
        logger.debug("Processing %d reads" % len(assignment_storage))
        for read_assignment in assignment_storage:
            if read_assignment is None:
                continue
            aggregator.read_stat_counter.add(read_assignment.assignment_type)
            aggregator.global_printer.add_read_info(read_assignment)
            aggregator.global_counter.add_read_info(read_assignment)

        if construct_models:
            model_constructor = GraphBasedModelConstructor(gene_info, current_chr_record, args,
                                                           aggregator.transcript_model_global_counter)
            model_constructor.process(assignment_storage)
            if args.check_canonical:
                io_support.add_canonical_info(model_constructor.transcript_model_storage, gene_info)
            tmp_gff_printer.dump(model_constructor.gene_info, model_constructor.transcript_model_storage)
            tmp_gff_printer.dump_read_assignments(model_constructor)
            for m in model_constructor.transcript_model_storage:
                if m.transcript_type != TranscriptModelType.known:
                    novel_model_storage.append(m)
            for a in model_constructor.transcript2transcript:
                sqanti_t2t_printer.add_read_info(a)
            for t in model_constructor.transcript_model_storage:
                transcript_stat_counter.add(t.transcript_type)

    aggregator.global_counter.dump()
    aggregator.read_stat_counter.dump(read_stat_file)
    if construct_models:
        if gffutils_db:
            all_models, gene_info = create_extended_storage(gffutils_db, chr_id, current_chr_record, novel_model_storage)
            if args.check_canonical:
                io_support.add_canonical_info(all_models, gene_info)
            tmp_extended_gff_printer.dump(gene_info, all_models)
        aggregator.transcript_model_global_counter.dump()
        transcript_stat_counter.dump(transcript_stat_file)
    logger.info("Finished processing chromosome " + chr_id)
    open(lock_file, "w").close()

    return aggregator.read_stat_counter, transcript_stat_counter


class ReadAssignmentAggregator:
    def __init__(self, args, sample, read_groups, gffutils_db=None, chr_id=None, gzipped=False):
        self.args = args
        self.read_groups = read_groups
        self.common_header = "# Command line: " + args._cmd_line + "\n# IsoQuant version: " + args._version + "\n"
        self.io_support = IOSupport(self.args)

        self.gene_set = set()
        self.transcript_set = set()
        if gffutils_db and chr_id:
            self.gene_set = set(get_all_chromosome_genes(gffutils_db, chr_id))
            self.transcript_set = set(get_all_chromosome_transcripts(gffutils_db, chr_id))

        self.read_stat_counter = EnumStats()
        self.corrected_bed_printer = BEDPrinter(sample.out_corrected_bed,
                                                self.args,
                                                print_corrected=True,
                                                gzipped=gzipped)
        printer_list = [self.corrected_bed_printer]
        if self.args.genedb:
            self.basic_printer = BasicTSVAssignmentPrinter(sample.out_assigned_tsv, self.args, self.io_support,
                                                           additional_header=self.common_header, gzipped=gzipped)
            printer_list.append(self.basic_printer)
        if self.args.sqanti_output:
            self.t2t_sqanti_printer = SqantiTSVPrinter(sample.out_t2t_tsv, self.args, self.io_support)
        self.global_printer = ReadAssignmentCompositePrinter(printer_list)

        self.global_counter = CompositeCounter([])
        if self.args.genedb:
            self.gene_counter = create_gene_counter(sample.out_gene_counts_tsv,
                                                    self.args.gene_quantification,
                                                    complete_feature_list=self.gene_set,
                                                    ignore_read_groups=True, output_zeroes=True)
            self.transcript_counter = create_transcript_counter(sample.out_transcript_counts_tsv,
                                                                self.args.transcript_quantification,
                                                                complete_feature_list=self.transcript_set,
                                                                ignore_read_groups=True, output_zeroes=True)
            self.global_counter.add_counters([self.gene_counter, self.transcript_counter])

        self.transcript_model_global_counter = CompositeCounter([])
        if not self.args.no_model_construction:
            self.transcript_model_counter = create_transcript_counter(sample.out_transcript_model_counts_tsv,
                                                                      self.args.transcript_quantification,
                                                                      ignore_read_groups=True, output_zeroes=False)
            self.transcript_model_global_counter.add_counters([self.transcript_model_counter])

        if self.args.count_exons and self.args.genedb:
            self.exon_counter = ExonCounter(sample.out_exon_counts_tsv, ignore_read_groups=True)
            self.intron_counter = IntronCounter(sample.out_intron_counts_tsv, ignore_read_groups=True)
            self.global_counter.add_counters([self.exon_counter, self.intron_counter])

        if self.args.read_group and self.args.genedb:
            self.gene_grouped_counter = create_gene_counter(sample.out_gene_grouped_counts_tsv,
                                                            self.args.gene_quantification,
                                                            complete_feature_list=self.gene_set,
                                                            read_groups=self.read_groups)
            self.transcript_grouped_counter = create_transcript_counter(sample.out_transcript_grouped_counts_tsv,
                                                                        self.args.transcript_quantification,
                                                                        complete_feature_list=self.transcript_set,
                                                                        read_groups=self.read_groups)
            self.global_counter.add_counters([self.gene_grouped_counter, self.transcript_grouped_counter])

            if self.args.count_exons:
                self.exon_grouped_counter = ExonCounter(sample.out_exon_grouped_counts_tsv)
                self.intron_grouped_counter = IntronCounter(sample.out_intron_grouped_counts_tsv)
                self.global_counter.add_counters([self.exon_grouped_counter, self.intron_grouped_counter])

        if self.args.read_group and not self.args.no_model_construction:
            self.transcript_model_grouped_counter = create_transcript_counter(
                sample.out_transcript_model_grouped_counts_tsv,
                self.args.transcript_quantification,
                read_groups=self.read_groups, output_zeroes=False)
            self.transcript_model_global_counter.add_counters([self.transcript_model_grouped_counter])

    def finalize_aggregators(self, sample):
        if self.args.genedb:
            logger.info("Gene counts are stored in " + self.gene_counter.output_counts_file_name)
            logger.info("Transcript counts are stored in " + self.transcript_counter.output_counts_file_name)
            logger.info("Read assignments are stored in " + self.basic_printer.output_file_name +
                        (".gz" if self.args.gzipped else ""))
        self.read_stat_counter.print_start("Read assignment statistics")


# Class for processing all samples against gene database
class DatasetProcessor:
    def __init__(self, args):
        self.args = args
        self.args.gunzipped_reference = None
        self.common_header = "# Command line: " + args._cmd_line + "\n# IsoQuant version: " + args._version + "\n"
        self.io_support = IOSupport(self.args)
        self.all_read_groups = set()
        self.alignment_stat_counter = EnumStats()

        if args.genedb:
            logger.info("Loading gene database from " + self.args.genedb)
            self.gffutils_db = gffutils.FeatureDB(self.args.genedb)
        else:
            self.gffutils_db = None

        if self.args.needs_reference:
            logger.info("Loading reference genome from %s" % self.args.reference)
            ref_dir = os.path.dirname(self.args.reference)
            ref_file_name = os.path.basename(self.args.reference)
            ref_name, outer_ext = os.path.splitext(ref_file_name)

            # make symlink for pyfaidx index
            args.fai_file_name = self.args.reference + ".fai"
            if not os.path.exists(args.fai_file_name) and not os.access(ref_dir, os.W_OK):
                # index does not exist near the reference and reference folder is not writable
                # store index in the output folder in this case
                args.fai_file_name = os.path.join(args.output, ref_file_name  + ".fai")

            low_ext = outer_ext.lower()
            if low_ext in ['.gz', '.gzip', '.bgz']:
                try:
                    self.reference_record_dict = Fasta(self.args.reference, indexname=args.fai_file_name)
                except UnsupportedCompressionFormat:
                    gunzipped_reference = os.path.join(args.output, ref_name)
                    if not os.path.exists(gunzipped_reference) or not self.args.resume:
                        with open(gunzipped_reference, "w") as outf:
                            shutil.copyfileobj(gzip.open(self.args.reference, "rt"), outf)
                        logger.info("Loading uncompressed reference from " + gunzipped_reference)
                    self.args.reference = gunzipped_reference
                    self.reference_record_dict = Fasta(self.args.reference, indexname=args.fai_file_name)
            else:
                self.reference_record_dict = Fasta(self.args.reference, indexname=args.fai_file_name)
        else:
            self.reference_record_dict = None

    def __del__(self):
        if not self.args.keep_tmp and self.args.gunzipped_reference:
            if os.path.exists(self.args.gunzipped_reference):
                os.remove(self.args.gunzipped_reference)

    def process_all_samples(self, input_data):
        logger.info("Processing " + proper_plural_form("experiment", len(input_data.samples)))
        for sample in input_data.samples:
            self.process_sample(sample)
        logger.info("Processed " + proper_plural_form("experiment", len(input_data.samples)))

    # Run through all genes in db and count stats according to alignments given in bamfile_name
    def process_sample(self, sample):
        logger.info("Processing experiment " + sample.prefix)
        logger.info("Experiment has " + proper_plural_form("BAM file", len(sample.file_list)) + ": " + ", ".join(
            map(lambda x: x[0], sample.file_list)))
        self.args.use_technical_replicas = self.args.read_group == "file_name" and len(sample.file_list) > 1

        self.all_read_groups = set()
        if self.args.resume and os.path.exists(sample.read_group_file + "_lock"):
            logger.info("Read group table was split during the previous run, existing files will be used")
        else:
            fname = read_group_lock_filename(sample)
            if os.path.exists(fname):
                os.remove(fname)
            prepare_read_groups(self.args, sample)
            open(fname, "w").close()

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

        total_alignments, polya_found, self.all_read_groups = self.load_read_info(saves_file)

        polya_fraction = polya_found / total_alignments if total_alignments > 0 else 0.0
        logger.info("Total alignments used for analysis: %d, polyA tail detected in %d (%.1f%%)" %
                    (total_alignments, polya_found, polya_fraction * 100.0))
        if (polya_fraction < self.args.low_polya_percentage_threshold and
                self.args.polya_requirement_strategy != PolyAUsageStrategies.never):
            logger.warning("PolyA percentage is suspiciously low. IsoQuant expects non-polya-trimmed reads. "
                           "If you aim to construct transcript models, consider using --polya_requirement option.")

        self.args.requires_polya_for_construction = set_polya_requirement_strategy(
            polya_fraction >= self.args.polya_percentage_threshold,
            self.args.polya_requirement_strategy)
        self.args.require_monointronic_polya = set_polya_requirement_strategy(
            # do not require polyA tails for mono-intronic only if the data is reliable and polyA percentage is low
            self.args.require_monointronic_polya or self.args.requires_polya_for_construction,
            self.args.polya_requirement_strategy)
        self.args.require_monoexonic_polya = set_polya_requirement_strategy(
            # do not require polyA tails for mono-intronic only if the data is reliable and polyA percentage is low
            self.args.require_monoexonic_polya or self.args.requires_polya_for_construction,
            self.args.polya_requirement_strategy)

        self.process_assigned_reads(sample, saves_file)
        if not self.args.read_assignments and not self.args.keep_tmp:
            for f in glob.glob(saves_file + "_*"):
                os.remove(f)
            for f in glob.glob(sample.read_group_file + "*"):
                os.remove(f)
        logger.info("Processed experiment " + sample.prefix)

    def collect_reads(self, sample):
        logger.info('Collecting read alignments')
        info_file = sample.out_raw_file + "_info"
        lock_file = sample.out_raw_file + "_lock"

        if os.path.exists(lock_file):
            if self.args.resume:
                logger.info("Collected reads detected, will not process")
                return
            else:
                os.remove(lock_file)

        chr_ids = sorted(
            self.reference_record_dict.keys(),
            key=lambda x: len(self.reference_record_dict[x]),
            reverse=True,
        )
        if not self.args.resume:
            clean_locks(chr_ids, sample.out_raw_file, reads_collected_lock_file_name)
            clean_locks(chr_ids, sample.out_raw_file, reads_processed_lock_file_name)

        all_read_groups = set()
        multimapped_reads = defaultdict(list)

        read_gen = (
            collect_reads_in_parallel,
            itertools.repeat(sample),
            chr_ids,
            itertools.repeat(self.args),
        )

        if self.args.threads > 1:
            with ProcessPoolExecutor(max_workers=self.args.threads) as proc:
                results = proc.map(*read_gen, chunksize=1)

                for storage, read_groups, alignment_stats in results:
                    for basic_read_assignment in storage:
                        multimapped_reads[basic_read_assignment.read_id].append(basic_read_assignment)
                    all_read_groups.update(read_groups)
                    self.alignment_stat_counter.merge(alignment_stats)
        else:
            for storage, read_groups, alignment_stats in map(*read_gen):
                for basic_read_assignment in storage:
                    multimapped_reads[basic_read_assignment.read_id].append(basic_read_assignment)
                all_read_groups.update(read_groups)
                self.alignment_stat_counter.merge(alignment_stats)

        for bam_file in list(map(lambda x: x[0], sample.file_list)):
            bam = pysam.AlignmentFile(bam_file, "rb", require_index=True)
            self.alignment_stat_counter.add(AlignmentType.unaligned, bam.unmapped)
        self.alignment_stat_counter.print_start("Alignments collected, overall alignment statistics:")

        logger.info("Resolving multimappers")
        multimap_resolver = MultimapResolver(self.args.multimap_strategy)
        multimap_dumper = {}
        for chr_id in chr_ids:
            multimap_dumper[chr_id] = open(sample.out_raw_file + "_multimappers_" + chr_id, "wb")
        total_assignments = 0
        polya_assignments = 0

        for assignment_list in multimapped_reads.values():
            if len(assignment_list) > 1:
                multimap_resolver.resolve(assignment_list)
                resolved_lists = defaultdict(list)
                for a in assignment_list:
                    resolved_lists[a.chr_id].append(a)
                for chr_id in resolved_lists.keys():
                    write_list(resolved_lists[chr_id], multimap_dumper[chr_id], BasicReadAssignment.serialize)

            for a in assignment_list:
                if a.assignment_type != ReadAssignmentType.suspended:
                    total_assignments += 1
                    if a.polyA_found:
                        polya_assignments += 1

        for chr_id in chr_ids:
            write_int(TERMINATION_INT, multimap_dumper[chr_id])
            multimap_dumper[chr_id].close()

        info_dumper = open(info_file, "wb")
        write_int(total_assignments, info_dumper)
        write_int(polya_assignments, info_dumper)
        write_list(list(all_read_groups), info_dumper, write_string)
        info_dumper.close()
        open(lock_file, "w").close()

        if total_assignments == 0:
            logger.warning("No reads were assigned to isoforms, check your input files")
        else:
            logger.info('Finishing read assignment, total assignments %d, polyA percentage %.1f' %
                        (total_assignments, 100 * polya_assignments / total_assignments))

    def process_assigned_reads(self, sample, dump_filename):
        chr_ids = sorted(
            self.reference_record_dict.keys(),
            key=lambda x: len(self.reference_record_dict[x]),
            reverse=True
        )
        logger.info("Processing assigned reads " + sample.prefix)
        logger.info("Transcript models construction is turned %s" %
                    ("off" if self.args.no_model_construction else "on"))

        # set up aggregators and outputs
        aggregator = ReadAssignmentAggregator(self.args, sample, self.all_read_groups, gzipped=self.args.gzipped)
        transcript_stat_counter = EnumStats()

        gff_printer = VoidTranscriptPrinter()
        extended_gff_printer = VoidTranscriptPrinter()
        if not self.args.no_model_construction:
            logger.info("Transcript construction options:")
            logger.info("  Novel monoexonic transcripts will be reported: %s"
                        % ("yes" if self.args.report_novel_unspliced else "no"))
            logger.info("  PolyA tails are required for multi-exon transcripts to be reported: %s"
                        % ("yes" if self.args.requires_polya_for_construction else "no"))
            logger.info("  PolyA tails are required for 2-exon transcripts to be reported: %s"
                        % ("yes" if self.args.require_monointronic_polya else "no"))
            logger.info("  PolyA tails are required for known monoexon transcripts to be reported: %s"
                        % ("yes" if self.args.require_monoexonic_polya else "no"))
            logger.info("  PolyA tails are required for novel monoexon transcripts to be reported: %s" % "yes")
            logger.info("  Splice site reporting level: %s" % self.args.report_canonical_strategy.name)

            gff_printer = GFFPrinter(
                sample.out_dir, sample.prefix, header=self.common_header, gzipped=self.args.gzipped
            )
            if self.args.genedb:
                extended_gff_printer = GFFPrinter(
                    sample.out_dir, sample.prefix,
                    gtf_suffix=".extended_annotation.gtf", output_r2t=False,
                    header=self.common_header
                )

        model_gen = (
            construct_models_in_parallel,
            (SampleData(sample.file_list, f"{sample.prefix}_{chr_id}", sample.out_dir, sample.readable_names_dict, sample.illumina_bam) for chr_id in chr_ids),
            chr_ids,
            itertools.repeat(dump_filename),
            itertools.repeat(self.args),
            itertools.repeat(self.all_read_groups),
        )

        if self.args.threads > 1:
            with ProcessPoolExecutor(max_workers=self.args.threads) as proc:
                results = proc.map(*model_gen, chunksize=1)

                for read_stat_counter, tsc in results:
                    for k, v in read_stat_counter.stats_dict.items():
                        aggregator.read_stat_counter.stats_dict[k] += v

                    if not self.args.no_model_construction:
                        for k, v in tsc.stats_dict.items():
                            transcript_stat_counter.stats_dict[k] += v
        else:
            for read_stat_counter, tsc in map(*model_gen):
                for k, v in read_stat_counter.stats_dict.items():
                    aggregator.read_stat_counter.stats_dict[k] += v

                if not self.args.no_model_construction:
                    for k, v in tsc.stats_dict.items():
                        transcript_stat_counter.stats_dict[k] += v

        if not self.args.no_model_construction:
            self.merge_transcript_models(sample.prefix, aggregator, chr_ids, gff_printer)
            logger.info("Transcript model file " + gff_printer.model_fname)
            if self.args.genedb:
                merge_files(extended_gff_printer.model_fname, sample.prefix, chr_ids,
                            extended_gff_printer.out_gff, copy_header=False)
                logger.info("Extended annotation is saved to " + extended_gff_printer.model_fname)
            transcript_stat_counter.print_start("Transcript model statistics")

        self.merge_assignments(sample, aggregator, chr_ids)
        if self.args.sqanti_output:
            merge_files(sample.out_t2t_tsv, sample.prefix, chr_ids, open(sample.out_t2t_tsv, "w"), copy_header=False)

        aggregator.finalize_aggregators(sample)

    def load_read_info(self, dump_filename):
        info_loader = open(dump_filename + "_info", "rb")
        total_assignments = read_int(info_loader)
        polya_assignments = read_int(info_loader)
        all_read_groups = set(read_list(info_loader, read_string))
        info_loader.close()
        return total_assignments, polya_assignments, all_read_groups

    def merge_assignments(self, sample, aggregator, chr_ids):
        if self.args.genedb:
            merge_files(sample.out_assigned_tsv, sample.prefix, chr_ids,
                        aggregator.basic_printer.output_file, copy_header=False)
        merge_files(sample.out_corrected_bed, sample.prefix, chr_ids,
                    aggregator.corrected_bed_printer.output_file, copy_header=False)

        for counter in aggregator.global_counter.counters:
            unaligned = self.alignment_stat_counter.stats_dict[AlignmentType.unaligned]
            merge_counts(counter, sample.prefix, chr_ids, unaligned)
            counter.convert_counts_to_tpm(self.args.normalization_method)

    def merge_transcript_models(self, label, aggregator, chr_ids, gff_printer):
        merge_files(gff_printer.model_fname, label, chr_ids, gff_printer.out_gff, copy_header=False)
        merge_files(gff_printer.r2t_fname, label, chr_ids, gff_printer.out_r2t, copy_header=False)
        for counter in aggregator.transcript_model_global_counter.counters:
            unaligned = self.alignment_stat_counter.stats_dict[AlignmentType.unaligned]
            merge_counts(counter, label, chr_ids, unaligned)
            counter.convert_counts_to_tpm(self.args.normalization_method)
