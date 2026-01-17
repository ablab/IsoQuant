############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import gc
import glob
import gzip
import itertools
import logging
import multiprocessing
import shutil
import sys
from enum import Enum, unique
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import gffutils
import pysam
from pyfaidx import Fasta

from .modes import IsoQuantMode
from .common import proper_plural_form, large_output_enabled
from .error_codes import IsoQuantExitCode
from .serialization import *
from .stats import EnumStats
from .file_utils import merge_files, merge_counts
from .input_data_storage import SampleData
from .alignment_processor import AlignmentType
from .read_groups import prepare_read_groups, get_grouping_strategy_names
from .assignment_io import IOSupport
from .processed_read_manager import ProcessedReadsManagerHighMemory, ProcessedReadsManagerNoSecondary, ProcessedReadsManagerNormalMemory
from .id_policy import SimpleIDDistributor, FeatureIdStorage
from .file_naming import *
from .transcript_printer import GFFPrinter, VoidTranscriptPrinter
from .barcode_calling.umi_filtering import create_transcript_info_dict
from .table_splitter import split_read_table_parallel
from .assignment_aggregator import ReadAssignmentAggregator
from .parallel_workers import (
    collect_reads_in_parallel,
    construct_models_in_parallel,
    filter_umis_in_parallel,
)

logger = logging.getLogger('IsoQuant')


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


# Class for processing all samples against gene database
class DatasetProcessor:
    def __init__(self, args):
        self.args = args
        self.input_data = args.input_data
        self.args.gunzipped_reference = None
        self.common_header = "# Command line: " + args._cmd_line + "\n# IsoQuant version: " + args._version + "\n"
        self.io_support = IOSupport(self.args)
        self.all_read_groups = []  # Will be initialized per sample as list of sets
        self.grouping_strategy_names = get_grouping_strategy_names(self.args)
        self.alignment_stat_counter = EnumStats()
        self.transcript_type_dict = {}

        if args.genedb:
            logger.info("Loading gene database from " + self.args.genedb)
            self.gffutils_db = gffutils.FeatureDB(self.args.genedb)
            # TODO remove
            if self.args.mode.needs_pcr_deduplication():
                self.transcript_type_dict = create_transcript_info_dict(self.args.genedb)
        else:
            self.gffutils_db = None

        if self.args.needs_reference:
            logger.info("Loading reference genome from %s" % self.args.reference)
            self.reference_record_dict = Fasta(self.args.reference, indexname=args.fai_file_name)
        else:
            self.reference_record_dict = None
        self.chr_ids = []

    def __del__(self):
        pass

    def clean_up(self):
        if not self.args.keep_tmp and self.args.gunzipped_reference:
            if os.path.exists(self.args.gunzipped_reference):
                os.remove(self.args.gunzipped_reference)

        for sample in self.input_data.samples:
            if not self.args.read_assignments and not self.args.keep_tmp:
                for f in glob.glob(sample.out_raw_file + "_*"):
                    os.remove(f)
                for f in glob.glob(sample.read_group_file + "*"):
                    os.remove(f)
                if self.args.mode.needs_pcr_deduplication():
                    os.remove(umi_filtered_global_lock_file_name(sample.out_umi_filtered_done))

    def process_all_samples(self, input_data):
        logger.info("Processing " + proper_plural_form("experiment", len(input_data.samples)))
        logger.info("Secondary alignments will%s be used" % ("" if self.args.use_secondary else " not"))
        for sample in input_data.samples:
            self.process_sample(sample)
        self.clean_up()
        logger.info("Processed " + proper_plural_form("experiment", len(self.input_data.samples)))

    # Run through all genes in db and count stats according to alignments given in bamfile_name
    def process_sample(self, sample):
        logger.info("Processing experiment " + sample.prefix)
        logger.info("Experiment has " + proper_plural_form("BAM file", len(sample.file_list)) + ": " + ", ".join(
            map(lambda x: x[0], sample.file_list)))
        self.chr_ids = self.get_chromosome_ids(sample)

        logger.info("Total number of chromosomes to be processed %d: %s " %
                    (len(self.chr_ids), ", ".join(map(lambda x: str(x), sorted(self.chr_ids)))))

        # Check if file_name grouping is enabled for this sample (for technical replicas)
        sample.use_technical_replicas = (self.args.use_replicas and
                                         len(sample.file_list) > 1 and
                                         self.args.read_group is not None and
                                         "file_name" in self.args.read_group)
        if sample.use_technical_replicas:
            logger.info("Technical replicas filtering enabled: novel transcripts must be confirmed by at least 2 files")
        elif len(sample.file_list) > 1 and not self.args.use_replicas:
            logger.info("Technical replicas filtering disabled by --use_replicas false")

        # Initialize all_read_groups as list of sets (one per grouping strategy)
        self.all_read_groups = [set() for _ in self.grouping_strategy_names]
        fname = read_group_lock_filename(sample)
        if self.args.resume and os.path.exists(fname):
            logger.info("Read group table was split during the previous run, existing files will be used")
        else:
            if os.path.exists(fname):
                os.remove(fname)
            prepare_read_groups(self.args, sample)
            open(fname, "w").close()

        if self.args.mode.needs_pcr_deduplication():
            if self.args.barcoded_reads:
                sample.barcoded_reads = self.args.barcoded_reads

            split_barcodes_dict = {}
            for chr_id in self.get_chr_list():
                split_barcodes_dict[chr_id] = sample.barcodes_split_reads + "_" + chr_id
            barcode_split_done = split_barcodes_lock_filename(sample)
            if self.args.resume and os.path.exists(barcode_split_done):
                logger.info("Barcode table was split during the previous run, existing files will be used")
            else:
                if os.path.exists(barcode_split_done):
                    os.remove(barcode_split_done)
                self.split_read_barcode_table(sample, split_barcodes_dict)
                open(barcode_split_done, "w").close()

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

            if self.args.mode.needs_pcr_deduplication():
                self.filter_umis(sample)

        total_assignments, polya_found, self.all_read_groups = self.load_read_info(saves_file)

        if self.args.mode.needs_pcr_deduplication():
            # move clean-up somewhere else
            for bc_split_file in split_barcodes_dict.values():
                os.remove(bc_split_file)
            os.remove(barcode_split_done)

        polya_fraction = polya_found / total_assignments if total_assignments > 0 else 0.0
        logger.info("Total assignments used for analysis: %d, polyA tail detected in %d (%.1f%%)" %
                    (total_assignments, polya_found, polya_fraction * 100.0))
        if (polya_fraction < self.args.low_polya_percentage_threshold and
                self.args.polya_requirement_strategy != PolyAUsageStrategies.never):
            logger.warning("PolyA percentage is suspiciously low. IsoQuant expects non-polya-trimmed reads. "
                           "If you aim to construct transcript models, consider using --polya_trimmed and --polya_requirement options.")

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
        logger.info("Processed experiment " + sample.prefix)

    def keep_only_defined_chromosomes(self, chr_set: set):
        if self.args.process_only_chr:
            chr_set.intersection_update(self.args.process_only_chr)
        elif self.args.discard_chr:
            chr_set.difference_update(self.args.discard_chr)

        return chr_set

    def get_chromosome_ids(self, sample):
        genome_chromosomes = set(self.reference_record_dict.keys())
        genome_chromosomes = self.keep_only_defined_chromosomes(genome_chromosomes)

        bam_chromosomes = set()
        for bam_file in list(map(lambda x: x[0], sample.file_list)):
            bam = pysam.AlignmentFile(bam_file, "rb", require_index=True)
            bam_chromosomes.update(bam.references)
        bam_chromosomes = self.keep_only_defined_chromosomes(bam_chromosomes)

        bam_genome_overlap = genome_chromosomes.intersection(bam_chromosomes)
        if len(bam_genome_overlap) != len(genome_chromosomes) or len(bam_genome_overlap) != len(bam_chromosomes):
            if len(bam_genome_overlap) == 0:
                logger.critical("Chromosomes in the BAM file(s) have different names than chromosomes in the reference"
                                " genome. Make sure that the same genome was used to generate your BAM file(s).")
                sys.exit(IsoQuantExitCode.CHROMOSOME_MISMATCH)
            else:
                logger.warning("Chromosome list from the reference genome is not the same as the chromosome list from"
                               " the BAM file(s). Make sure that the same genome was used to generate the BAM file(s).")
                logger.warning("Only %d overlapping chromosomes will be processed." % len(bam_genome_overlap))

        if not self.args.genedb:
            return list(sorted(
                bam_genome_overlap,
                key=lambda x: len(self.reference_record_dict[x]),
                reverse=True,
            ))

        gene_annotation_chromosomes = set()
        gffutils_db = gffutils.FeatureDB(self.args.genedb)
        for feature in gffutils_db.all_features():
            gene_annotation_chromosomes.add(feature.seqid)
        gene_annotation_chromosomes = self.keep_only_defined_chromosomes(gene_annotation_chromosomes)

        common_overlap = gene_annotation_chromosomes.intersection(bam_genome_overlap)
        if len(common_overlap) != len(gene_annotation_chromosomes):
            if len(common_overlap) == 0:
                logger.critical("Chromosomes in the gene annotation have different names than chromosomes in the "
                                "reference genome or BAM file(s). Please, check the input data.")
                sys.exit(IsoQuantExitCode.CHROMOSOME_MISMATCH)
            else:
                logger.warning("Chromosome list from the gene annotation is not the same as the chromosome list from"
                               " the reference genomes or BAM file(s). Please, check you input data.")
                logger.warning("Only %d overlapping chromosomes will be processed." % len(common_overlap))

        return list(sorted(
            common_overlap,
            key=lambda x: len(self.reference_record_dict[x]),
            reverse=True,
        ))

    def get_chr_list(self):
        return self.chr_ids

    def collect_reads(self, sample):
        logger.info('Collecting read alignments')
        chr_ids = self.get_chr_list()
        info_file = sample.get_info_file()
        lock_file = sample.get_collection_lock_file()

        if os.path.exists(lock_file):
            if self.args.resume:
                logger.info("Collected reads detected, will not process")
                return
            else:
                os.remove(lock_file)

        if not self.args.resume:
            clean_locks(chr_ids, sample.out_raw_file, reads_collected_lock_file_name)
            clean_locks(chr_ids, sample.out_raw_file, reads_processed_lock_file_name)

        if not self.args.use_secondary:
            processed_read_manager_type = ProcessedReadsManagerNoSecondary
        elif self.args.high_memory:
            processed_read_manager_type = ProcessedReadsManagerHighMemory
        else:
            processed_read_manager_type = ProcessedReadsManagerNormalMemory

        read_gen = (
            collect_reads_in_parallel,
            itertools.repeat(sample),
            chr_ids,
            itertools.repeat(chr_ids),
            itertools.repeat(self.args),
            itertools.repeat(processed_read_manager_type)
        )

        # Initialize all_read_groups as list of sets (one per grouping strategy)
        all_read_groups = [set() for _ in self.grouping_strategy_names]

        if self.args.threads > 1:
            # Clean up parent memory before spawning workers
            gc.collect()
            mp_context = multiprocessing.get_context('fork')
            with ProcessPoolExecutor(max_workers=self.args.threads, mp_context=mp_context) as proc:
                results = proc.map(*read_gen, chunksize=1)
        else:
            results = map(*read_gen)

        sample_procesed_read_manager = processed_read_manager_type(sample, self.args.multimap_strategy, chr_ids, self.args.genedb)
        logger.info("Counting multimapped reads")
        for chr_id, read_groups, alignment_stats, processed_reads in results:
            logger.info("Counting reads from %s" % chr_id)
            # read_groups can be either a single set or a list of sets
            if isinstance(read_groups, list):
                # MultiReadGrouper returns list of sets
                for i, group_set in enumerate(read_groups):
                    all_read_groups[i].update(group_set)
            else:
                # Single grouper returns a set
                all_read_groups[0].update(read_groups)
            self.alignment_stat_counter.merge(alignment_stats)
            sample_procesed_read_manager.merge(processed_reads, chr_id)

        logger.info("Resolving multimappers")
        total_assignments, polya_assignments = sample_procesed_read_manager.resolve()
        logger.info("Multimappers resolved")

        for bam_file in list(map(lambda x: x[0], sample.file_list)):
            bam = pysam.AlignmentFile(bam_file, "rb", require_index=True)
            self.alignment_stat_counter.add(AlignmentType.unaligned, bam.unmapped)
        self.alignment_stat_counter.print_start("Alignments collected, overall alignment statistics:")

        info_dumper = open(info_file, "wb")
        write_int(total_assignments, info_dumper)
        write_int(polya_assignments, info_dumper)
        # Save all_read_groups as list of lists (convert sets to lists)
        write_int(len(all_read_groups), info_dumper)  # Number of grouping strategies
        for group_set in all_read_groups:
            write_list(list(group_set), info_dumper, write_string)
        info_dumper.close()
        open(lock_file, "w").close()

        if total_assignments == 0:
            logger.warning("No reads were assigned to isoforms, check your input files")
        else:
            logger.info('Finishing read assignment, total assignments %d, polyA percentage %.1f' %
                        (total_assignments, 100 * polya_assignments / total_assignments))

    def process_assigned_reads(self, sample, dump_filename):
        chr_ids = self.get_chr_list()
        logger.info("Processing assigned reads " + sample.prefix)
        logger.info("Transcript models construction is turned %s" %
                    ("off" if self.args.no_model_construction else "on"))

        # set up aggregators and outputs
        aggregator = ReadAssignmentAggregator(self.args, sample, self.all_read_groups, gzipped=self.args.gzipped,
                                             grouping_strategy_names=self.grouping_strategy_names)
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

            # GFF printers below only serve for creating the main output files,
            # not intended for dumping transcript models directly
            exon_id_storage = FeatureIdStorage(SimpleIDDistributor())
            gff_printer = GFFPrinter(
                sample.out_dir, sample.prefix, exon_id_storage,
                output_r2t=large_output_enabled(self.args, "read2transcripts"),
                header=self.common_header,
                gzipped=self.args.gzipped
            )
            if self.args.genedb:
                extended_gff_printer = GFFPrinter(
                    sample.out_dir, sample.prefix, exon_id_storage,
                    gtf_suffix=".extended_annotation.gtf", output_r2t=False,
                    header=self.common_header
                )

        model_gen = (
            construct_models_in_parallel,
            itertools.repeat(sample),
            chr_ids,
            itertools.repeat(chr_ids),
            itertools.repeat(dump_filename),
            itertools.repeat(self.args),
            itertools.repeat(self.all_read_groups),
        )

        if self.args.threads > 1:
            # Clean up parent memory before spawning workers
            gc.collect()
            mp_context = multiprocessing.get_context('fork')
            with ProcessPoolExecutor(max_workers=self.args.threads, mp_context=mp_context) as proc:
                results = proc.map(*model_gen, chunksize=1)
        else:
            results = map(*model_gen)

        for read_stat_counter, tsc in results:
            for k, v in read_stat_counter.stats_dict.items():
                aggregator.read_stat_counter.stats_dict[k] += v

            if not self.args.no_model_construction:
                for k, v in tsc.stats_dict.items():
                    transcript_stat_counter.stats_dict[k] += v

        self.merge_assignments(sample, aggregator, chr_ids)
        if self.args.sqanti_output:
            merge_files(sample.out_t2t_tsv, sample.prefix, chr_ids, aggregator.t2t_sqanti_printer.output_file, copy_header=False)
        self.finalize(aggregator)

        if not self.args.no_model_construction:
            transcript_stat_counter.print_start("Transcript model statistics")
            self.merge_transcript_models(sample.prefix, aggregator, chr_ids, gff_printer)
            logger.info("Transcript model file " + gff_printer.model_fname)
            if self.args.genedb:
                merge_files(extended_gff_printer.model_fname, sample.prefix, chr_ids,
                            extended_gff_printer.out_gff, copy_header=False)
                logger.info("Extended annotation is saved to " + extended_gff_printer.model_fname)

            logger.info("Counts for generated transcript models are saves to: " +
                        aggregator.transcript_model_counter.output_counts_file_name)
            if self.args.read_group:
                for counter in aggregator.transcript_model_global_counter.counters:
                    if counter.ignore_read_groups:
                        continue
                    logger.info("Grouped counts for discovered transcript models are saves to: " +
                                counter.output_counts_file_name)
            aggregator.transcript_model_global_counter.finalize(self.args)
            aggregator.gene_model_global_counter.finalize(self.args)

    def finalize(self, aggregator):
        if aggregator.basic_printer:
            logger.info("Read assignments are stored in " + aggregator.basic_printer.output_file_name +
                        (".gz" if self.args.gzipped else ""))
        if self.args.genedb:
            aggregator.read_stat_counter.print_start("Read assignment statistics")

            logger.info("Gene counts are stored in " + aggregator.gene_counter.output_counts_file_name)
            logger.info("Transcript counts are stored in " + aggregator.transcript_counter.output_counts_file_name)
            if self.args.read_group:
                for counter in aggregator.global_counter.counters:
                    if counter.ignore_read_groups:
                        continue
                    logger.info("Grouped counts are saves to: " + counter.output_counts_file_name)
            logger.info("Counts can be converted to other formats using src/convert_grouped_counts.py")
            aggregator.global_counter.finalize(self.args)

    def filter_umis(self, sample):
        umi_filtering_done = umi_filtered_global_lock_file_name(sample.out_umi_filtered_done)
        if os.path.exists(umi_filtering_done):
            if self.args.resume:
                logger.info("UMI filtering detecting, skipping")
                return
            os.remove(umi_filtering_done)

        # edit distances for UMI filtering, first one will be used for counts
        umi_ed_dict = {IsoQuantMode.bulk: [],
                       IsoQuantMode.tenX_v3: [3],
                       IsoQuantMode.visium_5prime: [3],
                       IsoQuantMode.curio: [3],
                       IsoQuantMode.visium_hd: [4],
                       IsoQuantMode.stereoseq: [4],
                       IsoQuantMode.stereoseq_nosplit: [4]}

        for i, edit_distance in enumerate(umi_ed_dict[self.args.mode]):
            logger.info("Filtering PCR duplicates with edit distance %d" % edit_distance)
            umi_ed_filtering_done = umi_filtered_lock_file_name(sample.out_umi_filtered_done, "", edit_distance)
            if os.path.exists(umi_ed_filtering_done):
                if self.args.resume:
                    logger.info("Filtering was done previously, skipping edit distance %d" % edit_distance)
                    return
                os.remove(umi_ed_filtering_done)

            output_prefix = umi_output_prefix(sample.out_umi_filtered, edit_distance)
            logger.info("Results will be saved to %s" % output_prefix)
            output_filtered_reads = i == 0

            umi_gen = (
                filter_umis_in_parallel,
                itertools.repeat(sample),
                self.get_chr_list(),
                itertools.repeat(self.get_chr_list()),
                itertools.repeat(self.args),
                itertools.repeat(edit_distance),
                itertools.repeat(output_filtered_reads),
            )

            if self.args.threads > 1:
                # Clean up parent memory before spawning workers
                gc.collect()
                mp_context = multiprocessing.get_context('fork')
                with ProcessPoolExecutor(max_workers=self.args.threads, mp_context=mp_context) as proc:
                    results = proc.map(*umi_gen, chunksize=1)
            else:
                results = map(*umi_gen)

            stat_dict = defaultdict(int)
            files_to_remove = []
            save_allinfo = large_output_enabled(self.args, "allinfo")
            if save_allinfo:
                allinfo_fname = output_prefix + ".allinfo"
                if self.args.gzipped:
                    allinfo_fname += ".gz"
                    allinfo_outf = gzip.open(allinfo_fname, "wt")
                else:
                    allinfo_outf = open(allinfo_fname, "w")

            for all_info_file_name, stats_output_file_name, umi_filter_done in results:
                if save_allinfo:
                    shutil.copyfileobj(open(all_info_file_name, "r"), allinfo_outf)
                for l in open(stats_output_file_name, "r"):
                    v = l.strip().split("\t")
                    if len(v) != 2:
                        continue
                    stat_dict[v[0]] += int(v[1])
                files_to_remove.append(all_info_file_name)
                files_to_remove.append(stats_output_file_name)
                files_to_remove.append(umi_filter_done)

            if save_allinfo:
                allinfo_outf.close()

            logger.info("PCR duplicates filtered with edit distance %d, filtering stats:" % edit_distance)
            with open(output_prefix + ".stats.tsv", "w") as outf:
                for k, v in stat_dict.items():
                    logger.info("  %s: %d" % (k, v))
                    outf.write("%s\t%d\n" % (k, v))

            open(umi_ed_filtering_done, "w").close()
            for f in files_to_remove:
                os.remove(f)

        open(umi_filtering_done, "w").close()

    def split_read_barcode_table(self, sample, split_barcodes_file_names):
        logger.info("Splitting read barcode table")
        # TODO: untrusted UMIs and third party format, both can be done by passing parsing function instead of columns
        split_read_table_parallel(sample, sample.barcoded_reads, split_barcodes_file_names, self.args.threads,
                                  read_column=0, group_columns=(1, 2, 3, 4), delim='\t')
        logger.info("Read barcode table was split")


    @staticmethod
    def load_read_info(dump_filename):
        info_loader = open(info_file_name(dump_filename), "rb")
        total_assignments = read_int(info_loader)
        polya_assignments = read_int(info_loader)
        # Load all_read_groups as list of sets
        num_strategies = read_int(info_loader)
        all_read_groups = []
        for _ in range(num_strategies):
            group_list = read_list(info_loader, read_string)
            all_read_groups.append(set(group_list))
        info_loader.close()
        return total_assignments, polya_assignments, all_read_groups

    def merge_assignments(self, sample, aggregator, chr_ids):
        if self.args.genedb and aggregator.basic_printer:
            merge_files(sample.out_assigned_tsv, sample.prefix, chr_ids,
                        aggregator.basic_printer.output_file, copy_header=False)
        if aggregator.corrected_bed_printer:
            merge_files(sample.out_corrected_bed, sample.prefix, chr_ids,
                        aggregator.corrected_bed_printer.output_file, copy_header=False)

        for counter in aggregator.global_counter.counters:
            unaligned = self.alignment_stat_counter.stats_dict[AlignmentType.unaligned]
            merge_counts(counter, sample.prefix, chr_ids, unaligned)
            # counter.convert_counts_to_tpm(self.args.normalization_method)

    def merge_transcript_models(self, label, aggregator, chr_ids, gff_printer):
        merge_files(gff_printer.model_fname, label, chr_ids, gff_printer.out_gff, copy_header=False)
        if gff_printer.output_r2t:
            merge_files(gff_printer.r2t_fname, label, chr_ids, gff_printer.out_r2t, copy_header=False)
        for counter in aggregator.transcript_model_global_counter.counters:
            unaligned = self.alignment_stat_counter.stats_dict[AlignmentType.unaligned]
            merge_counts(counter, label, chr_ids, unaligned)
        for counter in aggregator.gene_model_global_counter.counters:
            merge_counts(counter, label, chr_ids)
