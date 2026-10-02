############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import gc
import glob
import itertools
import logging
import multiprocessing
import os
import shutil
import sys
from enum import Enum, unique
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from typing import Dict, List, Optional

import gffutils
import pysam
from pyfaidx import Fasta

from .modes import IsoQuantMode
from .common import proper_plural_form, large_output_enabled, setup_worker_logging, _get_log_params
from isoquant_lib.utils.error_codes import IsoQuantExitCode
from isoquant_lib.utils.serialization import (
    read_int,
    read_list,
    read_string,
    write_int,
    write_list,
    write_string,
)
from isoquant_lib.utils.stats import EnumStats
from isoquant_lib.utils.file_utils import (merge_files, merge_counts, merge_file_list, gzip_file_in_place,
                                          load_usable_fragments, usable_fragment_list,
                                          open_text_write, read_stats_tsv,
                                          resolve_optionally_gzipped)
from isoquant_lib.utils.checkpoints import (CheckpointStore, enum_stats_to_payload, mark_stage_done, marker_name,
                                           payload_to_enum_stats, remove_files, run_checkpoint_dir, run_stage)
from isoquant_lib.utils.bam_utils import (PLACEHOLDERS, collect_unmapped_read_ids,
                                         load_barcode_umi_tags, merge_bam_files,
                                         references_with_alignments, write_unmapped_bam)
from .alignment.alignment_processor import AlignmentType
from .report.run_summary import RunSummary
from .report.html_report import render_html
from .terminal_prediction.external_polya import normalize_polya_reads
from .assignment.read_groups import prepare_read_groups, get_grouping_strategy_names
from .assignment.assignment_io import (IOSupport, BasicTSVAssignmentPrinter, BEDPrinter, ReadInfoPrinter,
                                      SqantiTSVPrinter)
from .terminal_prediction.terminal_counter import TRAINING_SUFFIX
from .processed_read_manager import ProcessedReadsManagerHighMemory, ProcessedReadsManagerNoSecondary, ProcessedReadsManagerNormalMemory
from isoquant_lib.utils.id_policy import SimpleIDDistributor, FeatureIdStorage
from isoquant_lib.utils.file_naming import (
    allinfo_file_name,
    allinfo_stats_file_name,
    info_file_name,
    read_stat_file_name,
    saves_file_name,
    transcript_stat_file_name,
    umi_barcode2barcode_prefix,
    umi_output_prefix,
    tagged_bam_fragment_name,
)
from isoquant_lib.model_construction.transcript_printer import GFFPrinter
from .barcode_calling.umi_filtering import create_transcript_info_dict
from isoquant_lib.utils.table_splitter import split_read_table_parallel
from .assignment.assignment_aggregator import ReadAssignmentAggregator
from isoquant_lib.utils.input_data_storage import SampleData
from isoquant_lib.utils.string_pools import setup_string_pools
from .parallel_workers import (
    collect_reads_in_parallel,
    construct_models_in_parallel,
    filter_umis_in_parallel,
    write_deduplicated_bam_in_parallel,
    write_tagged_bam_in_parallel,
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


# Edit distances for UMI filtering, per mode. The first one is the round whose
# survivors are counted, and the one the run summary reports.
UMI_EDIT_DISTANCES = {IsoQuantMode.bulk: [],
                      IsoQuantMode.tenX_v3: [3],
                      IsoQuantMode.tenX_v2: [3],
                      IsoQuantMode.visium_5prime: [3],
                      IsoQuantMode.curio: [3],
                      IsoQuantMode.visium_hd: [4],
                      IsoQuantMode.stereoseq: [4],
                      IsoQuantMode.custom_sc: [4]}


# Class for processing all samples against gene database
class DatasetProcessor:
    def __init__(self, args):
        self.args = args
        self.input_data = args.input_data
        self.common_header = "# Command line: " + args._cmd_line + "\n# IsoQuant version: " + args._version + "\n"
        self.io_support = IOSupport(self.args)
        self.all_read_groups = []  # Will be initialized per sample as list of sets
        self.grouping_strategy_names = get_grouping_strategy_names(self.args)
        self.alignment_stat_counter = EnumStats()
        self.transcript_type_dict = {}
        # Per-sample QC summary, written as SAMPLE.summary_stats.json / SAMPLE.summary.html
        self.run_summary: Optional[RunSummary] = None
        # each experiment derives its polyA requirements from the options as given,
        # never from the flags the previous experiment left in args
        self.original_require_monointronic_polya = args.require_monointronic_polya
        self.original_require_monoexonic_polya = args.require_monoexonic_polya

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

    def sample_store(self, sample: SampleData) -> CheckpointStore:
        return CheckpointStore(sample.checkpoint_dir, self.args.resume)

    def process_all_samples(self, input_data):
        logger.info("Processing " + proper_plural_form("experiment", len(input_data.samples)))
        logger.info("Secondary alignments will%s be used" % ("" if self.args.use_secondary else " not"))
        run_store = CheckpointStore(run_checkpoint_dir(self.args.output), self.args.resume)
        for sample in input_data.samples:
            sample_marker = marker_name("sample", sample.prefix)
            if run_store.is_done(sample_marker):
                logger.info("Experiment %s was completed during the previous run, skipping" % sample.prefix)
            else:
                self.process_sample(sample)
                mark_stage_done(run_store, sample_marker)
            # after the sample marker only: a crash before it must find every intermediate
            self.remove_sample_intermediates(sample)
        logger.info("Processed " + proper_plural_form("experiment", len(self.input_data.samples)))

    def remove_sample_intermediates(self, sample: SampleData) -> None:
        """Delete what the stages of a completed sample shared; by glob, so it also works
        for a sample skipped on resume. Idempotent."""
        if self.args.keep_tmp:
            return
        patterns = [glob.escape(sample.barcodes_split_reads) + "_*",
                    glob.escape(sample.polya_split_reads) + "_*",
                    glob.escape(sample.polya_reads_normalized)]
        if not self.args.read_assignments:
            # saves of another run (--read_assignments) are not ours to delete
            patterns += [glob.escape(sample.out_raw_file) + "_*",
                         glob.escape(sample.read_group_file) + "*"]
        for pattern in patterns:
            remove_files(glob.glob(pattern))
        shutil.rmtree(sample.chr_fragment_dir, ignore_errors=True)

    # Run through all genes in db and count stats according to alignments given in bamfile_name
    def process_sample(self, sample):
        """Run the stages of one experiment; see .claude/RESUME_CHECKPOINTS.md.

        Prologue, derive and post-construct steps always run: they rebuild the in-memory
        state later steps need, whichever stages the previous run already finished.
        """
        logger.info("Processing experiment " + sample.prefix)
        store = self.sample_store(sample)
        keep_tmp = self.args.keep_tmp

        # ---- prologue (always runs)
        self.run_summary = RunSummary(sample.prefix, isoquant_version=self.args._version,
                                      command_line=self.args._cmd_line, mode=self.args.mode.name)
        # Alignment statistics are per experiment: without the reset a second sample
        # would report (and store as __not_aligned) the counts of the previous one too.
        self.alignment_stat_counter = EnumStats()
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
        if not self.args.no_model_construction:
            if sample.use_technical_replicas:
                logger.info("Technical replicas filtering enabled: novel transcripts must be confirmed by at least 2 files")
            elif len(sample.file_list) > 1 and not self.args.use_replicas:
                logger.info("Technical replicas filtering disabled by --use_replicas false")

        # Initialize all_read_groups as list of sets (one per grouping strategy)
        self.all_read_groups = [set() for _ in self.grouping_strategy_names]

        uses_barcode_table = self.args.mode.needs_pcr_deduplication() and not getattr(self.args, 'barcoded_bam', False)
        split_barcodes_dict = {}
        tagged_bam_references = None
        if uses_barcode_table:
            if self.args.barcoded_reads:
                sample.barcoded_reads = self.args.barcoded_reads
            # a tagged BAM copies every reference, including the unplaced scaffolds IsoQuant
            # does not analyse, so those reads need a barcode table of their own to be tagged.
            # The same list drives the split and the copy, so every fragment finds its table.
            if large_output_enabled(self.args, "tagged_bam"):
                tagged_bam_references = references_with_alignments([f[0] for f in sample.file_list])
            split_chr_ids = self.get_chr_list() if tagged_bam_references is None else tagged_bam_references
            split_barcodes_dict = {chr_id: sample.get_barcodes_split_file(chr_id) for chr_id in split_chr_ids}

        # ---- input preparation stages
        run_stage(store, "read_groups_split", run=lambda: prepare_read_groups(self.args, sample))
        run_stage(store, "polya_split", run=lambda: self.split_external_polya(sample),
                  cleanup=lambda: remove_files([sample.polya_reads_normalized]),
                  enabled=self.args.polya_trimmed.uses_read_table() and not self.args.read_assignments,
                  keep_tmp=keep_tmp)
        run_stage(store, "barcode_split", run=lambda: self.split_read_barcode_table(sample, split_barcodes_dict),
                  enabled=uses_barcode_table)
        # nothing to tag with otherwise; the user was warned about that at startup
        run_stage(store, "tagged_bam", run=lambda: self.write_tagged_bam(sample, tagged_bam_references) and None,
                  enabled=uses_barcode_table and tagged_bam_references is not None)

        # ---- read collection and UMI filtering
        if self.args.read_assignments:
            saves_file = self.args.read_assignments[0]
            logger.info('Using read assignments from {}*'.format(saves_file))
        else:
            run_stage(store, "collect", run=lambda: self.collect_reads(sample, store),
                      restore=self.restore_collection_stats)
            saves_file = sample.out_raw_file
            logger.info('Read assignments files saved to {}*. '.
                        format(sample.out_raw_file))
            if not keep_tmp:
                logger.info("To keep these intermediate files for debug purposes use --keep_tmp flag")

            if self.args.mode.needs_pcr_deduplication():
                self.filter_umis(sample, store)
                # the survivors files live under out_raw_file and are deleted with the sample intermediates
                run_stage(store, "dedup_bam", run=lambda: self.write_deduplicated_bam(sample) and None,
                          enabled=large_output_enabled(self.args, "deduplicated_bam"))

        # ---- derive (always runs)
        total_assignments, polya_found, self.all_read_groups = self.load_read_info(saves_file)

        # Warn if large barcode count with per-barcode grouping
        if "barcode" in self.grouping_strategy_names:
            barcode_idx = self.grouping_strategy_names.index("barcode")
            barcode_count = len(self.all_read_groups[barcode_idx])
            if barcode_count > 10000:
                logger.warning("Large number of barcodes detected (%d). Per-barcode grouping may consume "
                               "substantial RAM. Consider using '--barcode2spot' to group barcodes by "
                               "cell type or spatial region.", barcode_count)

        polya_fraction = polya_found / total_assignments if total_assignments > 0 else 0.0
        logger.info("Total assignments used for analysis: %d, polyA tail detected in %d (%.1f%%)" %
                    (total_assignments, polya_found, polya_fraction * 100.0))
        self.run_summary.set_polya_stats(total_assignments, polya_found)
        if (polya_fraction < self.args.low_polya_percentage_threshold and
                self.args.polya_requirement_strategy != PolyAUsageStrategies.never):
            logger.warning("PolyA percentage is suspiciously low. IsoQuant expects non-polya-trimmed reads. "
                           "If you aim to construct transcript models, consider using --polya_trimmed and --polya_requirement options.")

        self.args.requires_polya_for_construction = set_polya_requirement_strategy(
            polya_fraction >= self.args.polya_percentage_threshold,
            self.args.polya_requirement_strategy)
        # derived from the options as given, not from the previous experiment's flags
        self.args.require_monointronic_polya = set_polya_requirement_strategy(
            # do not require polyA tails for mono-intronic only if the data is reliable and polyA percentage is low
            self.original_require_monointronic_polya or self.args.requires_polya_for_construction,
            self.args.polya_requirement_strategy)
        self.args.require_monoexonic_polya = set_polya_requirement_strategy(
            # do not require polyA tails for mono-intronic only if the data is reliable and polyA percentage is low
            self.original_require_monoexonic_polya or self.args.requires_polya_for_construction,
            self.args.polya_requirement_strategy)

        # ---- transcript model construction and quantification
        run_stage(store, "construct", run=lambda: self.construct_models(sample, saves_file, store))
        self.report_construction_stats(saves_file)
        self.merge_outputs(sample, store)

        # ---- finish
        self.write_run_summary(sample)
        self.compress_barcode_tables(sample)
        logger.info("Processed experiment " + sample.prefix)

    def write_run_summary(self, sample: SampleData) -> None:
        """Write the QC summary of this experiment. Never fails the run: the numbers
        are a report, the analysis outputs are already on disk at this point."""
        if getattr(self.args, "no_report", False) or self.run_summary is None:
            return
        try:
            # The counts come from the first filtering round, so that is the one the
            # summary describes: several .ED<N>.stats.tsv can be left in the directory.
            edit_distances = UMI_EDIT_DISTANCES.get(self.args.mode) or []
            self.run_summary.collect_output_files(
                sample, self.grouping_strategy_names,
                umi_edit_distance=edit_distances[0] if edit_distances else None)
            json_file = sample.out_summary_json
            html_file = sample.out_summary_html
            self.run_summary.write_json(json_file)
            render_html(self.run_summary, html_file)
            logger.info("Run summary is stored in %s (and %s)" % (html_file, json_file))
        except Exception as e:
            logger.warning("Could not write the run summary: %s" % e)
            logger.debug("Run summary error", exc_info=True)

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

    def collect_reads(self, sample, store: CheckpointStore) -> Dict:
        """The collect stage: per-chromosome alignment processing, multimapper resolution.

        Returns the alignment statistics as the stage payload; a resumed run restores them
        with restore_collection_stats, everything else is read back from the _info file.
        """
        logger.info('Collecting read alignments')
        chr_ids = self.get_chr_list()
        info_file = sample.get_info_file()

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
            itertools.repeat(processed_read_manager_type),
            itertools.repeat(store),
            itertools.repeat("collect"),
        )

        # Initialize all_read_groups as list of sets (one per grouping strategy)
        all_read_groups = [set() for _ in self.grouping_strategy_names]

        if self.args.threads > 1:
            # Clean up parent memory before spawning workers
            gc.collect()
            mp_context = multiprocessing.get_context('fork')
            log_file, log_level = _get_log_params()
            with ProcessPoolExecutor(max_workers=self.args.threads, mp_context=mp_context,
                                     initializer=setup_worker_logging,
                                     initargs=(log_file, log_level)) as proc:
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

        bam_files = list(map(lambda x: x[0], sample.file_list))
        for bam_file in bam_files:
            bam = pysam.AlignmentFile(bam_file, "rb", require_index=True)
            self.alignment_stat_counter.add(AlignmentType.unaligned, bam.unmapped)
        skipped_references = set(references_with_alignments(bam_files)) - set(chr_ids)
        if skipped_references:
            logger.info("%d reference(s) carrying alignments were not processed, "
                        "the run summary reports reads on processed references"
                        % len(skipped_references))
        self.alignment_stat_counter.print_start("Alignments collected, overall alignment statistics:")
        # Primary alignments are only counted for the chromosomes that were processed,
        # while bam.unmapped covers the whole file: their sum is the number of input
        # reads only when nothing was left out (--process_only_chr / --discard_chr, or
        # an annotation and a genome that do not cover the same references).
        covers_all_reads = not skipped_references
        self.run_summary.set_alignment_stats(self.alignment_stat_counter.stats_dict,
                                             covers_all_reads=covers_all_reads)

        with open(info_file, "wb") as info_dumper:
            write_int(total_assignments, info_dumper)
            write_int(polya_assignments, info_dumper)
            # Save all_read_groups as list of lists (convert sets to lists)
            write_int(len(all_read_groups), info_dumper)  # Number of grouping strategies
            for group_set in all_read_groups:
                write_list(list(group_set), info_dumper, write_string)

        if total_assignments == 0:
            logger.warning("No reads were assigned to isoforms, check your input files")
        else:
            logger.info('Finishing read assignment, total assignments %d, polyA percentage %.1f' %
                        (total_assignments, 100 * polya_assignments / total_assignments))
        return {"alignment_stats": enum_stats_to_payload(self.alignment_stat_counter.stats_dict),
                "covers_all_reads": covers_all_reads}

    def restore_collection_stats(self, payload: Optional[Dict]) -> None:
        """Re-apply what a skipped collect stage would have left in memory."""
        if not payload:
            return
        for element, count in payload_to_enum_stats(payload["alignment_stats"], AlignmentType).items():
            self.alignment_stat_counter.add(element, count)
        self.alignment_stat_counter.print_start("Alignment statistics from the previous run:")
        self.run_summary.set_alignment_stats(self.alignment_stat_counter.stats_dict,
                                             covers_all_reads=payload["covers_all_reads"])

    def construct_models(self, sample, saves_file: str, store: CheckpointStore) -> None:
        """The construct stage: per-chromosome assignment processing and model construction."""
        chr_ids = self.get_chr_list()
        logger.info("Processing assigned reads " + sample.prefix)
        logger.info("Transcript models construction is turned %s" %
                    ("off" if self.args.no_model_construction else "on"))
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

        self.map_over_chromosomes(construct_models_in_parallel, sample, chr_ids, saves_file, self.args,
                                  store, "construct")

    def report_construction_stats(self, saves_file: str) -> None:
        """Post-construct step (always runs): sum the per-chromosome statistics files."""
        read_stat_counter = EnumStats()
        transcript_stat_counter = EnumStats()
        for chr_id in self.get_chr_list():
            chr_dump_file = saves_file_name(saves_file, chr_id)
            read_stat_counter.merge(EnumStats(read_stat_file_name(chr_dump_file)))
            if not self.args.no_model_construction:
                transcript_stat_counter.merge(EnumStats(transcript_stat_file_name(chr_dump_file)))

        if self.args.genedb:
            read_stat_counter.print_start("Read assignment statistics")
            self.run_summary.set_assignment_stats(read_stat_counter.stats_dict)
        if not self.args.no_model_construction:
            transcript_stat_counter.print_start("Transcript model statistics")
            self.run_summary.set_transcript_model_stats(transcript_stat_counter.stats_dict)

    def merge_outputs(self, sample, store: CheckpointStore) -> None:
        """Merge the per-chromosome fragments into the final outputs, one stage per output.

        Every unit opens (and truncates) only its own output, and the aggregator is built
        with truncate_outputs=False, so rerunning a unit cannot empty an output finished
        before. Fragments are removed by each unit after its marker. A parent marker skips
        everything, including building the aggregator (unit names depend on its string pools).
        """
        if store.is_done("merge"):
            logger.info("merge: done in the previous run, skipping")
            return
        chr_ids = self.get_chr_list()
        label = sample.prefix
        fragment_dir = sample.chr_fragment_dir
        keep_tmp = self.args.keep_tmp
        construct = not self.args.no_model_construction

        # Build string pools for the merge phase (to convert group names <-> IDs)
        # This must match the pools used by parallel workers
        string_pools = None
        if self.args.read_group:
            string_pools = setup_string_pools(self.args, sample, chr_ids, chr_id=None,
                                              load_barcode_pool=False, load_tsv_pools=False)
        aggregator = ReadAssignmentAggregator(self.args, sample, string_pools, gzipped=self.args.gzipped,
                                              grouping_strategy_names=self.grouping_strategy_names,
                                              truncate_outputs=False)

        def fragments_of(file_name: str) -> List[str]:
            return merge_file_list(file_name, label, chr_ids, fragment_dir)

        def text_unit(name: str, enabled: bool, final_name: str, make_printer, header_lines: Optional[int],
                      message: str) -> None:
            def run() -> None:
                printer = make_printer()
                merge_files(final_name, label, chr_ids, printer.output_file, copy_header=False,
                            header_lines=header_lines, remove_inputs=False, fragment_dir=fragment_dir)
                printer.close()
                logger.info(message + printer.output_file_name)
            run_stage(store, marker_name("merge", name), run=run,
                      cleanup=lambda: remove_files(fragments_of(final_name)),
                      enabled=enabled, keep_tmp=keep_tmp)

        genedb = bool(self.args.genedb)
        text_unit("read_info", genedb and large_output_enabled(self.args, "read_info"), sample.out_read_info_tsv,
                  lambda: ReadInfoPrinter(sample.out_read_info_tsv, self.args, self.io_support,
                                          additional_header=aggregator.common_header, gzipped=self.args.gzipped),
                  3, "Read info is stored in ")
        text_unit("assignments", genedb and large_output_enabled(self.args, "read_assignments"), sample.out_assigned_tsv,
                  lambda: BasicTSVAssignmentPrinter(sample.out_assigned_tsv, self.args, self.io_support,
                                                    additional_header=aggregator.common_header,
                                                    gzipped=self.args.gzipped),
                  3, "Read assignments are stored in ")
        text_unit("corrected_bed", large_output_enabled(self.args, "corrected_bed"), sample.out_corrected_bed,
                  lambda: BEDPrinter(sample.out_corrected_bed, self.args, print_corrected=True,
                                     gzipped=self.args.gzipped),
                  None, "Corrected read alignments are stored in ")
        text_unit("sqanti_t2t", bool(self.args.sqanti_output), sample.out_t2t_tsv,
                  lambda: SqantiTSVPrinter(sample.out_t2t_tsv, self.args, self.io_support),
                  1, "SQANTI-like novel vs known comparison is stored in ")

        self.training_unit(store, "polya_training", getattr(self.args, "collect_polya_training", None),
                           sample.out_polya_prediction_tsv + TRAINING_SUFFIX, label, chr_ids, fragment_dir)
        tss_csv = getattr(self.args, "collect_tss_training", None) if self.args.fl_data else None
        self.training_unit(store, "tss_training", tss_csv,
                           sample.out_tss_prediction_tsv + TRAINING_SUFFIX, label, chr_ids, fragment_dir)

        unaligned = self.alignment_stat_counter.stats_dict.get(AlignmentType.unaligned, 0)
        if genedb:
            if self.args.run_quantification:
                logger.info("Gene counts are stored in " + aggregator.gene_counter.output_counts_file_name)
                logger.info("Transcript counts are stored in " + aggregator.transcript_counter.output_counts_file_name)
            if self.args.read_group:
                for counter in aggregator.global_counter.counters:
                    if counter.ignore_read_groups:
                        continue
                    logger.info("Grouped counts are saves to: " + counter.output_counts_file_name)
            logger.info("Counts can be converted to other formats using isoquant_lib/quantification/convert_grouped_counts.py")
        for counter in aggregator.global_counter.counters:
            self.counter_units(store, counter, label, chr_ids, fragment_dir, unaligned, finalize=genedb)

        if construct:
            exon_id_storage = FeatureIdStorage(SimpleIDDistributor())
            gff_name = os.path.join(sample.out_dir, sample.prefix + ".transcript_models.gtf")
            extended_gff_name = os.path.join(sample.out_dir, sample.prefix + ".extended_annotation.gtf")

            def gtf_unit(name: str, enabled: bool, final_name: str, gtf_suffix: str, message: str) -> None:
                def run() -> None:
                    gff_printer = GFFPrinter(sample.out_dir, sample.prefix, exon_id_storage,
                                             gtf_suffix=gtf_suffix, header=self.common_header)
                    merge_files(gff_printer.model_fname, label, chr_ids, gff_printer.out_gff, copy_header=False,
                                remove_inputs=False, fragment_dir=fragment_dir)
                    gff_printer.close()
                    logger.info(message + gff_printer.model_fname)
                run_stage(store, marker_name("merge", name), run=run,
                          cleanup=lambda: remove_files(fragments_of(final_name)),
                          enabled=enabled, keep_tmp=keep_tmp)

            gtf_unit("gtf", True, gff_name, ".transcript_models.gtf", "Transcript model file ")
            gtf_unit("extended_gtf", genedb, extended_gff_name, ".extended_annotation.gtf",
                     "Extended annotation is saved to ")
            # Per-chr model read_info files carry common_header (2) + column header
            # (1) = 3 header lines, matching the main read_info merge. The base
            # (ungzipped) path names the per-chr files; the final handle may be gzip.
            text_unit("model_reads", large_output_enabled(self.args, "read2transcripts"),
                      sample.out_transcript_model_reads_tsv,
                      lambda: ReadInfoPrinter(sample.out_transcript_model_reads_tsv, self.args, IOSupport(self.args),
                                              additional_header=self.common_header, gzipped=self.args.gzipped),
                      3, "Reads to transcript models are stored in ")

            logger.info("Counts for generated transcript models are saves to: " +
                        aggregator.transcript_model_counter.output_counts_file_name)
            if self.args.read_group:
                for counter in aggregator.transcript_model_global_counter.counters:
                    if counter.ignore_read_groups:
                        continue
                    logger.info("Grouped counts for discovered transcript models are saves to: " +
                                counter.output_counts_file_name)
            for counter in aggregator.transcript_model_global_counter.counters:
                self.counter_units(store, counter, label, chr_ids, fragment_dir, unaligned, finalize=True)
            for counter in aggregator.gene_model_global_counter.counters:
                self.counter_units(store, counter, label, chr_ids, fragment_dir, 0, finalize=True)

        if not keep_tmp:
            self._cleanup_per_chr_temp_files(sample, chr_ids)
        mark_stage_done(store, "merge")

    def counter_units(self, store: CheckpointStore, counter, label: str, chr_ids: List[str],
                      fragment_dir: str, unaligned: int, finalize: bool) -> None:
        """Two stages per counter: a mechanical merge of the fragments, then the conversion
        (TPM, matrix, MTX, loom), which depends on the options and the number of groups."""
        name = os.path.basename(counter.output_counts_file_name)
        keep_tmp = self.args.keep_tmp

        def merge_run() -> None:
            open(counter.output_file, "w").close()
            merge_counts(counter, label, chr_ids, unaligned, remove_inputs=False, fragment_dir=fragment_dir)

        def merge_cleanup() -> None:
            fragments = merge_file_list(counter.output_counts_file_name, label, chr_ids, fragment_dir)
            if counter.output_stats_file_name:
                fragments += merge_file_list(counter.output_stats_file_name, label, chr_ids, fragment_dir)
            remove_files(fragments)

        def finalize_run() -> None:
            load_usable_fragments(counter, label, chr_ids, fragment_dir)
            counter.finalize(self.args)

        run_stage(store, marker_name("merge", "counter", name), run=merge_run, cleanup=merge_cleanup,
                  keep_tmp=keep_tmp)
        run_stage(store, marker_name("merge", "counter_finalize", name), run=finalize_run,
                  cleanup=lambda: remove_files(usable_fragment_list(counter, label, chr_ids, fragment_dir)),
                  enabled=finalize, keep_tmp=keep_tmp)

    def training_unit(self, store: CheckpointStore, name: str, csv_path: Optional[str], fragment_base: str,
                      label: str, chr_ids: List[str], fragment_dir: str) -> None:
        """Developer training-data collection: concatenate the per-chromosome training
        fragments (written by TerminalCounter._dump_training_features) into the CSV."""
        def run() -> None:
            with open(csv_path, "w") as fh:
                merge_files(fragment_base, label, chr_ids, fh, header_lines=1,
                            remove_inputs=False, fragment_dir=fragment_dir)
            logger.info("%s features written to %s" % (name, csv_path))
        run_stage(store, marker_name("merge", name), run=run,
                  cleanup=lambda: remove_files(merge_file_list(fragment_base, label, chr_ids, fragment_dir)),
                  enabled=bool(csv_path), keep_tmp=self.args.keep_tmp)

    def _cleanup_per_chr_temp_files(self, sample, chr_ids: list):
        """Remove any remaining per-chromosome temporary files after merging.

        Each merge unit removes the fragments it consumed, but files from disabled output
        types (e.g. corrected_bed, read2transcripts not in --large_output) or from
        previous runs with different options may remain. This catch-all cleanup
        removes them.
        """
        for chr_id in chr_ids:
            chr_prefix = sample.get_chr_prefix(chr_id)
            pattern = os.path.join(glob.escape(sample.chr_fragment_dir), glob.escape(chr_prefix) + ".*")
            remove_files(glob.glob(pattern))

    def filter_umis(self, sample, store: CheckpointStore) -> None:
        """UMI deduplication: one stage per edit distance, then per barcode2barcode column."""
        keep_tmp = self.args.keep_tmp
        for i, edit_distance in enumerate(UMI_EDIT_DISTANCES[self.args.mode]):
            stage_name = marker_name("umi", "ED%d" % edit_distance)
            output_prefix = umi_output_prefix(sample.out_umi_filtered, edit_distance)
            run_stage(store, stage_name,
                      run=partial(self.filter_umis_round, sample, store, stage_name, edit_distance,
                                  output_prefix, i == 0, None, None,
                                  "Filtering PCR duplicates with edit distance %d" % edit_distance),
                      cleanup=partial(self.remove_umi_fragments, sample.out_umi_filtered_tmp, edit_distance),
                      keep_tmp=keep_tmp)

        # Barcode2barcode spot-based UMI dedup rounds
        if not getattr(self.args, 'barcode2barcode', None):
            return
        from .assignment.read_groups import parse_barcode2spot_spec, load_barcode2barcode_mapping
        filename, barcode_col, spot_cols = parse_barcode2spot_spec(self.args.barcode2barcode)
        mapping_cache = {}

        def barcode_remap(col_idx: int) -> Dict[str, str]:
            if "full" not in mapping_cache:
                mapping_cache["full"] = load_barcode2barcode_mapping(filename, barcode_col, spot_cols)
            return {bc: spots[col_idx] for bc, spots in mapping_cache["full"].items()
                    if col_idx < len(spots) and spots[col_idx]}

        for col_idx in range(len(spot_cols)):
            bc2bc_prefix = umi_barcode2barcode_prefix(sample.out_umi_filtered_tmp, col_idx)
            for edit_distance in UMI_EDIT_DISTANCES[self.args.mode]:
                stage_name = marker_name("umi_bc2bc", "col%d" % col_idx, "ED%d" % edit_distance)
                output_prefix = umi_output_prefix(
                    sample.out_umi_filtered + ".barcode_barcode_col%d" % col_idx, edit_distance)

                def run(col_idx=col_idx, edit_distance=edit_distance, stage_name=stage_name,
                        output_prefix=output_prefix, bc2bc_prefix=bc2bc_prefix) -> None:
                    self.filter_umis_round(sample, store, stage_name, edit_distance, output_prefix, False,
                                           barcode_remap(col_idx), bc2bc_prefix,
                                           "Filtering UMIs for barcode2barcode column %d with edit distance %d"
                                           % (col_idx, edit_distance))
                run_stage(store, stage_name, run=run,
                          cleanup=partial(self.remove_umi_fragments, bc2bc_prefix, edit_distance),
                          keep_tmp=keep_tmp)

    def filter_umis_round(self, sample, store: CheckpointStore, stage_name: str, edit_distance: int,
                          output_prefix: str, output_filtered_reads: bool, barcode_remap: Optional[Dict[str, str]],
                          tmp_prefix: Optional[str], message: str) -> None:
        logger.info(message)
        logger.info("Results will be saved to %s" % output_prefix)
        results = self.map_over_chromosomes(filter_umis_in_parallel, sample, self.get_chr_list(), self.args,
                                            edit_distance, store, stage_name, output_filtered_reads,
                                            barcode_remap, tmp_prefix)

        stat_dict = defaultdict(int)
        if large_output_enabled(self.args, "allinfo"):
            allinfo_fname = output_prefix + ".allinfo"
            if self.args.gzipped:
                allinfo_fname += ".gz"
            with open_text_write(allinfo_fname) as allinfo_outf:
                for all_info_file_name, _ in results:
                    with open(all_info_file_name, "r") as allinfo_inf:
                        shutil.copyfileobj(allinfo_inf, allinfo_outf)
        for _, stats_output_file_name in results:
            read_stats_tsv(stats_output_file_name, stat_dict)

        logger.info("PCR duplicates filtered with edit distance %d, filtering stats:" % edit_distance)
        with open(output_prefix + ".stats.tsv", "w") as outf:
            for k, v in stat_dict.items():
                logger.info("  %s: %d" % (k, v))
                outf.write("%s\t%d\n" % (k, v))

    def remove_umi_fragments(self, tmp_prefix: str, edit_distance: int) -> None:
        fragments = []
        for chr_id in self.get_chr_list():
            fragments.append(allinfo_file_name(tmp_prefix, chr_id, edit_distance))
            fragments.append(allinfo_stats_file_name(tmp_prefix, chr_id, edit_distance))
        remove_files(fragments)

    def map_over_chromosomes(self, worker, sample, *extra_args, chr_ids=None):
        """Run worker(sample, chr_id, *extra_args) for every chromosome, in parallel.

        Defaults to the chromosomes this run analysed; pass chr_ids to cover others.
        """
        gen = (worker, itertools.repeat(sample),
               self.get_chr_list() if chr_ids is None else chr_ids,
               *(itertools.repeat(a) for a in extra_args))
        if self.args.threads > 1:
            gc.collect()
            mp_context = multiprocessing.get_context('fork')
            log_file, log_level = _get_log_params()
            with ProcessPoolExecutor(max_workers=self.args.threads, mp_context=mp_context,
                                     initializer=setup_worker_logging,
                                     initargs=(log_file, log_level)) as proc:
                return list(proc.map(*gen, chunksize=1))
        return list(map(*gen))

    def write_tagged_bam(self, sample, references):
        """A copy of the input BAM(s) with barcode and UMI tags, keeping every alignment.

        Purely a side output: the barcode table split it reads from happens regardless, and
        nothing downstream looks at the result. references is the list the barcode table was
        split over, so each fragment has a table to read its tags from.
        """
        logger.info("Writing tagged BAM")
        bam_files = [f[0] for f in sample.file_list]
        fragments = self.map_over_chromosomes(write_tagged_bam_in_parallel, sample, self.args,
                                              chr_ids=references)

        # fetch() by chromosome never returns unmapped reads, so they need a pass of their own.
        # They belong to no chromosome and so appear in no split table, but a barcode is called
        # from the read sequence and does not need an alignment -- look theirs up directly.
        unmapped_fragment = tagged_bam_fragment_name(sample.out_raw_file, "unmapped")
        unmapped_ids = collect_unmapped_read_ids(bam_files)
        unmapped_tags = {}
        for table in sample.barcoded_reads:
            unmapped_tags.update(load_barcode_umi_tags(resolve_optionally_gzipped(table),
                                                       read_ids=unmapped_ids))
        if write_unmapped_bam(bam_files, unmapped_fragment, unmapped_tags,
                              self.args.barcode_tag, self.args.umi_tag):
            # a table row exists for every read, but its barcode may be the uncalled placeholder
            barcoded = sum(1 for tags in unmapped_tags.values() if tags[0] not in PLACEHOLDERS)
            logger.info("Copied %d unmapped reads, %d of them barcoded"
                        % (len(unmapped_ids), barcoded))
            fragments.append(unmapped_fragment)
        else:
            os.remove(unmapped_fragment)

        merged = merge_bam_files(sample.out_tagged_bam, fragments, self.args.threads)
        if merged is None:
            logger.warning("No alignments to tag, %s was not written" % sample.out_tagged_bam)
        else:
            logger.info("Tagged BAM saved to %s" % merged)
        return merged

    def write_deduplicated_bam(self, sample):
        """A primary-only BAM of the reads that survived UMI filtering, carrying their tags."""
        logger.info("Writing deduplicated BAM")
        fragments = self.map_over_chromosomes(write_deduplicated_bam_in_parallel, sample, self.args)
        merged = merge_bam_files(sample.out_deduplicated_bam, fragments, self.args.threads)
        if merged is None:
            logger.warning("No reads survived UMI filtering, %s was not written"
                           % sample.out_deduplicated_bam)
        else:
            logger.info("Deduplicated BAM saved to %s" % merged)
        return merged

    def compress_barcode_tables(self, sample):
        """Gzip the barcoded read tables once nothing needs them any more.

        They stay plain for the whole run because split_read_barcode_table reads them, and
        only tables this run produced are touched -- files passed via --barcoded_reads belong
        to the user and are left alone.
        """
        if not self.args.gzipped:
            return
        produced = sorted(glob.glob(sample.barcodes_tsv + "_*.tsv"))
        for table in produced:
            gzipped = gzip_file_in_place(table)
            logger.info("Compressed barcode table to %s" % gzipped)

    def split_external_polya(self, sample) -> None:
        """Normalize and split the --polya_trimmed list:/flnc: table per chromosome."""
        split_polya_dict = {chr_id: sample.get_polya_split_file(chr_id) for chr_id in self.get_chr_list()}
        logger.info("Splitting external polyA table")
        normalize_polya_reads(self.args.polya_trimmed_file, self.args.polya_trimmed.name,
                              sample.polya_reads_normalized)
        split_read_table_parallel(sample, [sample.polya_reads_normalized], split_polya_dict, self.args.threads,
                                  read_column=0, group_columns=(1,), delim='\t')
        logger.info("External polyA table was split")

    def split_read_barcode_table(self, sample, split_barcodes_file_names):
        logger.info("Splitting read barcode table")
        # Supports both IsoQuant's 6-column format and third-party 3-column format
        # Only columns 0, 1, 2 (read_id, barcode, umi) are required and preserved
        barcode_tables = [resolve_optionally_gzipped(f) for f in sample.barcoded_reads]
        split_read_table_parallel(sample, barcode_tables, split_barcodes_file_names, self.args.threads,
                                  read_column=0, group_columns=(1, 2), delim='\t')
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
