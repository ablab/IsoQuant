############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import glob
import itertools
import logging
import shutil
from enum import Enum, unique
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import gffutils
import pysam
from pyfaidx import Fasta

from .modes import IsoQuantMode
from .common import proper_plural_form
from .serialization import *
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
from .read_groups import (
    create_read_grouper,
    prepare_read_groups,
    load_table,
    get_grouping_strategy_names
)
from .assignment_io import (
    IOSupport,
    BEDPrinter,
    ReadAssignmentCompositePrinter,
    SqantiTSVPrinter,
    BasicTSVAssignmentPrinter,
    TmpFileAssignmentPrinter,
)
from .processed_read_manager import ProcessedReadsManagerHighMemory, ProcessedReadsManagerNoSecondary, ProcessedReadsManagerNormalMemory
from .id_policy import SimpleIDDistributor, ExcludingIdDistributor, FeatureIdStorage
from .file_naming import *
from .transcript_printer import GFFPrinter, VoidTranscriptPrinter, create_extended_storage
from .graph_based_model_construction import GraphBasedModelConstructor
from .gene_info import TranscriptModelType, get_all_chromosome_genes, get_all_chromosome_transcripts
from .assignment_loader import create_assignment_loader, BasicReadAssignmentLoader, load_genedb
from .barcode_calling.umi_filtering import create_transcript_info_dict, UMIFilter
from .table_splitter import split_read_table_parallel
from .string_pools import StringPoolManager

logger = logging.getLogger('IsoQuant')


def setup_string_pools(args, sample, chr_ids, chr_id=None, gffutils_db=None,
                       load_barcode_pool=False, load_tsv_pools=False, read_group_file_prefix=None):
    """
    Set up string pools for memory optimization during parallel processing.

    Args:
        args: Command-line arguments
        sample: SampleData object
        chr_ids: List of all chromosome IDs to process (for chromosome pool)
        chr_id: Current chromosome ID (required if loading per-chromosome pools)
        gffutils_db: Pre-loaded gffutils database (if None, will load from args.genedb)
        load_barcode_pool: Whether to load per-chromosome barcode pool
        load_tsv_pools: Whether to load per-chromosome TSV read group pools
        read_group_file_prefix: Override for read_group_file path (used when sample has modified prefix)

    Returns:
        StringPoolManager with pools configured
    """
    from .read_groups import get_grouping_pool_types

    string_pools = StringPoolManager()

    # Build chromosome pool from the definitive list of chromosomes to process
    string_pools.build_chromosome_pool(chr_ids)

    # Build gene/transcript pools from annotation (if available)
    if gffutils_db is None and args.genedb:
        gffutils_db = load_genedb(args.genedb)
    if gffutils_db:
        string_pools.build_from_gffutils(gffutils_db)

    # Set up read group pool type mapping
    pool_types = get_grouping_pool_types(args)
    for spec_idx, pool_type in pool_types.items():
        string_pools.set_group_spec_pool_type(spec_idx, pool_type)

    # Build file_name pool if needed
    if any(pt == 'file_name' for pt in pool_types.values()):
        string_pools.build_file_name_pool(sample)

    # Build barcode_spot pool if needed
    if any(pt == 'barcode_spot' for pt in pool_types.values()):
        if args.barcode2spot:
            string_pools.build_barcode_spot_pool(args.barcode2spot)

    # Load per-chromosome pools if requested
    if chr_id:
        if load_barcode_pool and sample.barcodes_split_reads:
            barcode_file = sample.barcodes_split_reads + "_" + chr_id
            string_pools.load_barcode_pool(barcode_file)

        # Use override if provided, otherwise use sample's read_group_file
        read_group_file = read_group_file_prefix if read_group_file_prefix else sample.read_group_file
        if load_tsv_pools and read_group_file:
            import re
            # Build mapping from tsv spec_index to (col_index, delimiter) from pool_types
            # pool_type format: 'tsv:SPEC_INDEX:COL_INDEX:DELIMITER'
            tsv_pool_info = {}  # spec_index -> list of (col_index, delimiter, pool_key)
            for pool_type in pool_types.values():
                if pool_type.startswith('tsv:'):
                    parts = pool_type.split(':')
                    if len(parts) >= 4:
                        tsv_spec_idx = int(parts[1])
                        col_idx = int(parts[2])
                        delimiter = parts[3]
                        # pool_key is spec_idx:col_idx to uniquely identify multi-column pools
                        pool_key = f"{tsv_spec_idx}:{col_idx}"
                        if tsv_spec_idx not in tsv_pool_info:
                            tsv_pool_info[tsv_spec_idx] = []
                        tsv_pool_info[tsv_spec_idx].append((col_idx, delimiter, pool_key))

            base_pattern = read_group_file + "_spec*_" + chr_id
            spec_files = sorted(glob.glob(base_pattern))
            for spec_file in spec_files:
                match = re.search(r'_spec(\d+)_', spec_file)
                if match:
                    spec_index = int(match.group(1))
                    # Load pool for each column in this spec file
                    for col_idx, delimiter, pool_key in tsv_pool_info.get(spec_index, [(1, '\t', f"{spec_index}:1")]):
                        string_pools.load_read_group_tsv_pool(spec_file, pool_key, col_idx, delimiter)

    return string_pools


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


def collect_reads_in_parallel(sample, chr_id, chr_ids, args, processed_read_manager_type):
    current_chr_record = Fasta(args.reference, indexname=args.fai_file_name)[chr_id]
    if args.high_memory:
        current_chr_record = str(current_chr_record)
    read_grouper = create_read_grouper(args, sample, chr_id)
    lock_file = reads_collected_lock_file_name(sample.out_raw_file, chr_id)
    save_file = saves_file_name(sample.out_raw_file, chr_id)
    group_file = save_file + "_groups"
    bamstat_file = save_file + "_bamstat"
    processed_reads_manager = processed_read_manager_type(sample, args.multimap_strategy, chr_ids, args.genedb)

    if os.path.exists(lock_file) and args.resume:
        logger.info("Detected processed reads for " + chr_id)
        if os.path.exists(group_file) and os.path.exists(save_file):
            # Load read_groups from file
            with open(group_file) as f:
                lines = [line.strip() for line in f]
                num_strategies = int(lines[0])
                if isinstance(read_grouper.read_groups, list):
                    # MultiReadGrouper: list of sets
                    read_grouper.read_groups = []
                    for i in range(num_strategies):
                        groups = set(lines[i + 1].split(";")) if lines[i + 1] else set()
                        read_grouper.read_groups.append(groups)
                else:
                    # Single grouper: single set
                    groups_str = lines[1] if len(lines) > 1 else ""
                    read_grouper.read_groups = set(groups_str.split(";")) if groups_str else set()
            alignment_stat_counter = EnumStats(bamstat_file)

            # Build string pools for loading serialized data
            gffutils_db = gffutils.FeatureDB(args.genedb) if args.genedb else None
            string_pools = setup_string_pools(args, sample, chr_ids, chr_id, gffutils_db,
                                              load_barcode_pool=True, load_tsv_pools=True)

            # Load dynamic pools (for read groups from BAM tags/read IDs)
            dynamic_pools_file = dynamic_pools_file_name(sample.out_raw_file, chr_id)
            if os.path.exists(dynamic_pools_file):
                logger.debug(f"Loading dynamic pools from {dynamic_pools_file}")
                with open(dynamic_pools_file, 'rb') as f:
                    string_pools.deserialize_dynamic_pools(f)

            loader = BasicReadAssignmentLoader(save_file, string_pools)
            while loader.has_next():
                for read_assignment in loader.get_next():
                    if read_assignment is None: continue
                    processed_reads_manager.load_read(read_assignment)
            logger.info("Loaded data for " + chr_id)
            return chr_id, read_grouper.read_groups, alignment_stat_counter, processed_reads_manager
        else:
            logger.warning("Something is wrong with save files for %s, will process from scratch " % chr_id)
            if not os.path.exists(group_file):
                logger.warning("%s does not exist" % group_file)
            if not os.path.exists(save_file):
                logger.warning("%s does not exist" % save_file)
            os.remove(lock_file)

    if os.path.exists(lock_file):
        os.remove(lock_file)

    tmp_printer = TmpFileAssignmentPrinter(save_file, args)
    bam_files = list(map(lambda x: x[0], sample.file_list))
    bam_file_pairs = [(pysam.AlignmentFile(bam, "rb", require_index=True), bam) for bam in bam_files]
    gffutils_db = gffutils.FeatureDB(args.genedb) if args.genedb else None
    illumina_bam = sample.illumina_bam

    # Build string pools for memory optimization
    string_pools = setup_string_pools(args, sample, chr_ids, chr_id, gffutils_db,
                                      load_barcode_pool=True, load_tsv_pools=True)

    # Load barcode dict for this chromosome if available (for backward compatibility during transition)
    barcode_dict = {}
    if sample.barcodes_split_reads:
        barcode_file = sample.barcodes_split_reads + "_" + chr_id
        if os.path.exists(barcode_file):
            logger.debug(f"Loading barcodes from {barcode_file}")
            for line in open(barcode_file):
                if line.startswith("#"):
                    continue
                # Use split('\t') to preserve empty fields (empty UMI = consecutive tabs)
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    # Empty strings become None (not added to pool)
                    barcode = parts[1] if parts[1] else None
                    umi = parts[2] if parts[2] else None
                    barcode_dict[parts[0]] = (barcode, umi)
            logger.debug("Loaded %d barcodes" % len(barcode_dict))

    logger.info("Processing chromosome " + chr_id)
    alignment_collector = \
        AlignmentCollector(chr_id, bam_file_pairs, args, illumina_bam, gffutils_db, current_chr_record, read_grouper,
                           barcode_dict, args.max_coverage_small_chr, args.max_coverage_normal_chr, string_pools)

    for gene_info, assignment_storage in alignment_collector.process():
        tmp_printer.add_gene_info(gene_info)
        for read_assignment in assignment_storage:
            tmp_printer.add_read_info(read_assignment)
            processed_reads_manager.add_read(read_assignment)
    with open(group_file, "w") as group_dump:
        # Save read_groups (can be set or list of sets)
        if isinstance(read_grouper.read_groups, list):
            # MultiReadGrouper: list of sets
            group_dump.write("%d\n" % len(read_grouper.read_groups))  # Number of strategies
            for group_set in read_grouper.read_groups:
                group_dump.write("%s\n" % ";".join(sorted(group_set)))  # Semicolon-separated groups
        else:
            # Single grouper: single set
            group_dump.write("1\n")  # One strategy
            group_dump.write("%s\n" % ";".join(sorted(read_grouper.read_groups)))
    alignment_collector.alignment_stat_counter.dump(bamstat_file)

    for bam in bam_file_pairs:
        bam[0].close()

    tmp_printer.close()

    # Save dynamic pools if they have data (for read groups from BAM tags/read IDs)
    if string_pools.has_dynamic_pools():
        dynamic_pools_file = dynamic_pools_file_name(sample.out_raw_file, chr_id)
        logger.debug(f"Saving dynamic pools to {dynamic_pools_file}")
        with open(dynamic_pools_file, 'wb') as f:
            string_pools.serialize_dynamic_pools(f)

    processed_reads_manager.finalize(chr_id)
    logger.info("Finished processing chromosome " + chr_id)
    open(lock_file, "w").close()

    return chr_id, read_grouper.read_groups, alignment_collector.alignment_stat_counter, processed_reads_manager


def construct_models_in_parallel(sample, chr_id, chr_ids, saves_prefix, args, read_groups):
    logger.info("Processing chromosome " + chr_id)
    use_filtered_reads = args.mode.needs_pcr_deduplication()

    # Derive read_group_file prefix from saves_prefix (replace .save with .read_group)
    read_group_file_prefix = saves_prefix.replace('.save', '.read_group') if saves_prefix.endswith('.save') else None

    # Build string pools for memory optimization
    string_pools = setup_string_pools(args, sample, chr_ids, chr_id,
                                      load_barcode_pool=False, load_tsv_pools=True,
                                      read_group_file_prefix=read_group_file_prefix)

    # Load dynamic pools (for read groups from BAM tags/read IDs)
    dynamic_pools_file = dynamic_pools_file_name(saves_prefix, chr_id)
    if os.path.exists(dynamic_pools_file):
        logger.debug(f"Loading dynamic pools from {dynamic_pools_file}")
        with open(dynamic_pools_file, 'rb') as f:
            string_pools.deserialize_dynamic_pools(f)

    loader = create_assignment_loader(chr_id, saves_prefix, args.genedb, args.reference, args.fai_file_name, string_pools, use_filtered_reads)

    chr_dump_file = saves_file_name(saves_prefix, chr_id)
    lock_file = reads_processed_lock_file_name(saves_prefix, chr_id)
    read_stat_file = "{}_read_stat".format(chr_dump_file)
    transcript_stat_file = "{}_transcript_stat".format(chr_dump_file)
    construct_models = not args.no_model_construction

    if os.path.exists(lock_file):
        if args.resume:
            logger.info("Processed assignments from chromosome " + chr_id + " detected")
            read_stat = EnumStats(read_stat_file)
            transcript_stat = EnumStats(transcript_stat_file) if construct_models else EnumStats()
            return read_stat, transcript_stat
        os.remove(lock_file)

    grouping_strategy_names = get_grouping_strategy_names(args)
    aggregator = ReadAssignmentAggregator(args, sample, read_groups, loader.genedb, chr_id,
                                         grouping_strategy_names=grouping_strategy_names)

    transcript_stat_counter = EnumStats()
    io_support = IOSupport(args)
    transcript_id_distributor = ExcludingIdDistributor(loader.genedb, chr_id)
    exon_id_storage = FeatureIdStorage(SimpleIDDistributor(), loader.genedb, chr_id, "exon")

    if construct_models:
        tmp_gff_printer = GFFPrinter(sample.out_dir, sample.prefix, exon_id_storage,
                                     output_r2t=not args.no_large_files,
                                     check_canonical=args.check_canonical)
    else:
        tmp_gff_printer = VoidTranscriptPrinter()
    if construct_models and args.genedb:
        tmp_extended_gff_printer = GFFPrinter(sample.out_dir, sample.prefix, exon_id_storage,
                                              gtf_suffix=".extended_annotation.gtf",
                                              output_r2t=False, check_canonical=args.check_canonical)
    else:
        tmp_extended_gff_printer = VoidTranscriptPrinter()

    sqanti_t2t_printer = SqantiTSVPrinter(sample.out_t2t_tsv, args, IOSupport(args)) \
        if args.sqanti_output else VoidTranscriptPrinter()
    novel_model_storage = []

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
            strategy_names = aggregator.grouping_strategy_names if hasattr(aggregator, 'grouping_strategy_names') else []
            model_constructor = GraphBasedModelConstructor(gene_info, loader.chr_record, args,
                                                           aggregator.transcript_model_global_counter,
                                                           aggregator.gene_model_global_counter,
                                                           transcript_id_distributor,
                                                           grouping_strategy_names=strategy_names,
                                                           use_technical_replicas=sample.use_technical_replicas,
                                                           string_pools=string_pools)
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
        if loader.genedb:
            all_models, gene_info = create_extended_storage(loader.genedb, chr_id, loader.chr_record, novel_model_storage)
            if args.check_canonical:
                io_support.add_canonical_info(all_models, gene_info)
            tmp_extended_gff_printer.dump(gene_info, all_models)
        aggregator.transcript_model_global_counter.dump()
        aggregator.gene_model_global_counter.dump()
        transcript_stat_counter.dump(transcript_stat_file)
    logger.info("Finished processing chromosome " + chr_id)
    open(lock_file, "w").close()

    return aggregator.read_stat_counter, transcript_stat_counter


def filter_umis_in_parallel(sample, chr_id, chr_ids, args, edit_distance, output_filtered_reads=False):
    transcript_type_dict = create_transcript_info_dict(args.genedb, [chr_id])
    umi_filtered_done = umi_filtered_lock_file_name(sample.out_umi_filtered_done, chr_id, edit_distance)
    all_info_file_name = allinfo_file_name(sample.out_umi_filtered_done, chr_id, edit_distance)
    stats_output_file_name = allinfo_stats_file_name(sample.out_umi_filtered_done, chr_id, edit_distance)

    if os.path.exists(umi_filtered_done):
        if args.resume:
            return all_info_file_name, stats_output_file_name, umi_filtered_done
        os.remove(umi_filtered_done)

    logger.info("Filtering PCR duplicates for chromosome " + chr_id)
    barcode_feature_table = {}
    if args.barcode2spot:
        for barcode2spot_file in args.barcode2spot:
            barcode_feature_table.update(load_table(barcode2spot_file, 0, 1, '\t'))

    # Build string pools for memory optimization
    # Barcode pool is critical: ReadAssignments were serialized with barcode_id values
    # that reference the barcode pool from read collection
    string_pools = setup_string_pools(args, sample, chr_ids, chr_id,
                                      load_barcode_pool=True, load_tsv_pools=False)

    # Load dynamic pools (for read groups from BAM tags/read IDs)
    dynamic_pools_file = dynamic_pools_file_name(sample.out_raw_file, chr_id)
    if os.path.exists(dynamic_pools_file):
        logger.debug(f"Loading dynamic pools from {dynamic_pools_file}")
        with open(dynamic_pools_file, 'rb') as f:
            string_pools.deserialize_dynamic_pools(f)

    umi_filter = UMIFilter(args.umi_length, edit_distance)
    filtered_reads = filtered_reads_file_name(sample.out_raw_file, chr_id) if output_filtered_reads else None
    umi_filter.process_single_chr(chr_id, sample.out_raw_file,
                                  transcript_type_dict,
                                  barcode_feature_table,
                                  all_info_file_name,
                                  filtered_reads,
                                  stats_output_file_name,
                                  string_pools)
    open(umi_filtered_done, "w").close()
    logger.info("PCR duplicates filtered for chromosome " + chr_id)

    return all_info_file_name, stats_output_file_name, umi_filtered_done


class ReadAssignmentAggregator:
    def __init__(self, args, sample, read_groups, gffutils_db=None, chr_id=None, gzipped=False, grouping_strategy_names=None):
        self.args = args
        self.read_groups = read_groups
        self.grouping_strategy_names = grouping_strategy_names if grouping_strategy_names else ["default"]
        self.common_header = "# Command line: " + args._cmd_line + "\n# IsoQuant version: " + args._version + "\n"
        self.io_support = IOSupport(self.args)

        self.gene_set = set()
        self.transcript_set = set()
        if gffutils_db and chr_id:
            self.gene_set = set(get_all_chromosome_genes(gffutils_db, chr_id))
            self.transcript_set = set(get_all_chromosome_transcripts(gffutils_db, chr_id))

        self.read_stat_counter = EnumStats()

        printer_list = []
        self.corrected_bed_printer = None
        if not self.args.no_large_files:
            self.corrected_bed_printer = BEDPrinter(sample.out_corrected_bed,
                                                    self.args,
                                                    print_corrected=True,
                                                    gzipped=gzipped)
            printer_list.append(self.corrected_bed_printer)
        self.basic_printer = None
        if self.args.genedb and not self.args.no_large_files:
            self.basic_printer = BasicTSVAssignmentPrinter(sample.out_assigned_tsv, self.args, self.io_support,
                                                           additional_header=self.common_header, gzipped=gzipped)
            sample.out_assigned_tsv_result = self.basic_printer.output_file_name
            printer_list.append(self.basic_printer)
        self.t2t_sqanti_printer = VoidTranscriptPrinter()
        if self.args.sqanti_output:
            self.t2t_sqanti_printer = SqantiTSVPrinter(sample.out_t2t_tsv, self.args, self.io_support)
        self.global_printer = ReadAssignmentCompositePrinter(printer_list)

        self.global_counter = CompositeCounter()
        if self.args.genedb:
            self.gene_counter = create_gene_counter(sample.out_gene_counts_tsv,
                                                    self.args.gene_quantification,
                                                    complete_feature_list=self.gene_set)
            self.transcript_counter = create_transcript_counter(sample.out_transcript_counts_tsv,
                                                                self.args.transcript_quantification,
                                                                complete_feature_list=self.transcript_set)
            self.global_counter.add_counters([self.gene_counter, self.transcript_counter])

        self.transcript_model_global_counter = CompositeCounter()
        self.gene_model_global_counter = CompositeCounter()
        if not self.args.no_model_construction:
            self.transcript_model_counter = create_transcript_counter(sample.out_transcript_model_counts_tsv,
                                                                      self.args.transcript_quantification)
            self.gene_model_counter = create_gene_counter(sample.out_gene_model_counts_tsv,
                                                          self.args.gene_quantification)

            self.transcript_model_global_counter.add_counter(self.transcript_model_counter)
            self.gene_model_global_counter.add_counter(self.gene_model_counter)

        if self.args.count_exons and self.args.genedb:
            self.exon_counter = ExonCounter(sample.out_exon_counts_tsv, ignore_read_groups=True)
            self.intron_counter = IntronCounter(sample.out_intron_counts_tsv, ignore_read_groups=True)
            self.global_counter.add_counters([self.exon_counter, self.intron_counter])

        if self.args.read_group and self.args.genedb:
            for group_idx, strategy_name in enumerate(self.grouping_strategy_names):
                # Add strategy name as suffix to output file
                gene_out_file = f"{sample.out_gene_grouped_counts_tsv}_{strategy_name}"
                transcript_out_file = f"{sample.out_transcript_grouped_counts_tsv}_{strategy_name}"

                gene_counter = create_gene_counter(gene_out_file,
                                                   self.args.gene_quantification,
                                                   complete_feature_list=self.gene_set,
                                                   read_groups=self.read_groups[group_idx],
                                                   group_index=group_idx)
                transcript_counter = create_transcript_counter(transcript_out_file,
                                                              self.args.transcript_quantification,
                                                              complete_feature_list=self.transcript_set,
                                                              read_groups=self.read_groups[group_idx],
                                                              group_index=group_idx)

                self.global_counter.add_counters([gene_counter, transcript_counter])

                if self.args.count_exons:
                    exon_out_file = f"{sample.out_exon_grouped_counts_tsv}_{strategy_name}"
                    intron_out_file = f"{sample.out_intron_grouped_counts_tsv}_{strategy_name}"
                    exon_counter = ExonCounter(exon_out_file, group_index=group_idx)
                    intron_counter = IntronCounter(intron_out_file, group_index=group_idx)
                    self.global_counter.add_counters([exon_counter, intron_counter])

        if self.args.read_group and not self.args.no_model_construction:
            for group_idx, strategy_name in enumerate(self.grouping_strategy_names):
                transcript_model_out_file = f"{sample.out_transcript_model_grouped_counts_tsv}_{strategy_name}"
                gene_model_out_file = f"{sample.out_gene_model_grouped_counts_tsv}_{strategy_name}"

                transcript_model_counter = create_transcript_counter(
                    transcript_model_out_file,
                    self.args.transcript_quantification,
                    read_groups=self.read_groups[group_idx],
                    group_index=group_idx)
                gene_model_counter = create_gene_counter(
                    gene_model_out_file,
                    self.args.gene_quantification,
                    read_groups=self.read_groups[group_idx],
                    group_index=group_idx)

                self.transcript_model_global_counter.add_counter(transcript_model_counter)
                self.gene_model_global_counter.add_counter(gene_model_counter)


# Class for processing all samples against gene database
class DatasetProcessor:
    def __init__(self, args):
        self.args = args
        self.input_data = args.input_data
        self.args.gunzipped_reference = None
        self.common_header = "# Command line: " + args._cmd_line + "\n# IsoQuant version: " + args._version + "\n"
        self.io_support = IOSupport(self.args)
        self.all_read_groups = []  # Will be initialized per sample as list of sets
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
        sample.use_technical_replicas = (len(sample.file_list) > 1 and
                                         self.args.read_group is not None and
                                         "file_name" in self.args.read_group)

        # Initialize all_read_groups as list of sets (one per grouping strategy)
        grouping_strategy_names = get_grouping_strategy_names(self.args)
        self.all_read_groups = [set() for _ in grouping_strategy_names]
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
                exit(-1)
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
                exit(-1)
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
        info_file = sample.out_raw_file + "_info"
        lock_file = sample.out_raw_file + "_lock"

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
        grouping_strategy_names = get_grouping_strategy_names(self.args)
        all_read_groups = [set() for _ in grouping_strategy_names]

        if self.args.threads > 1:
            with ProcessPoolExecutor(max_workers=self.args.threads) as proc:
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
        grouping_strategy_names = get_grouping_strategy_names(self.args)
        aggregator = ReadAssignmentAggregator(self.args, sample, self.all_read_groups, gzipped=self.args.gzipped,
                                             grouping_strategy_names=grouping_strategy_names)
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
                output_r2t=not self.args.no_large_files,
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
            (SampleData(sample.file_list, f"{sample.prefix}_{chr_id}", sample.out_dir, sample.readable_names_dict, sample.illumina_bam, sample.barcoded_reads) for chr_id in chr_ids),
            chr_ids,
            itertools.repeat(chr_ids),
            itertools.repeat(dump_filename),
            itertools.repeat(self.args),
            itertools.repeat(self.all_read_groups),
        )

        if self.args.threads > 1:
            with ProcessPoolExecutor(max_workers=self.args.threads) as proc:
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

            output_prefix = sample.out_umi_filtered + (".ALL" if edit_distance < 0 else ".ED%d" % edit_distance)
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
                with ProcessPoolExecutor(max_workers=self.args.threads) as proc:
                    results = proc.map(*umi_gen, chunksize=1)
            else:
                results = map(*umi_gen)

            stat_dict = defaultdict(int)
            files_to_remove = []
            with open(output_prefix + ".allinfo", "w") as outf:
                for all_info_file_name, stats_output_file_name, umi_filter_done in results:
                    shutil.copyfileobj(open(all_info_file_name, "r"), outf)
                    for l in open(stats_output_file_name, "r"):
                        v = l.strip().split("\t")
                        if len(v) != 2:
                            continue
                        stat_dict[v[0]] += int(v[1])
                    files_to_remove.append(all_info_file_name)
                    files_to_remove.append(stats_output_file_name)
                    files_to_remove.append(umi_filter_done)

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
        info_loader = open(dump_filename + "_info", "rb")
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
