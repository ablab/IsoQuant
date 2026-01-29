############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import os

import gffutils
import pysam
from pyfaidx import Fasta

from .stats import EnumStats
from .alignment_processor import AlignmentCollector
from .assignment_io import IOSupport, TmpFileAssignmentPrinter, SqantiTSVPrinter
from .assignment_loader import create_assignment_loader, BasicReadAssignmentLoader
from .barcode_calling.umi_filtering import create_transcript_info_dict, UMIFilter
from .file_naming import (
    read_groups_file_name,
    bamstat_file_name,
    read_group_file_from_saves,
    dynamic_pools_file_name,
    saves_file_name,
    reads_processed_lock_file_name,
    read_stat_file_name,
    transcript_stat_file_name,
    umi_filtered_lock_file_name,
    allinfo_file_name,
    allinfo_stats_file_name,
)
from .gene_info import TranscriptModelType
from .graph_based_model_construction import GraphBasedModelConstructor
from .id_policy import SimpleIDDistributor, ExcludingIdDistributor, FeatureIdStorage
from .read_groups import create_read_grouper, get_grouping_strategy_names, load_table
from .transcript_printer import GFFPrinter, VoidTranscriptPrinter, create_extended_storage
from .assignment_aggregator import ReadAssignmentAggregator
from .string_pools import setup_string_pools
from .common import large_output_enabled

logger = logging.getLogger('IsoQuant')


# Helper functions to reduce code duplication

def load_dynamic_pools(string_pools, dynamic_pools_file):
    """Load dynamic pools from file if it exists."""
    if os.path.exists(dynamic_pools_file):
        logger.debug(f"Loading dynamic pools from {dynamic_pools_file}")
        with open(dynamic_pools_file, 'rb') as f:
            string_pools.deserialize_dynamic_pools(f)


def save_dynamic_pools(string_pools, dynamic_pools_file):
    """Save dynamic pools to file if they have data."""
    if string_pools.has_dynamic_pools():
        logger.debug(f"Saving dynamic pools to {dynamic_pools_file}")
        with open(dynamic_pools_file, 'wb') as f:
            string_pools.serialize_dynamic_pools(f)


def save_read_groups(group_file, read_groups):
    """Save read_groups to file (handles both single set and list of sets)."""
    with open(group_file, "w") as group_dump:
        if isinstance(read_groups, list):
            # MultiReadGrouper: list of sets
            group_dump.write("%d\n" % len(read_groups))
            for group_set in read_groups:
                group_dump.write("%s\n" % ";".join(sorted(group_set)))
        else:
            # Single grouper: single set
            group_dump.write("1\n")
            group_dump.write("%s\n" % ";".join(sorted(read_groups)))


def load_read_groups(group_file, read_grouper):
    """Load read_groups from file into read_grouper."""
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


def load_barcode_dict(sample, chr_id):
    """Load barcode dictionary for a chromosome."""
    barcode_dict = {}
    if sample.barcodes_split_reads:
        barcode_file = sample.get_barcodes_split_file(chr_id)
        if os.path.exists(barcode_file):
            logger.debug(f"Loading barcodes from {barcode_file}")
            for line in open(barcode_file):
                if line.startswith("#"):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    barcode = parts[1] if parts[1] else None
                    umi = parts[2] if parts[2] else None
                    barcode_dict[parts[0]] = (barcode, umi)
            logger.debug("Loaded %d barcodes" % len(barcode_dict))
    return barcode_dict


def collect_reads_in_parallel(sample, chr_id, chr_ids, args, processed_read_manager_type):
    current_chr_record = Fasta(args.reference, indexname=args.fai_file_name)[chr_id]
    if args.high_memory:
        current_chr_record = str(current_chr_record)
    read_grouper = create_read_grouper(args, sample, chr_id)
    lock_file = sample.get_collected_lock_file(chr_id)
    save_file = sample.get_save_file(chr_id)
    group_file = read_groups_file_name(save_file)
    bamstat_file = bamstat_file_name(save_file)
    processed_reads_manager = processed_read_manager_type(sample, args.multimap_strategy, chr_ids, args.genedb)

    if os.path.exists(lock_file) and args.resume:
        logger.info("Detected processed reads for " + chr_id)
        if os.path.exists(group_file) and os.path.exists(save_file):
            load_read_groups(group_file, read_grouper)
            alignment_stat_counter = EnumStats(bamstat_file)

            # Build string pools for loading serialized data
            gffutils_db = gffutils.FeatureDB(args.genedb) if args.genedb else None
            string_pools = setup_string_pools(args, sample, chr_ids, chr_id, gffutils_db,
                                              load_barcode_pool=True, load_tsv_pools=True)

            # Load dynamic pools (for read groups from BAM tags/read IDs)
            load_dynamic_pools(string_pools, sample.get_dynamic_pools_file(chr_id))

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

    # Load barcode dict for this chromosome if available
    barcode_dict = load_barcode_dict(sample, chr_id)

    logger.info("Processing chromosome " + chr_id)
    alignment_collector = \
        AlignmentCollector(chr_id, bam_file_pairs, args, illumina_bam, gffutils_db, current_chr_record, read_grouper,
                           barcode_dict, args.max_coverage_small_chr, args.max_coverage_normal_chr, string_pools)

    for gene_info, assignment_storage in alignment_collector.process():
        tmp_printer.add_gene_info(gene_info)
        for read_assignment in assignment_storage:
            tmp_printer.add_read_info(read_assignment)
            processed_reads_manager.add_read(read_assignment)
    save_read_groups(group_file, read_grouper.read_groups)
    alignment_collector.alignment_stat_counter.dump(bamstat_file)

    for bam in bam_file_pairs:
        bam[0].close()

    tmp_printer.close()

    # Save dynamic pools if they have data (for read groups from BAM tags/read IDs)
    save_dynamic_pools(string_pools, sample.get_dynamic_pools_file(chr_id))

    processed_reads_manager.finalize(chr_id)
    logger.info("Finished processing chromosome " + chr_id)
    open(lock_file, "w").close()

    return chr_id, read_grouper.read_groups, alignment_collector.alignment_stat_counter, processed_reads_manager


def construct_models_in_parallel(sample, chr_id, chr_ids, saves_prefix, args, read_groups):
    from .read_groups import get_grouping_pool_types
    logger.info("Processing chromosome " + chr_id)
    use_filtered_reads = args.mode.needs_pcr_deduplication()

    # Derive read_group_file prefix from saves_prefix
    read_group_file_prefix = read_group_file_from_saves(saves_prefix)

    # Check if barcode pool is needed for any grouper
    pool_types = get_grouping_pool_types(args)
    needs_barcode_pool = any(pt == 'barcode' for pt in pool_types.values())

    # Build string pools for memory optimization
    string_pools = setup_string_pools(args, sample, chr_ids, chr_id,
                                      load_barcode_pool=needs_barcode_pool, load_tsv_pools=True,
                                      read_group_file_prefix=read_group_file_prefix)

    # Load dynamic pools (for read groups from BAM tags/read IDs)
    load_dynamic_pools(string_pools, dynamic_pools_file_name(saves_prefix, chr_id))

    loader = create_assignment_loader(chr_id, saves_prefix, args.genedb, args.reference, args.fai_file_name, string_pools, use_filtered_reads)

    chr_dump_file = saves_file_name(saves_prefix, chr_id)
    lock_file = reads_processed_lock_file_name(saves_prefix, chr_id)
    chr_read_stat_file = read_stat_file_name(chr_dump_file)
    chr_transcript_stat_file = transcript_stat_file_name(chr_dump_file)
    construct_models = not args.no_model_construction

    if os.path.exists(lock_file):
        if args.resume:
            logger.info("Processed assignments from chromosome " + chr_id + " detected")
            read_stat = EnumStats(chr_read_stat_file)
            transcript_stat = EnumStats(chr_transcript_stat_file) if construct_models else EnumStats()
            return read_stat, transcript_stat
        os.remove(lock_file)

    grouping_strategy_names = get_grouping_strategy_names(args)
    aggregator = ReadAssignmentAggregator(args, sample, string_pools, loader.genedb, chr_id,
                                         grouping_strategy_names=grouping_strategy_names)

    transcript_stat_counter = EnumStats()
    io_support = IOSupport(args)
    transcript_id_distributor = ExcludingIdDistributor(loader.genedb, chr_id)
    exon_id_storage = FeatureIdStorage(SimpleIDDistributor(), loader.genedb, chr_id, "exon")

    # Use chromosome-specific prefix for output files
    chr_prefix = sample.get_chr_prefix(chr_id)

    if construct_models:
        tmp_gff_printer = GFFPrinter(sample.out_dir, chr_prefix, exon_id_storage,
                                     output_r2t=large_output_enabled(args, "read2transcripts"),
                                     check_canonical=args.check_canonical)
    else:
        tmp_gff_printer = VoidTranscriptPrinter()
    if construct_models and args.genedb:
        tmp_extended_gff_printer = GFFPrinter(sample.out_dir, chr_prefix, exon_id_storage,
                                              gtf_suffix=".extended_annotation.gtf",
                                              output_r2t=False, check_canonical=args.check_canonical)
    else:
        tmp_extended_gff_printer = VoidTranscriptPrinter()

    sqanti_t2t_printer = SqantiTSVPrinter(sample.get_t2t_tsv_file(chr_id), args, IOSupport(args)) \
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
    aggregator.read_stat_counter.dump(chr_read_stat_file)
    if construct_models:
        if loader.genedb:
            all_models, gene_info = create_extended_storage(loader.genedb, chr_id, loader.chr_record, novel_model_storage)
            if args.check_canonical:
                io_support.add_canonical_info(all_models, gene_info)
            tmp_extended_gff_printer.dump(gene_info, all_models)
        aggregator.transcript_model_global_counter.dump()
        aggregator.gene_model_global_counter.dump()
        transcript_stat_counter.dump(chr_transcript_stat_file)
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
        from .read_groups import parse_barcode2spot_spec
        filename, barcode_col, spot_cols = parse_barcode2spot_spec(args.barcode2spot)
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) > barcode_col:
                    barcode = parts[barcode_col]
                    # Extract only the specified spot columns
                    cell_types = [parts[col] for col in spot_cols if col < len(parts)]
                    barcode_feature_table[barcode] = cell_types

    # Build string pools for memory optimization
    # Barcode pool is critical: ReadAssignments were serialized with barcode_id values
    # that reference the barcode pool from read collection
    string_pools = setup_string_pools(args, sample, chr_ids, chr_id,
                                      load_barcode_pool=True, load_tsv_pools=False)

    # Load dynamic pools (for read groups from BAM tags/read IDs)
    load_dynamic_pools(string_pools, sample.get_dynamic_pools_file(chr_id))

    umi_filter = UMIFilter(args.umi_length, edit_distance)
    filtered_reads = sample.get_filtered_reads_file(chr_id) if output_filtered_reads else None
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
