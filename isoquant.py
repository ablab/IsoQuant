#!/usr/bin/env python3
#
# ############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################
# Cap per-process math-library thread pools to 1 BEFORE importing numpy /
# scipy / xgboost (OpenBLAS/OpenMP read these at import time). IsoQuant already
# parallelizes at the chromosome/process level, so letting each forked worker
# spin up an all-cores BLAS/OpenMP team only oversubscribes the CPU: massive
# cpu_time inflation (idle OpenMP spin-wait) with no wall-clock gain.
# setdefault() leaves any user-provided override in place.
import os as _os
from random import choices

for _thread_env in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    _os.environ.setdefault(_thread_env, "1")
############################################################################
import argparse
import glob
import json
import logging
import os.path
import pickle
import shutil
import sys
import time
import gzip
from collections import namedtuple
from io import StringIO
from traceback import print_exc
from concurrent.futures import ProcessPoolExecutor
import concurrent.futures

import pysam
import gffutils
import pyfaidx

from isoquant_lib.utils.error_codes import IsoQuantExitCode
from isoquant_lib.modes import IsoQuantMode, ISOQUANT_MODES, AnalysisType, ANALYSIS_ALIASES, ANALYSIS_CHOICES
from isoquant_lib.gtf2db import convert_gtf_to_db
from isoquant_lib.utils.read_mapper import (
    DATA_TYPE_ALIASES,
    SUPPORTED_STRANDEDNESS,
    SUPPORTED_ALIGNERS,
    ASSEMBLY,
    PACBIO_CCS_DATA,
    NANOPORE_DATA,
    DataSetReadMapper
)
from isoquant_lib.alignment.alignment_processor import PolyATrimmed
from isoquant_lib.dataset_processor import DatasetProcessor, PolyAUsageStrategies
from isoquant_lib.model_construction.model_construction import StrandnessReportingLevel
from isoquant_lib.assignment.long_read_assigner import AmbiguityResolvingMethod
from isoquant_lib.quantification.long_read_counter import COUNTING_STRATEGIES, CountingStrategy, GroupedOutputFormat, NormalizationMethod
from isoquant_lib.utils.input_data_storage import InputDataStorage, InputDataType
from isoquant_lib.assignment.multimap_resolver import MultimapResolvingStrategy
from isoquant_lib.utils.stats import combine_counts
from isoquant_lib.barcode_calling import process_single_thread, process_in_parallel, get_umi_length
from isoquant_lib.common import setup_worker_logging, _get_log_params


logger = logging.getLogger('IsoQuant')

# Large output file types for --large_output option
LARGE_OUTPUT_TYPES = ["read_info", "read_assignments", "corrected_bed", "read2transcripts", "allinfo", "none"]


def bool_str(s):
    s = s.lower()
    if s not in {'false', 'true', '0', '1'}:
        raise ValueError('Not a valid boolean string')
    return s == 'true' or s == '1'


def parse_args(cmd_args=None, namespace=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    ref_args_group = parser.add_argument_group('Reference data')
    input_args_group = parser.add_argument_group('Input data')
    output_args_group = parser.add_argument_group('Output naming')
    pipeline_args_group = parser.add_argument_group('Pipeline options')
    algo_args_group = parser.add_argument_group('Algorithm settings')
    output_setup_args_group = parser.add_argument_group('Output configuration')
    align_args_group = parser.add_argument_group('Aligner settings')
    filer_args_group = parser.add_argument_group('Read filtering options')
    sc_args_group = parser.add_argument_group('Single-cell/spatial-related options:')

    show_full_help = '--full_help' in cmd_args

    def add_option_to_group(opt_group, *args, **kwargs):  # show command only with --full-help
        opt_group.add_argument(*args, **kwargs)

    def add_additional_option_to_group(opt_group, *args, **kwargs):  # show command only with --full-help
        if not show_full_help:
            kwargs['help'] = argparse.SUPPRESS
        opt_group.add_argument(*args, **kwargs)

    def add_hidden_option(*args, **kwargs):  # show command only with --full-help
        kwargs['help'] = argparse.SUPPRESS
        parser.add_argument(*args, **kwargs)

    parser.add_argument("--full_help", action='help', help="show full list of options")

    parser.add_argument("--test", action=TestMode, nargs=0, help="run IsoQuant on toy dataset")

    add_hidden_option('--debug', action='store_true', default=False,
                      help='debug log output.')

    # REFERENCE
    add_option_to_group(ref_args_group,"--reference", "-r",
                        help="reference genome in FASTA format (can be gzipped)",
                        type=str)
    add_option_to_group(ref_args_group,"--genedb", "-g", help="gene database in gffutils DB format or GTF/GFF "
                                                              "format (optional)", type=str)
    add_option_to_group(ref_args_group,'--complete_genedb', action='store_true', default=False,
                        help="use this flag if gene annotation contains transcript and gene features (e.g. for official annotations)")
    add_additional_option_to_group(ref_args_group, "--discard_chr", nargs="+", help="chromosome IDs to ignore",
                                   type=str, default=[])
    add_additional_option_to_group(ref_args_group, "--process_only_chr", nargs="+", help="chromosome IDs to process",
                                   type=str, default=None)
    add_additional_option_to_group(ref_args_group, "--index", help="genome index for specified aligner (optional)",
                                   type=str)

    # INPUT READS
    add_option_to_group(input_args_group, "--data_type", "-d", type=str, choices=DATA_TYPE_ALIASES.keys(),
                        help="type of data to process")
    input_args = input_args_group.add_mutually_exclusive_group()
    input_args.add_argument('--bam', nargs='+', type=str,
                            help='sorted and indexed BAM file(s), each file will be treated as a separate sample')
    input_args.add_argument('--fastq', nargs='+', type=str,
                            help='input FASTQ/FASTA file(s) with reads, each file will be treated as a separate sample')
    input_args.add_argument('--unmapped_bam', nargs='+', type=str,
                            help='unmapped BAM file(s), each file will be treated as a separate sample')
    input_args.add_argument('--yaml', type=str, help='yaml file containing all input files, one entry per sample')

    add_option_to_group(input_args_group, '--illumina_bam', nargs='+', type=str,
                                  help='sorted and indexed file(s) with Illumina reads from the same sample')

    add_option_to_group(input_args_group, "--read_group", nargs='+', type=str,
                        help="read grouping rules (space-separated); supported values: "
                             "tag:TAG (BAM tag), "
                             "file:FILE:READ_COLUMNS:GROUP_COLUMNS:DELIM (TSV file), "
                             "read_id:DELIM (read ID suffix), "
                             "file_name (filename), "
                             "barcode_spot (barcodes to spots/cell via --barcode2spot), "
                             "barcode_barcode (barcodes to spots via --barcode2barcode), "
                             "barcode (barcode--barcoded_reads), "
                             "no_auto (only explicitly stated rules), "
                             "none (disable grouping)")

    add_additional_option_to_group(input_args_group, "--read_assignments", nargs='+', type=str,
                                   help="reuse read assignments (binary format)", default=None)

    # INPUT PROPERTIES
    add_option_to_group(input_args_group, '--stranded',  type=str, choices=SUPPORTED_STRANDEDNESS,
                        help="reads strandness type [none]", default="none")
    add_option_to_group(input_args_group, '--polya_trimmed', default=PolyATrimmed.none.name, type=str,
                        choices=[e.name for e in PolyATrimmed],
                        help="define reads which had polyA tail trimmed [%s]" % PolyATrimmed.none.name)
    add_option_to_group(input_args_group, '--fl_data', action='store_true', default=False,
                        help="reads represent FL transcripts; both ends of the read are considered to be reliable")

    # OUTPUT
    add_option_to_group(output_args_group, "--output", "-o",
                        help="output folder, will be created automatically [isoquant_output]",
                        type=str, default="isoquant_output")
    add_option_to_group(output_args_group,'--prefix', '-p', type=str,
                        help='experiment name; to be used for folder and file naming [OUT]',
                        default="OUT")
    add_option_to_group(output_args_group,'--labels', '-l', nargs='+', type=str,
                        help='sample/replica labels to be used as column names; input file names are used '
                             'if not set; must be equal to the number of input files given via --fastq/--bam')

    # PIPELINE STEPS
    add_option_to_group(pipeline_args_group, "--analysis", nargs='+', type=str,
                        choices=ANALYSIS_CHOICES, default=None,
                        metavar="ANALYSIS",
                        help="analyses to run (space-separated); supported values: "
                             "quantification/quant, transcript_discovery/td, "
                             "exon_quantification/ex_quant, fusion [auto]")

    add_option_to_group(pipeline_args_group, "--threads", "-t", help="number of threads [16]", type=int,
                        default="16")

    # deprecated
    add_additional_option_to_group(pipeline_args_group, "--count_exons",
                                   help="deprecated: use --analysis exon_quantification",
                                   action='store_true', default=False)
    add_additional_option_to_group(pipeline_args_group, "--count_intron_retentions",
                                   help="deprecated: use --analysis exon_quantification",
                                   action='store_true', default=False)
    add_additional_option_to_group(pipeline_args_group, "--no_model_construction", action="store_true",
                                   default=False, help="deprecated: omit transcript_discovery from --analysis")

    resume_args = pipeline_args_group.add_mutually_exclusive_group()
    resume_args.add_argument("--resume", action="store_true", default=False,
                             help="resume failed run, specify output folder, input options are not allowed")
    resume_args.add_argument("--force", action="store_true", default=False,
                             help="force to overwrite the previous run")

    add_additional_option_to_group(pipeline_args_group, '--clean_start', action='store_true', default=False,
                                   help='Do not use previously generated index, feature db or alignments.')

    add_additional_option_to_group(pipeline_args_group, "--run_aligner_only", action="store_true", default=False,
                                   help="align reads to reference without running further analysis")
    add_additional_option_to_group(pipeline_args_group, "--no_gtf_check", help="do not perform GTF checks",
                                   dest="gtf_check",
                                   action='store_false', default=True)
    add_additional_option_to_group(pipeline_args_group, "--high_memory",
                                   help="increase RAM consumption (store alignment and the genome in RAM)",
                                   action='store_true', default=False)
    add_additional_option_to_group(pipeline_args_group, "--keep_tmp", help="do not remove temporary files "
                                                                           "in the end", action='store_true',
                                   default=False)

    # OUTPUT SETUP
    add_option_to_group(output_setup_args_group, '--check_canonical', action='store_true', default=False,
                                     help="report whether splice junctions are canonical")
    add_option_to_group(output_setup_args_group, "--sqanti_output", help="produce SQANTI-like TSV output",
                                     action='store_true', default=False)
    add_option_to_group(output_setup_args_group, "--bam_tags",
                                         help="comma separated list of BAM tags to be imported to read_assignments.tsv",
                                         type=str)
    add_option_to_group(output_setup_args_group, "--large_output", nargs='+', type=str,
                                         choices=LARGE_OUTPUT_TYPES,
                                         default=["read_info", "read2transcripts"],
                                         help="large output files to generate [read_info read2transcripts]")

    add_additional_option_to_group(output_setup_args_group, "--genedb_output", help="output folder for converted gene "
                                                                                    "database, will be created automatically "
                                                                                    " (same as output by default)", type=str)

    add_additional_option_to_group(output_setup_args_group, "--no_gzip", help="do not gzip large output files",
                                   dest="gzipped", action='store_false', default=True)
    add_additional_option_to_group(output_setup_args_group, "--normalization_method",
                                   type=str, choices=[e.name for e in NormalizationMethod],
                                   help="TPM normalization method",
                                   default=NormalizationMethod.simple.name)
    add_additional_option_to_group(output_setup_args_group, "--counts_format", type=str, nargs='+',
                                   choices=[e.name for e in GroupedOutputFormat],
                                   help="output format for grouped counts",
                                   default=[GroupedOutputFormat.default.name])
    add_additional_option_to_group(output_setup_args_group, "--emit_read_ids", action='store_true', default=False,
                                   help="add read id columns to exon splice-site counts output")
    add_additional_option_to_group(output_setup_args_group, "--old_exon_count_format", action='store_true', default=False,
                                   help="output old exon inclusion/exclusion counts (deprecated)")

    # ALIGNER
    add_additional_option_to_group(align_args_group, "--aligner", choices=SUPPORTED_ALIGNERS,
                                   help="use this alignmer method, can be [minimap2]", type=str)
    add_additional_option_to_group(align_args_group, "--no_junc_bed", action="store_true", default=False,
                                   help="do NOT use annotation for read mapping")
    add_additional_option_to_group(align_args_group, "--junc_bed_file", type=str,
                                   help="annotation in BED format produced by minimap's paftools.js gff2bed "
                                        "(will be created automatically if not given)")
    add_additional_option_to_group(align_args_group, "--indexing_options", type=str,
                                   help="additional options that will be passed to the aligner indexer")
    add_additional_option_to_group(align_args_group, "--mapping_options", type=str,
                                   help="additional options that will be passed to the aligner")


    # READ FILTERING
    add_additional_option_to_group(filer_args_group, "--use_secondary",
                                   help="use secondary alignments (slower processing)",
                                   action='store_true', default=False)
    add_additional_option_to_group(filer_args_group, "--min_mapq",
                                   help="ignore alignments with MAPQ < this (including secondary alignments) [None]", type=int)
    add_additional_option_to_group(filer_args_group, "--inconsistent_mapq_cutoff",
                                   help="ignore inconsistent alignments with MAPQ < cutoff (works only with the reference annotation) [5]",
                                   type=int, default=5)
    add_additional_option_to_group(filer_args_group, "--simple_alignments_mapq_cutoff",
                                   help="ignore alignments with 1 or 2 exons and MAPQ < cutoff "
                                        "(works only in annotation-free mode) [1]", type=int, default=1)
    add_additional_option_to_group(filer_args_group, "--max_coverage_small_chr",
                                   help="process only a fraction of reads for high-coverage loci on small chromosomes (e.g. MT), "
                                        "improves running time and RAM [1000000]",
                                   type=int, default=1000000)
    add_additional_option_to_group(filer_args_group, "--max_coverage_normal_chr",
                                   help="process only a fraction of reads for high-coverage loci on usual chromosomes, "
                                        "improves running time and RAM [-1]",
                                   type=int, default=-1)

    # SC ARGUMENTS
    add_option_to_group(sc_args_group, "--mode", "-m", type=str, choices=ISOQUANT_MODES,
                                   help="IsoQuant mode [%s]" % IsoQuantMode.bulk.name, default=IsoQuantMode.bulk.name)
    add_option_to_group(sc_args_group, '--barcode_whitelist', type=str, nargs='+',
                                   help='file(s) with barcode whitelist(s) for barcode calling')
    add_option_to_group(sc_args_group, "--barcoded_reads", type=str, nargs='+',
                                   help='TSV file(s) with barcoded reads')
    add_option_to_group(sc_args_group, "--barcoded_bam", action='store_true', default=False,
                                   help='extract barcodes and UMIs from BAM tags (CB/UB by default)')
    add_option_to_group(sc_args_group, "--barcode2spot", type=str,
                                   help='TSV file mapping barcode to cell type / spot id.')
    add_option_to_group(sc_args_group, "--molecule", type=str,
                                   help='molecule definition file (MDF) for custom_sc mode: '
                                        'defines molecule structure for universal barcode extraction')
    add_additional_option_to_group(sc_args_group, "--barcode_tag", type=str, default="CB",
                                   help='BAM tag for cell barcode [CB]')
    add_additional_option_to_group(sc_args_group, "--umi_tag", type=str, default="UB",
                                   help='BAM tag for UMI [UB]')
    add_additional_option_to_group(sc_args_group, "--strip_barcode_suffix", action='store_true', default=False,
                                   help='remove suffix after dash from barcodes extracted from BAM tag')
    add_additional_option_to_group(sc_args_group, "--barcode2barcode", type=str,
                                   help='TSV file mapping barcode to spot IDs for UMI deduplication')

    # ALGORITHM
    add_additional_option_to_group(algo_args_group, "--report_novel_unspliced", "-u", type=bool_str,
                                   help="report novel monoexonic transcripts (true/false), "
                                        "default: false for ONT, true for other data types")
    add_additional_option_to_group(algo_args_group, "--report_canonical",  type=str,
                                   choices=[e.name for e in StrandnessReportingLevel],
                                   help="reporting level for novel transcripts based on canonical splice sites "
                                        "[%s]" % StrandnessReportingLevel.auto.name,
                                   default=StrandnessReportingLevel.auto.name)
    add_additional_option_to_group(algo_args_group, "--polya_requirement", type=str,
                                   choices=[e.name for e in PolyAUsageStrategies],
                                   help="require polyA tails to be present when reporting transcripts; "
                                        "default: auto (requires polyA only when polyA percentage is >= 70%%)",
                                   default=PolyAUsageStrategies.auto.name)
    # Alternative-polyA/TSS isoform discovery is default ON when an annotation is
    # given (--genedb). --novel_apa (default off) extends it to novel
    # (non-reference) transcripts, not only known ones.
    add_additional_option_to_group(algo_args_group, "--novel_apa", action="store_true", default=False,
                                   help="generate novel isoforms with identical intron chain but distinct polyA sites")

    # Splice-site correction shifts novel transcript-model junctions onto canonical
    # motifs using clustered read deletions near junctions; default ON.
    add_additional_option_to_group(algo_args_group, "--no_splice_site_correction", action="store_true",
                                   default=False,
                                   help="disable canonical splice-site correction for novel transcript models")

    add_additional_option_to_group(algo_args_group, "--transcript_quantification", choices=COUNTING_STRATEGIES,
                                   help="transcript quantification strategy [%s]" % CountingStrategy.unique_only.name,
                                   type=str, default=CountingStrategy.unique_only.name)
    add_additional_option_to_group(algo_args_group, "--gene_quantification", choices=COUNTING_STRATEGIES,
                                   help="gene quantification strategy [%s]" % CountingStrategy.unique_splicing_consistent.name,
                                   type=str, default=CountingStrategy.unique_splicing_consistent.name)

    add_additional_option_to_group(algo_args_group, "--matching_strategy",
                                   choices=["exact", "precise", "default", "loose"],
                                   help="read-to-isoform matching strategy from the most strict to least",
                                   type=str, default=None)
    add_additional_option_to_group(algo_args_group, "--splice_correction_strategy",
                                   choices=["none", "default_pacbio", "default_ont",
                                            "conservative_ont", "all", "assembly"],
                                   help="read alignment correction strategy to use", type=str, default=None)
    add_additional_option_to_group(algo_args_group, "--model_construction_strategy",
                                   choices=["reliable", "default_pacbio", "sensitive_pacbio", "pacbio_all",
                                            "fl_pacbio", "default_ont", "sensitive_ont", "all", "assembly"],
                                   help="transcript model construction strategy to use", type=str, default=None)
    add_additional_option_to_group(algo_args_group, "--delta", type=int, default=None,
                                   help="delta for inexact splice junction comparison [auto]")
    add_additional_option_to_group(algo_args_group, "--use_replicas", type=bool_str, default=True,
                                   help="require novel transcripts to be confirmed by multiple files "
                                        "when file_name grouping is used [true]")

    # REST
    add_hidden_option("--graph_clustering_distance", type=int, default=None,
                      help="intron graph clustering distance, "
                           "splice junctions less that this number of bp apart will not be differentiated")
    add_hidden_option("--cage", help="bed file with CAGE peaks", type=str, default=None)
    add_hidden_option("--cage-shift", type=int, default=50, help="interval before read start to look for CAGE peak")

    # PolyA / TSS training-data collection (developer-only).
    # When either flag is given to a normal IsoQuant run, the matching
    # terminal-position counter switches from inference to dumping per-peak
    # features (the eight FEATURE_COLUMNS plus `chromosome` and a `true_peak`
    # 0/1 label computed against the annotated transcript end) to the CSV
    # path supplied as the flag value. The XGBoost model is not consulted;
    # rows are accumulated, not filtered. Feed the CSV to
    # `misc/train_polya_tss_model.py` to fit a fresh model.
    #   python isoquant.py … --genedb GENEDB --collect_polya_training peaks.csv
    #   python misc/train_polya_tss_model.py --features peaks.csv \
    #       --output isoquant_lib/data/model_polya.json
    # For TSS, add --fl_data and --collect_tss_training tss_peaks.csv.
    # These flags are intentionally hidden from --help/--full_help; the run
    # emits a clearly marked warning when they are used.
    add_hidden_option("--collect_polya_training", type=str, default=None,
                      help="developer: dump per-peak features + true_peak label for polyA training to this CSV path.")
    add_hidden_option("--collect_tss_training", type=str, default=None,
                      help="developer: dump per-peak features + true_peak label for TSS training to this CSV path.")

    isoquant_version = "4.0.0"
    try:
        with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "VERSION")) as version_f:
            isoquant_version = version_f.readline().strip()
    except FileNotFoundError:
        try:
            from importlib.metadata import version as _get_version
            isoquant_version = _get_version("isoquant")
        except Exception:
            pass
    parser.add_argument('--version', '-v', action='version', version='IsoQuant ' + isoquant_version,
                        help="show IsoQuant version and exit")

    args = parser.parse_args(cmd_args, namespace)

    if args.resume:
        resume_parser = argparse.ArgumentParser(add_help=False)
        resume_parser.add_argument("--resume", action="store_true", default=False,
                                   help="resume failed run, specify only output folder, "
                                        "input options are not allowed")
        resume_parser.add_argument("--output", "-o",
                                   help="output folder, will be created automatically [isoquant_output]",
                                   type=str, required=True)
        resume_parser.add_argument('--debug', action='store_true', default=argparse.SUPPRESS,
                                   help='debug log output.')
        resume_parser.add_argument("--threads", "-t", help="number of threads to use",
                                   type=int, default=argparse.SUPPRESS)
        resume_parser.add_argument("--high_memory",
                                   help="increase RAM consumption (store alignment and the genome in RAM)",
                                   action='store_true', default=False)
        resume_parser.add_argument("--keep_tmp", help="do not remove temporary files in the end",
                                   action='store_true', default=argparse.SUPPRESS)

        args, unknown_args = resume_parser.parse_known_args(cmd_args)
        if unknown_args:
            logger.error("You cannot specify options other than --output/--threads/--debug/--high_memory "
                         "with --resume option")
            parser.print_usage()
            sys.exit(IsoQuantExitCode.INCOMPATIBLE_OPTIONS)

    args._cmd_line = " ".join(sys.argv)
    args._version = isoquant_version

    args.output_exists = os.path.exists(args.output)
    if not args.output_exists:
        os.makedirs(args.output)

    return args, parser


def get_bam_files_from_samples(input_data) -> list:
    """Extract all BAM file paths from input_data.samples.

    Returns a list of BAM file paths (long-read and Illumina).
    """
    bam_files: list[str] = []
    for sample in input_data.samples:
        for lib in sample.file_list:
            for in_file in lib:
                bam_files.append(in_file)
        if getattr(sample, "illumina_bam", None):
            bam_files.extend(sample.illumina_bam)
    return [f for f in bam_files if os.path.isfile(f)]


def run_fusion_detection_on_samples(fd, samples: list) -> dict:
    """Run fusion detection per sample using a shared FusionDetector.

    The report is written alongside the other per-sample outputs as
    ``<sample.out_dir>/<sample.prefix>.fusions.tsv``.

    Args:
        fd: FusionDetector instance (initialized with first BAM, will be reused)
        samples: List of SampleData objects to process

    Returns:
        Dictionary with summary: {"total": int, "successful": int, "failed": int, "skipped": list}
    """
    summary = {"total": 0, "successful": 0, "failed": 0, "skipped": []}
    for sample in samples:
        sample_bams = [f for lib in sample.file_list for f in lib if os.path.isfile(f)]
        if not sample_bams:
            continue
        summary["total"] += 1
        out_fname = os.path.join(sample.out_dir, sample.prefix + ".fusions.tsv")
        try:
            fd.clear_state()  # Accumulate all of this sample's BAMs into one report
            for bam_path in sample_bams:
                logger.info("Running fusion detection on %s" % bam_path)
                fd.bam_path = bam_path  # Switch BAM
                fd.detect_fusions()
            fd.report(output_path=out_fname)
            logger.info("Fusion candidates for sample %s written to %s" % (sample.prefix, out_fname))
            summary["successful"] += 1
        except Exception as e:
            logger.error("Fusion detection failed for sample %s: %s" % (sample.prefix, str(e)))
            logger.debug("Traceback:", exc_info=True)
            summary["failed"] += 1
            summary["skipped"].append(sample.prefix)
    return summary


def check_and_load_args(args, parser):
    args.param_file = os.path.join(args.output, ".params")
    if args.resume:
        if not os.path.exists(args.output) or not os.path.exists(args.param_file):
            # logger is not defined yet
            logger.error("Previous run config was not detected, cannot resume. "
                         "Check that output folder is correctly specified.")
            sys.exit(IsoQuantExitCode.RESUME_CONFIG_NOT_FOUND)
        args = load_previous_run(args)
    elif args.output_exists:
        if os.path.exists(args.param_file):
            if args.force:
                logger.warning("Output folder already contains a previous run, will be overwritten.")
            else:
                logger.warning("Output folder already contains a previous run, some files may be overwritten. "
                               "Use --resume to resume a failed run. Use --force to avoid this message.")
                logger.warning("Press Ctrl+C to interrupt the run now.")
                delay = 9
                for i in range(delay):
                    countdown = delay - i
                    sys.stdout.write("Resuming the run in %d second%s\r" % (countdown, "s" if countdown > 1 else ""))
                    time.sleep(1)
                logger.info("Overwriting the previous run")
                time.sleep(1)
        else:
            logger.warning("Output folder already exists, some files may be overwritten.")

    if args.genedb_output is None:
        args.genedb_output = args.output
    elif not os.path.exists(args.genedb_output):
        os.makedirs(args.genedb_output)
    if not args.genedb:
        args.genedb_filename = None
    elif args.genedb.lower().endswith("db"):
        args.genedb_filename = args.genedb
    else:
        args.genedb_filename = os.path.join(args.genedb_output, os.path.splitext(os.path.basename(args.genedb))[0] + ".db")

    if not check_input_params(args):
        parser.print_usage()
        sys.exit(IsoQuantExitCode.INVALID_PARAMETER)

    # Validate --read_group none
    if args.read_group:
        if "none" in args.read_group and len(args.read_group) > 1:
            logger.error("--read_group 'none' cannot be combined with other values")
            sys.exit(IsoQuantExitCode.INVALID_PARAMETER)

    # Validate --large_output values
    if args.large_output:
        if "none" in args.large_output and len(args.large_output) > 1:
            logger.error("--large_output 'none' cannot be combined with other values")
            sys.exit(IsoQuantExitCode.INVALID_PARAMETER)
        for val in args.large_output:
            if val not in LARGE_OUTPUT_TYPES:
                logger.error("Invalid --large_output value: %s. Valid values: %s" % (val, ", ".join(LARGE_OUTPUT_TYPES)))
                sys.exit(IsoQuantExitCode.INVALID_PARAMETER)

    save_params(args)
    return args


def load_previous_run(args):
    logger.info("Loading parameters from the previous run")
    logger.error("Only --output/--threads/--debug/--high_memory are compatible with --resume option")
    unpickler = pickle.Unpickler(open(args.param_file, "rb"), fix_imports=False)
    loaded_args = unpickler.load()

    for option in args.__dict__:
        loaded_args.__dict__[option] = args.__dict__[option]

    if loaded_args.debug:
        logger.setLevel(logging.DEBUG)
        logger.handlers[0].setLevel(logging.DEBUG)

    return loaded_args


def save_params(args):
    for file_opt in ["genedb", "reference", "index", "bam", "fastq", "junc_bed_file",
                     "cage", "genedb_output", "read_assignments"]:
        if file_opt in args.__dict__ and args.__dict__[file_opt]:
            if isinstance(args.__dict__[file_opt], list):
                args.__dict__[file_opt] = list(map(os.path.abspath, args.__dict__[file_opt]))
            else:
                args.__dict__[file_opt] = os.path.abspath(args.__dict__[file_opt])

    if "read_group" in args.__dict__ and args.__dict__["read_group"]:
        updated_specs = []
        for spec in args.read_group:
            vals = spec.split(":")
            if len(vals) > 1 and vals[0] == 'file':
                vals[1] = os.path.abspath(vals[1])
                updated_specs.append(":".join(vals))
            else:
                updated_specs.append(spec)
        args.read_group = updated_specs

    pickler = pickle.Pickler(open(args.param_file, "wb"),  -1)
    pickler.dump(args)
    pass


# Translate the --analysis option (and the deprecated stage flags) into the
# internal boolean flags consumed across the pipeline. Requires args.mode and
# args.genedb to already be resolved. Warns but never aborts on infeasible
# combinations.
def resolve_analyses(args):
    is_single_cell = args.mode.needs_barcode_calling()

    # capture legacy stage flags before they are overwritten below
    legacy_count_exons = args.count_exons
    legacy_intron_retentions = args.count_intron_retentions
    legacy_no_model = args.no_model_construction

    if args.analysis is not None:
        analyses = {ANALYSIS_ALIASES[a] for a in args.analysis}
    elif is_single_cell:
        analyses = {AnalysisType.quantification}
    elif args.genedb:
        analyses = {AnalysisType.quantification, AnalysisType.transcript_discovery}
    else:
        analyses = {AnalysisType.transcript_discovery}

    # apply deprecated stage flags additively, on top of the chosen analyses
    if legacy_count_exons:
        logger.warning("--count_exons is deprecated, use --analysis exon_quantification")
        analyses.add(AnalysisType.exon_quantification)
    if legacy_intron_retentions:
        logger.warning("--count_intron_retentions is deprecated, use --analysis exon_quantification")
    if legacy_no_model:
        logger.warning("--no_model_construction is deprecated, omit transcript_discovery from --analysis")
        analyses.discard(AnalysisType.transcript_discovery)

    # feasibility: annotation-dependent analyses cannot run without a gene database
    if not args.genedb:
        for analysis in (AnalysisType.quantification, AnalysisType.exon_quantification, AnalysisType.fusion):
            if analysis in analyses:
                logger.warning("%s requires a reference annotation (--genedb) and will not be run" % analysis.name)
                analyses.discard(analysis)
        # without an annotation transcript discovery is the only available analysis
        if not analyses and not legacy_no_model:
            logger.info("Without a reference annotation only transcript discovery is available; "
                        "running transcript_discovery")
            analyses.add(AnalysisType.transcript_discovery)

    # feasibility: model construction in single-cell/spatial modes is limited
    if is_single_cell and AnalysisType.transcript_discovery in analyses:
        logger.warning("transcript_discovery in single-cell/spatial mode runs after UMI deduplication: "
                       "it will not yield novel genes and may be incomplete; "
                       "consider a pseudo-bulk run for transcript discovery")

    # canonical internal flags (names reused across the pipeline)
    args.run_quantification = AnalysisType.quantification in analyses
    args.count_exons = AnalysisType.exon_quantification in analyses or legacy_count_exons
    args.count_intron_retentions = AnalysisType.exon_quantification in analyses or legacy_intron_retentions
    args.fusion = AnalysisType.fusion in analyses
    args.no_model_construction = AnalysisType.transcript_discovery not in analyses
    # polyA/TSS site prediction is part of quantification and is also needed for
    # model construction; gated by the annotation (TSS additionally needs --fl_data)
    args.predict_terminal_sites = bool(args.genedb) and (args.run_quantification or not args.no_model_construction)
    args.analyses = analyses


# Check user's params
def _validate_data_type_and_input(args):
    if not args.reference:
        logger.error("Reference genome was not provided")
        return False
    if not args.data_type:
        logger.error("Data type is not provided, choose one of " + " ".join(DATA_TYPE_ALIASES.keys()))
        return False
    elif args.data_type not in DATA_TYPE_ALIASES.keys():
        logger.error("Unsupported data type " + args.data_type + ", choose one of: " + " ".join(DATA_TYPE_ALIASES.keys()))
        return False
    args.data_type = DATA_TYPE_ALIASES[args.data_type]

    if not any([args.fastq, args.bam, args.unmapped_bam, args.read_assignments, args.yaml]):
        logger.error("No input data was provided")
        return False

    if args.yaml and args.illumina_bam:
        logger.error("When providing a yaml file it should include all input files, including the illumina bam file.")
        return False

    args.input_data = InputDataStorage(args)
    return True


def _validate_alignment_options(args):
    if args.aligner is not None and args.aligner not in SUPPORTED_ALIGNERS:
        logger.error(" Unsupported aligner " + args.aligner + ", choose one of: " + " ".join(SUPPORTED_ALIGNERS))
        return False

    if args.run_aligner_only and not args.input_data.input_type.needs_mapping():
        logger.error("Data type %s cannot be mapped and thus incompatible with --run_aligner_only option." % args.input_data.input_type.name)
        return False
    if args.stranded not in SUPPORTED_STRANDEDNESS:
        logger.error("Unsupported strandness " + args.stranded + ", choose one of: " + " ".join(SUPPORTED_STRANDEDNESS))
        return False
    return True


def _apply_stage_warnings(args):
    if not args.genedb and args.sqanti_output:
        args.sqanti_output = False
        logger.warning("--sqanti_output option has no effect without gene annotation")

    if args.no_model_construction and args.sqanti_output:
        args.sqanti_output = False
        logger.warning("--sqanti_output option has no effect without model construction")

    if args.process_only_chr and args.discard_chr:
        args.discard_chr = []
        logger.warning("--discard_chr has not effect when --process_only_chr is set and will be ignored")


def _dedup_read_group_specs(args):
    if "read_group" in args.__dict__ and args.__dict__["read_group"]:
        updated_specs = []
        spec_set = set()
        for spec in args.read_group:
            if spec in spec_set:
                logger.warning("Read group %s is set twice, which has no effect, duplicate will be ignored" % spec)
                continue
            updated_specs.append(spec)
            spec_set.add(spec)
        args.read_group = updated_specs


def _validate_barcode_calling(args):
    args.umi_length = 0
    if not args.mode.needs_barcode_calling():
        return
    barcode_sources = sum([bool(args.barcode_whitelist), bool(args.barcoded_reads), bool(args.barcoded_bam)])
    if barcode_sources > 1:
        logger.critical("Options --barcode_whitelist, --barcoded_reads, and --barcoded_bam are mutually exclusive")
        sys.exit(IsoQuantExitCode.INVALID_PARAMETER)
    if args.mode == IsoQuantMode.custom_sc:
        if not any([args.molecule, args.barcoded_reads, args.barcoded_bam]):
            logger.critical("custom_sc mode requires --molecule, --barcoded_reads, or --barcoded_bam")
            sys.exit(IsoQuantExitCode.BARCODE_WHITELIST_MISSING)
    elif not any([args.barcode_whitelist, args.barcoded_reads, args.barcoded_bam]):
        logger.critical("You have chosen single-cell/spatial mode %s, please specify barcode whitelist, "
                        "file with barcoded reads, or --barcoded_bam" % args.mode.name)
        sys.exit(IsoQuantExitCode.BARCODE_WHITELIST_MISSING)
    if args.barcoded_bam:
        args.umi_length = _detect_umi_length_from_bam(args.input_data.samples[0].file_list[0][0], args.umi_tag)
    else:
        args.umi_length = get_umi_length(args.mode)


def check_input_params(args):
    if not _validate_data_type_and_input(args):
        return False
    if not _validate_alignment_options(args):
        return False

    if not isinstance(args.mode, IsoQuantMode):
        args.mode = IsoQuantMode[args.mode]

    # translate --analysis (and the deprecated stage flags) into internal booleans
    resolve_analyses(args)

    _apply_stage_warnings(args)
    _dedup_read_group_specs(args)
    _validate_barcode_calling(args)

    check_input_files(args)
    return True


def _detect_umi_length_from_bam(bam_path: str, umi_tag: str) -> int:
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.has_tag(umi_tag):
                return len(read.get_tag(umi_tag))
    return 0


def check_bam_file(bam_path: str, check_index: bool = True):
    """Check BAM file exists and optionally has index."""
    if not os.path.isfile(bam_path):
        logger.critical("BAM file " + bam_path + " does not exist")
        sys.exit(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)
    if check_index:
        bamfile_in = pysam.AlignmentFile(bam_path, "rb")
        if not bamfile_in.has_index():
            logger.critical("BAM file " + bam_path + " is not indexed, run samtools sort and samtools index")
            sys.exit(IsoQuantExitCode.BAM_NOT_INDEXED)
        bamfile_in.close()


def check_file_exists(file_path: str, description: str):
    """Check that a file exists, exit with error if not."""
    if not os.path.isfile(file_path):
        logger.critical(f"{description} {file_path} does not exist")
        sys.exit(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)


def extract_read_group_file_path(spec: str):
    """
    Extract file path from read_group spec if it's a file-based spec.

    Returns file path for 'file:path:...' specs, None otherwise.
    """
    parts = spec.split(":")
    if len(parts) >= 2 and parts[0] == 'file':
        return parts[1]
    return None


def _check_input_read_file(args, in_file):
    """Validate a single input read file (save prefix / BAM / FASTQ)."""
    if args.input_data.input_type == InputDataType.save:
        saves = glob.glob(in_file + "*")
        if not saves:
            logger.critical("Input files " + in_file + "* do not exist")
        return
    if not os.path.isfile(in_file):
        logger.critical("Input file " + in_file + " does not exist")
        sys.exit(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)
    if args.input_data.input_type == InputDataType.bam:
        check_bam_file(in_file, check_index=True)


def _check_reference_and_reads(args):
    # Check reference genome
    if args.reference and not os.path.isfile(args.reference):
        logger.critical("Reference genome " + args.reference + " does not exist")
        sys.exit(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)

    # Check input reads (BAM/FASTQ/save files)
    for sample in args.input_data.samples:
        for lib in sample.file_list:
            for in_file in lib:
                _check_input_read_file(args, in_file)

        # Check Illumina BAM files
        if sample.illumina_bam is not None:
            for illumina in sample.illumina_bam:
                check_bam_file(illumina, check_index=True)


def _check_barcode_input_files(args):
    # Check barcoded reads files (from args, not sample - sample.barcoded_reads is set later)
    if hasattr(args, 'barcoded_reads') and args.barcoded_reads:
        if isinstance(args.barcoded_reads, list):
            for bc_file in args.barcoded_reads:
                check_file_exists(bc_file, "Barcoded reads file")
        else:
            check_file_exists(args.barcoded_reads, "Barcoded reads file")

    # Check molecule definition file
    if hasattr(args, 'molecule') and args.molecule:
        check_file_exists(args.molecule, "Molecule definition file")

    # Check barcode whitelist files
    if hasattr(args, 'barcode_whitelist') and args.barcode_whitelist:
        for wl_file in args.barcode_whitelist:
            check_file_exists(wl_file, "Barcode whitelist file")


def _check_barcode_mapping_files(args):
    from isoquant_lib.assignment.read_groups import parse_barcode2spot_spec

    # Check barcode2spot file (parse spec to extract filename)
    if hasattr(args, 'barcode2spot') and args.barcode2spot:
        bc2spot_file, _, _ = parse_barcode2spot_spec(args.barcode2spot)
        check_file_exists(bc2spot_file, "Barcode to spot mapping file")

    # Check barcode2barcode file (parse spec to extract filename)
    if hasattr(args, 'barcode2barcode') and args.barcode2barcode:
        bc2bc_file, _, _ = parse_barcode2spot_spec(args.barcode2barcode)
        check_file_exists(bc2bc_file, "Barcode to barcode mapping file")


def _check_read_group_files(args):
    # Check read_group file specs
    if hasattr(args, 'read_group') and args.read_group:
        for spec in args.read_group:
            file_path = extract_read_group_file_path(spec)
            if file_path:
                check_file_exists(file_path, "Read group file")


def _check_annotation_files(args):
    # Check junction BED file
    if hasattr(args, 'junc_bed_file') and args.junc_bed_file:
        check_file_exists(args.junc_bed_file, "Junction BED file")

    # Check CAGE file (currently not supported)
    if args.cage is not None:
        logger.critical("CAGE data is not supported yet")
        sys.exit(IsoQuantExitCode.INVALID_PARAMETER)
        if not os.path.isfile(args.cage):
            logger.critical("Bed file with CAGE peaks " + args.cage + " does not exist")
            sys.exit(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)

    # Check gene database
    if args.genedb is not None:
        if not os.path.isfile(args.genedb):
            logger.critical("Gene database " + args.genedb + " does not exist")
            sys.exit(IsoQuantExitCode.GENE_DB_NOT_FOUND)
    else:
        args.no_junc_bed = True

    # Check read assignments
    if args.read_assignments is not None:
        for r in args.read_assignments:
            if not glob.glob(r + "*"):
                logger.critical("No files found with prefix " + str(r))
                sys.exit(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)


def check_input_files(args):
    _check_reference_and_reads(args)
    _check_barcode_input_files(args)
    _check_barcode_mapping_files(args)
    _check_read_group_files(args)
    _check_annotation_files(args)


def create_output_dirs(args):
    for sample in args.input_data.samples:
        sample_dir = sample.out_dir
        if os.path.exists(sample_dir):
            if not args.resume:
                logger.warning(sample_dir + " folder already exists, some files may be overwritten")
        else:
            os.makedirs(sample_dir)
        sample_aux_dir = sample.aux_dir
        if os.path.exists(sample_aux_dir):
            if not args.resume:
                logger.warning(sample_aux_dir + " folder already exists, some files may be overwritten")
        else:
            os.makedirs(sample_aux_dir)


def set_logger(args):
    output_level = logging.DEBUG if args.__dict__.get('debug') else logging.INFO
    log_file = os.path.join(args.output, "isoquant.log")
    if os.path.exists(log_file):
        old_log_file = os.path.join(args.output, "isoquant.log.old")
        with open(old_log_file, "a") as olf:
            olf.write("\n")
            shutil.copyfileobj(open(log_file, "r"), olf)
    with open(log_file, "w") as f:
        f.write("Command line: " + args._cmd_line + '\n')
    setup_worker_logging(log_file, output_level)
    logger.info("Running IsoQuant version " + args._version)


def set_data_dependent_options(args):
    matching_strategies = {ASSEMBLY: "precise", PACBIO_CCS_DATA: "precise", NANOPORE_DATA: "default"}
    if args.matching_strategy is None:
        args.matching_strategy = matching_strategies[args.data_type]

    model_construction_strategies = {ASSEMBLY: "assembly", PACBIO_CCS_DATA: "default_pacbio", NANOPORE_DATA: "default_ont"}
    if args.model_construction_strategy is None:
        args.model_construction_strategy = model_construction_strategies[args.data_type]
        if args.fl_data and args.model_construction_strategy == "default_pacbio":
            args.model_construction_strategy = "fl_pacbio"

    splice_correction_strategies = {ASSEMBLY: "assembly", PACBIO_CCS_DATA: "default_pacbio", NANOPORE_DATA: "default_ont"}
    if args.splice_correction_strategy is None:
        args.splice_correction_strategy = splice_correction_strategies[args.data_type]

    args.resolve_ambiguous = 'monoexon_and_fsm' if args.fl_data else 'default'
    args.requires_polya_for_construction = False

    # Handle --read_group none: disable all auto-added groupings
    if args.read_group is not None and "none" in args.read_group:
        args.read_group = None
        return

    # Automatically add file_name grouping when multiple files are present
    if args.input_data.has_replicas():
        if args.read_group is None:
            # No read grouping specified, use file_name
            args.read_group = ["file_name"]
        else:
            # Read grouping specified, ensure file_name is included
            if "file_name" not in args.read_group and "no_auto" not in args.read_group:
                args.read_group.append("file_name")

    # Automatically add barcode_spot grouping when --barcode2spot is set
    if hasattr(args, 'barcode2spot') and args.barcode2spot:
        if args.read_group is None:
            args.read_group = ["barcode_spot"]
        elif "barcode_spot" not in args.read_group and "no_auto" not in args.read_group:
            args.read_group.append("barcode_spot")

    if hasattr(args, 'barcode2barcode') and args.barcode2barcode:
        # Automatically add barcode_barcode grouping when --barcode2barcode is set
        if args.read_group is None:
            args.read_group = ["barcode_barcode"]
        elif "barcode_barcode" not in args.read_group and "no_auto" not in args.read_group:
            args.read_group.append("barcode_barcode")

        # Automatically add allinfo output
        if args.large_output is None:
            args.large_output = ["allinfo"]
        elif "none" in args.large_output:
            pass
        elif "allinfo" not in args.large_output:
            args.large_output.append("allinfo")

    # In SC modes, auto-add barcode grouping if no barcode-related grouping is set
    if args.mode.needs_barcode_calling():
        barcode_groupings = {"barcode", "barcode_spot", "barcode_barcode", "no_auto"}
        has_barcode_grouping = (args.read_group is not None and
                                any(rg in barcode_groupings for rg in args.read_group))
        if not has_barcode_grouping:
            if args.read_group is None:
                args.read_group = ["barcode"]
            else:
                args.read_group.append("barcode")
            logger.info("Single-cell/spatial mode: automatically adding '--read_group barcode'. "
                        "Use '--read_group none' or `--read_group no_auto` to disable.")


def set_matching_options(args):
    MatchingStrategy = namedtuple('MatchingStrategy',
                                  ('delta', 'max_intron_shift', 'max_missed_exon_len', 'max_fake_terminal_exon_len',
                                   'max_suspicious_intron_abs_len', 'max_suspicious_intron_rel_len',
                                   'resolve_ambiguous', 'correct_minor_errors'))

    strategies = {
        'exact':   MatchingStrategy(0, 0, 0, 0, 0, 0.0, 'monoexon_only', False),
        'precise': MatchingStrategy(4, 30, 50, 20, 0, 0.0, 'monoexon_and_fsm', True),
        'default': MatchingStrategy(6, 60, 100, 40, 60, 1.0, 'monoexon_and_fsm', True),
        'loose':   MatchingStrategy(12, 60, 100, 40, 60, 1.0, 'all',  True),
    }

    strategy = strategies[args.matching_strategy]

    if args.delta is None:
        args.delta = strategy.delta
    elif args.delta < 0:
        logger.error("--delta can not be negative")
        sys.exit(IsoQuantExitCode.INVALID_PARAMETER)
    args.minor_exon_extension = 50
    args.major_exon_extension = 300
    args.max_intron_shift = strategy.max_intron_shift
    args.max_missed_exon_len = strategy.max_missed_exon_len
    args.max_fake_terminal_exon_len = strategy.max_fake_terminal_exon_len
    # short introns that are actually long deletions, fix minimaps logic
    args.max_suspicious_intron_abs_len = strategy.max_suspicious_intron_abs_len
    args.max_suspicious_intron_rel_len = strategy.max_suspicious_intron_rel_len
    args.min_abs_exon_overlap = 10
    args.min_rel_exon_overlap = 0.2
    args.micro_intron_length = 50
    args.max_intron_abs_diff = min(30, args.max_intron_shift)
    args.max_intron_rel_diff = 0.2
    args.apa_delta = args.minor_exon_extension
    args.minimal_exon_overlap = 5
    args.minimal_intron_absence_overlap = 20
    args.polya_window = 16
    args.polya_fraction = 0.75
    if args.resolve_ambiguous == 'default':
        args.resolve_ambiguous = strategy.resolve_ambiguous
    if args.resolve_ambiguous not in AmbiguityResolvingMethod.__dict__:
        logger.error("Incorrect resolving ambiguity method: " + args.resolve_ambiguous + ", default will be used")
        args.resolve_ambiguous = strategy.resolve_ambiguous
    args.resolve_ambiguous = AmbiguityResolvingMethod[args.resolve_ambiguous]
    args.correct_minor_errors = strategy.correct_minor_errors

    updated_strategy = MatchingStrategy(args.delta, args.max_intron_shift, args.max_missed_exon_len,
                                        args.max_fake_terminal_exon_len,
                                        args.max_suspicious_intron_abs_len, args.max_suspicious_intron_rel_len,
                                        args.resolve_ambiguous, args.correct_minor_errors)
    logger.debug('Using %s strategy. Updated strategy: %s.' % (args.matching_strategy, updated_strategy))


def set_splice_correction_options(args):
    SplicSiteCorrectionStrategy = namedtuple('SplicSiteCorrectionStrategy',
                                             ('fuzzy_junctions', 'intron_shifts', 'skipped_exons',
                                              'terminal_exons', 'fake_terminal_exons', 'microintron_retention'))
    strategies = {
        'none': SplicSiteCorrectionStrategy(False, False, False, False, False, False),
        'default_pacbio': SplicSiteCorrectionStrategy(True, False, True, False, False, True),
        'conservative_ont': SplicSiteCorrectionStrategy(True, False, True, False, False, False),
        'default_ont': SplicSiteCorrectionStrategy(True, False, True, False, True, True),
        'all': SplicSiteCorrectionStrategy(True, True, True, True, True, True),
        'assembly': SplicSiteCorrectionStrategy(False, False, True, False, False, False)
    }
    strategy = strategies[args.splice_correction_strategy]
    args.correct_fuzzy_junctions = strategy.fuzzy_junctions
    args.correct_intron_shifts = strategy.intron_shifts
    args.correct_skipped_exons = strategy.skipped_exons
    args.correct_terminal_exons = strategy.terminal_exons
    args.correct_fake_terminal_exons = strategy.fake_terminal_exons
    args.correct_microintron_retention = strategy.microintron_retention


def set_model_construction_options(args):
    ModelConstructionStrategy = namedtuple('ModelConstructionStrategy',
                                           ('min_novel_intron_count',
                                            'graph_clustering_ratio', 'graph_clustering_distance',
                                            'min_novel_isolated_intron_abs', 'min_novel_isolated_intron_rel',
                                            'terminal_position_abs', 'terminal_position_rel',
                                            'terminal_internal_position_rel',
                                            'min_known_count', 'min_nonfl_count',
                                            'min_novel_count', 'min_novel_count_rel',
                                            'min_mono_count_rel', 'singleton_adjacent_cov',
                                            'fl_only', 'novel_monoexonic',
                                            'require_monointronic_polya', 'require_monoexonic_polya',
                                            'report_canonical'))
    strategies = {
        'reliable':        ModelConstructionStrategy(2, 0.5, 20,  5, 0.05,  1, 0.1,  0.1,  2, 4, 8, 0.05, 0.05, 50,
                                                     True, False, True, True, StrandnessReportingLevel.only_canonical),
        'default_pacbio':  ModelConstructionStrategy(1, 0.5, 10,  2, 0.02,  1, 0.05,  0.05,  1, 2, 2, 0.02, 0.005, 100,
                                                     False, True, False, True, StrandnessReportingLevel.only_canonical),
        'sensitive_pacbio': ModelConstructionStrategy(1, 0.5, 5,   2, 0.005,  1, 0.01,  0.02,  1, 2, 2, 0.005, 0.001, 100,
                                                      False, True, False, False, StrandnessReportingLevel.only_stranded),
        # Like sensitive_pacbio but admits single-read novel isoforms (min_novel_count=1)
        # with a lower relative floor. Very sensitive discovery for clean full-length
        # PacBio reads (e.g. Iso-Seq FLNC); recovers 5'-truncated genes at the cost of
        # some novel precision, so not recommended for noisy data.
        'pacbio_all':      ModelConstructionStrategy(1, 0.5, 5,   2, 0.005,  1, 0.01,  0.02,  1, 2, 1, 0.002, 0.001, 100,
                                                     False, True, False, False, StrandnessReportingLevel.only_stranded),
        'default_ont':     ModelConstructionStrategy(1, 0.5, 20,  3, 0.02,  1, 0.05,  0.05,  1, 3, 3, 0.02, 0.02, 10,
                                                     False, False, True, True, StrandnessReportingLevel.only_canonical),
        'sensitive_ont':   ModelConstructionStrategy(1, 0.5, 20,  3, 0.005,  1, 0.01,  0.02,  1, 2, 3, 0.005, 0.005, 10,
                                                     False, True, False, False, StrandnessReportingLevel.only_stranded),
        'fl_pacbio':       ModelConstructionStrategy(1, 0.5, 10,  2, 0.02,  1, 0.05,  0.01,  1, 2, 3, 0.02, 0.005, 100,
                                                     True, True, False, False, StrandnessReportingLevel.only_canonical),
        'all':             ModelConstructionStrategy(0, 0.3, 5,   1, 0.002,  1, 0.01, 0.01, 1, 1, 1, 0.002, 0.001, 500,
                                                     False, True, False, False, StrandnessReportingLevel.all),
        'assembly':        ModelConstructionStrategy(0, 0.3, 5,   1, 0.05,  1, 0.01, 0.02,  1, 1, 1, 0.05, 0.01, 50,
                                                     False, True, False, False, StrandnessReportingLevel.only_stranded)
    }
    strategy = strategies[args.model_construction_strategy]

    args.min_novel_intron_count = strategy.min_novel_intron_count
    args.graph_clustering_ratio = strategy.graph_clustering_ratio
    if args.graph_clustering_distance is None:
        args.graph_clustering_distance = strategy.graph_clustering_distance
    elif args.graph_clustering_distance < 0:
        logger.error("--graph_clustering_distance can not be negative")
        sys.exit(IsoQuantExitCode.INVALID_PARAMETER)
    args.min_novel_isolated_intron_abs = strategy.min_novel_isolated_intron_abs
    args.min_novel_isolated_intron_rel = strategy.min_novel_isolated_intron_rel
    args.terminal_position_abs = strategy.terminal_position_abs
    args.terminal_position_rel = strategy.terminal_position_rel
    args.terminal_internal_position_rel = strategy.terminal_internal_position_rel

    args.min_known_count = strategy.min_known_count
    args.min_nonfl_count = strategy.min_nonfl_count
    args.min_novel_count = strategy.min_novel_count
    args.min_mono_count_rel = strategy.min_mono_count_rel
    args.min_novel_count_rel = strategy.min_novel_count_rel
    args.singleton_adjacent_cov = strategy.singleton_adjacent_cov
    args.fl_only = strategy.fl_only
    args.min_mono_exon_coverage = 0.75

    if args.report_novel_unspliced is None:
        args.report_novel_unspliced = strategy.novel_monoexonic

    if not args.report_novel_unspliced and not args.no_model_construction:
        logger.info("Novel unspliced transcripts will not be reported, "
                    "set --report_novel_unspliced true to discover them")

    args.require_monointronic_polya = strategy.require_monointronic_polya
    args.require_monoexonic_polya = strategy.require_monoexonic_polya
    args.polya_requirement_strategy = PolyAUsageStrategies[args.polya_requirement]
    args.polya_trimmed = PolyATrimmed[args.polya_trimmed]
    args.report_canonical_strategy = StrandnessReportingLevel[args.report_canonical]
    if args.report_canonical_strategy == StrandnessReportingLevel.auto:
        args.report_canonical_strategy = strategy.report_canonical


def set_configs_directory(args):
    config_dir = os.path.join(os.environ['HOME'], '.config', 'IsoQuant')
    os.makedirs(config_dir, exist_ok=True)

    args.db_config_path = os.path.join(config_dir, 'db_config.json')
    args.index_config_path = os.path.join(config_dir, 'index_config.json')
    args.bed_config_path = os.path.join(config_dir, 'bed_config.json')
    args.alignment_config_path = os.path.join(config_dir, 'alignment_config.json')
    for config_path in (args.db_config_path, args.index_config_path, args.bed_config_path, args.alignment_config_path):
        if not os.path.exists(config_path):
            with open(config_path, 'w') as f_out:
                json.dump({}, f_out, indent=2)


def set_additional_params(args):
    set_configs_directory(args)
    set_data_dependent_options(args)
    set_matching_options(args)
    set_model_construction_options(args)
    set_splice_correction_options(args)

    args.print_additional_info = True
    args.indel_near_splice_site_dist = 10
    args.upstream_region_len = 20

    args.multimap_strategy = "take_best"
    multimap_strategies = {}
    for e in MultimapResolvingStrategy:
        multimap_strategies[e.name] = e.value
    args.multimap_strategy = MultimapResolvingStrategy(multimap_strategies[args.multimap_strategy])

    args.needs_reference = True
    if args.needs_reference and not args.reference:
        logger.warning("Reference genome is not provided! This may affect quality of the results!")
        args.needs_reference = False

    args.simple_models_mapq_cutoff = 30
    args.polya_percentage_threshold = 0.7
    args.low_polya_percentage_threshold = 0.1

    if args.bam_tags:
        args.bam_tags = args.bam_tags.split(",")
    else:
        args.bam_tags = []
    args.original_annotation = None


def prepare_reference_genome(args):
    if not args.needs_reference:
        return
    logger.info("Reading reference genome from %s" % args.reference)
    ref_dir = os.path.dirname(args.reference)
    ref_file_name = os.path.basename(args.reference)
    ref_name, outer_ext = os.path.splitext(ref_file_name)

    # make symlink for pyfaidx index
    args.fai_file_name = args.reference + ".fai"
    if not os.path.exists(args.fai_file_name) and not os.access(ref_dir, os.W_OK):
        # index does not exist near the reference and reference folder is not writable
        # store index in the output folder in this case
        args.fai_file_name = os.path.join(args.output, ref_file_name + ".fai")

    low_ext = outer_ext.lower()
    if low_ext in ['.gz', '.gzip', '.bgz']:
        gunzipped_reference = os.path.join(args.output, ref_name)
        if not os.path.exists(gunzipped_reference) or not args.resume:
            logger.info("Decompressing reference to " + str(gunzipped_reference))
            with open(gunzipped_reference, "w") as outf:
                shutil.copyfileobj(gzip.open(args.reference, "rt"), outf)
        args.reference = gunzipped_reference


class BarcodeCallingArgs:
    def __init__(self, input, barcode_whitelist, mode, output, out_fasta, tmp_dir, threads,
                 molecule: str = None):
        self.input = input  # Can be a single file (str) or list of files
        self.barcodes = barcode_whitelist
        self.mode = mode
        self.output_tsv = output  # Can be a single filename (str) or list of filenames
        self.out_fasta = out_fasta  # Can be a single filename (str), list of filenames, or None
        self.tmp_dir = tmp_dir
        self.threads = threads
        self.molecule = molecule


def call_barcodes(args):
    if args.barcoded_bam:
        logger.info("Barcodes will be extracted from BAM tags (%s/%s)" % (args.barcode_tag, args.umi_tag))
        return
    if args.barcoded_reads:
        # TODO barcoded files via YAML
        args.input_data.samples[0].barcoded_reads = args.barcoded_reads
        return
    for sample in args.input_data.samples:
        # Collect all input files for this sample
        input_files = [files[0] for files in sample.file_list]
        output_barcodes_list = [sample.barcodes_tsv + "_%d.tsv" % i for i in range(len(input_files))]
        barcodes_done_list = [sample.barcodes_done + "_%d.tsv" % i for i in range(len(input_files))]

        output_fasta_list = None
        new_reads = []
        if args.mode.produces_new_fasta():
            output_fasta_list = [sample.split_reads_fasta + "_%d.fa" % i for i in range(len(input_files))]
            new_reads = [[fasta] for fasta in output_fasta_list]

        # Check if all files were already processed during resume
        all_done = all(os.path.exists(done) for done in barcodes_done_list)
        if all_done and args.resume:
            logger.info("Barcodes were called during the previous run, skipping")
            sample.barcoded_reads.extend(output_barcodes_list)
            if args.mode.produces_new_fasta():
                sample.file_list = new_reads
            continue

        # Remove existing done markers
        for barcodes_done in barcodes_done_list:
            if os.path.exists(barcodes_done):
                os.remove(barcodes_done)

            bc_threads = 1 if args.mode.enforces_single_thread() else args.threads
            bc_args = BarcodeCallingArgs(input_files, args.barcode_whitelist, args.mode,
                                         output_barcodes_list, output_fasta_list, sample.aux_dir, bc_threads,
                                         molecule=getattr(args, 'molecule', None))
            # Launching barcode calling in a separate process has the following reason:
            # Read chunks are not cleared by the GC in the end of barcode calling, leaving the main
            # IsoQuant process to consume ~2,5 GB even when barcode calling is done.
            # Once 16 child processes are created later, IsoQuant instantly takes threads x 2,5 GB for nothing.
            log_file, log_level = _get_log_params()
            with ProcessPoolExecutor(max_workers=1,
                                     initializer=setup_worker_logging,
                                     initargs=(log_file, log_level)) as proc:
                logger.info("Detecting barcodes for %d file(s)" % len(input_files))
                if bc_threads == 1:
                    future_res = proc.submit(process_single_thread, bc_args)
                else:
                    future_res = proc.submit(process_in_parallel, bc_args)

            concurrent.futures.wait([future_res],  return_when=concurrent.futures.ALL_COMPLETED)
            if future_res.exception() is not None:
                raise future_res.exception()

        # Mark all files as done and add to barcoded_reads
        for i, (input_file, output_barcodes, barcodes_done) in enumerate(zip(input_files, output_barcodes_list, barcodes_done_list)):
            sample.barcoded_reads.append(output_barcodes)
            open(barcodes_done, "w").close()
            logger.info("Processed %s, barcodes are stored in %s" % (input_file, output_barcodes))

        if args.mode.produces_new_fasta():
            logger.info("Reads were split during barcode calling")
            logger.info("The following files will be used instead of original reads %s " % ", ".join(map(lambda x: x[0], new_reads)))
            sample.file_list = new_reads


def run_pipeline(args):
    logger.info(" === IsoQuant pipeline started === ")
    logger.info("Python version: %s" % sys.version)
    logger.info("gffutils version: %s" % gffutils.__version__)
    logger.info("pysam version: %s" % pysam.__version__)
    logger.info("pyfaidx version: %s" % pyfaidx.__version__)

    if args.mode.needs_barcode_calling():
        # call barcodes
        call_barcodes(args)

    # gunzip refernece genome if needed
    prepare_reference_genome(args)

    # convert GTF/GFF if needed
    if args.genedb and not args.genedb.lower().endswith('db'):
        args.original_annotation = args.genedb
        args.genedb = convert_gtf_to_db(args)

    # map reads if fastqs are provided
    if args.input_data.input_type.needs_mapping():
        # substitute input reads with bams
        dataset_mapper = DataSetReadMapper(args)
        args.index = dataset_mapper.index_fname
        args.input_data = dataset_mapper.map_reads(args)

    if args.run_aligner_only:
        logger.info("Isoform assignment step is skipped because --run-aligner-only option was used")
        return

    # run isoform assignment
    dataset_processor = DatasetProcessor(args)
    dataset_processor.process_all_samples(args.input_data)

    # aggregate counts for all samples
    if len(args.input_data.samples) > 1 and args.genedb and args.run_quantification:
        combine_counts(args.input_data, args.output)

    # Run fusion detection after isoform detection when fusion is enabled
    if getattr(args, "fusion", False):
        logger.info(" === Isoform detection completed, starting fusion detection === ")
        bam_files = get_bam_files_from_samples(args.input_data)
        if not args.genedb:
            logger.warning("Fusion detection requires --genedb; skipping")
        elif not bam_files:
            logger.warning("No BAM files available for fusion detection; skipping")
        else:
            try:
                from isoquant_lib.fusion.fusion_detector import FusionDetector
                fd = FusionDetector(bam_files[0], args.genedb, reference_fasta=args.reference)
                summary = run_fusion_detection_on_samples(fd, args.input_data.samples)
                logger.info("Fusion detection summary: %d total, %d successful, %d failed" %
                            (summary["total"], summary["successful"], summary["failed"]))
                logger.info(" === Fusion detection finished === ")
            except Exception as e:
                logger.warning("Fusion detection encountered an error and was skipped: %s" % str(e))
                logger.debug("Traceback:", exc_info=True)

    logger.info(" === IsoQuant pipeline finished === ")


# Test mode is triggered by --test option
class TestMode(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        self.out_dir = 'isoquant_test'
        if os.path.exists(self.out_dir):
            shutil.rmtree(self.out_dir)
        source_dir = os.path.dirname(os.path.realpath(__file__))
        test_data_dir = os.path.join(source_dir, 'isoquant_lib', 'test_data')
        if not os.path.isdir(test_data_dir):
            # pip-installed: find test data via src package
            import isoquant_lib
            test_data_dir = os.path.join(os.path.dirname(os.path.realpath(isoquant_lib.__file__)), 'test_data')
        if not os.path.isdir(test_data_dir):
            sys.stderr.write("ERROR: Test data not found. Cannot run in test mode.\n")
            sys.exit(1)
        options = ['--output', self.out_dir, '--threads', '2',
                   '--fastq', os.path.join(test_data_dir, 'chr9.4M.ont.sim.fq.gz'),
                   '--reference', os.path.join(test_data_dir, 'chr9.4M.fa.gz'),
                   '--genedb', os.path.join(test_data_dir, 'chr9.4M.gtf.gz'),
                   '--clean_start', '--data_type', 'nanopore', '--complete_genedb', '--force', '-p', 'TEST_DATA']
        print('=== Running in test mode === ')
        print("Running IsoQuant in test mode with the following options:")
        print(' '.join(options))
        print('Any other option is ignored ')
        main(options)
        if self._check_log():
            logger.info(' === TEST PASSED CORRECTLY === ')
        else:
            logger.error(' === TEST FAILED ===')
            sys.exit(IsoQuantExitCode.TEST_FAILED)
        parser.exit()

    def _check_log(self):
        with open(os.path.join(self.out_dir, 'isoquant.log'), 'r') as f:
            log = f.read()

        correct_results = ['total assignments 4', 'polyA tail detected in 2', 'unique: 1', 'known: 2', 'Processed 1 experiment']
        return all([result in log for result in correct_results])


def main(cmd_args):
    args, parser = parse_args(cmd_args)
    if not cmd_args:
        parser.print_usage()
        sys.exit(IsoQuantExitCode.SUCCESS)
    set_logger(args)
    if getattr(args, "collect_polya_training", None) or getattr(args, "collect_tss_training", None):
        logger.warning("=" * 78)
        logger.warning("DEVELOPER MODE: polyA/TSS training-data collection enabled.")
        logger.warning("Per-peak features will be written to the supplied CSV path INSTEAD of")
        logger.warning("running the production XGBoost peak filter. Predictions in")
        logger.warning("*.polyA_prediction.tsv / *.TSS_prediction.tsv will be empty.")
        logger.warning("This mode is intended only for retraining the shipped model; see")
        logger.warning("misc/train_polya_tss_model.py and .claude/POLYA_TSS_TRAINING.md.")
        logger.warning("=" * 78)
    args = check_and_load_args(args, parser)
    create_output_dirs(args)
    set_additional_params(args)
    run_pipeline(args)


def main_entry():
    """Entry point for console_scripts (pip install)."""
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except KeyboardInterrupt:
        raise
    except Exception:
        if logger.handlers:
            strout = StringIO()
            print_exc(file=strout)
            s = strout.getvalue()
            if s:
                logger.critical("IsoQuant failed with the following error, please, submit this issue to "
                                "https://github.com/ablab/IsoQuant/issues\n" + s)
            else:
                print_exc()
        else:
            sys.stderr.write("IsoQuant failed with the following error, please, submit this issue to "
                             "https://github.com/ablab/IsoQuant/issues\n")
            print_exc()
        sys.exit(IsoQuantExitCode.UNCAUGHT_EXCEPTION)


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except KeyboardInterrupt:
        raise
    except Exception:
        if logger.handlers:
            strout = StringIO()
            print_exc(file=strout)
            s = strout.getvalue()
            if s:
                logger.critical("IsoQuant failed with the following error, please, submit this issue to "
                                "https://github.com/ablab/IsoQuant/issues\n" + s)
            else:
                print_exc()
        else:
            sys.stderr.write("IsoQuant failed with the following error, please, submit this issue to "
                             "https://github.com/ablab/IsoQuant/issues\n")
            print_exc()
        sys.exit(IsoQuantExitCode.UNCAUGHT_EXCEPTION)
