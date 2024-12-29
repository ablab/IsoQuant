#!/usr/bin/env python3
#
# ############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
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
from collections import namedtuple
from io import StringIO
from traceback import print_exc
from concurrent.futures import ProcessPoolExecutor
import concurrent.futures

import pysam
import gffutils
import pyfaidx

from src.gtf2db import convert_gtf_to_db
from src.read_mapper import (
    DATA_TYPE_ALIASES,
    SUPPORTED_STRANDEDNESS,
    SUPPORTED_ALIGNERS,
    ASSEMBLY,
    PACBIO_CCS_DATA,
    NANOPORE_DATA,
    DataSetReadMapper
)
from src.dataset_processor import DatasetProcessor, PolyAUsageStrategies,  ISOQUANT_MODES, IsoQuantMode
from src.graph_based_model_construction import StrandnessReportingLevel
from src.long_read_assigner import AmbiguityResolvingMethod
from src.long_read_counter import COUNTING_STRATEGIES, CountingStrategy, NormalizationMethod, GroupedOutputFormat
from src.input_data_storage import InputDataStorage
from src.multimap_resolver import MultimapResolvingStrategy
from src.stats import combine_counts
from detect_barcodes import process_single_thread, process_in_parallel
from src.barcode_calling.umi_filtering import UMIFilter, create_transcript_info_dict, load_barcodes


logger = logging.getLogger('IsoQuant')


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
    sc_args_group = parser.add_argument_group('Single-cell/spatial-related options:')

    other_options = parser.add_argument_group("Additional options:")
    show_full_help = '--full_help' in cmd_args

    def add_additional_option(*args, **kwargs):  # show command only with --full-help
        if not show_full_help:
            kwargs['help'] = argparse.SUPPRESS
        other_options.add_argument(*args, **kwargs)

    def add_additional_option_to_group(opt_group, *args, **kwargs):  # show command only with --full-help
        if not show_full_help:
            kwargs['help'] = argparse.SUPPRESS
        opt_group.add_argument(*args, **kwargs)

    def add_hidden_option(*args, **kwargs):  # show command only with --full-help
        kwargs['help'] = argparse.SUPPRESS
        parser.add_argument(*args, **kwargs)

    parser.add_argument("--full_help", action='help', help="show full list of options")
    add_hidden_option('--debug', action='store_true', default=False,
                      help='Debug log output.')

    output_args_group.add_argument("--output", "-o", help="output folder, will be created automatically "
                                                          "[default=isoquant_output]",
                                   type=str, default="isoquant_output")
    output_args_group.add_argument('--prefix', '-p', type=str,
                                   help='experiment name; to be used for folder and file naming; default is OUT',
                                   default="OUT")
    output_args_group.add_argument('--labels', '-l', nargs='+', type=str,
                                   help='sample/replica labels to be used as column names; input file names are used '
                                        'if not set; must be equal to the number of input files given via --fastq/--bam')
    # REFERENCE
    ref_args_group.add_argument("--reference", "-r", help="reference genome in FASTA format (can be gzipped)",
                                type=str)
    ref_args_group.add_argument("--genedb", "-g", help="gene database in gffutils DB format or GTF/GFF "
                                                       "format (optional)", type=str)
    ref_args_group.add_argument('--complete_genedb', action='store_true', default=False,
                                help="use this flag if gene annotation contains transcript and gene metafeatures, "
                                     "e.g. with official annotations, such as GENCODE; "
                                     "speeds up gene database conversion")
    add_additional_option_to_group(ref_args_group, "--index", help="genome index for specified aligner (optional)",
                                   type=str)

    # INPUT READS

    input_args = input_args_group.add_mutually_exclusive_group()
    input_args.add_argument('--bam', nargs='+', type=str,
                            help='sorted and indexed BAM file(s), each file will be treated as a separate sample')
    input_args.add_argument('--fastq', nargs='+', type=str,
                            help='input FASTQ file(s), each file will be treated as a separate sample; '
                                 'reference genome should be provided when using reads as input')
    add_additional_option_to_group(input_args,'--bam_list', type=str,
                                   help='text file with list of BAM files, one file per line, '
                                        'leave empty line between samples')
    add_additional_option_to_group(input_args,'--fastq_list', type=str,
                                   help='text file with list of FASTQ files, one file per line, '
                                        'leave empty line between samples')
    input_args.add_argument('--yaml', type=str, help='yaml file containing all input files, one entry per sample'
                                                     ', check readme for format info')

    input_args_group.add_argument('--illumina_bam', nargs='+', type=str,
                                  help='sorted and indexed file(s) with Illumina reads from the same sample')

    input_args_group.add_argument("--read_group", help="a way to group feature counts (no grouping by default): "
                                             "by BAM file tag (tag:TAG); "
                                             "using additional file (file:FILE:READ_COL:GROUP_COL:DELIM); "
                                             "using read id (read_id:DELIM); "
                                             "by original file name (file_name)", type=str)

    # INPUT PROPERTIES
    input_args_group.add_argument("--data_type", "-d", type=str, choices=DATA_TYPE_ALIASES.keys(),
                        help="type of data to process, supported types are: " + ", ".join(DATA_TYPE_ALIASES.keys()))
    input_args_group.add_argument('--stranded',  type=str, help="reads strandness type, supported values are: " +
                        ", ".join(SUPPORTED_STRANDEDNESS), default="none")
    input_args_group.add_argument('--fl_data', action='store_true', default=False,
                        help="reads represent FL transcripts; both ends of the read are considered to be reliable")

    # SC ARGUMENTS
    sc_args_group.add_argument("--mode", "-m", type=str, choices=ISOQUANT_MODES,
                               help="IsoQuant modes: " + ", ".join(ISOQUANT_MODES) +
                                    "; default:%s" % IsoQuantMode.bulk.name, default=IsoQuantMode.bulk.name)
    sc_args_group.add_argument('--barcode_whitelist', type=str,
                               help='file with barcode whitelist for barcode calling')
    sc_args_group.add_argument("--barcoded_reads", type=str, nargs='+',
                               help='file with barcoded reads; barcodes will be called automatically if not provided')
    sc_args_group.add_argument("--barcode_column", type=str,
                               help='column with barcodes in barcoded_reads file, default=1; read id column is 0',
                               default=1)


    # ALGORITHM
    add_additional_option_to_group(algo_args_group, "--report_novel_unspliced", "-u", type=bool_str,
                                   help="report novel monoexonic transcripts (true/false), "
                                        "default: false for ONT, true for other data types")
    add_additional_option_to_group(algo_args_group, "--report_canonical",  type=str,
                                   choices=[e.name for e in StrandnessReportingLevel],
                                   help="reporting level for novel transcripts based on canonical splice sites;"
                                        " default: " + StrandnessReportingLevel.auto.name,
                                   default=StrandnessReportingLevel.only_stranded.name)
    add_additional_option_to_group(algo_args_group, "--polya_requirement", type=str,
                                   choices=[e.name for e in PolyAUsageStrategies],
                                   help="require polyA tails to be present when reporting transcripts; "
                                        "default: auto (requires polyA only when polyA percentage is >= 70%%)",
                                   default=PolyAUsageStrategies.auto.name)

    add_additional_option_to_group(algo_args_group, "--transcript_quantification", choices=COUNTING_STRATEGIES,
                                   help="transcript quantification strategy", type=str,
                                   default=CountingStrategy.unique_only.name)
    add_additional_option_to_group(algo_args_group, "--gene_quantification", choices=COUNTING_STRATEGIES,
                                   help="gene quantification strategy", type=str,
                                   default=CountingStrategy.unique_splicing_consistent.name)

    add_additional_option_to_group(algo_args_group, "--matching_strategy",
                                   choices=["exact", "precise", "default", "loose"],
                                   help="read-to-isoform matching strategy from the most strict to least",
                                   type=str, default=None)
    add_additional_option_to_group(algo_args_group, "--splice_correction_strategy",
                                   choices=["none", "default_pacbio", "default_ont",
                                            "conservative_ont", "all", "assembly"],
                                   help="read alignment correction strategy to use", type=str, default=None)
    add_additional_option_to_group(algo_args_group, "--model_construction_strategy",
                                   choices=["reliable", "default_pacbio", "sensitive_pacbio", "fl_pacbio",
                                            "default_ont", "sensitive_ont", "all", "assembly"],
                                   help="transcript model construction strategy to use", type=str, default=None)

    # OUTPUT PROPERTIES
    pipeline_args_group.add_argument("--threads", "-t", help="number of threads to use", type=int,
                                     default="16")
    pipeline_args_group.add_argument('--check_canonical', action='store_true', default=False,
                                     help="report whether splice junctions are canonical")
    pipeline_args_group.add_argument("--sqanti_output", help="produce SQANTI-like TSV output",
                                     action='store_true', default=False)
    pipeline_args_group.add_argument("--count_exons", help="perform exon and intron counting",
                                     action='store_true', default=False)
    add_additional_option_to_group(pipeline_args_group,"--bam_tags",
                                   help="comma separated list of BAM tags to be imported to read_assignments.tsv",
                                   type=str)

    # PIPELINE STEPS
    resume_args = pipeline_args_group.add_mutually_exclusive_group()
    resume_args.add_argument("--resume", action="store_true", default=False,
                             help="resume failed run, specify output folder, input options are not allowed")
    resume_args.add_argument("--force", action="store_true", default=False,
                             help="force to overwrite the previous run")
    add_additional_option_to_group(pipeline_args_group, '--clean_start', action='store_true', default=False,
                                   help='Do not use previously generated index, feature db or alignments.')

    add_additional_option_to_group(pipeline_args_group, "--no_model_construction", action="store_true",
                                   default=False, help="run only read assignment and quantification")
    add_additional_option_to_group(pipeline_args_group, "--run_aligner_only", action="store_true", default=False,
                                   help="align reads to reference without running further analysis")

    # ADDITIONAL
    add_additional_option("--delta", type=int, default=None,
                          help="delta for inexact splice junction comparison, chosen automatically based on data type")
    add_hidden_option("--graph_clustering_distance", type=int, default=None,
                      help="intron graph clustering distance, "
                           "splice junctions less that this number of bp apart will not be differentiated")
    add_additional_option("--no_gzip", help="do not gzip large output files", dest="gzipped",
                          action='store_false', default=True)
    add_additional_option("--no_gtf_check", help="do not perform GTF checks", dest="gtf_check",
                          action='store_false', default=True)
    add_additional_option("--high_memory", help="increase RAM consumption (store alignment and the genome in RAM)",
                          action='store_true', default=False)
    add_additional_option("--no_junc_bed", action="store_true", default=False,
                          help="do NOT use annotation for read mapping")
    add_additional_option("--junc_bed_file", type=str,
                          help="annotation in BED format produced by minimap's paftools.js gff2bed "
                               "(will be created automatically if not given)")
    add_additional_option("--no_secondary", help="ignore secondary alignments (not recommended)", action='store_true',
                          default=False)
    add_additional_option("--min_mapq", help="ignore alignments with MAPQ < this"
                                             "(also filters out secondary alignments, default: None)", type=int)
    add_additional_option("--inconsistent_mapq_cutoff", help="ignore inconsistent alignments with MAPQ < this "
                                                             "(works only with the reference annotation, default=5)",
                          type=int, default=5)
    add_additional_option("--simple_alignments_mapq_cutoff", help="ignore alignments with 1 or 2 exons and "
                                                                  "MAPQ < this (works only in annotation-free mode, "
                                                                  "default=1)", type=int, default=1)
    add_additional_option("--normalization_method", type=str, choices=[e.name for e in NormalizationMethod],
                          help="TPM normalization method: simple - conventional normalization using all counted reads;"
                               "usable_reads - includes all assigned reads.",
                          default=NormalizationMethod.simple.name)
    add_additional_option("--counts_format", type=str, choices=[e.name for e in GroupedOutputFormat],
                          help="output format for grouped counts",
                          default=GroupedOutputFormat.both.name)

    add_additional_option_to_group(pipeline_args_group, "--keep_tmp", help="do not remove temporary files "
                                                                           "in the end", action='store_true',
                                   default=False)
    add_additional_option_to_group(input_args_group, "--read_assignments", nargs='+', type=str,
                                   help="reuse read assignments (binary format)", default=None)
    add_hidden_option("--aligner", help="force to use this alignment method, can be " + ", ".join(SUPPORTED_ALIGNERS)
                                        + "; chosen based on data type if not set", type=str)
    add_additional_option_to_group(output_args_group, "--genedb_output", help="output folder for converted gene "
                                                                              "database, will be created automatically "
                                                                              " (same as output by default)", type=str)
    add_hidden_option("--cage", help="bed file with CAGE peaks", type=str, default=None)
    add_hidden_option("--cage-shift", type=int, default=50, help="interval before read start to look for CAGE peak")
    parser.add_argument("--test", action=TestMode, nargs=0, help="run IsoQuant on toy dataset")

    isoquant_version = "3.4.0"
    try:
        with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "VERSION")) as version_f:
            isoquant_version = version_f.readline().strip()
    except FileNotFoundError:
        pass
    parser.add_argument('--version', '-v', action='version', version='IsoQuant ' + isoquant_version)

    args = parser.parse_args(cmd_args, namespace)

    if args.resume:
        resume_parser = argparse.ArgumentParser(add_help=False)
        resume_parser.add_argument("--resume", action="store_true", default=False,
                                   help="resume failed run, specify only output folder, "
                                        "input options are not allowed")
        resume_parser.add_argument("--output", "-o",
                                   help="output folder, will be created automatically [default=isoquant_output]",
                                   type=str, required=True)
        resume_parser.add_argument('--debug', action='store_true', default=argparse.SUPPRESS,
                                   help='Debug log output.')
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
            exit(-2)

    args._cmd_line = " ".join(sys.argv)
    args._version = isoquant_version

    args.output_exists = os.path.exists(args.output)
    if not args.output_exists:
        os.makedirs(args.output)

    return args, parser


def check_and_load_args(args, parser):
    args.param_file = os.path.join(args.output, ".params")
    if args.resume:
        if not os.path.exists(args.output) or not os.path.exists(args.param_file):
            # logger is not defined yet
            logger.error("Previous run config was not detected, cannot resume. "
                         "Check that output folder is correctly specified.")
            exit(-3)
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
        args.genedb_filename = os.path.join(args.output, os.path.splitext(os.path.basename(args.genedb))[0] + ".db")

    if not check_input_params(args):
        parser.print_usage()
        exit(-1)

    save_params(args)
    return args


def load_previous_run(args):
    logger.info("Loading parameters of the previous run, all arguments will be ignored")
    unpickler = pickle.Unpickler(open(args.param_file, "rb"), fix_imports=False)
    loaded_args = unpickler.load()

    for option in args.__dict__:
        loaded_args.__dict__[option] = args.__dict__[option]

    if loaded_args.debug:
        logger.setLevel(logging.DEBUG)
        logger.handlers[0].setLevel(logging.DEBUG)

    return loaded_args


def save_params(args):
    for file_opt in ["genedb", "reference", "index", "bam", "fastq", "bam_list", "fastq_list", "junc_bed_file",
                     "cage", "genedb_output", "read_assignments"]:
        if file_opt in args.__dict__ and args.__dict__[file_opt]:
            if isinstance(args.__dict__[file_opt], list):
                args.__dict__[file_opt] = list(map(os.path.abspath, args.__dict__[file_opt]))
            else:
                args.__dict__[file_opt] = os.path.abspath(args.__dict__[file_opt])

    if "read_group" in args.__dict__ and args.__dict__["read_group"]:
        vals = args.read_group.split(":")
        if len(vals) > 1 and vals[0] == 'file':
            vals[1] = os.path.abspath(vals[1])
            args.read_group = ":".join(vals)

    pickler = pickle.Pickler(open(args.param_file, "wb"),  -1)
    pickler.dump(args)
    pass


# Check user's params
def check_input_params(args):
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

    if not args.fastq and not args.fastq_list and not args.bam and not args.bam_list and not args.read_assignments and not args.yaml:
        logger.error("No input data was provided")
        return False

    if args.yaml and args.illumina_bam:
        logger.error("When providing a yaml file it should include all input files, including the illumina bam file.")
        return False

    if args.illumina_bam and (args.fastq_list or args.bam_list):
        logger.error("Unsupported combination of list of input files and Illumina bam file."
                     "To combine multiple experiments with short read correction please use yaml input.")
        return False

    args.input_data = InputDataStorage(args)
    if args.aligner is not None and args.aligner not in SUPPORTED_ALIGNERS:
        logger.error(" Unsupported aligner " + args.aligner + ", choose one of: " + " ".join(SUPPORTED_ALIGNERS))
        return False

    if args.run_aligner_only and args.input_data.input_type == "bam":
        logger.error("Do not use BAM files with --run_aligner_only option.")
        return False
    if args.stranded not in SUPPORTED_STRANDEDNESS:
        logger.error("Unsupported strandness " + args.stranded + ", choose one of: " + " ".join(SUPPORTED_STRANDEDNESS))
        return False

    if not args.genedb:
        if args.count_exons:
            logger.warning("--count_exons option has no effect without gene annotation")
        if args.sqanti_output:
            args.sqanti_output = False
            logger.warning("--sqanti_output option has no effect without gene annotation")
        if args.no_model_construction:
            logger.warning("Setting --no_model_construction without providing a gene "
                           "annotation will not produce any meaningful results")

    if args.no_model_construction and args.sqanti_output:
        args.sqanti_output = False
        logger.warning("--sqanti_output option has no effect without model construction")

    if not isinstance(args.mode, IsoQuantMode):
        args.mode = IsoQuantMode[args.mode]
    if args.mode in [IsoQuantMode.double, IsoQuantMode.tenX]:
        if not args.barcode_whitelist and not args.barcoded_reads:
            logger.critical("You have chosen single-cell mode %s, please specify barcode whitelist or file with "
                            "barcoded reads" % args.mode.name)
            exit(-3)

    check_input_files(args)
    return True


def check_input_files(args):
    for sample in args.input_data.samples:
        for lib in sample.file_list:
            for in_file in lib:
                if args.input_data.input_type == "save":
                    saves = glob.glob(in_file + "*")
                    if not saves:
                        logger.critical("Input files " + in_file + "* do not exist")
                    continue
                if not os.path.isfile(in_file):
                    logger.critical("Input file " + in_file + " does not exist")
                    exit(-1)
                if args.input_data.input_type == "bam":
                    bamfile_in = pysam.AlignmentFile(in_file, "rb")
                    if not bamfile_in.has_index():
                        logger.critical("BAM file " + in_file + " is not indexed, run samtools sort and samtools index")
                        exit(-1)
                    bamfile_in.close()
        if sample.illumina_bam is not None:
            for illumina in sample.illumina_bam:
                bamfile_in = pysam.AlignmentFile(illumina, "rb")
                if not bamfile_in.has_index():
                    logger.critical("BAM file " + illumina + " is not indexed, run samtools sort and samtools index")
                    exit(-1)
                bamfile_in.close()

    if args.cage is not None:
        logger.critical("CAGE data is not supported yet")
        exit(-1)
        if not os.path.isfile(args.cage):
            logger.critical("Bed file with CAGE peaks " + args.cage + " does not exist")
            exit(-1)

    if args.genedb is not None:
        if not os.path.isfile(args.genedb):
            logger.critical("Gene database " + args.genedb + " does not exist")
            exit(-1)
    else:
        args.no_junc_bed = True

    if args.read_assignments is not None:
        for r in args.read_assignments:
            if not glob.glob(r + "*"):
                logger.critical("No files found with prefix " + str(r))
                exit(-1)


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


def set_logger(args, logger_instance):
    if "debug" not in args.__dict__ or not args.debug:
        output_level = logging.INFO
    else:
        output_level = logging.DEBUG

    logger_instance.setLevel(output_level)
    log_file = os.path.join(args.output, "isoquant.log")
    if os.path.exists(log_file):
        old_log_file = os.path.join(args.output, "isoquant.log.old")
        with open(old_log_file, "a") as olf:
            olf.write("\n")
            shutil.copyfileobj(open(log_file, "r"), olf)

    f = open(log_file, "w")
    f.write("Command line: " + args._cmd_line + '\n')
    f.close()
    fh = logging.FileHandler(log_file)
    fh.set_name("isoquant_file_log")
    fh.setLevel(output_level)
    ch = logging.StreamHandler(sys.stdout)
    ch.set_name("isoquant_screen_log")
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    if all(fh.get_name() != h.get_name() for h in logger_instance.handlers):
        logger_instance.addHandler(fh)
    if all(ch.get_name() != h.get_name() for h in logger_instance.handlers):
        logger_instance.addHandler(ch)

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
    if args.read_group is None and args.input_data.has_replicas():
        args.read_group = "file_name"
    args.use_technical_replicas = args.read_group == "file_name"


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
        exit(-3)
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
        'sensitive_pacbio':ModelConstructionStrategy(1, 0.5, 5,   2, 0.005,  1, 0.01,  0.02,  1, 2, 2, 0.005, 0.001, 100,
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
        exit(-3)
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
                json.dump({}, f_out)


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


class BarcodeCallingArgs:
    def __init__(self, input, barcode_whitelist, mode, output, tmp_dir, threads):
        self.input = input
        self.barcodes = barcode_whitelist
        self.mode = mode
        self.output = output
        self.tmp_dir = tmp_dir
        self.threads = threads
        self.min_score = None


def call_barcodes(args):
    if not args.barcoded_reads:
        sample = args.input_data.samples[0]
        for i, files in enumerate(sample.file_list):
            output_barcodes = sample.barcodes_tsv + "_%d.tsv" % i
            if args.resume and os.path.exists(output_barcodes):
                # FIXME could be incomplete barcode calling run
                logger.info("Barcodes were called during the previous run, skipping")
            else:
                bc_args = BarcodeCallingArgs(files[0], args.barcode_whitelist, args.mode.name,
                                             output_barcodes, sample.aux_dir, args.threads)
                # Launching barcode calling in a separate process has the following reason:
                # Read chunks are not cleared by the GC in the end of barcode calling, leaving the main
                # IsoQuant process to consume ~2,5 GB even when barcode calling is done.
                # Once 16 child processes are created later, IsoQuant instantly takes threads x 2,5 GB for nothing.
                with ProcessPoolExecutor(max_workers=1) as proc:
                    logger.info("Detecting barcodes")
                    if args.threads == 1:
                        future_res = proc.submit(process_single_thread, bc_args)
                    else:
                        future_res = proc.submit(process_in_parallel, bc_args)

                concurrent.futures.wait([future_res],  return_when=concurrent.futures.ALL_COMPLETED)
                if future_res.exception() is not None:
                    raise future_res.exception()

            args.input_data.samples[0].barcoded_reads.append(output_barcodes)
    else:
        args.input_data.samples[0].barcoded_reads = args.barcoded_reads


def filter_umis(args):
    if args.barcoded_reads:
        args.input_data.samples[0].barcoded_reads = args.barcoded_reads

    if args.genedb:
        transcript_type_dict = create_transcript_info_dict(args.genedb)
    else:
        transcript_type_dict = {}

    barcode_umi_dict = load_barcodes(args.input_data.samples[0].barcoded_reads, True)
    for d in [2, -1]:
        logger.info("== Filtering by UMIs with edit distance %d ==" % d)
        output_prefix = args.input_data.samples[0].out_umi_filtered + (".ALL" if d < 0 else "ED%d" % d)
        logger.info("Results will be saved to %s" % output_prefix)
        umi_filter = UMIFilter(barcode_umi_dict, d)
        umi_filter.process(args.input_data.samples[0].out_assigned_tsv, output_prefix, transcript_type_dict)
        logger.info("== Done filtering by UMIs with edit distance %d ==" % d)


def run_pipeline(args):
    logger.info(" === IsoQuant pipeline started === ")
    logger.info("gffutils version: %s" % gffutils.__version__)
    logger.info("pysam version: %s" % pysam.__version__)
    logger.info("pyfaidx version: %s" % pyfaidx.__version__)
    if args.mode in [IsoQuantMode.double, IsoQuantMode.tenX]:
        # call barcodes
        call_barcodes(args)

    # convert GTF/GFF if needed
    if args.genedb and not args.genedb.lower().endswith('db'):
        args.genedb = convert_gtf_to_db(args)

    # map reads if fastqs are provided
    if args.input_data.input_type == "fastq":
        # substitute input reads with bams
        dataset_mapper = DataSetReadMapper(args)
        args.index = dataset_mapper.index_fname
        args.input_data = dataset_mapper.map_reads(args)

    if args.run_aligner_only:
        logger.info("Isoform assignment step is skipped because --run-aligner-only option was used")
    else:
        # run isoform assignment
        dataset_processor = DatasetProcessor(args)
        dataset_processor.process_all_samples(args.input_data)

        # aggregate counts for all samples
        if len(args.input_data.samples) > 1 and args.genedb:
            combine_counts(args.input_data, args.output)

    if args.mode in [IsoQuantMode.double, IsoQuantMode.tenX]:
        filter_umis(args)

    logger.info(" === IsoQuant pipeline finished === ")


# Test mode is triggered by --test option
class TestMode(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        out_dir = 'isoquant_test'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        source_dir = os.path.dirname(os.path.realpath(__file__))
        options = ['--output', out_dir, '--threads', '2',
                   '--fastq', os.path.join(source_dir, 'tests/simple_data/chr9.4M.ont.sim.fq.gz'),
                   '--reference', os.path.join(source_dir, 'tests/simple_data/chr9.4M.fa.gz'),
                   '--genedb', os.path.join(source_dir, 'tests/simple_data/chr9.4M.gtf.gz'),
                   '--clean_start', '--data_type', 'nanopore', '--complete_genedb', '--force', '-p', 'TEST_DATA']
        print('=== Running in test mode === ')
        print('Any other option is ignored ')
        main(options)
        if self._check_log():
            logger.info(' === TEST PASSED CORRECTLY === ')
        else:
            logger.error(' === TEST FAILED ===')
            exit(-1)
        parser.exit()

    @staticmethod
    def _check_log():
        with open('isoquant_test/isoquant.log', 'r') as f:
            log = f.read()

        correct_results = ['total assignments 4', 'polyA tail detected in 2', 'unique: 1', 'known: 2', 'Processed 1 experiment']
        return all([result in log for result in correct_results])


def main(cmd_args):
    args, parser = parse_args(cmd_args)
    if not cmd_args:
        parser.print_usage()
        exit(0)
    set_logger(args, logger)
    args = check_and_load_args(args, parser)
    create_output_dirs(args)
    set_additional_params(args)
    run_pipeline(args)


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except KeyboardInterrupt:
        raise
    except:
        if logger.handlers:
            strout = StringIO()
            print_exc(file=strout)
            s = strout.getvalue()
            if s:
                logger.critical("IsoQuant failed with the following error, please, submit this issue to "
                                "https://github.com/ablab/IsoQuant/issues" + s)
            else:
                print_exc()
        else:
            sys.stderr.write("IsoQuant failed with the following error, please, submit this issue to "
                             "https://github.com/ablab/IsoQuant/issues")
            print_exc()
        sys.exit(-1)
