#!/usr/bin/env python3
#
# ############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################
import glob
import logging
import sys
from shutil import rmtree
from io import StringIO

from src.gtf2db import *
from src.read_mapper import *
from src.dataset_processor import *
from src.long_read_counter import COUNTING_STRATEGIES

logger = logging.getLogger('IsoQuant')


def parse_args(args=None, namespace=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    show_full_help = '--full_help' in sys.argv

    def add_additional_option(*args, **kwargs):  # show command only with --full-help
        if not show_full_help:
            kwargs['help'] = argparse.SUPPRESS
        parser.add_argument(*args, **kwargs)

    parser.add_argument("--output", "-o", help="output folder, will be created automatically [default=isoquant_output]",
                        type=str, default="isoquant_output")
    parser.add_argument('--debug', action='store_true', default=False,
                        help='Debug log output.')
    # REFERENCE
    parser.add_argument("--genedb", "-g", help="gene database in gffutils DB format or GTF/GFF format", type=str)

    parser.add_argument('--complete_genedb', action='store_true', default=False,
                        help="use this flag if gene annotation contains transcript and gene metafeatures, "
                             "e.g. with official annotations, such as GENCODE; "
                             "speeds up gene database conversion")
    parser.add_argument("--reference", "-r", help="reference genome in FASTA format (can be gzipped)",
                        type=str)
    parser.add_argument("--index", help="genome index for specified aligner (optional)", type=str)
    parser.add_argument('--clean_start', action='store_true', default=False,
                        help='Do not use previously generated index, feature db or alignments.')
    # INPUT READS
    input_args = parser.add_mutually_exclusive_group()
    input_args.add_argument('--bam', nargs='+', type=str,
                            help='sorted and indexed BAM file(s), each file will be treated as a separate sample')
    input_args.add_argument('--fastq', nargs='+', type=str,
                            help='input FASTQ file(s), each file will be treated as a separate sample; '
                                 'reference genome should be provided when using reads as input')
    input_args.add_argument('--bam_list', type=str, help='text file with list of BAM files, one file per line'
                                                         ', leave empty line between samples')
    input_args.add_argument('--fastq_list', type=str, help='text file with list of FASTQ files, one file per line'
                                                           ', leave empty line between samples')

    parser.add_argument('--labels', '-l', nargs='+', type=str,
                        help='sample names to be used; input file names are used if not set')
    parser.add_argument("--read_group", help="a way to group feature counts (no grouping by default): "
                                             "by BAM file tag (tag:TAG); "
                                             "using additional file (file:FILE:READ_COL:GROUP_COL:DELIM); "
                                             "using read id (read_id:DELIM); "
                                             "by original file name (file_name)", type=str)

    # INPUT PROPERTIES
    parser.add_argument("--data_type", "-d", type=str, choices=DATA_TYPE_ALIASES.keys(),
                        help="type of data to process, supported types are: " + ", ".join(DATA_TYPE_ALIASES.keys()))
    parser.add_argument('--stranded',  type=str, help="reads strandness type, supported values are: " +
                        ", ".join(SUPPORTED_STRANDEDNESS), default="none")
    parser.add_argument('--fl_data', action='store_true', default=False,
                        help="reads represent FL transcripts; both ends of the read are considered to be reliable")

    # ALGORITHM
    add_additional_option("--delta", type=int, default=None,
                          help="delta for inexact splice junction comparison, chosen automatically based on data type")
    parser.add_argument("--matching_strategy", choices=["exact", "precise", "default", "loose"],
                        help="read-to-isoform matching strategy from the most strict to least", type=str, default=None)
    parser.add_argument("--splice_correction_strategy", choices=["none", "default_pacbio", "default_ont", "conservative_ont", "all", "assembly"],
                        help="read alignment correction strategy to use",
                        type=str, default=None)
    parser.add_argument("--model_construction_strategy", choices=["reliable", "default_pacbio", "sensitive_pacbio", "fl_pacbio",
                                                                  "default_ont", "sensitive_ont", "all", "assembly"],
                        help="transcript model construction strategy to use",
                        type=str, default=None)

    parser.add_argument("--transcript_quantification", choices=COUNTING_STRATEGIES,
                        help="transcript quantification strategy", type=str, default="with_ambiguous")
    parser.add_argument("--gene_quantification", choices=COUNTING_STRATEGIES,
                        help="gene quantification strategy", type=str, default="with_inconsistent")

    # OUTPUT PROPERTIES
    parser.add_argument("--full_help", action='help', help="show full list of options")
    parser.add_argument("--test", action=TestMode, nargs=0, help="run IsoQuant on toy dataset")
    parser.add_argument("--threads", "-t", help="number of threads to use", type=int, default="16")

    parser.add_argument('--check_canonical', action='store_true', default=False,
                        help="report whether splice junctions are canonical")
    parser.add_argument("--sqanti_output", help="produce SQANTI-like TSV output (requires more time)",
                        action='store_true', default=False)
    parser.add_argument("--count_exons", help="perform exon and intron counting", action='store_true', default=False)

    # PIPELINE STEPS
    resume_args = parser.add_mutually_exclusive_group()
    resume_args.add_argument("--resume", action="store_true", default=False,
                             help="resume failed run, specify output folder, input options are not allowed")
    resume_args.add_argument("--force", action="store_true", default=False,
                             help="force to overwrite the previous run")

    parser.add_argument("--no_model_construction", action="store_true", default=False,
                          help="run only read assignment and quantification")
    parser.add_argument("--run_aligner_only", action="store_true", default=False,
                          help="align reads to reference without running further analysis")

    # ADDITIONAL
    add_additional_option("--low_memory", help="decrease RAM consumption (leads to slower processing)",
                          action='store_true', default=False)
    add_additional_option("--no_junc_bed", action="store_true", default=False,
                          help="do NOT use annotation for read mapping")
    add_additional_option("--junc_bed_file", type=str,
                          help="annotation in BED format produced by minimap's paftools.js gff2bed "
                               "(will be created automatically if not given)")
    add_additional_option("--no_secondary", help="ignore secondary alignments (not recommended)", action='store_true',
                          default=False)
    add_additional_option("--keep_tmp", help="do not remove temporary files in the end", action='store_true',
                          default=False)
    add_additional_option("--read_assignments", nargs='+', type=str,
                          help="reuse read assignments (binary format) to construct transcript models",
                          default=None)
    add_additional_option("--aligner", help="force to use this alignment method, can be " + ", ".join(SUPPORTED_ALIGNERS) +
                                            "; chosen based on data type if not set", type=str)
    add_additional_option("--genedb_output", help="output folder for converted gene database,"
                                                  " will be created automatically (same as output by default)", type=str)
    add_additional_option("--cage", help="bed file with CAGE peaks", type=str, default=None)
    add_additional_option("--cage-shift", type=int, default=50, help="interval before read start to look for CAGE peak")

    isoquant_version = "3.0.0"
    try:
        with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "VERSION")) as version_f:
            isoquant_version = version_f.readline().strip()
    except FileNotFoundError:
        pass
    parser.add_argument('--version', '-v', action='version', version='IsoQuant ' + isoquant_version)

    args = parser.parse_args(args, namespace)

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
        resume_parser.add_argument("--threads", "-t", help="number of threads to use", type=int, default=argparse.SUPPRESS)
        resume_parser.add_argument("--low_memory", help="decrease RAM consumption (leads to slower processing)",
                                   action='store_true', default=argparse.SUPPRESS)
        resume_parser.add_argument("--keep_tmp", help="do not remove temporary files in the end", action='store_true',
                                   default=argparse.SUPPRESS)
        args, unknown_args = resume_parser.parse_known_args(sys.argv[1:])
        if unknown_args:
            logger.error("You cannot specify options other than --output/--threads/--debug/--low_memorty "
                         "with --resume option")
            exit(-2)

    args._cmd_line = " ".join(sys.argv)
    args._version = isoquant_version

    args.output_exists = os.path.exists(args.output)
    if not args.output_exists:
        os.makedirs(args.output)

    return args


def check_and_load_args(args):
    args.param_file = os.path.join(args.output, ".params")
    if args.resume:
        if not os.path.exists(args.output) or not os.path.exists(args.param_file):
            # logger is not defined yet
            logger.error("Previous run config was not detected, cannot resume. "
                         "Check that output folder is correctly specified.")
            exit(-1)
        args = load_previous_run(args)
    elif args.output_exists:
        if os.path.exists(args.param_file):
            if args.force:
                logger.warning("Output folder already contains previous run. Will be overwritten.")
            else:
                logger.error("Output folder already contains previous run. Set --force to overwrite or use "
                             "--resume to continue unfinished run.")
                exit(-2)
        else:
            logger.warning("Output folder already exists, some files may be overwritten.")

    args.gtf = args.genedb
    if args.genedb_output is None:
        args.genedb_output = args.output
    elif not os.path.exists(args.genedb_output):
        os.makedirs(args.genedb_output)

    if not check_input_params(args):
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

    if not args.fastq and not args.fastq_list and not args.bam and not args.bam_list:
        logger.error("No input data was provided")
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

    check_input_files(args)
    return True


def check_input_files(args):
    for sample in args.input_data.samples:
        for lib in sample.file_list:
            for in_file in lib:
                if args.input_data.input_type == "save":
                    saves = glob.glob(in_file + "*")
                    if not saves:
                        print("ERROR! Input files " + in_file + "* do not exist")
                    continue
                if not os.path.isfile(in_file):
                    print("ERROR! Input file " + in_file + " does not exist")
                    exit(-1)
                if args.input_data.input_type == "bam":
                    # TODO: sort and index file if needed
                    bamfile_in = pysam.AlignmentFile(in_file, "rb")
                    if not bamfile_in.has_index():
                        print("ERROR! BAM file " + in_file + " is not indexed, run samtools sort and samtools index")
                        exit(-1)
                    bamfile_in.close()

    if args.cage is not None:
        logger.critical("CAGE data is not supported yet")
        exit(-1)
        if not os.path.isfile(args.cage):
            print("ERROR! Bed file with CAGE peaks " + args.cage + " does not exist")
            exit(-1)

    if args.genedb is not None:
        if not os.path.isfile(args.genedb):
            print("ERROR! Gene database " + args.genedb + " does not exist")
            exit(-1)
    else:
        args.no_junc_bed = True

    if args.read_assignments is not None:
        for r in args.read_assignments:
            if not glob.glob(r + "*"):
                print("No files found with prefix " + str(r))
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
    file_mode = "a" if args.resume else "w"
    f = open(log_file, file_mode)
    f.write("Command line: " + args._cmd_line + '\n')
    f.close()
    fh = logging.FileHandler(log_file)
    fh.setLevel(output_level)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger_instance.addHandler(fh)
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
    args.needs_polya_for_construction = False
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

    args.delta = args.delta if args.delta is not None else strategy.delta
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
                                            'fl_only'))
    strategies = {
        'reliable':        ModelConstructionStrategy(2, 0.5, 20,  5, 0.05,  1, 0.1,  0.1,  2, 4, 8, 0.05, 0.05, 50, True),
        'default_pacbio':  ModelConstructionStrategy(1, 0.5, 10,  2, 0.02,  1, 0.05,  0.05,  1, 2, 2, 0.02, 0.005, 100, False),
        'sensitive_pacbio':ModelConstructionStrategy(1, 0.5, 10,  2, 0.005,  1, 0.01,  0.02,  1, 2, 2, 0.005, 0.001, 100, False),
        'default_ont':     ModelConstructionStrategy(1, 0.5, 20,  3, 0.02,  1, 0.05,  0.05,  1, 3, 3, 0.02, 0.02, 10, False),
        'sensitive_ont':   ModelConstructionStrategy(1, 0.5, 20,  3, 0.005,  1, 0.01,  0.02,  1, 2, 3, 0.005, 0.005, 10, False),
        'fl_pacbio':       ModelConstructionStrategy(1, 0.5, 10,  2, 0.02,  1, 0.05,  0.01,  1, 2, 3, 0.02, 0.005, 100, True),
        'all':             ModelConstructionStrategy(0, 0.3, 10,  1, 0.002,  1, 0.01, 0.01, 1, 1, 1, 0.002, 0.001, 500, False),
        'assembly':        ModelConstructionStrategy(0, 0.3, 10,  1, 0.05,  1, 0.01, 0.02,  1, 1, 1, 0.05, 0.01, 50, False)
    }
    strategy = strategies[args.model_construction_strategy]

    args.min_novel_intron_count = strategy.min_novel_intron_count
    args.graph_clustering_ratio = strategy.graph_clustering_ratio
    args.graph_clustering_distance = strategy.graph_clustering_distance
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
    args.no_polya = False

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

    args.multi_intron_mapping_quality_cutoff = 5
    args.mono_mapping_quality_cutoff = 1
    args.simple_models_mapq_cutoff = 30


def run_pipeline(args):
    logger.info(" === IsoQuant pipeline started === ")

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
        if len(args.input_data.samples) > 1:
            combine_counts(args.input_data, args.output)

    logger.info(" === IsoQuant pipeline finished === ")



# Test mode is triggered by --test option
class TestMode(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        source_dir = os.path.dirname(os.path.realpath(__file__))
        options = ['--output', 'isoquant_test', '--threads', '2',
                   '--fastq', os.path.join(source_dir, 'tests/simple_data/chr9.4M.ont.sim.fq.gz'),
                   '--reference', os.path.join(source_dir, 'tests/simple_data/chr9.4M.fa.gz'),
                   '--genedb', os.path.join(source_dir, 'tests/simple_data/chr9.4M.gtf.gz'),
                   '--clean_start', '--data_type', 'nanopore', '--complete_genedb', '--force']
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

        correct_results = ['total assignments 4', 'inconsistent: 1', 'unique: 1', 'known: 2', 'Processed 1 sample']
        return all([result in log for result in correct_results])


def main(args):
    args = parse_args(args)
    set_logger(args, logger)
    args = check_and_load_args(args)
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
        strout = StringIO()
        if logger.handlers:
            print_exc(file=strout)
            s = strout.read()
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
