#!/usr/bin/env python
#
# ############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import json
import logging
import argparse
from traceback import print_exc
from collections import namedtuple
from shutil import rmtree

import gffutils
import pysam
from Bio import SeqIO

from src.input_data_storage import *
from src.gtf2db import *
from src.read_mapper import *
from src.dataset_processor import *

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

    # REFERENCE
    parser.add_argument("--genedb", "-g", help="gene database in gffutils DB format or GTF/GFF format", type=str,
                        required='--run_aligner_only' not in sys.argv)
    parser.add_argument('--complete_genedb', action='store_true', default=False,
                        help="use this flag if gene annotation contains transcript and gene metafeatures, "
                             "e.g. with official annotations, such as GENCODE; "
                             "speeds up gene database conversion")
    parser.add_argument("--reference", "-r", help="reference genome in FASTA format, "
                                                  "should be provided to compute some additional stats and "
                                                  "when raw reads are used as an input", type=str)
    parser.add_argument("--index", help="genome index for specified aligner, "
                                        "should be provided only when raw reads are used as an input", type=str)
    parser.add_argument('--clean-start', action='store_true', default=False,
                        help='Do not use previously generated index, feature db or alignments.')
    # INPUT READS
    input_args = parser.add_mutually_exclusive_group(required=True)
    input_args.add_argument('--bam', nargs='+', type=str,
                            help='sorted and indexed BAM file(s), each file will be treated as a separate sample')
    input_args.add_argument('--fastq', nargs='+', type=str,
                            help='input FASTQ file(s), each file will be treated as a separate sample; '
                                 'reference genome should be provided when using raw reads')
    input_args.add_argument('--bam_list', type=str, help='text file with list of BAM files, one file per line'
                                                         ', leave empty line between samples')
    input_args.add_argument('--fastq_list', type=str, help='text file with list of FASTQ files, one file per line'
                                                           ', leave empty line between samples')
    parser.add_argument("--data_type", "-d", type=str, required=True, choices=DATATYPE_TO_ALIGNER.keys(),
                        help="type of data to process, supported types are: " + ", ".join(DATATYPE_TO_ALIGNER.keys()))
    parser.add_argument('--stranded',  type=str, help="reads strandness type, supported values are: " +
                        ", ".join(SUPPORTED_STRANDEDNESS), default="none")
    parser.add_argument('--has_polya', action='store_true', default=False,
                        help="set if reads were not polyA trimmed; polyA tails will be detected and further "
                             " required for transcript model construction")
    parser.add_argument('--fl_data', action='store_true', default=False,
                        help="reads represent FL transcripts; both ends of the read are considered to be reliable")

    # PIPELINE AND OUTPUT
    parser.add_argument("--full_help", action='help', help="show full list of options")
    parser.add_argument("--test", action=TestMode, nargs=0, help="run IsoQuant on toy dataset")
    parser.add_argument("--threads", "-t", help="number of threads to use", type=int, default="16")

    add_additional_option("--run_aligner_only", action="store_true",
                          help="align reads to reference without isoform assignment")
    parser.add_argument('--labels', '-l', nargs='+', type=str,
                        help='sample names to be used; input file names are used if not set')
    parser.add_argument("--read_group", help="a way to group feature counts (no grouping by default): "
                                             "by BAM file tag (tag:TAG), "
                                             "using additional file (file:FILE:READ_COL:GROUP_COL:DELIM), "
                                             "using read id (read_id:DELIM)", type=str)

    parser.add_argument("--sqanti_output", help="produce SQANTI-like TSV output (requires more time)",
                        action='store_true', default=False)
    parser.add_argument("--count_exons", help="perform exon and intron counting", action='store_true', default=False)
    add_additional_option("--use_secondary", help="do not ignore secondary alignments", action='store_true', default=False)

    # ADDITIONAL OPTIONS
    add_additional_option("--aligner", help="force to use this alignment method, can be " + ", ".join(SUPPORTED_ALIGNERS) +
                                            "; chosen based on data type if not set", type=str)
    #add_additional_option("--path_to_aligner", help="folder with the aligner, $PATH is used by default", type=str)
    add_additional_option("--keep_tmp", help="do not remove temporary files in the end", action='store_true',
                          default=False)
    add_additional_option("--cage", help="bed file with CAGE peaks", type=str, default=None)
    add_additional_option("--cage-shift", type=int, default=50, help="interval before read start to look for CAGE peak")

    # ALGORITHM
    parser.add_argument("--matching_strategy", choices=["exact", "precise", "default", "loose"],
                        help="matching strategy to use from most strict to least", type=str, default=None)
    add_additional_option("--delta", type=int, default=None,
                          help="delta for inexact splice junction comparison, chosen automatically based on data type")
    add_additional_option("--correct_minor_errors", type=bool, default=None,
                          help="do not treat alignment artefacts as modification events")
    add_additional_option("--max_intron_shift", type=int, default=None,
                          help="set maximum length for intron shift")
    add_additional_option("--max_missed_exon_len", type=int, default=None,
                          help="set maximum length for skipped exon")

    # TODO: add read-type presets, transcript model contruction presets and counting
    parser.add_argument("--model_construction_strategy", choices=["reliable", "default", "fl", "all", "assembly"],
                        help="transcritp model construnction strategy to use",
                        type=str, default=None)
    add_additional_option("--report_intron_retention", type=bool, default=None,
                          help="report intron retention events in transcript model files")
    add_additional_option("--collapse_subisoforms", type=bool, default=None,
                          help="collapse isoforms whose intron chain is a subsequence of other intron chain")
    add_additional_option("--min_ref_fsm_supporting_reads", type=int, default=None,
                          help="minimal number of FSM reads that support known isoform")
    add_additional_option("--min_ref_supporting_reads", type=int, default=None,
                          help="minimal number of reads that support known isoform")
    add_additional_option("--min_novel_fsm_supporting_reads", type=int, default=None,
                          help="minimal number of FSM reads that support novel isoform")
    add_additional_option("--min_novel_supporting_reads", type=int, default=None,
                          help="minimal number of reads that support novel isoform")
    add_additional_option("--min_reads_supporting_tsts", type=int, default=None,
                          help="minimal number of reads that support isoform terminal sites")

    args = parser.parse_args(args, namespace)

    if os.path.exists(args.output):
        # logger is not defined yet
        print("WARNING! Output folder already exists, some files may be overwritten")
    else:
        os.makedirs(args.output)

    if not check_params(args):
        parser.print_usage()
        exit(-1)
    return args


# Test mode is triggered by --test option
class TestMode(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        source_dir = os.path.dirname(os.path.realpath(__file__))
        options = ['--output', 'isoquant_test', '--threads', '1',
                   '--fastq', os.path.join(source_dir, 'tests/toy_data/MAPT.Mouse.ONT.simulated.fastq'),
                   '--reference', os.path.join(source_dir, 'tests/toy_data/MAPT.Mouse.reference.fasta'),
                   '--genedb', os.path.join(source_dir, 'tests/toy_data/MAPT.Mouse.genedb.gtf'),
                   '--cage', os.path.join(source_dir, 'tests/toy_data/MAPT.Mouse.CAGE.bed'),
                   '--clean-start',
                   '--data_type', 'nanopore', '--complete_genedb']
        print('=== Running in test mode === ')
        print('Any other option is ignored ')
        main(options)
        if self._check_log():
            logger.info(' === TEST PASSED CORRECTLY === ')
        else:
            logger.error(' === TEST FAILED ===')
        parser.exit()

    @staticmethod
    def _check_log():
        with open('isoquant_test/isoquant.log', 'r') as f:
            log = f.read()

        correct_results = ['empty: 15', 'unique: 117', 'known: 10', 'Processed 1 sample']
        return all([result in log for result in correct_results])


# Check user's params
def check_params(args):
    args.input_data = InputDataStorage(args)
    if args.input_data.input_type == "fastq" and args.reference is None and args.index is None:
        print("ERROR! Reference genome or index were not provided, raw reads cannot be processed")
        return False

    if args.aligner is not None and args.aligner not in SUPPORTED_ALIGNERS:
        print("ERROR! Unsupported aligner " + args.aligner + ", choose one of: " + " ".join(SUPPORTED_ALIGNERS))
        return False
    if args.data_type is None:
        print("ERROR! Data type is not provided, choose one of " + " ".join(DATATYPE_TO_ALIGNER.keys()))
        return False
    elif args.data_type not in DATATYPE_TO_ALIGNER.keys():
        print("ERROR! Unsupported data type " + args.data_type + ", choose one of: " + " ".join(DATATYPE_TO_ALIGNER.keys()))
        return False

    if args.run_aligner_only and args.input_data.input_type == "bam":
        print("ERROR! Do not use BAM files with --run_aligner_only option.")
        return False
    if args.stranded not in SUPPORTED_STRANDEDNESS:
        print("ERROR! Unsupported strandedness " + args.stranded + ", choose one of: " + " ".join(SUPPORTED_STRANDEDNESS))
        return False

    check_input_files(args)
    return True


def check_input_files(args):
    for sample in args.input_data.samples:
        for lib in sample.file_list:
            for in_file in lib:
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
        if not os.path.isfile(args.cage):
            print("ERROR! Bed file with CAGE peaks " + args.cage + " does not exist")
            exit(-1)

    if args.genedb is not None:
        if not os.path.isfile(args.genedb):
            print("ERROR! Gene database " + args.genedb + " does not exist")
            exit(-1)


def create_output_dirs(args):
    args.tmp_dir = os.path.join(args.output, "tmp")
    if os.path.exists(args.tmp_dir):
        logger.warning("Tmp folder already exists, some files may be overwritten")
    else:
        os.makedirs(args.tmp_dir)

    for sample in args.input_data.samples:
        sample_dir = sample.out_dir
        if os.path.exists(sample_dir):
            logger.warning(sample_dir + " folder already exists, some files may be overwritten")
        else:
            os.makedirs(sample_dir)


def set_logger(args, logger_instance):
    logger_instance.setLevel(logging.INFO)
    log_file = os.path.join(args.output, "isoquant.log")
    f = open(log_file, "w")
    f.write("CMD: " + ' '.join(sys.argv) + '\n')
    f.close()
    fh = logging.FileHandler(log_file)
    # FIXME
    fh.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger_instance.addHandler(fh)
    logger_instance.addHandler(ch)


def set_data_dependent_options(args):
    matching_strategies = {ASSEMBLY: "precise", PACBIO_CCS_DATA: "precise", PACBIO_DATA: "default", NANOPORE_DATA: "default"}
    if args.matching_strategy is None:
        args.matching_strategy = matching_strategies[args.data_type]

    model_construction_strategies = {ASSEMBLY: "assembly", PACBIO_CCS_DATA: "default", PACBIO_DATA: "default", NANOPORE_DATA: "default"}
    if args.model_construction_strategy is None:
        args.model_construction_strategy = model_construction_strategies[args.data_type]
        if args.fl_data and args.model_construction_strategy == "default":
            args.model_construction_strategy = "fl"

    args.resolve_ambiguous = ExonAmbiguityResolvingMethod.full_splice_matches_only if args.fl_data else None


def set_matching_options(args):
    MatchingStrategy = namedtuple('MatchingStrategy',
                                  ('delta', 'max_exon_extension', 'max_intron_shift', 'max_missed_exon_len',
                                   'allow_extra_terminal_introns', 'apa_delta',
                                   'resolve_ambiguous', 'correct_minor_errors'))

    strategies = {
        'exact':   MatchingStrategy(0,  5,   0,   0,   False, 20,  ExonAmbiguityResolvingMethod.mono_exonic_only, False),
        'precise': MatchingStrategy(3,  10,  30,  50,  False, 50, ExonAmbiguityResolvingMethod.mono_exonic_only, True),
        'default': MatchingStrategy(6,  60,  60,  200, False, 100, ExonAmbiguityResolvingMethod.mono_exonic_only, True),
        'loose':   MatchingStrategy(12, 100, 120, 300, True, 100, ExonAmbiguityResolvingMethod.all,  True),
    }

    strategy = strategies[args.matching_strategy]

    args.delta = args.delta or strategy.delta
    args.max_exon_extension = strategy.max_exon_extension
    args.max_intron_shift = args.max_intron_shift or strategy.max_intron_shift
    args.max_missed_exon_len = args.max_missed_exon_len or strategy.max_missed_exon_len
    args.allow_extra_terminal_introns = strategy.allow_extra_terminal_introns
    args.apa_delta = strategy.apa_delta
    args.resolve_ambiguous = args.resolve_ambiguous or strategy.resolve_ambiguous
    args.correct_minor_errors = \
        strategy.correct_minor_errors if args.correct_minor_errors is None else args.correct_minor_errors

    updated_strategy = MatchingStrategy(args.delta, args.max_exon_extension, args.max_intron_shift,
                                        args.max_missed_exon_len, args.allow_extra_terminal_introns,
                                        args.apa_delta, args.resolve_ambiguous, args.correct_minor_errors)
    logger.debug('Using %s strategy. Updated strategy: %s.' % (args.matching_strategy, updated_strategy))


def set_model_construction_options(args):
    ModelConstructionStrategy = namedtuple('ModelConstructionStrategy',
                                           ('min_ref_fsm_supporting_reads', 'min_ref_supporting_reads',
                                            'min_novel_fsm_supporting_reads', 'min_novel_supporting_reads',
                                            'report_intron_retention', 'max_dist_to_isoforms_tsts',
                                            'max_dist_to_novel_tsts', 'min_reads_supporting_tsts',
                                            'collapse_subisoforms', 'count_ambiguous'))

    strategies = {
        'reliable': ModelConstructionStrategy(2, 4, 5, 9, False, 30, 100, 8, True, True),
        'default': ModelConstructionStrategy(1, 2, 4, 7, True, 30, 100, 5, True, True),
        'fl': ModelConstructionStrategy(1, 2, 3, 3, True, 15, 30, 2, False, True),
        'all': ModelConstructionStrategy(1, 1, 2, 3, True, 30, 100, 3, False, True),
        'assembly': ModelConstructionStrategy(1, 1, 1, 1, True, 30, 60, 1, True, True)
    }

    strategy = strategies[args.model_construction_strategy]

    # transcript model construction
    args.min_ref_fsm_supporting_reads = args.min_ref_fsm_supporting_reads or strategy.min_ref_supporting_reads
    args.min_ref_supporting_reads = args.min_ref_supporting_reads or strategy.min_ref_supporting_reads
    args.min_novel_fsm_supporting_reads = args.min_novel_fsm_supporting_reads or strategy.min_novel_fsm_supporting_reads
    args.min_novel_supporting_reads = args.min_novel_supporting_reads or strategy.min_novel_supporting_reads
    args.report_intron_retention = \
        strategy.report_intron_retention if args.report_intron_retention is None else args.report_intron_retention
    args.max_dist_to_isoforms_tsts = strategy.max_dist_to_isoforms_tsts
    args.max_dist_to_novel_tsts = strategy.max_dist_to_novel_tsts
    args.min_reads_supporting_tsts = args.min_reads_supporting_tsts or strategy.min_reads_supporting_tsts
    args.collapse_subisoforms = \
        strategy.collapse_subisoforms if args.collapse_subisoforms is None else args.collapse_subisoforms
    args.count_ambiguous = strategy.count_ambiguous

    args.require_polyA = args.has_polya
    args.require_cage_peak = args.cage is not None

    updated_strategy = ModelConstructionStrategy(args.min_ref_fsm_supporting_reads, args.min_ref_supporting_reads,
                                                 args.min_novel_fsm_supporting_reads, args.min_novel_fsm_supporting_reads,
                                                 args.report_intron_retention, args.max_dist_to_isoforms_tsts,
                                                 args.max_dist_to_novel_tsts, args.min_reads_supporting_tsts,
                                                 args.require_polyA, args.require_cage_peak)
    logger.debug('Using %s strategy. Updated strategy: %s.' % (args.model_construction_strategy, updated_strategy))


def set_configs_directory(args):
    config_dir = os.path.join(os.environ['HOME'], '.config', 'IsoQuant')
    os.makedirs(config_dir, exist_ok=True)

    args.db_config_path = os.path.join(config_dir, 'db_config.json')
    args.index_config_path = os.path.join(config_dir, 'index_config.json')
    args.alignment_config_path = os.path.join(config_dir, 'alignment_config.json')
    for config_path in (args.db_config_path, args.index_config_path, args.alignment_config_path):
        if not os.path.exists(config_path):
            with open(config_path, 'w') as f_out:
                json.dump({}, f_out)


def set_additional_params(args):
    set_configs_directory(args)
    set_data_dependent_options(args)
    set_matching_options(args)
    set_model_construction_options(args)

    args.print_additional_info = True
    args.memory_efficient = False

    args.indel_near_splice_site_dist = 10
    args.upstream_region_len = 20

    # TODO move to options
    args.multimap_strategy = "take_best"
    multimap_strategies = {}
    for e in MultimapResolvingStrategy:
        multimap_strategies[e.name] = e.value
    args.multimap_strategy = MultimapResolvingStrategy(multimap_strategies[args.multimap_strategy])


def run_pipeline(args):
    logger.info(" === IsoQuant pipeline started === ")

    # map reads if fastqs are provided
    if args.input_data.input_type == "fastq":
        # substitute input reads with bams
        dataset_mapper = DataSetReadMapper(args)
        args.index = dataset_mapper.index_fname
        args.input_data = dataset_mapper.map_reads(args)

    if args.run_aligner_only:
        logger.info("Isoform assignment step is skipped because --run-aligner-only option was used")
    else:
        # convert GTF/GFF if needed
        if not args.genedb.endswith('db'):
            args.genedb = convert_gtf_to_db(args)
        # run isoform assignment
        dataset_processor = DatasetProcessor(args)
        dataset_processor.process_all_samples(args.input_data)

        # aggregate counts for all samples
        if len(args.input_data.samples) > 1:
            combine_counts(args.input_data, args.output)

    logger.info(" === IsoQuant pipeline finished === ")


def clean_up(args):
    if not args.keep_tmp:
        rmtree(args.tmp_dir)


def main(args):
    args = parse_args(args)
    set_logger(args, logger)
    create_output_dirs(args)
    set_additional_params(args)
    run_pipeline(args)
    clean_up(args)


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
