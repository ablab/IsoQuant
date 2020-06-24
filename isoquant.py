#!/usr/bin/env python
#
# ############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import logging
import argparse
from traceback import print_exc
from collections import namedtuple

import gffutils
import pysam
from Bio import SeqIO

from src.input_data_storage import *
from src.gtf2db import *
from src.map_input import *
from src.dataset_processor import *

logger = logging.getLogger('IsoQuant')


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    show_full_help = '--full-help' in sys.argv

    def add_additional_option(*args, **kwargs):  # show command only with --full-help
        if not show_full_help:
            kwargs['help'] = argparse.SUPPRESS
        parser.add_argument(*args, **kwargs)

    parser.add_argument("--full-help", action='help', help="show full list of options")

    input_args = parser.add_mutually_exclusive_group(required=True)
    input_args.add_argument('--bam', nargs='+', type=str, help='sorted and indexed BAM file(s), '
                                                            'each file will be treated as a separate sample')
    input_args.add_argument('--fastq', nargs='+', type=str, help='input FASTQ file(s), '
                                                             'each file will be treated as a separate sample'
                                                             'reference genome should be provided when using raw reads')
    input_args.add_argument('--bam_list', type=str, help='text file with list of BAM files, one file per line'
                                                     ', leave empty line between samples')
    input_args.add_argument('--fastq_list', type=str, help='text file with list of FASTQ files, one file per line'
                                                       ', leave empty line between samples')
    parser.add_argument("--data_type", "-d", type=str, required=True,
                        help="type of data to process, supported types are: " + " ".join(DATATYPE_TO_ALIGNER.keys()))
    parser.add_argument("--genedb", "-g", help="gene database in gffutils .db format or GTF/GFF format", type=str,
                        required='--run_aligner_only' not in sys.argv)

    parser.add_argument("--reference", help="reference genome in FASTA format, "
                                            "should be provided to compute some additional stats and"
                                            "when raw reads are used as an input", type=str)
    parser.add_argument("--index", help="genome index for specified aligner, "
                                        "should be provided only when raw reads are used as an input", type=str)

    parser.add_argument("--run_aligner_only", action="store_true", help="align reads to reference without isoform assignment")
    parser.add_argument("--threads", "-t", help="number of threads to use", type=int, default="16")
    parser.add_argument("--output", "-o", help="output folder, will be created automatically [default=isoquant_output]",
                        type=str, default="isoquant_output")
    parser.add_argument('--labels', nargs='+', type=str, help='sample names to be used')

    parser.add_argument("--keep_tmp", help="do not remove temporary files in the end", action='store_true', default=False)
    parser.add_argument("--prefix", help="prefix for output files", type=str, default="")
    parser.add_argument("--read_group", help="a way to groups feature counts (no grouping by default): "
                                             "by BAM file tag (tag:TAG),"
                                             "using information table (file:FILE:READ_COL:GROUP_COL:DELIM),"
                                             "using read id (read_id:DELIM)", type=str)
    parser.add_argument("--matching-strategy", choices=["exact", "precise", "default", "loose"],
                        help="matching strategy to use from most strict to least", type=str, default="default")
    parser.add_argument("--no-sqanti-output", help="produce SQANTI-like TSV output (saves time)", action='store_true', default=True)
    parser.add_argument("--count-exons", help="perform exont and intron counting", action='store_true', default=False)

    ## ADDITIONAL OPTIONS
    add_additional_option("--aligner", help="force to use this alignment method, can be " + ", ".join(SUPPORTED_ALIGNERS) +
                                            "chosen based on data type if not set", type=str)
    add_additional_option("--path_to_aligner", help="folder with the aligner, $PATH is used by default", type=str)
    add_additional_option("--delta", type=int, default=None,
                          help="delta for inexact splice junction comparison, chosen automatically based on data type")
    add_additional_option("--max-exon-extension", type=int, default=None, help="set maximum length for exon elongation")
    add_additional_option("--max-intron-shift", type=int, default=None, help="set maximum length for intron shift")
    add_additional_option("--max-missed-exon-len", type=int, default=None, help="set maximum length for skipped exon")

    args = parser.parse_args()

    if os.path.exists(args.output):
        # logger is not defined yet
        print("WARNING! Output folder already exists, some files may be overwritten")
    else:
        os.makedirs(args.output)

    if not check_params(args):
        parser.print_usage()
        exit(-1)
    return args


# Check user's params
def check_params(args):
    args.input_data = InputDataStorage(args)
    if args.input_data.input_type == "fastq" and args.reference is None and args.index is None:
        print("ERROR! Reference genome or index were not provided, raw reads cannot be processed")
        return False

    if args.aligner is not None and args.aligner not in SUPPORTED_ALIGNERS:
        print("ERROR! Unsupported aligner" + args.aligner + ", choose one of: " + " ".join(SUPPORTED_ALIGNERS))
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


def set_logger(args, logger_instnace):
    logger_instnace.setLevel(logging.DEBUG)
    log_file = os.path.join(args.output, "isoquant.log")
    f = open(log_file, "w")
    f.write("CMD: " + ' '.join(sys.argv) + '\n')
    f.close()
    fh = logging.FileHandler(log_file)
    #FIXME
    fh.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger_instnace.addHandler(fh)
    logger_instnace.addHandler(ch)


def set_additional_params(args):
    Strategy = namedtuple('Strategy',
                          ('delta', 'max_exon_extension', 'max_intron_shift', 'max_missed_exon_len',
                           'allow_extra_terminal_introns', 'resolve_ambiguous', 'correct_minor_errors'))

    strategies = {
        'exact':   Strategy(0,  5,   0,   0,   False, ExonAmbiguityResolvingMethod.mono_exonic_only, False),
        'precise': Strategy(3,  10,  30,  50,  False, ExonAmbiguityResolvingMethod.mono_exonic_only, True),
        'default': Strategy(6,  30,  60,  200, False, ExonAmbiguityResolvingMethod.mono_exonic_only, True),
        'loose':   Strategy(12, 100, 200, 300, True,  ExonAmbiguityResolvingMethod.all,  True),
    }

    strategy = strategies[args.matching_strategy]

    args.delta = args.delta or strategy.delta
    args.max_exon_extension = args.max_exon_extension or strategy.max_exon_extension
    args.max_intron_shift = args.max_intron_shift or strategy.max_intron_shift
    args.max_missed_exon_len = args.max_missed_exon_len or strategy.max_missed_exon_len
    args.allow_extra_terminal_introns = strategy.allow_extra_terminal_introns
    args.resolve_ambiguous = strategy.resolve_ambiguous
    args.correct_minor_errors = strategy.correct_minor_errors

    updated_strategy = Strategy(args.delta, args.max_exon_extension, args.max_intron_shift, args.max_missed_exon_len,
                                args.allow_extra_terminal_introns, args.resolve_ambiguous, args.correct_minor_errors)
    logger.debug(f'Using {args.matching_strategy} strategy. Updated strategy: {updated_strategy}.')

    # TODO proper options
    args.print_additional_info = True
    args.skip_secondary = True

    args.memory_efficient = False

    args.indel_near_splice_site_dist = 10
    args.upstream_region_len = 20
    args.apa_delta = 20

    # trascript model construction
    args.min_ref_fsm_supporting_reads = 1
    args.min_ref_supporting_reads = 2
    args.min_novel_fsm_supporting_reads = 2
    args.min_novel_supporting_reads = 3
    args.report_intron_retention = True
    args.max_dist_to_isoforms_tsts = 20
    args.max_dist_to_novel_tsts = 50
    args.min_reads_supporting_tsts = 3
    args.require_polyA = True
    args.require_cage_peak = False # TODO
    args.collapse_subisoforms = True
    args.count_ambiguous = True

    # TODO move to options
    args.matching_stategy = "take_best"
    multimap_strategies = {}
    for e in MultimapResolvingStrategy:
        multimap_strategies[e.name] = e.value
    args.matching_stategy = MultimapResolvingStrategy(multimap_strategies[args.matching_stategy])


def run_pipeline(args):
    logger.info(" === IsoQuant pipeline started === ")
    # convert GTF/GFF in needed
    if not args.run_aligner_only and not args.genedb.endswith('db'):
        genedb_filename = args.genedb
        args.genedb = os.path.join(args.output, os.path.splitext(os.path.basename(genedb_filename))[0] + ".db")
        logger.info("Converting gene annotation file to .db format (takes a while)...")
        gtf2db(genedb_filename, args.genedb)
        logger.info("Gene database written to " + args.genedb)
        logger.info("Provide this database next time to avoid excessive conversion")

    # map reads if fastqs are provided
    if args.input_data.input_type == "fastq":
        # substitute input reads with bams
        dataset_mapper = DataSetMapper(args)
        args.index = dataset_mapper.index_path
        args.input_data = dataset_mapper.map_reads(args)

    if args.run_aligner_only:
        logger.info("Isoform assignment step is skipped because --run-aligner-only option was used")
    else:
        # run isoform assignment
        dataset_processor = DatasetProcessor(args)
        dataset_processor.process_all_samples(args.input_data)
    logger.info(" === IsoQuant pipeline finished === ")

def clean_up(args):
    #TODO
    pass


def main():
    args = parse_args()
    set_logger(args, logger)
    create_output_dirs(args)
    set_additional_params(args)
    run_pipeline(args)
    clean_up(args)


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
