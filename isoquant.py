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
    parser.add_argument('--bam', nargs='+', type=str, help='sorted and indexed BAM file(s), '
                                                            'each file will be treated as a separate sample')
    parser.add_argument('--fastq', nargs='+', type=str, help='input FASTQ file(s), '
                                                             'each file will be treated as a separate sample'
                                                             'reference genome should be provided when using raw reads')
    parser.add_argument('--bam_list', type=str, help='text file with list of BAM files, one file per line '
                                                     '(two in case of paired end reads), leave empty line between samples')
    parser.add_argument('--fastq_list', type=str, help='text file with list of FASTQ files, one file per line '
                                                       '(two in case of paired end reads), leave empty line between samples')
    parser.add_argument('--labels', nargs='+', type=str, help='sample names to be used')
    parser.add_argument("--data_type", "-d", type=str, help="type of data to process, supported types are: "
                                                            "assembly, raw_long_reads, hq_long_reads, barcoded_se_reads, barcoded_pe_reads")

    parser.add_argument("--genedb", "-g", help="gene database in gffutils .db format", type=str)
    parser.add_argument("--gtf", help="gene database in GTF/GFF format", type=str)
    parser.add_argument("--reference", help="reference genome in FASTA format, "
                                            "should be provided only when raw reads are used as an input", type=str)
    parser.add_argument("--index", help="genome index for specified aligner, "
                                        "should be provided only when raw reads are used as an input", type=str)

    parser.add_argument("--output", "-o", help="output folder, will be created automatically [./]", type=str, default="./")
    parser.add_argument("--keep_tmp", help="do not remove temporary files in the end", action='store_true', default=False)
    parser.add_argument("--prefix", help="prefix for output files", type=str, default="")
    parser.add_argument("--threads", "-t", help="number of threads to use", type=int, default="16")
    parser.add_argument("--read_info", help="text file with tab-separated information about input reads, according to "
                                            "which counts are groupped, e.g. cell type, barcode, etc.", type=str)

    ## ADDITIONAL OPTIONS
    add_additional_option("--aligner", help="force to use this alignment method, can be minimap2, star, gmap, hisat2; "
                                            "chosen based on data type if not set", type=str)
    add_additional_option("--path_to_aligner", help="folder with the aligner, $PATH is used by default", type=str)
    add_additional_option("--delta", help="delta for inexact splice junction comparison, "
                                          "chosen automatically based on data type", type=int, default="6")

    args = parser.parse_args()

    if args.output is None:
        logger.critical("Output folder was not specified")
        exit(-1)

    if os.path.exists(args.output):
        logger.warning("Output folder already exists, some files may be overwritten")
    else:
        os.makedirs(args.output)

    return args


# Check user's params
def check_params(args):
    input = list(map(lambda x: x is not None, [args.bam, args.fastq, args.bam_list, args.fastq_list]))
    if input.count(True) == 0:
        logger.critical("No input data was provided")
        exit(-1)
    elif input.count(True) > 1:
        logger.critical("Input data was provided using more than one option")
        exit(-1)

    args.input_data = InputDataStorage(args)
    if args.input_data.input_type == "fastq" and args.reference is None and args.index is None:
        logger.critical("Reference genome or index were not provided, raw reads cannot be processed")
        exit(-1)

    if args.genedb is None and args.gtf is None:
        logger.critical("Gene database was not provided")
        exit(-1)
    if args.genedb is not None and args.gtf is not None:
        logger.critical("Gene database was provided twice")
        exit(-1)

    if args.aligner is not None and args.aligner not in SUPPORTED_ALIGNERS:
        logger.critical("Unsupported aligner" + args.aligner + ", choose one of " + " ".join(SUPPORTED_ALIGNERS))
        exit(-1)
    if args.data_type is None:
        logger.critical("Data type is not provided")
        exit(-1)
    elif args.data_type not in DATATYPE_TO_ALIGNER.keys():
        logger.critical("Unsupported data type " + args.data_type + ", choose one of " + " ".join(DATATYPE_TO_ALIGNER.keys()))
        exit(-1)

    check_input_files(args)


def check_input_files(args):
    for sample in args.input_data.samples:
        for lib in sample.file_list:
            for in_file in lib:
                if not os.path.isfile(in_file):
                    raise Exception("Input file " + in_file + " does not exist")
                if args.input_data.input_type == "bam":
                    if not os.path.isfile(in_file):
                        raise Exception("BAM file " + in_file + " does not exist")
                    bamfile_in = pysam.AlignmentFile(in_file, "rb")
                    if not bamfile_in.has_index:
                        raise Exception("BAM file " + in_file + " is not indexed, run samtools index")
                    bamfile_in.close()

    if args.genedb is not None:
        if not os.path.isfile(args.genedb):
            logger.critical("Gene database " + args.genedb + " does not exist")
            raise Exception("Gene database " + args.genedb + " does not exist")
    if args.gtf is not None:
        if not os.path.isfile(args.gtf):
            logger.critical("Gene database " + args.gtf + " does not exist")
            raise Exception("Gene database " + args.gtf + " does not exist")


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
    f.close()
    fh = logging.FileHandler(log_file)
    #FIXME
    fh.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger_instnace.addHandler(fh)
    logger_instnace.addHandler(ch)


def set_additional_params(args):
    #TODO proper options
    args.print_additional_info = True
    args.skip_secondary = True

    args.resolve_ambiguous = False
    args.max_exon_extension = 300
    args.allow_extra_terminal_introns = False
    args.max_intron_shift = 100
    args.max_missed_exon_len = 200
    args.correct_minor_errors = True


def run_pipeline(args):
    logger.info(" === IsoQuant pipeline started === ")
    # convert GTF/GFF in needed
    if args.genedb is None:
        args.genedb = os.path.join(args.output, os.path.splitext(os.path.basename(args.gtf))[0] + ".db")
        logger.info("Converting gene annotation file to .db format (takes a while)...")
        gtf2db(args.gtf, args.genedb)
        logger.info("Gene database written to " + args.genedb)
        logger.info("Provide this database next time to avoid excessive conversion")

    # map reads if fastqs are provided
    if args.input_data.input_type == "fastq":
        # substitute input reads with bams
        dataset_mapper = DataSetMapper(args)
        args.index = dataset_mapper.index_path
        args.input_data = dataset_mapper.map_reads(args)

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
    check_params(args)
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
