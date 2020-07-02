############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import logging
import argparse
from traceback import print_exc

from src.input_data_storage import *
from src.read_mapper import align_fasta, index_reference

logger = logging.getLogger('IsoQuant')


DATATYPE_TO_ALIGNER = {'assembly' : 'starlong', 'raw_long_reads' : 'minimap2', 'hq_long_reads' : 'starlong'}
                       #'barcoded_se_reads' : 'star', 'barcoded_pe_reads' : 'star'}

SUPPORTED_ALIGNERS = ['starlong', 'minimap2']


class DataSetMapper:
    def __init__(self, args):
        self.args = args
        self.aligner = self.choose_aligner()
        self.index_path = self.create_index(args)

    def choose_aligner(self,):
        if self.args.aligner is not None:
            return self.args.aligner
        else:
            return DATATYPE_TO_ALIGNER[self.args.data_type]

    def create_index(self, args):
        if args.index and os.path.exists(args.index):
            return args.index
        return index_reference(self.aligner, args)

    def map_input(self, input_data):
        # returns new InputDataStorage
        pass

    def map_sample(self, sample, output_folder):
        pass

    def map_reads(self, args):
        samples = []
        for sample in args.input_data.samples:
            bam_files = []
            for fastq_files in sample.file_list:
                bam_files.append([align_fasta(self.aligner, fastq_files, args)])
            samples.append(SampleData(bam_files, sample.label, sample.out_dir))
        args.input_data.samples = samples
        args.input_data.input_type = "bam"
        return args.input_data


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--fastq', nargs='+', type=str, help='input FASTQ file(s), '
                                                             'each file will be treated as a separate sample.'
                                                             ' Reference genome should be provided when using raw reads')
    parser.add_argument('--fastq_list', type=str, help='text file with list of FASTQ files, one file per line '
                                                       '(two in case of paired end reads), leave empty line between samples')
    parser.add_argument("--data_type", "-d", type=str,
                        help="type of data to process, supported types are: " + " ".join(DATATYPE_TO_ALIGNER.keys()))

    parser.add_argument("--reference", help="reference genome in FASTA format,"
                                            "should be provided only when raw reads are used as an input", type=str)
    parser.add_argument("--index", help="genome index for specified aligner,"
                                        "should be provided only when raw reads are used as an input", type=str)

    parser.add_argument("--output", "-o", help="output folder, will be created automatically", type=str)
    parser.add_argument("--threads", "-t", help="number of threads to use", type=int, default="8")
    parser.add_argument("--aligner", help="force to use this alignment method, can be minimap2, star, starlong, gmap, hisat2", type=str)
    parser.add_argument("--path_to_aligner", help="folder with the aligner, $PATH is used by default", type=str)

    args = parser.parse_args()

    if args.output is None:
        print("Output folder was not specified")
        exit(-1)

    if os.path.exists(args.output):
        print("Output folder already exists, some files may be overwritten")
    else:
        os.makedirs(args.output)

    return args


# Check user's params
def check_params(args):
    input = map(lambda x: x is not None, [args.fastq, fastq_list])
    if input.count(True) == 0:
        logger.critical("No input data was provided")
        exit(-1)
    elif input.count(True) > 1:
        logger.critical("Input data was provided using more than one option")
        exit(-1)

    args.input_data = InputDataStorage(args)
    if args.reference is None and args.index is None:
        logger.critical("Reference genome or index were not provided, raw reads cannot be processed")
        exit(-1)

    if args.aligner is not None and args.aligne not in SUPPORTED_ALIGNERS:
        logger.critical("Unsupported aligner" + args.aligner + ", choose one of " + " ".join(SUPPORTED_ALIGNERS))
        exit(-1)
    if args.data_type is None:
        logger.critical("Data type is not provided")
        exit(-1)
    elif args.data_type not in DATATYPE_TO_ALIGNER.keys():
        logger.critical("Unsupported data type " + args.data_type + ", choose one of " + " ".join(DATATYPE_TO_ALIGNER.keys()))
        exit(-1)


def set_logger(args, logger_instnace):
    logger_instnace.setLevel(logging.INFO)
    log_file = os.path.join(args.output, "isoquant_mapping.log")
    f = open(log_file, "w")
    f.write("CMD: " + ' '.join(sys.argv) + '\n')
    f.close()
    fh = logging.FileHandler(log_file)
    #FIXME
    fh.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger_instnace.addHandler(fh)
    logger_instnace.addHandler(ch)


def main():
    args = parse_args()
    set_logger(args, logger)
    check_params(args)


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
