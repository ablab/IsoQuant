#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import random
import sys
import argparse
import gzip
from traceback import print_exc
import itertools
import shutil
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict
from enum import Enum, unique

import pysam
from Bio import SeqIO
import logging

from src.sc.barcode_detection_methods import DoubleBarcodeDetector, TenXBarcodeDetector
from src.sc.reports import *


logger = logging.getLogger('IsoQuant')


@unique
class BarcodeCallingAlgorithm(Enum):
    tenX = 2
    curio = 3


ALGORITHMS = [BarcodeCallingAlgorithm.tenX.name, BarcodeCallingAlgorithm.curio.name]



READ_CHUNK_SIZE = 10000


class BarcodeCaller:
    def __init__(self, output_table, barcode_detector):
        self.barcode_detector = barcode_detector
        self.output_table = output_table
        self.output_file = open(output_table, "w")
        self.read_stat = ReadStats()

    def __del__(self):
        # logger.info("\n%s" % str(self.read_stat))
        stat_out = open(self.output_table + ".stats", "w")
        stat_out.write(str(self.read_stat))
        stat_out.close()
        self.output_file.close()

    def process(self, input_file):
        logger.info("Processing " + input_file)
        fname, outer_ext = os.path.splitext(os.path.basename(input_file))
        low_ext = outer_ext.lower()

        handle = input_file
        if low_ext in ['.gz', '.gzip']:
            handle = gzip.open(input_file, "rt")
            input_file = fname
            fname, outer_ext = os.path.splitext(os.path.basename(input_file))
            low_ext = outer_ext.lower()

        if low_ext in ['.fq', '.fastq']:
            self._process_fastx(SeqIO.parse(handle, "fastq"))
        elif low_ext in ['.fa', '.fasta']:
            self._process_fastx(SeqIO.parse(handle, "fasta"))
        elif low_ext in ['.bam', '.sam']:
            self._process_bam(pysam.AlignmentFile(input_file, "rb"))
        else:
            logger.error("Unknown file format " + input_file)
        logger.info("Finished " + input_file)

    def _process_fastx(self, read_handler):
        counter = 0
        for r in read_handler:
            if counter % 100 == 0:
                sys.stdout.write("Processed %d reads\r" % counter)
            counter += 1
            read_id = r.id
            seq = str(r.seq)
            self._process_read(read_id, seq)

    def _process_bam(self, read_handler):
        counter = 0
        for r in read_handler:
            if counter % 100 == 0:
                sys.stdout.write("Processed %d reads\r" % counter)
            counter += 1
            read_id = r.query_name
            seq = r.query_sequence
            self._process_read(read_id, seq)

    def _process_read(self, read_id, read_sequence):
        logger.debug("==== %s ====" % read_id)
        barcode_result = self.barcode_detector.find_barcode_umi(read_id, read_sequence)
        self.output_file.write("%s\n" % str(barcode_result))
        self.read_stat.add_read(barcode_result)

    def process_chunk(self, read_chunk):
        for read_id, seq in read_chunk:
            self._process_read(read_id, seq)


def fastx_file_chunk_reader(handler):
    current_chunk = []
    for r in handler:
        current_chunk.append((r.id, str(r.seq)))
        if len(current_chunk) >= READ_CHUNK_SIZE:
            yield current_chunk
            current_chunk = []
    yield current_chunk


def bam_file_chunk_reader(handler):
    current_chunk = []
    for r in handler:
        current_chunk.append((r.query_name, r.query_sequence))
        if len(current_chunk) >= READ_CHUNK_SIZE:
            yield current_chunk
            current_chunk = []
    yield current_chunk


def process_chunk(barcodes, read_chunk, output_file, num, min_score, method):
    output_file += "_" + str(num)
    barcode_detector = None
    if method == BarcodeCallingAlgorithm.tenX:
        barcode_detector = TenXBarcodeDetector()
    elif method == BarcodeCallingAlgorithm.curio:
        barcode_detector = DoubleBarcodeDetector(barcodes, min_score=min_score)
    else:
        logger.critical("Unknown barcode calling method")
    barcode_caller = BarcodeCaller(output_file, barcode_detector)
    barcode_caller.process_chunk(read_chunk)
    return output_file


def call_barcodes_single_thread(input_reads, barcodes, min_score, method, output):
    barcodes = load_barcodes(barcodes)
    barcode_detector = None
    if method == BarcodeCallingAlgorithm.tenX:
        barcode_detector = TenXBarcodeDetector()
    elif method == BarcodeCallingAlgorithm.curio:
        barcode_detector = DoubleBarcodeDetector(barcodes, min_score=min_score)
    else:
        logger.critical("Unknown barcode calling method")
    barcode_caller = BarcodeCaller(output, barcode_detector)
    barcode_caller.process(input_reads)


def call_barcodes_in_parallel(input_reads, barcodes, min_score, method, threads, aux, output):
    barcodes = load_barcodes(barcodes)
    logger.info("Loaded %d barcodes" % len(barcodes))

    input_file = input_reads
    logger.info("Processing " + input_file)
    fname, outer_ext = os.path.splitext(os.path.basename(input_file))
    low_ext = outer_ext.lower()

    handle = input_file
    if low_ext in ['.gz', '.gzip']:
        handle = gzip.open(input_file, "rt")
        input_file = fname
        fname, outer_ext = os.path.splitext(os.path.basename(input_file))
        low_ext = outer_ext.lower()

    if low_ext in ['.fq', '.fastq']:
        read_chunk_gen = fastx_file_chunk_reader(SeqIO.parse(handle, "fastq"))
    elif low_ext in ['.fa', '.fasta']:
        read_chunk_gen = fastx_file_chunk_reader(SeqIO.parse(handle, "fasta"))
    elif low_ext in ['.bam', '.sam']:
        read_chunk_gen = bam_file_chunk_reader(pysam.AlignmentFile(input_file, "rb"))
    else:
        logger.error("Unknown file format " + input_file)
        exit(-1)

    tmp_dir = os.path.join(aux, "barcode_calling_%x" % random.randint(0, 1 << 32))
    while os.path.exists(tmp_dir):
        tmp_dir = os.path.join(aux, "barcode_calling_%x" % random.randint(0, 1 << 32))
    os.makedirs(tmp_dir)

    barcode_calling_gen = (
        process_chunk,
        itertools.repeat(barcodes),
        read_chunk_gen,
        itertools.repeat(os.path.join(tmp_dir, "bc")),
        itertools.count(start=0, step=1),
        itertools.repeat(min_score),
        itertools.repeat(method)
    )

    with ProcessPoolExecutor(max_workers=threads) as proc:
        output_files = proc.map(*barcode_calling_gen, chunksize=1)
    outf = open(output, "w")
    stat_dict = defaultdict(int)
    for tmp_file in output_files:
        shutil.copyfileobj(open(tmp_file, "r"), outf)
        for l in open(tmp_file + ".stats", "r"):
            v = l.strip().split("\t")
            if len(v) != 2:
                continue
            stat_dict[v[0]] += int(v[1])

    logger.info(stat_dict)
    shutil.rmtree(tmp_dir)
    logger.info("Done")


def load_barcodes(inf):
    barcode_list = []
    for l in open(inf):
        barcode_list.append(l.strip().split()[0])
    return barcode_list


def set_logger(logger_instance):
    logger_instance.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def parse_args():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", required=True)
    parser.add_argument("--barcodes", "-b", type=str, help="barcode whitelist", required=True)
    # parser.add_argument("--umi", "-u", type=str, help="potential UMIs")
    parser.add_argument("--input", "-i", type=str, help="input reads in [gzipped] FASTA, FASTQ, BAM, SAM)", required=True)
    parser.add_argument("--threads", "-t", type=int, help="threads to use", default=8)
    parser.add_argument("--min_score", type=int, help="minimal barcode score", default=14)
    parser.add_argument("--method", "-m", type=str, choices=ALGORITHMS,
                        help="IsoQuant modes: " + ", ".join(ALGORITHMS) +
                             "; default:curio", default="curio")
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    set_logger(logger)
    algorithm = BarcodeCallingAlgorithm[args.method]
    if args.threads == 1:
        call_barcodes_single_thread(args.input, args.barcodes, args.min_score, algorithm, args.output)
    else:
        aux = os.path.dirname(args.output)
        call_barcodes_in_parallel(args.input, args.barcodes, args.min_score, algorithm, aux, args.output)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
