#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################
import collections
import os
import random
import sys
import argparse
import gzip
from traceback import print_exc
import itertools
import shutil
from concurrent.futures import ProcessPoolExecutor
import gc
from collections import defaultdict
import collections

import pysam
from Bio import SeqIO
import logging
from .barcode_calling.barcode_callers import (
    TenXBarcodeDetector,
    DoubleBarcodeDetector,
    IlluminaDoubleBarcodeDetector,
    BruteForceDoubleBarcodeDetector,
    ReadStats
)

logger = logging.getLogger('IsoQuant')


READ_CHUNK_SIZE = 100000
BARCODE_CALLING_MODES = {'tenX': TenXBarcodeDetector,
                         'double': DoubleBarcodeDetector,
                         'double_illumina': IlluminaDoubleBarcodeDetector,
                         'double_slow': BruteForceDoubleBarcodeDetector}


class BarcodeCaller:
    def __init__(self, output_table, barcode_detector):
        self.barcode_detector = barcode_detector
        self.output_table = output_table
        self.output_file = open(output_table, "w")
        self.output_file.write(barcode_detector.result_type().header() + "\n")
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


class FastaChunkReader():
    def __init__(self, handler, lock):
        self.handler = handler
        self.lock = lock

    def get_next(self):
        chunk = []
        while len(chunk) < READ_CHUNK_SIZE:
            with self.lock:
                r = next(self.handler, None)
            if not r: return chunk
            chunk.append((r.id, str(r.seq)))
        return chunk


class BamChunkReader():
    def __init__(self, handler, lock):
        self.handler = handler
        self.lock = lock

    def get_next(self):
        chunk = []
        while len(chunk) < READ_CHUNK_SIZE:
            with self.lock:
                r = next(self.handler, None)
            if not r: return chunk
            if r.is_secondary or r.is_supplementary:
                continue
            chunk.append((r.query_name, r.query_sequence))
        return chunk


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
        if r.is_secondary or r.is_supplementary:
            continue
        current_chunk.append((r.query_name, r.query_sequence))
        if len(current_chunk) >= READ_CHUNK_SIZE:
            yield current_chunk
            current_chunk = []
    yield current_chunk


def process_chunk(barcode_detector, read_processor, output_file, num):
    output_file += "_" + str(num)
    barcode_caller = BarcodeCaller(output_file, barcode_detector)
    barcode_caller.process_chunk(read_processor.get_next())
    return output_file


def process_single_thread(args):
    barcodes = load_barcodes(args.barcodes)
    logger.info("Loaded %d barcodes" % len(barcodes))
    logger.info("Processing " + args.input)
    barcode_detector = BARCODE_CALLING_MODES[args.mode](barcodes)
    if args.min_score:
        barcode_detector.min_score = args.min_score
    barcode_caller = BarcodeCaller(args.output, barcode_detector)
    barcode_caller.process(args.input)
    logger.info("Finished barcode calling")


def process_in_parallel(args):
    barcodes = load_barcodes(args.barcodes)
    logger.info("Loaded %d barcodes" % len(barcodes))

    input_file = args.input
    logger.info("Processing " + input_file)
    fname, outer_ext = os.path.splitext(os.path.basename(input_file))
    low_ext = outer_ext.lower()

    handle = input_file
    if low_ext in ['.gz', '.gzip']:
        handle = gzip.open(input_file, "rt")
        input_file = fname
        fname, outer_ext = os.path.splitext(os.path.basename(input_file))
        low_ext = outer_ext.lower()

    gc.disable()
    if low_ext in ['.fq', '.fastq']:
        read_chunk_gen = fastx_file_chunk_reader(SeqIO.parse(handle, "fastq"))
    elif low_ext in ['.fa', '.fasta']:
        read_chunk_gen = fastx_file_chunk_reader(SeqIO.parse(handle, "fasta"))
    elif low_ext in ['.bam', '.sam']:
        read_chunk_gen = bam_file_chunk_reader(pysam.AlignmentFile(input_file, "rb"))
    else:
        logger.error("Unknown file format " + input_file)
        exit(-1)

    tmp_dir = "barcode_calling_%x" % random.randint(0, 1 << 32)
    while os.path.exists(tmp_dir):
        tmp_dir = "barcode_calling_%x" % random.randint(0, 1 << 32)
    if args.tmp_dir:
        tmp_dir = os.path.join(args.tmp_dir, tmp_dir)
    os.makedirs(tmp_dir)

    barcode_detector = BARCODE_CALLING_MODES[args.mode](barcodes)
    if args.min_score:
        barcode_detector.min_score = args.min_score
    barcode_calling_gen = (
        process_chunk,
        itertools.repeat(barcode_detector),
        read_chunk_gen,
        itertools.repeat(os.path.join(tmp_dir, "bc")),
        itertools.count(start=0, step=1),
    )

    with ProcessPoolExecutor(max_workers=args.threads) as proc:
        output_files = proc.map(*barcode_calling_gen, chunksize=1)
    outf = open(args.output, "w")
    stat_dict = defaultdict(int)
    for tmp_file in output_files:
        shutil.copyfileobj(open(tmp_file, "r"), outf)
        for l in open(tmp_file + ".stats", "r"):
            v = l.strip().split("\t")
            if len(v) != 2:
                continue
            stat_dict[v[0]] += int(v[1])

    for k, v in stat_dict.items():
        logger.info("%s: %d" % (k, v))
    shutil.rmtree(tmp_dir)
    gc.enable()
    logger.info("Finished barcode calling")


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


def parse_args(sys_argv):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", required=True)
    parser.add_argument("--barcodes", "-b", type=str, help="barcode whitelist", required=True)
    # parser.add_argument("--umi", "-u", type=str, help="potential UMIs, detected de novo if not set")
    parser.add_argument("--mode", type=str, help="mode to be used", choices=BARCODE_CALLING_MODES.keys(),
                        default='double')
    parser.add_argument("--input", "-i", type=str, help="input reads in [gzipped] FASTA, FASTQ, BAM, SAM",
                        required=True)
    parser.add_argument("--threads", "-t", type=int, help="threads to use (16)", default=16)
    parser.add_argument("--tmp_dir", type=str, help="folder for temporary files")
    parser.add_argument("--min_score", type=int, help="minimal barcode score "
                                                      "(scoring system is +1, -1, -1, -1)", default=13)

    args = parser.parse_args(sys_argv)
    return args


def main(sys_argv):
    args = parse_args(sys_argv)
    set_logger(logger)
    if args.threads == 1:
        process_single_thread(args)
    else:
        process_in_parallel(args)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
