#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################
import concurrent.futures
import os
import random
import sys
import argparse
import gzip
from traceback import print_exc
import shutil
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict
# import h5py
import numpy

import pysam
from Bio import SeqIO
import logging
from src.barcode_calling.common import bit_to_str, str_to_2bit
from src.barcode_calling.barcode_callers import (
    TenXBarcodeDetector,
    DoubleBarcodeDetector,
    IlluminaDoubleBarcodeDetector,
    BruteForceDoubleBarcodeDetector,
    StereoBarcodeDetector,
    StereoBarcodeDetectorTSO,
    StereoBarcodeDetectorPC,
    ReadStats, DoubleBarcodeDetectionResult
)

logger = logging.getLogger('IsoQuant')


READ_CHUNK_SIZE = 100000
BARCODE_CALLING_MODES = {'tenX': TenXBarcodeDetector,
                         'double': DoubleBarcodeDetector,
                         'double_illumina': IlluminaDoubleBarcodeDetector,
                         'double_slow': BruteForceDoubleBarcodeDetector,
                         'stereo_tso': StereoBarcodeDetectorTSO,
                         'stereo_pc': StereoBarcodeDetectorPC}


class SimpleReadStorage:
    def __init__(self):
        self.read_ids = []
        self.sequences = []

    def add(self, read_id, seq):
        self.read_ids.append(read_id)
        self.sequences.append(seq)

    def clear(self):
        self.read_ids.clear()
        self.sequences.clear()

    def __len__(self):
        return len(self.read_ids)

    def __iter__(self):
        for i in range(len(self.read_ids)):
            yield self.read_ids[i], self.sequences[i]

    def __getstate__(self):
        return self.read_ids, self.sequences

    def __setstate__(self, state):
        self.read_ids = state[0]
        self.sequences = state[1]


class BarcodeCaller:
    def __init__(self, output_table, barcode_detector, header=False):
        self.barcode_detector = barcode_detector
        self.output_table = output_table
        self.output_file = open(output_table, "w")
        if header:
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

    # standard method, finds only a single barcode
    def _process_read_sinlge(self, read_id, read_sequence):
        logger.debug("==== %s ====" % read_id)
        barcode_result = self.barcode_detector.find_barcode_umi(read_id, read_sequence)
        self.output_file.write("%s\n" % str(barcode_result))
        if isinstance(barcode_result, list):
            barcode_result = barcode_result[0]
        self.read_stat.add_read(barcode_result)

    def _process_read(self, read_id, read_sequence):
        logger.debug("==== %s ====" % read_id)
        barcode_result = self.barcode_detector.find_barcode_umi(read_id, read_sequence)
        if isinstance(barcode_result, list):
            for r in barcode_result:
                self.output_file.write("%s\n" % str(r))
                self.read_stat.add_read(r)
            self.read_stat.additional_attributes_counts["Reads processed"] += 1
            if any(r.polyT != -1 for r in barcode_result):
                self.read_stat.additional_attributes_counts["Has polyT"] += 1
            if any(r.linker_start != -1 for r in barcode_result):
                self.read_stat.additional_attributes_counts["Has linker"] += 1
            if len(barcode_result) > 1 and isinstance(barcode_result[0], DoubleBarcodeDetectionResult):
                strands = set([r.strand for r in barcode_result])
                linkers = set([r.linker_start for r in barcode_result])
                if len(strands) == 1:
                    self.read_stat.additional_attributes_counts["Multiple polyAs (same strand)"] += 1
                    if len(linkers) > 1:
                        self.read_stat.additional_attributes_counts["Multiple linkers (same strand)"] += 1
                else:
                    self.read_stat.additional_attributes_counts["Multiple polyAs (opposite strand)"] += 1
                    if len(linkers) > 1:
                        self.read_stat.additional_attributes_counts["Multiple linkers (opposite strand)"] += 1

        else:
            self.output_file.write("%s\n" % str(barcode_result))
            self.read_stat.add_read(barcode_result)

    def process_chunk(self, read_chunk):
        for read_id, seq in read_chunk:
            self._process_read(read_id, seq)


def fastx_file_chunk_reader(handler):
    current_chunk = SimpleReadStorage()
    for r in handler:
        current_chunk.add(r.id, str(r.seq))
        if len(current_chunk) >= READ_CHUNK_SIZE:
            yield current_chunk
            current_chunk = SimpleReadStorage()
    yield current_chunk


def bam_file_chunk_reader(handler):
    current_chunk = SimpleReadStorage()
    for r in handler:
        if r.is_secondary or r.is_supplementary:
            continue
        current_chunk.add(r.query_name, r.query_sequence)
        if len(current_chunk) >= READ_CHUNK_SIZE:
            yield current_chunk
            current_chunk = SimpleReadStorage()
    yield current_chunk


def process_chunk(barcode_detector, read_chunk, output_file, num, min_score=None):
    output_file += "_" + str(num)
    if min_score:
        barcode_detector.min_score = min_score
    barcode_caller = BarcodeCaller(output_file, barcode_detector)
    barcode_caller.process_chunk(read_chunk)
    read_chunk.clear()
    return output_file


def process_single_thread(args):
    logger.info("Loading barcodes from %s" % args.barcodes)
    barcodes = load_barcodes_iter(args.barcodes)
    # logger.info("Loaded %d barcodes" % len(barcodes))
    logger.info("Preparing barcodes")
    barcode_detector = BARCODE_CALLING_MODES[args.mode](barcodes)
    if args.min_score:
        barcode_detector.min_score = args.min_score
    barcode_caller = BarcodeCaller(args.output, barcode_detector, header=True)
    logger.info("Processing " + args.input)
    barcode_caller.process(args.input)
    logger.info("Finished barcode calling")


def process_in_parallel(args):
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

    tmp_barcode_file = os.path.join(tmp_dir, "bc")
    count = 0
    future_results = []
    output_files = []

    # TODO substitute with a more elegant version
    # def lazy_multiprocessing_map(func, iterable, processes=None):
    #     with ProcessPoolExecutor(max_workers=processes) as executor:
    #         futures = []
    #         iterable = iter(iterable)
    #
    #         # Submit initial batch of tasks
    #         try:
    #             for _ in range(executor._max_workers):
    #                 item = next(iterable)
    #                 futures.append(executor.submit(func, item))
    #         except StopIteration:
    #             pass
    #
    #         while futures:
    #             # Wait for the next completed future
    #             for future in as_completed(futures):
    #                 futures.remove(future)
    #                 try:
    #                     yield future.result()
    #                 except Exception as e:
    #                     yield f"Error for input {future}: {e}"
    #                 break  # Process only one completed future at a time
    #
    #             # Submit the next item from the iterable
    #             try:
    #                 item = next(iterable)
    #                 futures.append(executor.submit(func, item))
    #             except StopIteration:
    #                 pass

    logger.info("Loading barcodes from %s" % args.barcodes)
    barcodes = load_barcodes_iter(args.barcodes)
    # logger.info("Loaded %d barcodes" % len(barcodes))
    barcode_detector = BARCODE_CALLING_MODES[args.mode](barcodes)
    logger.info("Barcode caller created")

    min_score = None
    if args.min_score:
        min_score = args.min_score

    with ProcessPoolExecutor(max_workers=args.threads) as proc:
        for chunk in read_chunk_gen:
            future_results.append(proc.submit(process_chunk, barcode_detector, chunk, tmp_barcode_file, count, min_score))
            count += 1
            if count >= args.threads:
                break

        reads_left = True
        while reads_left:
            completed_features, _ = concurrent.futures.wait(future_results, return_when=concurrent.futures.FIRST_COMPLETED)
            for c in completed_features:
                if c.exception() is not None:
                    raise c.exception()
                future_results.remove(c)
                output_files.append(c.result())
                if reads_left:
                    try:
                        chunk = next(read_chunk_gen)
                        future_results.append(proc.submit(process_chunk, barcode_detector, chunk, tmp_barcode_file, count, args.min_score))
                        count += 1
                    except StopIteration:
                        reads_left = False

        completed_features, _ = concurrent.futures.wait(future_results, return_when=concurrent.futures.ALL_COMPLETED)
        for c in completed_features:
            if c.exception() is not None:
                raise c.exception()
            output_files.append(c.result())

    outf = open(args.output, "w")
    header = BARCODE_CALLING_MODES[args.mode].result_type().header()
    outf.write(header + "\n")
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
    logger.info("Finished barcode calling")


def load_barcodes(inf):
    barcode_list = []
    if inf.endswith("h5") or inf.endswith("hdf5"):
        barcode_list = load_h5_barcodes_bit(inf)
    else:
        if inf.endswith("gz") or inf.endswith("gzip"):
            handle = gzip.open(inf, "rt")
        else:
            handle = open(inf, "r")
        for l in handle:
            barcode_list.append(l.strip().split()[0])
    return barcode_list


def load_h5_barcodes_bit(h5_file_path, dataset_name='bpMatrix_1'):
    barcode_list = []
    with h5py.File(h5_file_path, 'r') as h5_file:
        dataset = numpy.array(h5_file[dataset_name])
        for row in dataset:
            for col in row:
                barcode_list.append(bit_to_str(int(col[0])))
    return barcode_list


def load_barcodes_iter(inf):
    if inf.endswith("gz") or inf.endswith("gzip"):
        handle = gzip.open(inf, "rt")
    else:
        handle = open(inf, "r")
    for l in handle:
        yield l.strip().split()[0]


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
    parser.add_argument("--barcodes", "-b", type=str, help="barcode whitelist", required=False)
    # parser.add_argument("--umi", "-u", type=str, help="potential UMIs, detected de novo if not set")
    parser.add_argument("--mode", type=str, help="mode to be used", choices=BARCODE_CALLING_MODES.keys(),
                        default='double')
    parser.add_argument("--input", "-i", type=str, help="input reads in [gzipped] FASTA, FASTQ, BAM, SAM",
                        required=True)
    parser.add_argument("--threads", "-t", type=int, help="threads to use (16)", default=16)
    parser.add_argument("--tmp_dir", type=str, help="folder for temporary files")
    parser.add_argument("--min_score", type=int, help="minimal barcode score "
                                                      "(scoring system is +1, -1, -1, -1)", default=22)

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
