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

import pysam
from Bio import SeqIO
import logging

from kmer_index import KmerIndexer
from common import *
from reports import *


logger = logging.getLogger('IsoQuant')


READ_CHUNK_SIZE = 10000


class DoubleBarcodeDetector:
    LINKER = "TCTTCAGCGTTCCCGAGA"
    PCR_PRIMER = "TACACGACGCTCTTCCGATCT"
    LEFT_BC_LENGTH = 8
    RIGHT_BC_LENGTH = 6
    BC_LENGTH = LEFT_BC_LENGTH + RIGHT_BC_LENGTH
    MIN_BC_LEN = BC_LENGTH - 4
    UMI_LENGTH = 9
    NON_T_UMI_BASES = 2
    UMI_LEN_DELTA = 2
    SCORE_DIFF = 1

    TERMINAL_MATCH_DELTA = 1
    STRICT_TERMINAL_MATCH_DELTA = 0

    def __init__(self, joint_barcode_list, umi_list=None, min_score=14):
        self.pcr_primer_indexer = KmerIndexer([DoubleBarcodeDetector.PCR_PRIMER], kmer_size=6)
        self.linker_indexer = KmerIndexer([DoubleBarcodeDetector.LINKER], kmer_size=5)
        self.barcode_indexer = KmerIndexer(joint_barcode_list, kmer_size=5)
        self.umi_set = None
        if umi_list:
            self.umi_set =  set(umi_list)
            logger.debug("Loaded %d UMIs" % len(umi_list))
            self.umi_indexer = KmerIndexer(umi_list, kmer_size=5)
        self.min_score = min_score

    def find_barcode_umi(self, read_id, sequence):
        read_result = self._find_barcode_umi_fwd(read_id, sequence)
        if read_result.polyT != -1:
            read_result.set_strand("+")
        if read_result.is_valid():
            return read_result

        rev_seq = reverese_complement(sequence)
        read_rev_result = self._find_barcode_umi_fwd(read_id, rev_seq)
        if read_rev_result.polyT != -1:
            read_rev_result.set_strand("-")
        if read_rev_result.is_valid():
            return read_rev_result

        return read_result if read_result.more_informative_than(read_rev_result) else read_rev_result

    def _find_barcode_umi_fwd(self, read_id, sequence):
        # look for linker in the entire read
        linker_occurrences = self.linker_indexer.get_occurrences(sequence)
        linker_start, linker_end = detect_exact_positions(sequence, 0, len(sequence),
                                                          self.linker_indexer.k, self.LINKER,
                                                          linker_occurrences, min_score=14,
                                                          start_delta=self.TERMINAL_MATCH_DELTA,
                                                          end_delta=self.TERMINAL_MATCH_DELTA)

        if linker_start is None:
            return BarcodeDetectionResult(read_id)
        logger.debug("LINKER: %d-%d" % (linker_start, linker_end))
        primer_end = -1 # forget about primer

        # use relaxed parameters once the linker is found
        presumable_polyt_start = linker_end + self.RIGHT_BC_LENGTH + self.UMI_LENGTH
        search_start = presumable_polyt_start - 4
        search_end = min(len(sequence), presumable_polyt_start + 10)
        polyt_start = find_polyt_start(sequence[search_start:search_end], window_size=5, polya_fraction=1.0)
        if polyt_start != -1:
            polyt_start += search_start

        barcode_start = linker_start - self.LEFT_BC_LENGTH
        if barcode_start < 0:
            barcode_start = 0
        barcode_end = linker_end + self.RIGHT_BC_LENGTH + 1
        if barcode_end > len(sequence):
            barcode_end = len(sequence)

        potential_barcode = sequence[barcode_start:linker_start] + sequence[linker_end + 1:barcode_end]
        logger.debug("Barcode: %s" % (potential_barcode))
        if len(potential_barcode) < self.MIN_BC_LEN:
            return BarcodeDetectionResult(read_id, linker_start=linker_start, linker_end=linker_end)
        matching_barcodes = self.barcode_indexer.get_occurrences(potential_barcode)
        barcode, bc_score, bc_start, bc_end = \
            find_candidate_with_max_score_ssw(matching_barcodes, potential_barcode,
                                              min_score=len(potential_barcode) - 1, score_diff=self.SCORE_DIFF)

        if barcode is None:
            return BarcodeDetectionResult(read_id, polyt_start, -1, linker_start, linker_end)
        logger.debug("Found: %s %d-%d" % (barcode, bc_start, bc_end))
        # position of barcode end in the reference: end of potential barcode minus bases to the alignment end
        read_barcode_end = linker_end + self.RIGHT_BC_LENGTH + 1 - (len(potential_barcode) - bc_end - 1)

        potential_umi_start = read_barcode_end + 1
        potential_umi_end = polyt_start - 1
        if polyt_start != -1 or potential_umi_end - potential_umi_start <= 5:
            potential_umi_end = potential_umi_start + self.UMI_LENGTH - 1
        potential_umi = sequence[potential_umi_start:min(potential_umi_end + 1, len(sequence))]
        logger.debug("Potential UMI: %s" % potential_umi)

        umi = None
        good_umi = False
        if self.umi_set:
            matching_umis = self.umi_indexer.get_occurrences(potential_umi)
            umi, umi_score, umi_start, umi_end = \
                find_candidate_with_max_score_ssw(matching_umis, potential_umi, min_score=7)
            logger.debug("Found UMI %s %d-%d" % (umi, umi_start, umi_end))

        if not umi :
            umi = potential_umi
            if self.UMI_LENGTH - self.UMI_LEN_DELTA <= len(umi) <= self.UMI_LENGTH + self.UMI_LEN_DELTA and \
                    all(x != "T" for x in umi[-self.NON_T_UMI_BASES:]):
                good_umi = True

        if not umi:
            return BarcodeDetectionResult(read_id, polyt_start, primer_end, linker_start, linker_end, barcode, bc_score)
        return BarcodeDetectionResult(read_id, polyt_start, primer_end, linker_start, linker_end,
                                      barcode, bc_score, umi, good_umi)


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


def process_chunk(barcodes, read_chunk, output_file, num, min_score):
    output_file += "_" + str(num)
    barcode_detector = DoubleBarcodeDetector(barcodes, min_score=min_score)
    barcode_caller = BarcodeCaller(output_file, barcode_detector)
    barcode_caller.process_chunk(read_chunk)
    return output_file


def process_single_thread(args):
    barcodes = load_barcodes(args.barcodes)
    umis = load_barcodes(args.umi) if args.umi else None
    barcode_detector = DoubleBarcodeDetector(barcodes, umi_list=umis, min_score=args.min_score)
    barcode_caller = BarcodeCaller(args.output, barcode_detector)
    barcode_caller.process(args.input)


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
    os.makedirs(tmp_dir)

    barcode_calling_gen = (
        process_chunk,
        itertools.repeat(barcodes),
        read_chunk_gen,
        itertools.repeat(os.path.join(tmp_dir, "bc")),
        itertools.count(start=0, step=1),
        itertools.repeat(args.min_score)
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
    parser.add_argument("--umi", "-u", type=str, help="potential UMIs")
    parser.add_argument("--input", "-i", type=str, help="input reads in [gzipped] FASTA, FASTQ, BAM, SAM)", required=True)
    parser.add_argument("--threads", "-t", type=int, help="threads to use", default=8)
    parser.add_argument("--min_score", type=int, help="minimal barcode score", default=14)

    args = parser.parse_args()
    return args


def main():
    #set_logger(logger)
    args = parse_args()
    set_logger(logger)
    if args.threads == 1:
        process_single_thread(args)
    else:
        process_in_parallel(args)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
