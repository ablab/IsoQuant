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
import numpy
import pysam
from Bio import SeqIO, Seq, SeqRecord
import logging

from src.modes import IsoQuantMode
from src.barcode_calling.common import bit_to_str, reverese_complement
from src.barcode_calling.barcode_callers import (
    TenXBarcodeDetector,
    DoubleBarcodeDetector,
    SharedMemoryStereoBarcodeDetector,
    SharedMemoryStereoSplttingBarcodeDetector,
    ReadStats, 
    VisiumHDBarcodeDetector
)

logger = logging.getLogger('IsoQuant')


READ_CHUNK_SIZE = 100000

BARCODE_CALLING_MODES = {IsoQuantMode.tenX_v3: TenXBarcodeDetector,
                         IsoQuantMode.curio: DoubleBarcodeDetector,
                         IsoQuantMode.stereoseq_nosplit: SharedMemoryStereoBarcodeDetector,
                         IsoQuantMode.stereoseq: SharedMemoryStereoSplttingBarcodeDetector,
                         IsoQuantMode.visium_5prime: TenXBarcodeDetector,
                         IsoQuantMode.visium_hd: VisiumHDBarcodeDetector
                         }

BARCODE_FILES_REQUIRED = {IsoQuantMode.tenX_v3: [1],
                          IsoQuantMode.curio: [1, 2],
                          IsoQuantMode.stereoseq_nosplit: [1],
                          IsoQuantMode.stereoseq: [1],
                          IsoQuantMode.visium_5prime: [1],
                          IsoQuantMode.visium_hd: [2]
                          }


def stats_file_name(file_name):
    return file_name + ".stats"


def get_umi_length(isoquant_mode: IsoQuantMode):
    if isoquant_mode not in BARCODE_CALLING_MODES:
        return 0
    try:
        return BARCODE_CALLING_MODES[isoquant_mode].UMI_LEN
    except AttributeError:
        return 0


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
    def __init__(self, output_file_name, barcode_detector, header=False, output_sequences=None):
        self.barcode_detector = barcode_detector
        self.output_file_name = output_file_name
        self.output_file = open(self.output_file_name, "w")
        self.output_sequences = output_sequences
        self.output_sequences_file = None
        self.process_function = self._process_read_normal
        if self.output_sequences:
            self.output_sequences_file = open(self.output_sequences, "w")
            self.process_function = self._process_read_split
        if header:
            self.output_file.write(barcode_detector.result_type().header() + "\n")
        self.read_stat = ReadStats()

    def get_stats(self):
        return self.read_stat

    def dump_stats(self, file_name=None):
        if not file_name:
            file_name = stats_file_name(self.output_file_name)
        stat_out = open(file_name, "w")
        stat_out.write(str(self.read_stat))
        stat_out.close()

    def close(self):
        self.output_file.close()
        if self.output_sequences_file:
            self.output_sequences_file.close()

    def __del__(self):
        if not self.output_file.closed:
            self.output_file.close()
        if self.output_sequences_file and not self.output_sequences_file.closed:
            self.output_sequences_file.close()

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
            self._process_bam(pysam.AlignmentFile(input_file, "rb", check_sq=False))
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
            self.process_function(read_id, seq)

    def _process_bam(self, read_handler):
        counter = 0
        for r in read_handler:
            if counter % 100 == 0:
                sys.stdout.write("Processed %d reads\r" % counter)
            counter += 1
            read_id = r.query_name
            seq = r.query_sequence
            self.process_function(read_id, seq)

    # split read and find multiple barcodes
    def _process_read_split(self, read_id, read_sequence):
        logger.debug("==== %s ====" % read_id)
        barcode_result = self.barcode_detector.find_barcode_umi(read_id, read_sequence)

        seq_records = []
        require_tso = len(barcode_result.detected_patterns) > 1
        strands = set()
        for r in barcode_result.detected_patterns:
            self.read_stat.add_read(r)
            if not r.is_valid():
                self.output_file.write("%s\n" % str(r))
                continue

            read_segment_start = max(0, r.primer - 25, r.polyT - 75)
            read_segment_end = len(read_sequence) if r.tso5 == -1 else min(len(read_sequence), r.tso5 + 25)
            r.read_id = read_id + ("_%d_%d_%s" % (read_segment_start, read_segment_end, r.strand))
            if r.strand == "+":
                new_read_seq = read_sequence[read_segment_start:read_segment_end]
            else:
                new_read_seq = reverese_complement(read_sequence)[read_segment_start:read_segment_end]
            strands.add(r.strand)
            self.output_file.write("%s\n" % str(r))
            if self.output_sequences and (not require_tso or r.tso5 != -1):
                seq_records.append(SeqRecord.SeqRecord(seq=Seq.Seq(new_read_seq), id=r.read_id, description=""))

        self.read_stat.add_custom_stats("Splits", len(barcode_result.detected_patterns))
        # self.read_stat.add_custom_stats("Splits %d %s" % (len(barcode_result.detected_patterns), "".join(list(sorted(strands)))), 1)
        if self.output_sequences_file:
            SeqIO.write(seq_records, self.output_sequences_file, "fasta")

    def _process_read_normal(self, read_id, read_sequence):
        logger.debug("==== %s ====" % read_id)
        if read_sequence is None: return
        barcode_result = self.barcode_detector.find_barcode_umi(read_id, read_sequence)

        self.output_file.write("%s\n" % str(barcode_result))
        self.read_stat.add_read(barcode_result)

    def process_chunk(self, read_chunk):
        counter = 0
        for read_id, seq in read_chunk:
            self.process_function(read_id, seq)
            counter += 1
        return counter


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


def process_chunk(barcode_detector, read_chunk, output_file, num, out_fasta=None, min_score=None):
    output_file += "_" + str(num)
    if out_fasta:
        out_fasta += "_" + str(num)
    counter = 0
    if min_score:
        barcode_detector.min_score = min_score
    barcode_caller = BarcodeCaller(output_file, barcode_detector, output_sequences=out_fasta)
    counter += barcode_caller.process_chunk(read_chunk)
    read_chunk.clear()
    barcode_caller.dump_stats()
    barcode_caller.close()

    return output_file, out_fasta, counter


def prepare_barcodes(args):
    logger.info("Using barcodes from %s" % ", ".join(args.barcodes))
    barcode_files = len(args.barcodes)
    if barcode_files not in BARCODE_FILES_REQUIRED[args.mode]:
        logger.critical("Barcode calling mode %s requires %s files, %d provided" %
                        (args.mode.name, " or ".join([str(x) for x in BARCODE_FILES_REQUIRED[args.mode]]), barcode_files))
        exit(-3)
    barcodes = []
    for bc in args.barcodes:
        barcodes.append(load_barcodes(bc, needs_iterator=args.mode.needs_barcode_iterator()))

    if len(barcodes) == 1:
        barcodes = barcodes[0]
        if not args.mode.needs_barcode_iterator():
            logger.info("Loaded %d barcodes" % len(barcodes))
    else:
        if not args.mode.needs_barcode_iterator():
            for i, bc in enumerate(barcodes):
                logger.info("Loaded %d barcodes from %s" % (len(bc), args.barcodes[i]))
        barcodes = tuple(barcodes)
    return barcodes


def process_single_thread(args):
    logger.info("Preparing barcodes indices")
    barcodes = prepare_barcodes(args)
    barcode_detector = BARCODE_CALLING_MODES[args.mode](barcodes)
    if args.min_score:
        barcode_detector.min_score = args.min_score
    barcode_caller = BarcodeCaller(args.output_tsv, barcode_detector, header=True, output_sequences=args.out_fasta)
    barcode_caller.process(args.input)
    barcode_caller.dump_stats()
    for stat_line in barcode_caller.get_stats():
        logger.info("  " + stat_line)
    barcode_caller.close()
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
        read_chunk_gen = bam_file_chunk_reader(pysam.AlignmentFile(input_file, "rb", check_sq=False))
    else:
        logger.error("Unknown file format " + input_file)
        exit(-1)

    tmp_dir = "barcode_calling_%x" % random.randint(0, 1 << 32)
    while os.path.exists(tmp_dir):
        tmp_dir = "barcode_calling_%x" % random.randint(0, 1 << 32)
    if args.tmp_dir:
        tmp_dir = os.path.join(args.tmp_dir, tmp_dir)
    os.makedirs(tmp_dir)

    barcodes = prepare_barcodes(args)
    barcode_detector = BARCODE_CALLING_MODES[args.mode](barcodes)
    logger.info("Barcode caller created")

    min_score = None
    if args.min_score:
        min_score = args.min_score

    tmp_barcode_file = os.path.join(tmp_dir, "bc")
    tmp_fasta_file = os.path.join(tmp_dir, "subreads") if args.out_fasta else None
    chunk_counter = 0
    future_results = []
    output_files = []

    with ProcessPoolExecutor(max_workers=args.threads) as proc:
        for chunk in read_chunk_gen:
            future_results.append(proc.submit(process_chunk,
                                              barcode_detector,
                                              chunk,
                                              tmp_barcode_file,
                                              chunk_counter,
                                              tmp_fasta_file,
                                              min_score))
            chunk_counter += 1
            if chunk_counter >= args.threads:
                break

        reads_left = True
        read_counter = 0
        while future_results:
            completed_features, _ = concurrent.futures.wait(future_results,
                                                            return_when=concurrent.futures.FIRST_COMPLETED)
            for c in completed_features:
                if c.exception() is not None:
                    raise c.exception()
                res = c.result()
                out_file, out_fasta, read_count = res
                read_counter += read_count
                sys.stdout.write("Processed %d reads\r" % read_counter)
                output_files.append((out_file, out_fasta))
                future_results.remove(c)
                if reads_left:
                    try:
                        chunk = next(read_chunk_gen)
                        future_results.append(proc.submit(process_chunk,
                                                          barcode_detector,
                                                          chunk,
                                                          tmp_barcode_file,
                                                          chunk_counter,
                                                          tmp_fasta_file,
                                                          min_score))
                        chunk_counter += 1
                    except StopIteration:
                        reads_left = False

    with open(args.output_tsv, "w") as final_output_tsv:
        final_output_fasta = open(args.out_fasta, "w") if args.out_fasta else None
        header = BARCODE_CALLING_MODES[args.mode].result_type().header()
        final_output_tsv.write(header + "\n")
        stat_dict = defaultdict(int)
        for tmp_file, tmp_fasta in output_files:
            shutil.copyfileobj(open(tmp_file, "r"), final_output_tsv)
            if tmp_fasta and final_output_fasta:
                shutil.copyfileobj(open(tmp_fasta, "r"), final_output_fasta)
            for l in open(stats_file_name(tmp_file), "r"):
                v = l.strip().split("\t")
                if len(v) != 2:
                    continue
                stat_dict[v[0]] += int(v[1])

        if final_output_fasta is not None:
            final_output_fasta.close()

    with open(stats_file_name(args.output_tsv), "w") as out_stats:
        for k, v in stat_dict.items():
            logger.info("  %s: %d" % (k, v))
            out_stats.write("%s\t%d\n" % (k, v))
    shutil.rmtree(tmp_dir)
    logger.info("Finished barcode calling")


def load_barcodes(inf, needs_iterator=False):
    if inf.endswith("h5") or inf.endswith("hdf5"):
        return load_h5_barcodes_bit(inf)

    if inf.endswith("gz") or inf.endswith("gzip"):
        handle = gzip.open(inf, "rt")
    else:
        handle = open(inf, "r")

    barcode_iterator = iter(l.strip().split()[0] for l in handle)
    if needs_iterator:
        return barcode_iterator

    return [b for b in barcode_iterator]


def load_h5_barcodes_bit(h5_file_path, dataset_name='bpMatrix_1'):
    raise NotImplementedError()
    import h5py
    barcode_list = []
    with h5py.File(h5_file_path, 'r') as h5_file:
        dataset = numpy.array(h5_file[dataset_name])
        for row in dataset:
            for col in row:
                barcode_list.append(bit_to_str(int(col[0])))
    return barcode_list


def set_logger(logger_instance, args):
    logger_instance.setLevel(logging.INFO)
    if args.debug:
        logger_instance.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    if args.debug:
        ch.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def parse_args(sys_argv):
    def add_hidden_option(*args, **kwargs):  # show command only with --full-help
        kwargs['help'] = argparse.SUPPRESS
        parser.add_argument(*args, **kwargs)

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", required=True)
    parser.add_argument("--barcodes", "-b", nargs='+', type=str, help="barcode whitelist(s)", required=False)
    # parser.add_argument("--umi", "-u", type=str, help="potential UMIs, detected de novo if not set")
    parser.add_argument("--mode", type=str, help="mode to be used", choices=[x.name for x in BARCODE_CALLING_MODES.keys()],
                        default=IsoQuantMode.stereoseq.name)
    parser.add_argument("--input", "-i", type=str, help="input reads in [gzipped] FASTA, FASTQ, BAM, SAM",
                        required=True)
    parser.add_argument("--threads", "-t", type=int, help="threads to use (16)", default=16)
    parser.add_argument("--tmp_dir", type=str, help="folder for temporary files")
    parser.add_argument("--min_score", type=int, help="minimal barcode score "
                                                      "(scoring system is +1, -1, -1, -1)")
    add_hidden_option('--debug', action='store_true', default=False, help='Debug log output.')

    args = parser.parse_args(sys_argv)
    args.mode = IsoQuantMode[args.mode]
    args.out_fasta = None
    args.output_tsv = None
    return args


def check_args(args):
    if args.out_fasta is None and args.mode.produces_new_fasta():
        args.out_fasta = args.output + ".split_reads.fasta"
    if args.output_tsv is None:
        args.output_tsv = args.output + ".barcoded_reads.tsv"


def main(sys_argv):
    args = parse_args(sys_argv)
    set_logger(logger, args)
    check_args(args)
    out_dir = os.path.dirname(args.output)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    if args.threads == 1 or args.mode.enforces_single_thread():
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
