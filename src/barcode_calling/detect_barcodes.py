############################################################################
# Copyright (c) 2023-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Core barcode detection functionality.

This module contains the main barcode calling logic used by both
isoquant_detect_barcodes.py (standalone CLI) and isoquant.py (integrated pipeline).
"""

import concurrent.futures
import gc
import gzip
import logging
import multiprocessing
import os
import random
import shutil
import sys
from collections import defaultdict

import pysam
from Bio import SeqIO, Seq, SeqRecord

from ..modes import IsoQuantMode
from ..error_codes import IsoQuantExitCode
from ..common import setup_worker_logging, _get_log_params
from .common import reverese_complement, load_barcodes
from . import (
    TenXBarcodeDetector,
    TenXv2BarcodeDetector,
    CurioBarcodeDetector,
    SharedMemoryStereoBarcodeDetector,
    SharedMemoryStereoSplittingBarcodeDetector,
    ReadStats,
    VisiumHDBarcodeDetector,
    UniversalSingleMoleculeExtractor,
    MoleculeStructure
)

logger = logging.getLogger('IsoQuant')

READ_CHUNK_SIZE = 100000

BARCODE_CALLING_MODES = {
    IsoQuantMode.tenX_v3: TenXBarcodeDetector,
    IsoQuantMode.tenX_v2: TenXv2BarcodeDetector,
    IsoQuantMode.curio: CurioBarcodeDetector,
    IsoQuantMode.stereoseq_nosplit: SharedMemoryStereoBarcodeDetector,
    IsoQuantMode.stereoseq: SharedMemoryStereoSplittingBarcodeDetector,
    IsoQuantMode.visium_5prime: TenXBarcodeDetector,
    IsoQuantMode.visium_hd: VisiumHDBarcodeDetector,
    IsoQuantMode.custom_sc: UniversalSingleMoleculeExtractor
}

BARCODE_FILES_REQUIRED = {
    IsoQuantMode.tenX_v3: [1],
    IsoQuantMode.tenX_v2: [1],
    IsoQuantMode.curio: [1, 2],
    IsoQuantMode.stereoseq_nosplit: [1],
    IsoQuantMode.stereoseq: [1],
    IsoQuantMode.visium_5prime: [1],
    IsoQuantMode.visium_hd: [2],
    IsoQuantMode.custom_sc: [0]
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
            self.output_file.write(barcode_detector.header() + "\n")
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
        if self.output_sequences_file:
            SeqIO.write(seq_records, self.output_sequences_file, "fasta")

    def _process_read_normal(self, read_id, read_sequence):
        logger.debug("==== %s ====" % read_id)
        if read_sequence is None:
            return
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


def process_chunk(barcode_detector, read_chunk, output_file, num, out_fasta=None):
    output_file += "_" + str(num)
    if out_fasta:
        out_fasta += "_" + str(num)
    counter = 0

    barcode_caller = BarcodeCaller(output_file, barcode_detector, output_sequences=out_fasta)
    counter += barcode_caller.process_chunk(read_chunk)
    read_chunk.clear()
    barcode_caller.dump_stats()
    barcode_caller.close()

    return output_file, out_fasta, counter


def create_barcode_caller(args):
    logger.info("Creating barcode detector for mode %s" % args.mode.name)

    if args.mode == IsoQuantMode.custom_sc:
        if not args.molecule:
            logger.critical("Custom single-cell/spatial mode requires molecule description to be provided via --molecule")
            sys.exit(IsoQuantExitCode.INCOMPATIBLE_OPTIONS)
        if not os.path.isfile(args.molecule):
            logger.critical("Molecule file %s does not exist" % args.molecule)
            sys.exit(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)

        return BARCODE_CALLING_MODES[args.mode](MoleculeStructure(open(args.molecule)))

    logger.info("Using barcodes from %s" % ", ".join(args.barcodes))
    barcode_files = len(args.barcodes)
    if barcode_files not in BARCODE_FILES_REQUIRED[args.mode]:
        logger.critical("Barcode calling mode %s requires %s files, %d provided" %
                        (args.mode.name, " or ".join([str(x) for x in BARCODE_FILES_REQUIRED[args.mode]]), barcode_files))
        sys.exit(IsoQuantExitCode.BARCODE_FILE_COUNT_MISMATCH)

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

    barcode_detector = BARCODE_CALLING_MODES[args.mode](barcodes)

    return barcode_detector


def process_single_thread(args):
    barcode_detector = create_barcode_caller(args)

    # args.input, args.output_tsv are always lists
    # args.out_fasta is a list or None
    out_fastas = args.out_fasta if args.out_fasta else [None] * len(args.input)

    for idx, (input_file, output_tsv, out_fasta) in enumerate(zip(args.input, args.output_tsv, out_fastas)):
        if len(args.input) > 1:
            logger.info("Processing file %d/%d: %s" % (idx + 1, len(args.input), input_file))
        barcode_caller = BarcodeCaller(output_tsv, barcode_detector, header=True, output_sequences=out_fasta)
        barcode_caller.process(input_file)
        barcode_caller.dump_stats()
        for stat_line in barcode_caller.get_stats():
            logger.info("  " + stat_line)
        barcode_caller.close()
    logger.info("Finished barcode calling")


def _process_single_file_in_parallel(input_file, output_tsv, out_fasta, args, barcode_detector):
    """Process a single file in parallel (internal helper function)."""
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
        sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)

    tmp_dir = "barcode_calling_%x" % random.randint(0, 1 << 32)
    while os.path.exists(tmp_dir):
        tmp_dir = "barcode_calling_%x" % random.randint(0, 1 << 32)
    if args.tmp_dir:
        tmp_dir = os.path.join(args.tmp_dir, tmp_dir)
    else:
        tmp_dir = os.path.join(os.path.dirname(args.output), tmp_dir)
    os.makedirs(tmp_dir)

    tmp_barcode_file = os.path.join(tmp_dir, "bc")
    tmp_fasta_file = os.path.join(tmp_dir, "subreads") if out_fasta else None
    chunk_counter = 0
    future_results = []
    output_files = []

    # Clean up parent memory before spawning workers
    gc.collect()
    mp_context = multiprocessing.get_context('spawn')
    log_file, log_level = _get_log_params()
    # max_tasks_per_child requires Python 3.11+
    executor_kwargs = {
        'max_workers': args.threads,
        'mp_context': mp_context,
        'initializer': setup_worker_logging,
        'initargs': (log_file, log_level),
    }
    if sys.version_info >= (3, 11):
        executor_kwargs['max_tasks_per_child'] = 20
    with concurrent.futures.ProcessPoolExecutor(**executor_kwargs) as proc:
        for chunk in read_chunk_gen:
            future_results.append(proc.submit(process_chunk,
                                              barcode_detector,
                                              chunk,
                                              tmp_barcode_file,
                                              chunk_counter,
                                              tmp_fasta_file))
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
                tmp_out_file, tmp_out_fasta, read_count = res
                read_counter += read_count
                sys.stdout.write("Processed %d reads\r" % read_counter)
                output_files.append((tmp_out_file, tmp_out_fasta))
                future_results.remove(c)
                if reads_left:
                    try:
                        chunk = next(read_chunk_gen)
                        future_results.append(proc.submit(process_chunk,
                                                          barcode_detector,
                                                          chunk,
                                                          tmp_barcode_file,
                                                          chunk_counter,
                                                          tmp_fasta_file))
                        chunk_counter += 1
                    except StopIteration:
                        reads_left = False

    with open(output_tsv, "w") as final_output_tsv:
        final_output_fasta = open(out_fasta, "w") if out_fasta else None
        header = barcode_detector.header()
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

    with open(stats_file_name(output_tsv), "w") as out_stats:
        for k, v in stat_dict.items():
            logger.info("  %s: %d" % (k, v))
            out_stats.write("%s\t%d\n" % (k, v))
    shutil.rmtree(tmp_dir)


def process_in_parallel(args):
    """Process input files in parallel."""
    # args.input, args.output_tsv are always lists
    # args.out_fasta is a list or None
    out_fastas = args.out_fasta if args.out_fasta else [None] * len(args.input)
    barcode_detector = create_barcode_caller(args)

    for idx, (input_file, output_tsv, out_fasta) in enumerate(zip(args.input, args.output_tsv, out_fastas)):
        if len(args.input) > 1:
            logger.info("Processing file %d/%d: %s" % (idx + 1, len(args.input), input_file))
        _process_single_file_in_parallel(input_file, output_tsv, out_fasta, args, barcode_detector)

    logger.info("Finished barcode calling")
