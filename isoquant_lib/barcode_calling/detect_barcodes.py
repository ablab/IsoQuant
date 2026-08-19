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
from isoquant_lib.utils.error_codes import IsoQuantExitCode
from ..common import setup_worker_logging, _get_log_params
from .common import reverese_complement, load_barcodes
from .cell_selection import NOSEQ, CellBarcodeSelector, select_cell_barcodes
from . import (
    TenXBarcodeDetector,
    TenXv2BarcodeDetector,
    TenXSplittingBarcodeDetector,
    TenXv2SplittingBarcodeDetector,
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
    IsoQuantMode.stereoseq: SharedMemoryStereoBarcodeDetector,
    IsoQuantMode.visium_5prime: TenXBarcodeDetector,
    IsoQuantMode.visium_hd: VisiumHDBarcodeDetector,
    IsoQuantMode.custom_sc: UniversalSingleMoleculeExtractor
}

# Detectors that find several cDNA molecules per read, used when --split_molecules is on.
# Keys must match IsoQuantMode.supports_molecule_splitting().
SPLITTING_BARCODE_CALLING_MODES = {
    IsoQuantMode.tenX_v3: TenXSplittingBarcodeDetector,
    IsoQuantMode.tenX_v2: TenXv2SplittingBarcodeDetector,
    IsoQuantMode.stereoseq: SharedMemoryStereoSplittingBarcodeDetector,
    # same chemistry as 10x 3', so the same splitting detector applies
    IsoQuantMode.visium_5prime: TenXSplittingBarcodeDetector,
}

BARCODE_FILES_REQUIRED = {
    IsoQuantMode.tenX_v3: [1],
    IsoQuantMode.tenX_v2: [1],
    IsoQuantMode.curio: [1, 2],
    IsoQuantMode.stereoseq: [1],
    IsoQuantMode.visium_5prime: [1],
    IsoQuantMode.visium_hd: [2],
    IsoQuantMode.custom_sc: [0]
}


def barcode_detector_class(isoquant_mode: IsoQuantMode, split_molecules: bool):
    """Detector for a mode, splitting molecules or not."""
    if split_molecules:
        return SPLITTING_BARCODE_CALLING_MODES[isoquant_mode]
    return BARCODE_CALLING_MODES[isoquant_mode]


def stats_file_name(file_name):
    return file_name + ".stats"


def get_umi_length(isoquant_mode: IsoQuantMode):
    if isoquant_mode not in BARCODE_CALLING_MODES:
        return 0
    try:
        return BARCODE_CALLING_MODES[isoquant_mode].UMI_LEN
    except AttributeError:
        return 0


def get_barcode_length(isoquant_mode: IsoQuantMode):
    """Barcode length for modes that support graph-based correction."""
    if isoquant_mode not in BARCODE_CALLING_MODES:
        return 0
    try:
        return BARCODE_CALLING_MODES[isoquant_mode].BARCODE_LEN_10X
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
    def __init__(self, output_file_name, barcode_detector, header=False, output_sequences=None,
                 split_reads=False, barcode_counter=None):
        """
        split_reads must match the detector: a splitting detector returns one result per
        molecule and needs _process_read_split, whatever output_sequences says. The first
        pass of cell barcode detection splits but writes no FASTA, so the two cannot be
        inferred from each other.

        With barcode_counter set, detections are tallied into it instead of being written
        anywhere; output_file_name and output_sequences are then unused.
        """
        self.barcode_detector = barcode_detector
        self.barcode_counter = barcode_counter
        self.counting_only = barcode_counter is not None
        self.output_file_name = output_file_name
        self.output_file = None if self.counting_only else open(self.output_file_name, "w")
        self.output_sequences = None if self.counting_only else output_sequences
        self.output_sequences_file = None
        self.process_function = self._process_read_split if split_reads else self._process_read_normal
        if self.output_sequences:
            self.output_sequences_file = open(self.output_sequences, "w")
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

    def _emit(self, result):
        """Record one detection: tally its barcode, or write it out."""
        if self.counting_only:
            barcode = result.get_barcode()
            if barcode != NOSEQ:
                self.barcode_counter.add_barcode(barcode)
            return
        self.output_file.write("%s\n" % str(result))

    def close(self):
        if self.output_file:
            self.output_file.close()
        if self.output_sequences_file:
            self.output_sequences_file.close()

    def __del__(self):
        if self.output_file and not self.output_file.closed:
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
        valid_patterns = []
        for r in barcode_result.detected_patterns:
            self.read_stat.add_read(r)
            if not r.is_valid():
                self._emit(r)
                continue
            if self.counting_only:
                # the FASTA segment and its coordinates are output-only
                self._emit(r)
                valid_patterns.append(r)
                continue

            read_segment_start = r.get_fasta_segment_start()
            read_segment_end = r.get_fasta_segment_end(len(read_sequence))
            r.read_id = read_id + ("_%d_%d_%s" % (read_segment_start, read_segment_end, r.strand))
            if r.strand == "+":
                new_read_seq = read_sequence[read_segment_start:read_segment_end]
            else:
                new_read_seq = reverese_complement(read_sequence)[read_segment_start:read_segment_end]
            valid_patterns.append(r)
            self._emit(r)
            if self.output_sequences and (not require_tso or r.get_tso_position() != -1):
                seq_records.append(SeqRecord.SeqRecord(seq=Seq.Seq(new_read_seq), id=r.read_id, description=""))

        # Per-read split statistics
        self.read_stat.add_custom_stats("Input reads", 1)
        n_valid = len(valid_patterns)
        self.read_stat.add_custom_stats("Total cDNAs detected", n_valid)
        if n_valid == 0:
            self.read_stat.add_custom_stats("Reads with 0 cDNAs", 1)
        elif n_valid == 1:
            self.read_stat.add_custom_stats("Reads with 1 cDNA", 1)
        elif n_valid == 2:
            self.read_stat.add_custom_stats("Reads with 2 cDNAs", 1)
            orientation = valid_patterns[0].strand + valid_patterns[1].strand
            self.read_stat.add_custom_stats("2-cDNA orientation %s" % orientation, 1)
        else:
            self.read_stat.add_custom_stats("Reads with 3+ cDNAs", 1)

        tso_count = sum(1 for r in valid_patterns if r.get_tso_position() != -1)
        self.read_stat.add_custom_stats("TSO detected (valid cDNAs)", tso_count)

        if self.output_sequences_file:
            SeqIO.write(seq_records, self.output_sequences_file, "fasta")

    def _process_read_normal(self, read_id, read_sequence):
        logger.debug("==== %s ====" % read_id)
        if read_sequence is None:
            return
        barcode_result = self.barcode_detector.find_barcode_umi(read_id, read_sequence)

        self._emit(barcode_result)
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


# The detector is stored per worker process rather than passed to every task: it owns the
# whitelist index, which for a stock 10x whitelist is hundreds of megabytes. Shipping it with
# each 100k-read chunk would pickle it dozens of times per worker.
_WORKER_DETECTOR = {}


def setup_detector_worker(log_file, log_level, barcode_detector):
    setup_worker_logging(log_file, log_level)
    _WORKER_DETECTOR["detector"] = barcode_detector


def process_chunk(read_chunk, output_file, num, out_fasta=None, split_reads=False, barcode_detector=None):
    output_file += "_" + str(num)
    if out_fasta:
        out_fasta += "_" + str(num)
    counter = 0

    if barcode_detector is None:
        barcode_detector = _WORKER_DETECTOR["detector"]
    barcode_caller = BarcodeCaller(output_file, barcode_detector, output_sequences=out_fasta,
                                   split_reads=split_reads)
    counter += barcode_caller.process_chunk(read_chunk)
    read_chunk.clear()
    barcode_caller.dump_stats()
    barcode_caller.close()

    return output_file, out_fasta, counter


def count_chunk(read_chunk, split_reads=False, barcode_length=16, barcode_detector=None):
    """Worker: tally the barcodes in one chunk. Returns (counts, malformed, reads)."""
    if barcode_detector is None:
        barcode_detector = _WORKER_DETECTOR["detector"]
    selector = CellBarcodeSelector(barcode_length)
    barcode_caller = BarcodeCaller(None, barcode_detector, split_reads=split_reads,
                                   barcode_counter=selector)
    reads = barcode_caller.process_chunk(read_chunk)
    read_chunk.clear()
    # a chunk holds at most READ_CHUNK_SIZE reads, so this dict stays small enough to pickle
    return dict(selector.counts), selector.malformed_count, reads


def count_barcodes_in_reads(args, barcode_length):
    """Extract barcodes from every input file and return their counts.

    Nothing is written: the counts are all the cell barcode selection needs, and
    materialising a barcode per read would cost gigabytes of intermediate on a large run.
    """
    barcode_detector = create_barcode_caller(args)
    split_reads = bool(getattr(args, "split_molecules", False))
    selector = CellBarcodeSelector(barcode_length)
    single_thread = args.threads == 1 or args.mode.enforces_single_thread()

    for input_file in args.input:
        logger.info("Extracting barcodes from %s" % input_file)
        read_chunk_gen = open_read_chunks(input_file)
        if single_thread:
            barcode_caller = BarcodeCaller(None, barcode_detector, split_reads=split_reads,
                                           barcode_counter=selector)
            for chunk in read_chunk_gen:
                barcode_caller.process_chunk(chunk)
                chunk.clear()
            continue

        def submit(pool, chunk, num):
            return pool.submit(count_chunk, chunk, split_reads, barcode_length)

        def handle_result(result):
            counts, malformed, reads = result
            selector.merge_counts(counts, malformed)
            return reads

        run_chunks_in_parallel(read_chunk_gen, args, barcode_detector, submit, handle_result)

    logger.info("Extracted %d distinct barcodes" % len(selector.counts))
    return selector


def detect_cell_barcode_list(args, output_file, barcode_length, n_cells, n_cells_interval,
                             stats_file=None):
    """Count barcodes across the reads and write out the detected cell barcodes.

    Runs as one unit so the count table never leaves the process it was built in.
    """
    selector = count_barcodes_in_reads(args, barcode_length)
    return select_cell_barcodes(selector, output_file, args.barcodes,
                                n_cells=n_cells, n_cells_interval=n_cells_interval,
                                stats_file=stats_file)


def create_barcode_caller(args):
    split_molecules = bool(getattr(args, "split_molecules", False))
    logger.info("Creating barcode detector for mode %s%s" %
                (args.mode.name, ", splitting molecules" if split_molecules else ""))

    if not getattr(args, "whitelist_matching", True):
        # First pass of cell barcode detection: emit barcode windows verbatim so they can be
        # counted. No whitelist is needed here; it is only used afterwards, to filter the
        # candidate cell barcodes.
        if not args.mode.supports_cell_barcode_detection():
            logger.critical("Mode %s cannot extract barcodes without a whitelist" % args.mode.name)
            sys.exit(IsoQuantExitCode.INCOMPATIBLE_OPTIONS)
        return barcode_detector_class(args.mode, split_molecules)(None, whitelist_matching=False)

    if args.mode == IsoQuantMode.custom_sc:
        if not args.molecule:
            logger.critical("Custom single-cell/spatial mode requires molecule description to be provided via --molecule")
            sys.exit(IsoQuantExitCode.INCOMPATIBLE_OPTIONS)
        if not os.path.isfile(args.molecule):
            logger.critical("Molecule file %s does not exist" % args.molecule)
            sys.exit(IsoQuantExitCode.INPUT_FILE_NOT_FOUND)

        return barcode_detector_class(args.mode, split_molecules)(MoleculeStructure(open(args.molecule)))

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

    barcode_detector = barcode_detector_class(args.mode, split_molecules)(barcodes)

    return barcode_detector


def process_single_thread(args):
    barcode_detector = create_barcode_caller(args)
    split_reads = bool(getattr(args, "split_molecules", False))

    # args.input, args.output_tsv are always lists
    # args.out_fasta is a list or None
    out_fastas = args.out_fasta if args.out_fasta else [None] * len(args.input)

    for idx, (input_file, output_tsv, out_fasta) in enumerate(zip(args.input, args.output_tsv, out_fastas)):
        if len(args.input) > 1:
            logger.info("Processing file %d/%d: %s" % (idx + 1, len(args.input), input_file))
        barcode_caller = BarcodeCaller(output_tsv, barcode_detector, header=True, output_sequences=out_fasta,
                                       split_reads=split_reads)
        barcode_caller.process(input_file)
        barcode_caller.dump_stats()
        for stat_line in barcode_caller.get_stats():
            logger.info("  " + stat_line)
        barcode_caller.close()
    logger.info("Finished barcode calling")


def open_read_chunks(input_file):
    """Chunked reader for a [gzipped] FASTA/FASTQ or a BAM/SAM."""
    fname, outer_ext = os.path.splitext(os.path.basename(input_file))
    low_ext = outer_ext.lower()

    handle = input_file
    if low_ext in ['.gz', '.gzip']:
        handle = gzip.open(input_file, "rt")
        input_file = fname
        fname, outer_ext = os.path.splitext(os.path.basename(input_file))
        low_ext = outer_ext.lower()

    if low_ext in ['.fq', '.fastq']:
        return fastx_file_chunk_reader(SeqIO.parse(handle, "fastq"))
    if low_ext in ['.fa', '.fasta']:
        return fastx_file_chunk_reader(SeqIO.parse(handle, "fasta"))
    if low_ext in ['.bam', '.sam']:
        return bam_file_chunk_reader(pysam.AlignmentFile(input_file, "rb", check_sq=False))
    logger.error("Unknown file format " + input_file)
    sys.exit(IsoQuantExitCode.INVALID_FILE_FORMAT)


def run_chunks_in_parallel(read_chunk_gen, args, barcode_detector, submit, handle_result):
    """Feed read chunks to a worker pool, keeping args.threads tasks in flight.

    submit(pool, chunk, chunk_index) -> Future; handle_result(result) -> reads processed.
    The detector is shipped once per worker via the initializer rather than with every
    task: it owns the whitelist index, which is hundreds of megabytes for a stock whitelist.
    """
    # Clean up parent memory before spawning workers
    gc.collect()
    mp_context = multiprocessing.get_context('spawn')
    log_file, log_level = _get_log_params()
    # no max_tasks_per_child: recycling a worker would re-pickle the detector to its
    # replacement, and the detector is the expensive part
    with concurrent.futures.ProcessPoolExecutor(
            max_workers=args.threads,
            mp_context=mp_context,
            initializer=setup_detector_worker,
            initargs=(log_file, log_level, barcode_detector)) as proc:
        future_results = []
        chunk_counter = 0
        for chunk in read_chunk_gen:
            future_results.append(submit(proc, chunk, chunk_counter))
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
                read_counter += handle_result(c.result())
                sys.stdout.write("Processed %d reads\r" % read_counter)
                future_results.remove(c)
                if reads_left:
                    try:
                        future_results.append(submit(proc, next(read_chunk_gen), chunk_counter))
                        chunk_counter += 1
                    except StopIteration:
                        reads_left = False


def _process_single_file_in_parallel(input_file, output_tsv, out_fasta, args, barcode_detector):
    """Process a single file in parallel (internal helper function)."""
    split_reads = bool(getattr(args, "split_molecules", False))
    logger.info("Processing " + input_file)
    read_chunk_gen = open_read_chunks(input_file)

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
    output_files = []

    def submit(pool, chunk, num):
        return pool.submit(process_chunk, chunk, tmp_barcode_file, num, tmp_fasta_file, split_reads)

    def handle_result(result):
        tmp_out_file, tmp_out_fasta, read_count = result
        output_files.append((tmp_out_file, tmp_out_fasta))
        return read_count

    run_chunks_in_parallel(read_chunk_gen, args, barcode_detector, submit, handle_result)

    with open(output_tsv, "w") as final_output_tsv:
        final_output_fasta = open(out_fasta, "w") if out_fasta else None
        header = barcode_detector.header()
        final_output_tsv.write(header + "\n")
        stat_dict = defaultdict(int)
        for tmp_file, tmp_fasta in output_files:
            shutil.copyfileobj(open(tmp_file, "r"), final_output_tsv)
            if tmp_fasta and final_output_fasta:
                shutil.copyfileobj(open(tmp_fasta, "r"), final_output_fasta)
            for line in open(stats_file_name(tmp_file), "r"):
                v = line.strip().split("\t")
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
