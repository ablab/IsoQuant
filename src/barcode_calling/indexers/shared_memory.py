############################################################################
# Copyright (c) 2023-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Shared memory k-mer indexer for parallel processing of large barcode sets.

Uses multiprocessing.shared_memory to share index across worker processes.
Optimized with Numba JIT compilation when available.
"""

import gc
import logging
import math
import numpy
from collections import defaultdict
from multiprocessing import shared_memory
from typing import List, Tuple
from ..common import bit_to_str, str_to_2bit

logger = logging.getLogger('IsoQuant')

# Try to import numba for JIT compilation
try:
    from numba import njit
    NUMBA_AVAILABLE = True
except (ImportError, SystemError):
    # SystemError can occur with broken numba installations
    NUMBA_AVAILABLE = False


def _count_kmers_python(seqs, shifts, mask, kmer_counts, num_kmers_per_seq):
    """Pure Python fallback for k-mer counting."""
    for seq_idx in range(len(seqs)):
        seq = int(seqs[seq_idx])
        for k_pos in range(num_kmers_per_seq):
            kmer_idx = (seq >> int(shifts[k_pos])) & mask
            kmer_counts[kmer_idx] += 1


def _populate_index_python(seqs, shifts, mask, index_ranges, index, num_kmers_per_seq):
    """Pure Python fallback for index population."""
    total_kmers = len(index_ranges) - 1
    write_positions = numpy.zeros(total_kmers, dtype=numpy.uint64)

    for seq_idx in range(len(seqs)):
        seq = int(seqs[seq_idx])
        for k_pos in range(num_kmers_per_seq):
            kmer_idx = (seq >> int(shifts[k_pos])) & mask
            write_pos = int(index_ranges[kmer_idx] + write_positions[kmer_idx])
            index[write_pos] = seq_idx
            write_positions[kmer_idx] += 1


if NUMBA_AVAILABLE:
    @njit(cache=True)
    def _count_kmers_numba(seqs, shifts, mask, kmer_counts, num_kmers_per_seq):
        """Numba JIT-compiled k-mer counting."""
        for seq_idx in range(len(seqs)):
            seq = seqs[seq_idx]
            for k_pos in range(num_kmers_per_seq):
                kmer_idx = (seq >> shifts[k_pos]) & mask
                kmer_counts[kmer_idx] += 1

    @njit(cache=True)
    def _populate_index_numba(seqs, shifts, mask, index_ranges, index, num_kmers_per_seq):
        """Numba JIT-compiled index population - O(n) direct writes."""
        total_kmers = len(index_ranges) - 1
        write_positions = numpy.zeros(total_kmers, dtype=numpy.uint64)

        for seq_idx in range(len(seqs)):
            seq = seqs[seq_idx]
            for k_pos in range(num_kmers_per_seq):
                kmer_idx = (seq >> shifts[k_pos]) & mask
                write_pos = index_ranges[kmer_idx] + write_positions[kmer_idx]
                index[write_pos] = seq_idx
                write_positions[kmer_idx] += 1

    _count_kmers = _count_kmers_numba
    _populate_index = _populate_index_numba
else:
    _count_kmers = _count_kmers_python
    _populate_index = _populate_index_python


class SharedMemoryIndexInfo:
    """
    Metadata container for passing shared memory index info between processes.

    Stores shared memory block names and index parameters for reconstruction
    in worker processes.
    """

    def __init__(self, barcode_count: int, kmer_size: int, seq_len: int,
                 index_size: int, barcodes_sm_name: str, index_sm_name: str,
                 index_range_sm_name: str):
        self.barcode_count = barcode_count
        self.kmer_size = kmer_size
        self.seq_len = seq_len
        self.index_size = index_size
        self.barcodes_sm_name = barcodes_sm_name
        self.index_sm_name = index_sm_name
        self.index_range_sm_name = index_range_sm_name

    def __getstate__(self):
        return (self.barcode_count, self.kmer_size, self.seq_len,
                self.index_size, self.barcodes_sm_name,
                self.index_sm_name, self.index_range_sm_name)

    def __setstate__(self, state):
        self.barcode_count = state[0]
        self.kmer_size = state[1]
        self.seq_len = state[2]
        self.index_size = state[3]
        self.barcodes_sm_name = state[4]
        self.index_sm_name = state[5]
        self.index_range_sm_name = state[6]


class SharedMemoryArray2BitKmerIndexer:
    """
    Shared memory k-mer indexer for parallel processing of large barcode sets.

    Allocates index data in shared memory so worker processes can access
    the same index without copying. Uses Numba JIT compilation when available
    for fast index construction.

    Best for barcode sets > 1M where memory sharing is critical.
    """

    SEQ_DTYPE = numpy.uint64
    KMER_DTYPE = numpy.uint64
    INDEX_DTYPE = numpy.uint64

    def __init__(self, known_bin_seq: list, kmer_size: int = 12, seq_len: int = 25):
        """
        Initialize shared memory k-mer index.

        Args:
            known_bin_seq: Pre-encoded sequences as integers (use str_to_2bit)
            kmer_size: Length of k-mers
            seq_len: Length of sequences (all must be same length)
        """
        self.main_instance = True
        self.k = kmer_size
        total_kmers = int(math.pow(4, self.k))
        self.mask = (1 << (2 * self.k)) - 1
        self.seq_len = seq_len
        self.seq_mask = (1 << (2 * self.seq_len)) - 1
        self.total_sequences = len(known_bin_seq)

        # Handle empty case
        if self.total_sequences == 0:
            self.index_size = 0
            self.barcodes_shared_memory = None
            self.index_shared_memory = None
            self.index_ranges_shared_memory = None
            self.known_bin_seq = numpy.array([], dtype=self.SEQ_DTYPE)
            self.index = numpy.array([], dtype=self.KMER_DTYPE)
            self.index_ranges = numpy.zeros(total_kmers + 1, dtype=self.INDEX_DTYPE)
            return

        # Convert to numpy array for vectorized operations
        seqs = numpy.array(known_bin_seq, dtype=self.SEQ_DTYPE)

        # Precompute bit shifts for k-mer extraction
        num_kmers_per_seq = seq_len - kmer_size + 1
        shifts = numpy.array([(seq_len - kmer_size - i) * 2 for i in range(num_kmers_per_seq)],
                             dtype=numpy.uint64)

        logger.info("Building k-mer index for %d sequences (Numba: %s)",
                    self.total_sequences, "enabled" if NUMBA_AVAILABLE else "disabled")

        # Pass 1: Count k-mer occurrences (Numba JIT if available)
        kmer_counts = numpy.zeros(total_kmers, dtype=numpy.uint64)
        _count_kmers(seqs, shifts, self.mask, kmer_counts, num_kmers_per_seq)
        logger.info("K-mer counting complete")

        # Compute index size from counts
        self.index_size = int(numpy.sum(kmer_counts))

        # Allocate barcodes shared memory
        self.barcodes_shared_memory = shared_memory.SharedMemory(
            create=True, size=self.total_sequences * self.SEQ_DTYPE().nbytes)
        self.known_bin_seq = numpy.ndarray(
            shape=(self.total_sequences,), dtype=self.SEQ_DTYPE,
            buffer=self.barcodes_shared_memory.buf)
        self.known_bin_seq[:] = seqs[:]

        # Allocate index shared memory
        self.index_shared_memory = shared_memory.SharedMemory(
            create=True, size=self.index_size * self.KMER_DTYPE().nbytes)
        self.index = numpy.ndarray(
            shape=(self.index_size,), dtype=self.KMER_DTYPE,
            buffer=self.index_shared_memory.buf)

        # Allocate index_ranges shared memory
        self.index_ranges_shared_memory = shared_memory.SharedMemory(
            create=True, size=(total_kmers + 1) * self.INDEX_DTYPE().nbytes)
        self.index_ranges = numpy.ndarray(
            shape=(total_kmers + 1,), dtype=self.INDEX_DTYPE,
            buffer=self.index_ranges_shared_memory.buf)

        # Compute ranges from counts using cumsum
        self.index_ranges[0] = 0
        self.index_ranges[1:] = numpy.cumsum(kmer_counts)

        # Pass 2: Populate index using direct writes (Numba JIT if available)
        _populate_index(seqs, shifts, self.mask, self.index_ranges, self.index, num_kmers_per_seq)
        logger.info("K-mer index population complete")

        # Explicitly delete temporary arrays to free memory before forking
        del seqs, kmer_counts, shifts
        gc.collect()

    def __del__(self):
        barcodes_sm = getattr(self, 'barcodes_shared_memory', None)
        if barcodes_sm is not None:
            barcodes_sm.close()
            if getattr(self, 'main_instance', False):
                barcodes_sm.unlink()
        index_sm = getattr(self, 'index_shared_memory', None)
        if index_sm is not None:
            index_sm.close()
            if getattr(self, 'main_instance', False):
                index_sm.unlink()
        ranges_sm = getattr(self, 'index_ranges_shared_memory', None)
        if ranges_sm is not None:
            ranges_sm.close()
            if getattr(self, 'main_instance', False):
                ranges_sm.unlink()

    @classmethod
    def from_sharable_info(cls, shared_mem_index_info: SharedMemoryIndexInfo) -> 'SharedMemoryArray2BitKmerIndexer':
        """Reconstruct indexer in worker process from shared memory info."""
        kmer_index = cls.__new__(cls)
        kmer_index.main_instance = False
        kmer_index.k = shared_mem_index_info.kmer_size
        kmer_index.mask = (1 << (2 * kmer_index.k)) - 1
        kmer_index.seq_len = shared_mem_index_info.seq_len
        kmer_index.seq_mask = (1 << (2 * kmer_index.seq_len)) - 1
        kmer_index.index_size = shared_mem_index_info.index_size
        kmer_index.total_sequences = shared_mem_index_info.barcode_count
        total_kmers = int(math.pow(4, kmer_index.k))

        kmer_index.barcodes_shared_memory = shared_memory.SharedMemory(
            create=False, name=shared_mem_index_info.barcodes_sm_name)
        kmer_index.known_bin_seq = numpy.ndarray(
            shape=(kmer_index.total_sequences,),
            dtype=SharedMemoryArray2BitKmerIndexer.SEQ_DTYPE,
            buffer=kmer_index.barcodes_shared_memory.buf)

        kmer_index.index_shared_memory = shared_memory.SharedMemory(
            create=False, name=shared_mem_index_info.index_sm_name)
        kmer_index.index = numpy.ndarray(
            shape=(kmer_index.index_size,),
            dtype=kmer_index.KMER_DTYPE,
            buffer=kmer_index.index_shared_memory.buf)

        kmer_index.index_ranges_shared_memory = shared_memory.SharedMemory(
            create=False, name=shared_mem_index_info.index_range_sm_name)
        kmer_index.index_ranges = numpy.ndarray(
            shape=(total_kmers + 1,),
            dtype=kmer_index.INDEX_DTYPE,
            buffer=kmer_index.index_ranges_shared_memory.buf)

        return kmer_index

    def get_sharable_info(self) -> SharedMemoryIndexInfo:
        """Get info needed to reconstruct indexer in worker process."""
        return SharedMemoryIndexInfo(
            self.total_sequences, self.k, self.seq_len, self.index_size,
            self.barcodes_shared_memory.name, self.index_shared_memory.name,
            self.index_ranges_shared_memory.name)

    def _get_kmer_bin_indexes(self, bin_seq: int):
        """Extract k-mer indices from a 2-bit encoded sequence."""
        for i in range(self.seq_len - self.k + 1):
            yield (bin_seq >> ((self.seq_len - self.k - i) * 2)) & self.mask

    def empty(self) -> bool:
        """Check if index is empty."""
        return len(self.known_bin_seq) == 0

    def get_occurrences(self, sequence: str, max_hits: int = 0, min_kmers: int = 1,
                       hits_delta: int = 1, ignore_equal: bool = False) -> List[Tuple[str, int, List[int]]]:
        """
        Find indexed sequences with shared k-mers.

        Args:
            sequence: Query sequence (string, will be converted to 2-bit)
            max_hits: Maximum number of results (0 = unlimited)
            min_kmers: Minimum shared k-mers required
            hits_delta: Include results within this many k-mers of top hit
            ignore_equal: Skip exact matches

        Returns:
            List of (sequence_str, shared_kmer_count, kmer_positions) tuples
        """
        barcode_counts = defaultdict(int)
        barcode_positions = defaultdict(list)

        seq = str_to_2bit(sequence)
        for pos, kmer_idx in enumerate(self._get_kmer_bin_indexes(seq)):
            start_index = self.index_ranges[kmer_idx]
            end_index = self.index_ranges[kmer_idx + 1]
            for barcode_index in range(start_index, end_index):
                barcode = int(self.known_bin_seq[self.index[barcode_index]])
                barcode_counts[barcode] += 1
                barcode_positions[barcode].append(pos)

        result = []
        for barcode in barcode_counts.keys():
            count = barcode_counts[barcode]
            if count < min_kmers:
                continue
            if ignore_equal and barcode == seq:
                continue
            result.append((barcode, count, barcode_positions[barcode]))

        if not result:
            return []

        top_hits = max(result, key=lambda x: x[1])[1]
        result = filter(lambda x: x[1] >= top_hits - hits_delta, result)
        result = sorted(result, reverse=True, key=lambda x: x[1])

        if max_hits == 0:
            return [(bit_to_str(x[0], self.seq_len), x[1], x[2]) for x in result]
        return [(bit_to_str(x[0], self.seq_len), x[1], x[2]) for x in list(result)[:max_hits]]
