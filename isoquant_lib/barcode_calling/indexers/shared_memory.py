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


# Index entries are packed as little-endian uint40 (5 bytes), cutting the
# inverted-list footprint by 37.5% relative to uint64. This caps the supported
# barcode count at 2**40 - 1 (~1.1T), well above any realistic whitelist.
INDEX_BYTES_PER_ENTRY: int = 5
MAX_INDEXED_BARCODES: int = (1 << (INDEX_BYTES_PER_ENTRY * 8)) - 1


def _count_kmers_python(seqs, shifts, mask, kmer_counts, num_kmers_per_seq):
    """Pure Python fallback for k-mer counting."""
    for seq_idx in range(len(seqs)):
        seq = int(seqs[seq_idx])
        for k_pos in range(num_kmers_per_seq):
            kmer_idx = (seq >> int(shifts[k_pos])) & mask
            kmer_counts[kmer_idx] += 1


def _populate_index_python(seqs, shifts, mask, index_ranges, index_bytes, num_kmers_per_seq):
    """Pure Python fallback writing packed uint40 entries (little-endian)."""
    total_kmers = len(index_ranges) - 1
    write_positions = numpy.zeros(total_kmers, dtype=numpy.uint64)

    for seq_idx in range(len(seqs)):
        seq = int(seqs[seq_idx])
        for k_pos in range(num_kmers_per_seq):
            kmer_idx = (seq >> int(shifts[k_pos])) & mask
            write_pos = int(index_ranges[kmer_idx] + write_positions[kmer_idx])
            base = write_pos * 5
            index_bytes[base + 0] = seq_idx & 0xFF
            index_bytes[base + 1] = (seq_idx >> 8) & 0xFF
            index_bytes[base + 2] = (seq_idx >> 16) & 0xFF
            index_bytes[base + 3] = (seq_idx >> 24) & 0xFF
            index_bytes[base + 4] = (seq_idx >> 32) & 0xFF
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
    def _populate_index_numba(seqs, shifts, mask, index_ranges, index_bytes, num_kmers_per_seq):
        """Numba JIT-compiled index population: writes 5-byte little-endian entries."""
        total_kmers = len(index_ranges) - 1
        write_positions = numpy.zeros(total_kmers, dtype=numpy.uint64)

        for seq_idx in range(len(seqs)):
            seq = seqs[seq_idx]
            for k_pos in range(num_kmers_per_seq):
                kmer_idx = (seq >> shifts[k_pos]) & mask
                write_pos = index_ranges[kmer_idx] + write_positions[kmer_idx]
                base = write_pos * 5
                index_bytes[base + 0] = seq_idx & 0xFF
                index_bytes[base + 1] = (seq_idx >> 8) & 0xFF
                index_bytes[base + 2] = (seq_idx >> 16) & 0xFF
                index_bytes[base + 3] = (seq_idx >> 24) & 0xFF
                index_bytes[base + 4] = (seq_idx >> 32) & 0xFF
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
    INDEX_BYTE_DTYPE = numpy.uint8  # uint40 entries packed as 5 bytes each
    RANGES_DTYPE = numpy.uint64

    def __init__(self, known_bin_seq, kmer_size: int = 12, seq_len: int = 25):
        """
        Initialize shared memory k-mer index.

        Args:
            known_bin_seq: Pre-encoded sequences as integers (use str_to_2bit).
                Accepts a numpy.uint64 array (preferred — used without copy) or
                any sequence convertible by numpy.asarray.
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
            self.index = numpy.array([], dtype=self.INDEX_BYTE_DTYPE)
            self.index_ranges = numpy.zeros(total_kmers + 1, dtype=self.RANGES_DTYPE)
            return

        if self.total_sequences > MAX_INDEXED_BARCODES:
            raise ValueError(
                f"SharedMemoryArray2BitKmerIndexer supports at most "
                f"{MAX_INDEXED_BARCODES} sequences (uint40 packed entries); "
                f"got {self.total_sequences}"
            )

        # asarray reuses the buffer when input is already a uint64 numpy array,
        # avoiding a redundant 8-byte-per-entry copy on large whitelists.
        seqs = numpy.asarray(known_bin_seq, dtype=self.SEQ_DTYPE)

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

        # Compute index size (entries) from counts; bytes = entries * 5
        self.index_size = int(numpy.sum(kmer_counts))

        # Allocate barcodes shared memory and copy in once. After this we
        # discard the temporary `seqs` reference and use the shared-memory
        # view as the working array, so only one 8 GB-per-billion copy is
        # alive across the populate phase.
        self.barcodes_shared_memory = shared_memory.SharedMemory(
            create=True, size=self.total_sequences * self.SEQ_DTYPE().nbytes)
        self.known_bin_seq = numpy.ndarray(
            shape=(self.total_sequences,), dtype=self.SEQ_DTYPE,
            buffer=self.barcodes_shared_memory.buf)
        self.known_bin_seq[:] = seqs[:]
        del seqs

        # Allocate index shared memory as a packed byte buffer (uint40 entries)
        self.index_shared_memory = shared_memory.SharedMemory(
            create=True, size=self.index_size * INDEX_BYTES_PER_ENTRY)
        self.index = numpy.ndarray(
            shape=(self.index_size * INDEX_BYTES_PER_ENTRY,),
            dtype=self.INDEX_BYTE_DTYPE,
            buffer=self.index_shared_memory.buf)

        # Allocate index_ranges shared memory
        self.index_ranges_shared_memory = shared_memory.SharedMemory(
            create=True, size=(total_kmers + 1) * self.RANGES_DTYPE().nbytes)
        self.index_ranges = numpy.ndarray(
            shape=(total_kmers + 1,), dtype=self.RANGES_DTYPE,
            buffer=self.index_ranges_shared_memory.buf)

        # Compute ranges from counts using cumsum
        self.index_ranges[0] = 0
        self.index_ranges[1:] = numpy.cumsum(kmer_counts)

        # Pass 2: Populate index using direct writes (Numba JIT if available)
        _populate_index(self.known_bin_seq, shifts, self.mask,
                        self.index_ranges, self.index, num_kmers_per_seq)
        logger.info("K-mer index population complete")

        # Explicitly delete temporary arrays to free memory before forking
        del kmer_counts, shifts
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
            shape=(kmer_index.index_size * INDEX_BYTES_PER_ENTRY,),
            dtype=kmer_index.INDEX_BYTE_DTYPE,
            buffer=kmer_index.index_shared_memory.buf)

        kmer_index.index_ranges_shared_memory = shared_memory.SharedMemory(
            create=False, name=shared_mem_index_info.index_range_sm_name)
        kmer_index.index_ranges = numpy.ndarray(
            shape=(total_kmers + 1,),
            dtype=kmer_index.RANGES_DTYPE,
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
        index_buf = self.index
        for pos, kmer_idx in enumerate(self._get_kmer_bin_indexes(seq)):
            start_index = int(self.index_ranges[kmer_idx])
            end_index = int(self.index_ranges[kmer_idx + 1])
            for barcode_index in range(start_index, end_index):
                base = barcode_index * INDEX_BYTES_PER_ENTRY
                seq_idx = (int(index_buf[base])
                           | (int(index_buf[base + 1]) << 8)
                           | (int(index_buf[base + 2]) << 16)
                           | (int(index_buf[base + 3]) << 24)
                           | (int(index_buf[base + 4]) << 32))
                barcode = int(self.known_bin_seq[seq_idx])
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
