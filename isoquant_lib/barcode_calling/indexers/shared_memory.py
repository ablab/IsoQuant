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


# ---------------------------------------------------------------------------
# Sparse 3-anchor helpers (first / last / minimizer per sequence).
#
# The sparse indexer stores three anchors per barcode keyed in the same
# inverted-list layout (uint40 entries, 4^k bucket table). The minimizer is
# selected by argmin of splitmix64 over the *internal* window starts only:
# excluding p=0 and p=L-k guarantees the minimizer's k-mer is a different
# slice from the first/last anchors.
# ---------------------------------------------------------------------------

NUM_ANCHORS: int = 3
HAMMING_LOW_BIT_MASK: int = 0x5555_5555_5555_5555

# Anchor k-mer values are stored as uint32 in the transient (N, 3) buffer
# during index construction. 4**k must fit in uint32, i.e. k <= 16.
MAX_KMER_SIZE_FOR_UINT32_ANCHORS: int = 16


def _splitmix64_python(x: int) -> int:
    mask = 0xFFFFFFFFFFFFFFFF
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & mask
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & mask
    return (x ^ (x >> 31)) & mask


def _popcount64_python(x: int) -> int:
    x = x - ((x >> 1) & 0x5555555555555555)
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F
    return (x * 0x0101010101010101) >> 56


def _hamming_2bit_python(a: int, b: int) -> int:
    diff = a ^ b
    return _popcount64_python((diff | (diff >> 1)) & HAMMING_LOW_BIT_MASK)


def _compute_anchors_python(seqs, first_shift, last_shift, k_mask,
                            mid_shifts, anchors_out):
    n_mid = len(mid_shifts)
    sentinel = 0xFFFFFFFFFFFFFFFF
    for i in range(len(seqs)):
        seq = int(seqs[i])
        anchors_out[i, 0] = (seq >> int(first_shift)) & int(k_mask)
        anchors_out[i, 1] = (seq >> int(last_shift)) & int(k_mask)
        best_hash = sentinel
        best_kmer = 0
        for j in range(n_mid):
            kmer = (seq >> int(mid_shifts[j])) & int(k_mask)
            h = _splitmix64_python(kmer)
            if h < best_hash:
                best_hash = h
                best_kmer = kmer
        anchors_out[i, 2] = best_kmer


def _count_anchors_python(anchors, kmer_counts):
    for i in range(anchors.shape[0]):
        kmer_counts[anchors[i, 0]] += 1
        kmer_counts[anchors[i, 1]] += 1
        kmer_counts[anchors[i, 2]] += 1


def _populate_anchors_python(anchors, index_ranges, index_bytes):
    total_kmers = len(index_ranges) - 1
    write_positions = numpy.zeros(total_kmers, dtype=numpy.uint64)
    for i in range(anchors.shape[0]):
        for slot in range(NUM_ANCHORS):
            kmer_idx = int(anchors[i, slot])
            wp = int(index_ranges[kmer_idx] + write_positions[kmer_idx])
            base = wp * 5
            index_bytes[base + 0] = i & 0xFF
            index_bytes[base + 1] = (i >> 8) & 0xFF
            index_bytes[base + 2] = (i >> 16) & 0xFF
            index_bytes[base + 3] = (i >> 24) & 0xFF
            index_bytes[base + 4] = (i >> 32) & 0xFF
            write_positions[kmer_idx] += 1


if NUMBA_AVAILABLE:
    @njit(cache=True)
    def _splitmix64_numba(x):
        x = (x ^ (x >> numpy.uint64(30))) * numpy.uint64(0xBF58476D1CE4E5B9)
        x = (x ^ (x >> numpy.uint64(27))) * numpy.uint64(0x94D049BB133111EB)
        return x ^ (x >> numpy.uint64(31))

    @njit(cache=True)
    def _compute_anchors_numba(seqs, first_shift, last_shift, k_mask,
                               mid_shifts, anchors_out):
        n_mid = len(mid_shifts)
        sentinel = numpy.uint64(0xFFFFFFFFFFFFFFFF)
        for i in range(len(seqs)):
            seq = seqs[i]
            anchors_out[i, 0] = numpy.uint32((seq >> first_shift) & k_mask)
            anchors_out[i, 1] = numpy.uint32((seq >> last_shift) & k_mask)
            best_hash = sentinel
            best_kmer = numpy.uint64(0)
            for j in range(n_mid):
                kmer = (seq >> mid_shifts[j]) & k_mask
                h = _splitmix64_numba(kmer)
                if h < best_hash:
                    best_hash = h
                    best_kmer = kmer
            anchors_out[i, 2] = numpy.uint32(best_kmer)

    @njit(cache=True)
    def _count_anchors_numba(anchors, kmer_counts):
        for i in range(anchors.shape[0]):
            kmer_counts[anchors[i, 0]] += 1
            kmer_counts[anchors[i, 1]] += 1
            kmer_counts[anchors[i, 2]] += 1

    @njit(cache=True)
    def _populate_anchors_numba(anchors, index_ranges, index_bytes):
        total_kmers = len(index_ranges) - 1
        write_positions = numpy.zeros(total_kmers, dtype=numpy.uint64)
        for i in range(anchors.shape[0]):
            seq_idx = numpy.uint64(i)
            for slot in range(3):
                kmer_idx = anchors[i, slot]
                wp = index_ranges[kmer_idx] + write_positions[kmer_idx]
                base = wp * numpy.uint64(5)
                index_bytes[base + numpy.uint64(0)] = seq_idx & numpy.uint64(0xFF)
                index_bytes[base + numpy.uint64(1)] = (seq_idx >> numpy.uint64(8)) & numpy.uint64(0xFF)
                index_bytes[base + numpy.uint64(2)] = (seq_idx >> numpy.uint64(16)) & numpy.uint64(0xFF)
                index_bytes[base + numpy.uint64(3)] = (seq_idx >> numpy.uint64(24)) & numpy.uint64(0xFF)
                index_bytes[base + numpy.uint64(4)] = (seq_idx >> numpy.uint64(32)) & numpy.uint64(0xFF)
                write_positions[kmer_idx] += numpy.uint64(1)

    _compute_anchors = _compute_anchors_numba
    _count_anchors = _count_anchors_numba
    _populate_anchors = _populate_anchors_numba
else:
    _compute_anchors = _compute_anchors_python
    _count_anchors = _count_anchors_python
    _populate_anchors = _populate_anchors_python


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


class SharedMemorySparseAnchorIndexer:
    """
    Shared-memory inverted index that stores only 3 anchors per sequence:
    the first k-mer, the last k-mer, and one minimizer (splitmix64-hashed
    argmin) selected from the *internal* window starts (p in {1..L-k-1}).

    For 25-bp barcodes at k=14, this yields a 4x reduction of the inverted
    list compared to the dense :class:`SharedMemoryArray2BitKmerIndexer`,
    enabling 10B-scale whitelists. The error-tolerance ceiling for
    1-substitution reads is identical to dense (any base in the first/last
    overlap region appears in every 14-mer, so an error there breaks all
    anchors in either scheme); the minimizer's main wins are tier-3 ranking
    signal and rescue for 2-substitution reads when the errors split between
    front and back.

    :meth:`get_occurrences` returns candidates ordered as:
      Tier A (mask & 0b011 == 0b011): first AND last anchors both matched —
        all 25 positions are pinned to the query, so Hamming distance is 0.
      Tier B (popcount==2, not Tier A): one terminal anchor + minimizer,
        sorted by ascending Hamming distance.
      Tier C (popcount==1): single-anchor hits, sorted by ascending Hamming.
    The returned tuple shape ``(seq_str, score, positions)`` matches the
    dense indexer so existing SSW callers don't change. ``score`` is
    ``1000 * anchor_count - hamming`` so callers' "highest first" expectation
    is preserved.
    """

    SEQ_DTYPE = numpy.uint64
    INDEX_BYTE_DTYPE = numpy.uint8
    RANGES_DTYPE = numpy.uint64

    def __init__(self, known_bin_seq, kmer_size: int = 14, seq_len: int = 25):
        """
        Args:
            known_bin_seq: Pre-encoded sequences as integers (use str_to_2bit).
                Accepts a numpy.uint64 array (preferred — used without copy)
                or any sequence convertible by numpy.asarray.
            kmer_size: Length of k-mers (anchor length). Default 14.
            seq_len: Length of sequences. Default 25.
        """
        if seq_len < kmer_size + 2:
            raise ValueError(
                f"SharedMemorySparseAnchorIndexer requires seq_len >= "
                f"kmer_size + 2 (got seq_len={seq_len}, kmer_size={kmer_size}); "
                f"need at least one internal window for the minimizer."
            )
        if kmer_size > MAX_KMER_SIZE_FOR_UINT32_ANCHORS:
            raise ValueError(
                f"SharedMemorySparseAnchorIndexer stores anchor k-mer values "
                f"as uint32, so kmer_size must be <= "
                f"{MAX_KMER_SIZE_FOR_UINT32_ANCHORS} (got kmer_size={kmer_size})."
            )

        self.main_instance = True
        self.k = kmer_size
        self.seq_len = seq_len
        self.mask = (1 << (2 * self.k)) - 1
        self.seq_mask = (1 << (2 * self.seq_len)) - 1
        self.total_sequences = len(known_bin_seq)
        total_kmers = int(math.pow(4, self.k))

        # Cached shifts. The 2-bit encoding places string position 0 at the
        # most-significant bits, so the k-mer at start p is the slice of
        # length 2k anchored at LSB position (seq_len - k - p) * 2.
        self._first_shift = (seq_len - kmer_size) * 2
        self._last_shift = 0
        self._mid_shifts_py = tuple(
            (seq_len - kmer_size - p) * 2
            for p in range(1, seq_len - kmer_size)
        )

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
                f"SharedMemorySparseAnchorIndexer supports at most "
                f"{MAX_INDEXED_BARCODES} sequences (uint40 packed entries); "
                f"got {self.total_sequences}"
            )

        seqs = numpy.asarray(known_bin_seq, dtype=self.SEQ_DTYPE)

        first_shift = numpy.uint64(self._first_shift)
        last_shift = numpy.uint64(self._last_shift)
        k_mask = numpy.uint64(self.mask)
        mid_shifts = numpy.array(self._mid_shifts_py, dtype=numpy.uint64)

        logger.info(
            "Building sparse 3-anchor index for %d sequences (Numba: %s)",
            self.total_sequences, "enabled" if NUMBA_AVAILABLE else "disabled",
        )

        # Pass 1: extract the 3 anchors per sequence into a transient buffer.
        # uint32 is sufficient because 4**kmer_size fits in uint32 for
        # kmer_size <= 16 (validated above), and halves this buffer's size
        # vs. uint64 (e.g. 240 GB -> 120 GB at N=10B).
        anchors = numpy.empty(
            (self.total_sequences, NUM_ANCHORS), dtype=numpy.uint32,
        )
        _compute_anchors(seqs, first_shift, last_shift, k_mask, mid_shifts, anchors)
        logger.info("Anchor extraction complete")

        # Pass 2: count per-bucket occurrences for the cumulative offsets.
        kmer_counts = numpy.zeros(total_kmers, dtype=numpy.uint64)
        _count_anchors(anchors, kmer_counts)
        logger.info("Anchor counting complete")

        # Each barcode contributes exactly NUM_ANCHORS entries.
        self.index_size = int(NUM_ANCHORS * self.total_sequences)

        # Allocate barcodes shared memory and copy in once. After this the
        # temporary `seqs` reference is dropped — the shared-memory view is
        # the working array for the populate phase.
        self.barcodes_shared_memory = shared_memory.SharedMemory(
            create=True, size=self.total_sequences * self.SEQ_DTYPE().nbytes,
        )
        self.known_bin_seq = numpy.ndarray(
            shape=(self.total_sequences,), dtype=self.SEQ_DTYPE,
            buffer=self.barcodes_shared_memory.buf,
        )
        self.known_bin_seq[:] = seqs[:]
        del seqs

        self.index_shared_memory = shared_memory.SharedMemory(
            create=True, size=self.index_size * INDEX_BYTES_PER_ENTRY,
        )
        self.index = numpy.ndarray(
            shape=(self.index_size * INDEX_BYTES_PER_ENTRY,),
            dtype=self.INDEX_BYTE_DTYPE,
            buffer=self.index_shared_memory.buf,
        )

        self.index_ranges_shared_memory = shared_memory.SharedMemory(
            create=True, size=(total_kmers + 1) * self.RANGES_DTYPE().nbytes,
        )
        self.index_ranges = numpy.ndarray(
            shape=(total_kmers + 1,), dtype=self.RANGES_DTYPE,
            buffer=self.index_ranges_shared_memory.buf,
        )
        self.index_ranges[0] = 0
        self.index_ranges[1:] = numpy.cumsum(kmer_counts)

        # Pass 3: write packed uint40 entries into the inverted list.
        _populate_anchors(anchors, self.index_ranges, self.index)
        logger.info("Sparse anchor index population complete")

        del kmer_counts, anchors, mid_shifts
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
    def from_sharable_info(
        cls, shared_mem_index_info: SharedMemoryIndexInfo,
    ) -> 'SharedMemorySparseAnchorIndexer':
        kmer_index = cls.__new__(cls)
        kmer_index.main_instance = False
        kmer_index.k = shared_mem_index_info.kmer_size
        kmer_index.mask = (1 << (2 * kmer_index.k)) - 1
        kmer_index.seq_len = shared_mem_index_info.seq_len
        kmer_index.seq_mask = (1 << (2 * kmer_index.seq_len)) - 1
        kmer_index.index_size = shared_mem_index_info.index_size
        kmer_index.total_sequences = shared_mem_index_info.barcode_count
        total_kmers = int(math.pow(4, kmer_index.k))

        kmer_index._first_shift = (kmer_index.seq_len - kmer_index.k) * 2
        kmer_index._last_shift = 0
        kmer_index._mid_shifts_py = tuple(
            (kmer_index.seq_len - kmer_index.k - p) * 2
            for p in range(1, kmer_index.seq_len - kmer_index.k)
        )

        kmer_index.barcodes_shared_memory = shared_memory.SharedMemory(
            create=False, name=shared_mem_index_info.barcodes_sm_name)
        kmer_index.known_bin_seq = numpy.ndarray(
            shape=(kmer_index.total_sequences,),
            dtype=cls.SEQ_DTYPE,
            buffer=kmer_index.barcodes_shared_memory.buf)

        kmer_index.index_shared_memory = shared_memory.SharedMemory(
            create=False, name=shared_mem_index_info.index_sm_name)
        kmer_index.index = numpy.ndarray(
            shape=(kmer_index.index_size * INDEX_BYTES_PER_ENTRY,),
            dtype=cls.INDEX_BYTE_DTYPE,
            buffer=kmer_index.index_shared_memory.buf)

        kmer_index.index_ranges_shared_memory = shared_memory.SharedMemory(
            create=False, name=shared_mem_index_info.index_range_sm_name)
        kmer_index.index_ranges = numpy.ndarray(
            shape=(total_kmers + 1,), dtype=cls.RANGES_DTYPE,
            buffer=kmer_index.index_ranges_shared_memory.buf)
        return kmer_index

    def get_sharable_info(self) -> SharedMemoryIndexInfo:
        return SharedMemoryIndexInfo(
            self.total_sequences, self.k, self.seq_len, self.index_size,
            self.barcodes_shared_memory.name, self.index_shared_memory.name,
            self.index_ranges_shared_memory.name)

    def empty(self) -> bool:
        return len(self.known_bin_seq) == 0

    def _query_anchors(self, seq: int) -> Tuple[int, int, int]:
        """Compute (first14, last14, min14) for an encoded query."""
        first = (seq >> self._first_shift) & self.mask
        last = (seq >> self._last_shift) & self.mask
        best_hash = 0xFFFFFFFFFFFFFFFF
        best_kmer = 0
        for shift in self._mid_shifts_py:
            kmer = (seq >> shift) & self.mask
            h = _splitmix64_python(kmer)
            if h < best_hash:
                best_hash = h
                best_kmer = kmer
        return first, last, best_kmer

    def get_occurrences(self, sequence: str, max_hits: int = 0,
                        min_kmers: int = 1, hits_delta: int = 1,
                        ignore_equal: bool = False,
                        ) -> List[Tuple[str, int, List[int]]]:
        """
        Find indexed sequences that share at least one anchor with ``sequence``.

        Args:
            sequence: Query sequence (string, will be 2-bit encoded).
            max_hits: Maximum number of results (0 = unlimited).
            min_kmers: Minimum number of matching anchors required (1..3).
                With sparse indexing, the natural default is 1 — a single
                anchor hit is meaningful evidence at this k.
            hits_delta: Accepted but unused — kept for API compatibility
                with :class:`SharedMemoryArray2BitKmerIndexer`.
            ignore_equal: Skip exact matches.

        Returns:
            List of (sequence_str, score, positions) tuples. ``score`` is
            ``1000 * anchor_count - hamming_distance`` so the caller's
            descending-by-score sort produces tier-A → tier-B → tier-C
            order. ``positions`` carries the matched anchor offsets
            (0 = first, ``seq_len - k`` = last, ``-1`` = minimizer).
        """
        if len(self.known_bin_seq) == 0:
            return []

        seq = str_to_2bit(sequence)
        if seq > self.seq_mask:
            return []

        q_first, q_last, q_min = self._query_anchors(seq)

        slot_masks: dict[int, int] = {}
        index_buf = self.index
        for slot_bit, hash_val in (
                (0b001, q_first), (0b010, q_last), (0b100, q_min)):
            start_index = int(self.index_ranges[hash_val])
            end_index = int(self.index_ranges[hash_val + 1])
            for entry_idx in range(start_index, end_index):
                base = entry_idx * INDEX_BYTES_PER_ENTRY
                seq_idx = (int(index_buf[base])
                           | (int(index_buf[base + 1]) << 8)
                           | (int(index_buf[base + 2]) << 16)
                           | (int(index_buf[base + 3]) << 24)
                           | (int(index_buf[base + 4]) << 32))
                barcode = int(self.known_bin_seq[seq_idx])
                slot_masks[barcode] = slot_masks.get(barcode, 0) | slot_bit

        if not slot_masks:
            return []

        tier_a: List[Tuple[int, int, int]] = []
        tier_b: List[Tuple[int, int, int, int]] = []
        tier_c: List[Tuple[int, int, int, int]] = []
        first_pos = 0
        last_pos = self.seq_len - self.k
        for barcode, mask in slot_masks.items():
            if ignore_equal and barcode == seq:
                continue
            anchor_count = (mask & 1) + ((mask >> 1) & 1) + ((mask >> 2) & 1)
            if anchor_count < min_kmers:
                continue
            if (mask & 0b011) == 0b011:
                tier_a.append((barcode, anchor_count, mask))
            elif anchor_count == 2:
                ham = _hamming_2bit_python(seq, barcode)
                tier_b.append((barcode, anchor_count, mask, ham))
            else:
                ham = _hamming_2bit_python(seq, barcode)
                tier_c.append((barcode, anchor_count, mask, ham))

        # Tier A: prefer 3-anchor hits over 2-anchor first+last hits.
        tier_a.sort(key=lambda x: -x[1])
        tier_b.sort(key=lambda x: x[3])
        tier_c.sort(key=lambda x: x[3])

        result: List[Tuple[str, int, List[int]]] = []

        def _emit(barcode: int, anchor_count: int, mask: int, ham: int) -> None:
            score = 1000 * anchor_count - ham
            positions: List[int] = []
            if mask & 0b001:
                positions.append(first_pos)
            if mask & 0b010:
                positions.append(last_pos)
            if mask & 0b100:
                positions.append(-1)
            result.append((bit_to_str(barcode, self.seq_len), score, positions))

        for barcode, anchor_count, mask in tier_a:
            _emit(barcode, anchor_count, mask, 0)
            if max_hits and len(result) >= max_hits:
                return result
        for barcode, anchor_count, mask, ham in tier_b:
            _emit(barcode, anchor_count, mask, ham)
            if max_hits and len(result) >= max_hits:
                return result
        for barcode, anchor_count, mask, ham in tier_c:
            _emit(barcode, anchor_count, mask, ham)
            if max_hits and len(result) >= max_hits:
                return result

        return result
