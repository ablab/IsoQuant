############################################################################
# Copyright (c) 2023-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Memory-efficient k-mer indexers using 2-bit DNA encoding.

Dict2BitKmerIndexer: Dictionary with integer keys (~40% less memory than KmerIndexer)
Array2BitKmerIndexer: Flat array with range indexing for large barcode sets
"""

import math
from typing import List, Tuple, Iterable, DefaultDict
from collections import defaultdict

import numpy
from ..common import bit_to_str, str_to_2bit

# 2-bit encoded sequences, and offsets into the flat index. Matching what
# SharedMemoryArray2BitKmerIndexer uses, so the two stay interchangeable.
SEQ_DTYPE = numpy.uint64
INDEX_DTYPE = numpy.int64


class Dict2BitKmerIndexer:
    """
    Memory-efficient k-mer indexer using 2-bit encoded sequences and integer keys.

    Stores sequences as 2-bit encoded integers (8 bytes per 25bp barcode vs 25+ bytes for strings).
    Uses dictionary for k-mer index (vs 4^k array in Array2BitKmerIndexer).
    Best for medium to large barcode sets.
    """

    def __init__(self, known_bit_seqs: Iterable[int], kmer_size: int, seq_len: int):
        """
        Initialize k-mer index with 2-bit encoded sequences.

        Args:
            known_bit_seqs: Collection of 2-bit encoded sequences (use str_to_2bit or batch_str_to_2bit)
            kmer_size: Length of k-mers to use for indexing
            seq_len: Length of sequences in nucleotides
        """
        self.seq_list: List[int] = list(known_bit_seqs)
        self.k: int = kmer_size
        self.seq_len: int = seq_len
        self.mask: int = (1 << (2 * self.k)) - 1
        self.index: DefaultDict[int, List[int]] = defaultdict(list)
        self._index()

    def _get_kmer_indexes_from_bits(self, bin_seq: int) -> Iterable[int]:
        """Extract k-mer indices from a 2-bit encoded sequence."""
        # Convert numpy scalar to Python int for bitwise operations
        bin_seq = int(bin_seq)
        for i in range(self.seq_len - self.k + 1):
            yield (bin_seq >> ((self.seq_len - self.k - i) * 2)) & self.mask

    def _get_kmer_indexes_from_str(self, seq: str) -> Iterable[int]:
        """Generate binary-encoded k-mer indices from string using sliding window."""
        if len(seq) < self.k:
            return
        # Use same encoding as common.py: (ord(char) & 6) >> 1
        # A=0, C=1, T=2, G=3
        # Initialize first k-mer
        kmer_idx = 0
        for i in range(self.k):
            kmer_idx |= ((ord(seq[i]) & 6) >> 1) << ((self.k - i - 1) * 2)
        yield kmer_idx
        # Slide window
        for i in range(self.k, len(seq)):
            kmer_idx = ((kmer_idx << 2) & self.mask) | ((ord(seq[i]) & 6) >> 1)
            yield kmer_idx

    def _index(self) -> None:
        """Build k-mer index from all sequences."""
        for i, bin_seq in enumerate(self.seq_list):
            for kmer_idx in self._get_kmer_indexes_from_bits(bin_seq):
                self.index[kmer_idx].append(i)

    def empty(self) -> bool:
        """Check if index is empty."""
        return len(self.seq_list) == 0

    def get_occurrences(self, sequence: str, max_hits: int = 0, min_kmers: int = 1,
                        hits_delta: int = 1, ignore_equal: bool = False) -> List[Tuple[str, int, List[int]]]:
        """
        Find indexed sequences with shared k-mers.

        Args:
            sequence: Query sequence (string) to search
            max_hits: Maximum number of results (0 = unlimited)
            min_kmers: Minimum shared k-mers required
            hits_delta: Include results within this many k-mers of top hit
            ignore_equal: Skip exact matches

        Returns:
            List of (sequence_str, shared_kmer_count, kmer_positions) tuples
        """
        barcode_counts: DefaultDict[int, int] = defaultdict(int)
        barcode_positions: DefaultDict[int, List[int]] = defaultdict(list)

        query_bin = str_to_2bit(sequence)
        for pos, kmer_idx in enumerate(self._get_kmer_indexes_from_str(sequence)):
            for seq_index in self.index.get(kmer_idx, []):
                barcode_counts[seq_index] += 1
                barcode_positions[seq_index].append(pos)

        result = []
        for seq_index in barcode_counts.keys():
            count = barcode_counts[seq_index]
            if count < min_kmers:
                continue
            bin_seq = self.seq_list[seq_index]
            if ignore_equal and bin_seq == query_bin:
                continue
            result.append((bin_seq, count, barcode_positions[seq_index]))

        if not result:
            return []

        top_hits = max(result, key=lambda x: x[1])[1]
        result = filter(lambda x: x[1] >= top_hits - hits_delta, result)
        result = sorted(result, reverse=True, key=lambda x: x[1])

        # Convert 2-bit sequences back to strings
        if max_hits == 0:
            return [(bit_to_str(x[0], self.seq_len), x[1], x[2]) for x in result]
        return [(bit_to_str(x[0], self.seq_len), x[1], x[2]) for x in list(result)[:max_hits]]


class Array2BitKmerIndexer:
    """
    Memory-efficient k-mer indexer using 2-bit encoding for both k-mers and sequences.

    Stores sequences as integers (2 bits per nucleotide) to minimize memory.
    Uses flat array with range indexing for better cache performance.
    Best for large barcode sets (e.g., single-cell whitelists).
    """

    def __init__(self, known_bin_seq: Iterable[int], kmer_size: int, seq_len: int):
        """
        Initialize 2-bit k-mer index.

        Args:
            known_bin_seq: Pre-encoded sequences as integers (use str_to_2bit)
            kmer_size: Length of k-mers
            seq_len: Length of sequences (all must be same length)
        """
        self.k: int = kmer_size
        total_kmers = int(math.pow(4, kmer_size))
        self.mask: int = (1 << (2 * self.k)) - 1
        self.seq_len: int = seq_len
        self.seq_mask: int = (1 << (2 * self.seq_len)) - 1
        self.total_sequences: int = 0
        self._index(known_bin_seq, total_kmers)

    def _get_kmer_bin_indexes(self, bin_seq: int) -> Iterable[int]:
        """Extract k-mer indices from a 2-bit encoded sequence."""
        for i in range(self.seq_len - self.k + 1):
            yield (bin_seq >> ((self.seq_len - self.k - i) * 2)) & self.mask

    def _kmer_dtype(self):
        """Narrowest type holding a 2k-bit k-mer code; halves the sort's temporaries."""
        return numpy.uint32 if 2 * self.k <= 32 else numpy.uint64

    def _kmers_at(self, sequences: numpy.ndarray, offset: int) -> numpy.ndarray:
        """The k-mer every sequence carries at one offset, as a 2-bit array."""
        shift = numpy.uint64((self.seq_len - self.k - offset) * 2)
        kmers = numpy.bitwise_and(numpy.right_shift(sequences, shift), numpy.uint64(self.mask))
        return kmers.astype(self._kmer_dtype(), copy=False)

    def _index(self, known_bin_seq: Iterable[int], total_kmers: int) -> None:
        """Build a flat k-mer index from 2-bit encoded sequences.

        Counts each k-mer first and fills one flat array, rather than collecting entries
        into 4^k per-k-mer lists. Stays in numpy: a boxed code costs ~28 bytes extra.
        """
        if isinstance(known_bin_seq, numpy.ndarray):
            sequences = known_bin_seq.astype(SEQ_DTYPE, copy=False)
        else:
            sequences = numpy.fromiter(known_bin_seq, dtype=SEQ_DTYPE)
        self.total_sequences = len(sequences)
        n_offsets = max(0, self.seq_len - self.k + 1)

        # every (sequence, offset) pair, laid out sequence-major so that pair p belongs to
        # sequence p // n_offsets. Sorting these by k-mer groups the index and, being
        # stable, leaves each k-mer's entries in sequence-then-offset order.
        kmers = numpy.empty(self.total_sequences * n_offsets, dtype=self._kmer_dtype())
        for offset in range(n_offsets):
            kmers[offset::n_offsets] = self._kmers_at(sequences, offset)

        self.index_ranges = numpy.zeros(total_kmers + 1, dtype=INDEX_DTYPE)
        self.index_ranges[1:] = numpy.bincount(kmers, minlength=total_kmers)
        numpy.cumsum(self.index_ranges, out=self.index_ranges)

        order = numpy.argsort(kmers, kind="stable")
        del kmers
        # pair index -> sequence index, in place so the pair array is never held twice
        numpy.floor_divide(order, n_offsets, out=order)
        self.index = sequences[order]

    def empty(self) -> bool:
        """Check if index is empty."""
        return self.total_sequences == 0

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
        barcode_counts: DefaultDict[int, int] = defaultdict(int)
        barcode_positions: DefaultDict[int, List[int]] = defaultdict(list)

        seq = str_to_2bit(sequence)
        for pos, kmer_idx in enumerate(self._get_kmer_bin_indexes(seq)):
            start_index = self.index_ranges[kmer_idx]
            end_index = self.index_ranges[kmer_idx + 1]
            # tolist() unboxes the whole slice in C; leaving them as numpy scalars would
            # make every downstream shift a uint64/int64 mixed-type operation
            for barcode in self.index[start_index:end_index].tolist():
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

        # Convert 2-bit sequences back to strings
        if max_hits == 0:
            return [(bit_to_str(x[0], self.seq_len), x[1], x[2]) for x in result]
        return [(bit_to_str(x[0], self.seq_len), x[1], x[2]) for x in list(result)[:max_hits]]
