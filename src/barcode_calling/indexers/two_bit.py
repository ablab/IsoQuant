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
from typing import List, Tuple, Iterable, Dict, DefaultDict
from collections import defaultdict
from ..common import bit_to_str, str_to_2bit


class Dict2BitKmerIndexer:
    """
    Memory-efficient k-mer indexer using 2-bit encoded sequences and integer keys.

    Stores sequences as 2-bit encoded integers (8 bytes per 25bp barcode vs 25+ bytes for strings).
    Uses dictionary for k-mer index (vs 4^k array in Array2BitKmerIndexer).
    Best for medium to large barcode sets.
    """

    def __init__(self, known_bit_seqs: Iterable[int], kmer_size: int = 6, seq_len: int = 25):
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

    def __init__(self, known_bin_seq: Iterable[int], kmer_size: int = 12, seq_len: int = 25):
        """
        Initialize 2-bit k-mer index.

        Args:
            known_bin_seq: Pre-encoded sequences as integers (use str_to_2bit)
            kmer_size: Length of k-mers
            seq_len: Length of sequences (all must be same length)
        """
        self.k: int = kmer_size
        total_kmers = int(math.pow(4, kmer_size))
        tmp_index: List[List[int]] = [[] for _ in range(total_kmers)]
        self.mask: int = (1 << (2 * self.k)) - 1
        self.seq_len: int = seq_len
        self.seq_mask: int = (1 << (2 * self.seq_len)) - 1
        self.total_sequences: int = 0
        self._index(known_bin_seq, tmp_index)

        # Flatten index for better cache performance
        self.index: List[int] = []
        self.index_ranges: List[int] = [0]
        for l in tmp_index:
            self.index += l
            self.index_ranges.append(len(self.index))

    def _get_kmer_bin_indexes(self, bin_seq: int) -> Iterable[int]:
        """Extract k-mer indices from a 2-bit encoded sequence."""
        for i in range(self.seq_len - self.k + 1):
            yield (bin_seq >> ((self.seq_len - self.k - i) * 2)) & self.mask

    def _index(self, known_bin_seq: Iterable[int], tmp_index: List[List[int]]) -> None:
        """Build k-mer index from 2-bit encoded sequences."""
        for bin_seq in known_bin_seq:
            self.total_sequences += 1
            for kmer_idx in self._get_kmer_bin_indexes(bin_seq):
                tmp_index[kmer_idx].append(bin_seq)

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
            for barcode_index in range(start_index, end_index):
                barcode = self.index[barcode_index]
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
