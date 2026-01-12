############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
K-mer indexers for fast approximate string matching.

Used for barcode calling in single-cell data. Indexes known barcodes/UMIs
and allows fast lookup of similar sequences based on shared k-mers.
"""

import math
from typing import List, Tuple, Iterable, Dict, DefaultDict
from collections import defaultdict
from .common import bit_to_str, str_to_2bit


class KmerIndexer:
    """
    Basic k-mer indexer using dictionary-based storage.

    Indexes sequences by their k-mers for fast approximate matching.
    Best for small to medium barcode sets.
    """

    def __init__(self, known_strings: Iterable[str], kmer_size: int = 6):
        """
        Initialize k-mer index.

        Args:
            known_strings: Collection of reference sequences (barcodes/UMIs)
            kmer_size: Length of k-mers to use for indexing
        """
        self.seq_list: List[str] = list(known_strings)
        self.k: int = kmer_size
        self.index: DefaultDict[str, List[int]] = defaultdict(list)
        self._index()

    def _get_kmers(self, seq: str) -> Iterable[str]:
        """
        Generate all k-mers from a sequence using sliding window.

        Args:
            seq: Input sequence

        Yields:
            K-mer strings
        """
        if len(seq) < self.k:
            return
        kmer = seq[:self.k]
        yield kmer
        # Slide window by removing first char and adding next
        for i in range(self.k, len(seq)):
            kmer = kmer[1:] + seq[i]
            yield kmer

    def _index(self) -> None:
        """Build k-mer index from all sequences."""
        for i, barcode in enumerate(self.seq_list):
            for kmer in self._get_kmers(barcode):
                self.index[kmer].append(i)

    def append(self, barcode: str) -> None:
        """
        Add a new barcode to the index.

        Args:
            barcode: Sequence to add to index
        """
        self.seq_list.append(barcode)
        index = len(self.seq_list) - 1
        for kmer in self._get_kmers(barcode):
            self.index[kmer].append(index)

    def empty(self) -> bool:
        """Check if index is empty."""
        return len(self.seq_list) == 0

    def get_occurrences(self, sequence: str, max_hits: int = 0, min_kmers: int = 1,
                       hits_delta: int = 1, ignore_equal: bool = False) -> List[Tuple[str, int, List[int]]]:
        """
        Find indexed sequences with shared k-mers.

        Args:
            sequence: Query sequence to search
            max_hits: Maximum number of results (0 = unlimited)
            min_kmers: Minimum shared k-mers required
            hits_delta: Include results within this many k-mers of top hit
            ignore_equal: Skip exact matches

        Returns:
            List of (sequence, shared_kmer_count, kmer_positions) tuples,
            sorted by shared k-mer count (descending)
        """
        # Count shared k-mers for each indexed sequence
        barcode_counts: DefaultDict[int, int] = defaultdict(int)
        barcode_positions: DefaultDict[int, List[int]] = defaultdict(list)

        for pos, kmer in enumerate(self._get_kmers(sequence)):
            for i in self.index[kmer]:
                barcode_counts[i] += 1
                barcode_positions[i].append(pos)

        # Filter and build results
        result = []
        for i in barcode_counts.keys():
            count = barcode_counts[i]
            if count < min_kmers:
                continue
            if ignore_equal and self.seq_list[i] == sequence:
                continue
            result.append((self.seq_list[i], count, barcode_positions[i]))

        if not result:
            return []

        # Keep only top hits within hits_delta
        top_hits = max(result, key=lambda x: x[1])[1]
        result = filter(lambda x: x[1] >= top_hits - hits_delta, result)
        result = list(sorted(result, reverse=True, key=lambda x: x[1]))

        if max_hits == 0:
            return result
        return result[:max_hits]


class Dict2BitKmerIndexer:
    """
    Memory-efficient k-mer indexer using 2-bit encoded integer keys.

    Uses integer k-mer indices (2 bits per nucleotide) instead of string keys,
    reducing memory by ~40% compared to KmerIndexer while maintaining the same API.
    Best for medium to large barcode sets where ArrayKmerIndexer's O(4^k) memory
    would be prohibitive.
    """

    # Nucleotide to 2-bit encoding (A=00, C=01, G=10, T=11)
    NUCL2BIN: Dict[str, int] = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3}

    def __init__(self, known_strings: Iterable[str], kmer_size: int = 6):
        """
        Initialize k-mer index with integer keys.

        Args:
            known_strings: Collection of reference sequences (barcodes/UMIs)
            kmer_size: Length of k-mers to use for indexing
        """
        self.seq_list: List[str] = list(known_strings)
        self.k: int = kmer_size
        self.mask: int = (1 << (2 * self.k)) - 1  # Bit mask for k-mer
        self.index: DefaultDict[int, List[int]] = defaultdict(list)  # int keys!
        self._index()

    def _get_kmer_indexes(self, seq: str) -> Iterable[int]:
        """
        Generate binary-encoded k-mer indices using sliding window.

        Args:
            seq: Input sequence

        Yields:
            Integer k-mer indices
        """
        if len(seq) < self.k:
            return

        # Initialize first k-mer
        kmer_idx = 0
        for i in range(self.k):
            kmer_idx |= Dict2BitKmerIndexer.NUCL2BIN[seq[i]] << ((self.k - i - 1) * 2)
        yield kmer_idx

        # Slide window: shift left 2 bits, add new nucleotide
        for i in range(self.k, len(seq)):
            kmer_idx = ((kmer_idx << 2) & self.mask) | Dict2BitKmerIndexer.NUCL2BIN[seq[i]]
            yield kmer_idx

    def _index(self) -> None:
        """Build k-mer index from all sequences."""
        for i, seq in enumerate(self.seq_list):
            for kmer_idx in self._get_kmer_indexes(seq):
                self.index[kmer_idx].append(i)

    def append(self, seq: str) -> None:
        """
        Add a new sequence to the index.

        Args:
            seq: Sequence to add to index
        """
        i = len(self.seq_list)
        self.seq_list.append(seq)
        for kmer_idx in self._get_kmer_indexes(seq):
            self.index[kmer_idx].append(i)

    def empty(self) -> bool:
        """Check if index is empty."""
        return len(self.seq_list) == 0

    def get_occurrences(self, sequence: str, max_hits: int = 0, min_kmers: int = 1,
                       hits_delta: int = 1, ignore_equal: bool = False) -> List[Tuple[str, int, List[int]]]:
        """
        Find indexed sequences with shared k-mers.

        Args:
            sequence: Query sequence to search
            max_hits: Maximum number of results (0 = unlimited)
            min_kmers: Minimum shared k-mers required
            hits_delta: Include results within this many k-mers of top hit
            ignore_equal: Skip exact matches

        Returns:
            List of (sequence, shared_kmer_count, kmer_positions) tuples,
            sorted by shared k-mer count (descending)
        """
        barcode_counts: DefaultDict[int, int] = defaultdict(int)
        barcode_positions: DefaultDict[int, List[int]] = defaultdict(list)

        for pos, kmer_idx in enumerate(self._get_kmer_indexes(sequence)):
            for seq_index in self.index.get(kmer_idx, []):
                barcode_counts[seq_index] += 1
                barcode_positions[seq_index].append(pos)

        result = []
        for seq_index in barcode_counts.keys():
            count = barcode_counts[seq_index]
            if count < min_kmers:
                continue
            seq = self.seq_list[seq_index]
            if ignore_equal and seq == sequence:
                continue
            result.append((seq, count, barcode_positions[seq_index]))

        if not result:
            return []

        top_hits = max(result, key=lambda x: x[1])[1]
        result = filter(lambda x: x[1] >= top_hits - hits_delta, result)
        result = list(sorted(result, reverse=True, key=lambda x: x[1]))

        if max_hits == 0:
            return result
        return result[:max_hits]


class ArrayKmerIndexer:
    """
    Optimized k-mer indexer using array-based storage and binary encoding.

    Converts k-mers to integers (2 bits per nucleotide) for faster lookup.
    Memory usage: O(4^k) array entries. Best for k <= 8.
    """

    # Nucleotide to 2-bit encoding (A=00, C=01, G=10, T=11)
    NUCL2BIN: Dict[str, int] = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3}

    def __init__(self, known_strings: Iterable[str], kmer_size: int = 6):
        """
        Initialize array-based k-mer index.

        Args:
            known_strings: Collection of reference sequences
            kmer_size: Length of k-mers (recommended <= 8)
        """
        self.seq_list: List[str] = list(known_strings)
        self.k: int = kmer_size
        total_kmers = int(math.pow(4, kmer_size))  # 4^k possible k-mers
        self.index: List[List[int]] = [[] for _ in range(total_kmers)]
        self.mask: int = (1 << (2 * self.k)) - 1  # Bit mask for k-mer
        self._index()

    def _get_kmer_indexes(self, seq: str) -> Iterable[int]:
        """
        Generate binary-encoded k-mer indices using sliding window.

        Converts k-mers to integers (2 bits per nucleotide) for fast lookup.

        Args:
            seq: Input sequence

        Yields:
            Integer k-mer indices
        """
        if len(seq) < self.k:
            return

        # Initialize first k-mer
        kmer_idx = 0
        for i in range(self.k):
            # Encode nucleotide as 2 bits and shift into position
            kmer_idx |= ArrayKmerIndexer.NUCL2BIN[seq[i]] << ((self.k - i - 1) * 2)
        yield kmer_idx

        # Slide window: shift left 2 bits, add new nucleotide
        for i in range(self.k, len(seq)):
            kmer_idx <<= 2  # Shift left to drop leftmost nucleotide
            kmer_idx &= self.mask  # Keep only k nucleotides
            kmer_idx |= ArrayKmerIndexer.NUCL2BIN[seq[i]]  # Add rightmost nucleotide
            yield kmer_idx

    def _index(self) -> None:
        """Build k-mer index from all sequences."""
        for i, barcode in enumerate(self.seq_list):
            for kmer_idx in self._get_kmer_indexes(barcode):
                self.index[kmer_idx].append(i)

    def append(self, barcode: str) -> None:
        """Add a new barcode to the index."""
        self.seq_list.append(barcode)
        index = len(self.seq_list) - 1
        for kmer_idx in self._get_kmer_indexes(barcode):
            self.index[kmer_idx].append(index)

    def empty(self) -> bool:
        """Check if index is empty."""
        return len(self.seq_list) == 0

    def get_occurrences(self, sequence: str, max_hits: int = 0, min_kmers: int = 1,
                       hits_delta: int = 1, ignore_equal: bool = False) -> List[Tuple[str, int, List[int]]]:
        """Find indexed sequences with shared k-mers (same as KmerIndexer.get_occurrences)."""
        barcode_counts: DefaultDict[int, int] = defaultdict(int)
        barcode_positions: DefaultDict[int, List[int]] = defaultdict(list)

        for pos, kmer_idx in enumerate(self._get_kmer_indexes(sequence)):
            for i in self.index[kmer_idx]:
                barcode_counts[i] += 1
                barcode_positions[i].append(pos)

        result = []
        for i in barcode_counts.keys():
            count = barcode_counts[i]
            if count < min_kmers:
                continue
            if ignore_equal and self.seq_list[i] == sequence:
                continue
            result.append((self.seq_list[i], count, barcode_positions[i]))

        if not result:
            return []

        top_hits = max(result, key=lambda x: x[1])[1]
        result = filter(lambda x: x[1] >= top_hits - hits_delta, result)
        result = list(sorted(result, reverse=True, key=lambda x: x[1]))

        if max_hits == 0:
            return result
        return result[:max_hits]


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
            known_bin_seq: Pre-encoded sequences as integers
            kmer_size: Length of k-mers
            seq_len: Length of sequences (all must be same length)
        """
        self.k: int = kmer_size
        total_kmers = int(math.pow(4, kmer_size))
        tmp_index: List[List[int]] = [[] for _ in range(total_kmers)]
        self.mask: int = (1 << (2 * self.k)) - 1  # K-mer bit mask
        self.seq_len: int = seq_len
        self.seq_mask: int = (1 << (2 * self.seq_len)) - 1  # Sequence bit mask
        self.total_sequences: int = 0
        self._index(known_bin_seq, tmp_index)

        # Flatten index for better cache performance
        self.index: List[int] = []
        self.index_ranges: List[int] = [0]
        for l in tmp_index:
            self.index += l
            self.index_ranges.append(len(self.index))


    def _get_kmer_bin_indexes(self, bin_seq: int) -> Iterable[int]:
        """
        Extract k-mer indices from a 2-bit encoded sequence.

        Args:
            bin_seq: Sequence encoded as integer (2 bits per nucleotide)

        Yields:
            Integer k-mer indices
        """
        for i in range(self.seq_len - self.k + 1):
            # Extract k-mer by shifting and masking
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
            # Use flat index with ranges for cache-friendly access
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

        # Filter top hits
        top_hits = max(result, key=lambda x: x[1])[1]
        result = filter(lambda x: x[1] >= top_hits - hits_delta, result)
        result = sorted(result, reverse=True, key=lambda x: x[1])

        # Convert 2-bit sequences back to strings
        if max_hits == 0:
            return [(bit_to_str(x[0], self.seq_len), x[1], x[2]) for x in result]
        return [(bit_to_str(x[0], self.seq_len), x[1], x[2]) for x in list(result)[:max_hits]]