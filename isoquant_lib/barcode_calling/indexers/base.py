############################################################################
# Copyright (c) 2023-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Basic k-mer indexers for barcode calling.

KmerIndexer: Dictionary-based, best for small to medium sets
ArrayKmerIndexer: Array-based with 2-bit encoding, best for k <= 8
"""

import math
from typing import List, Tuple, Iterable, Dict, DefaultDict
from collections import defaultdict

from ..common import NUCL2BIN as SHARED_NUCL2BIN

# Returned for k-mers absent from an index. A shared empty tuple, so a miss costs no
# allocation -- lookups run once per base of every scanned region.
_NO_HITS: Tuple[int, ...] = ()


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
        # A plain dict, never a defaultdict: lookups scan whole cDNA regions, and a
        # subscript read on a defaultdict would insert an entry per foreign k-mer,
        # growing the index without bound for the life of the process.
        self.index: Dict[str, List[int]] = {}
        self._index()

    def _get_kmers(self, seq: str) -> Iterable[str]:
        """Generate all k-mers from a sequence using sliding window."""
        if len(seq) < self.k:
            return
        kmer = seq[:self.k]
        yield kmer
        for i in range(self.k, len(seq)):
            kmer = kmer[1:] + seq[i]
            yield kmer

    def _index(self) -> None:
        """Build k-mer index from all sequences."""
        for i, barcode in enumerate(self.seq_list):
            for kmer in self._get_kmers(barcode):
                self.index.setdefault(kmer, []).append(i)

    def append(self, barcode: str) -> None:
        """Add a new barcode to the index."""
        self.seq_list.append(barcode)
        index = len(self.seq_list) - 1
        for kmer in self._get_kmers(barcode):
            self.index.setdefault(kmer, []).append(index)

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

        for pos, kmer in enumerate(self._get_kmers(sequence)):
            for i in self.index.get(kmer, _NO_HITS):
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

    def _get_kmers_substr(self, seq, start, end):
        if end - start + 1 < self.k:
            return
        kmer = seq[start:start + self.k]
        yield kmer
        for i in range(start + self.k, min(len(seq), end + 1)):
            kmer = kmer[1:] + seq[i]
            yield kmer

    # [start, end]
    def get_occurrences_substr(self, sequence, start, end, max_hits=0, min_kmers=1, hits_delta=1) -> List[Tuple[str, int, List[int]]]:
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

        for pos, kmer in enumerate(self._get_kmers_substr(sequence, start, end)):
            for i in self.index.get(kmer, _NO_HITS):
                barcode_counts[i] += 1
                barcode_positions[i].append(pos)

        result = []
        for i in barcode_counts.keys():
            count = barcode_counts[i]
            if count < min_kmers:
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


class ArrayKmerIndexer:
    """
    Optimized k-mer indexer using array-based storage and binary encoding.

    Converts k-mers to integers (2 bits per nucleotide) for faster lookup.
    Memory usage: O(4^k) array entries. Best for k <= 8.
    """

    # the shared 2-bit encoding, so a code means the same thing here as everywhere else
    NUCL2BIN: Dict[str, int] = SHARED_NUCL2BIN

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
        self.mask: int = (1 << (2 * self.k)) - 1
        self._index()

    def _get_kmer_indexes(self, seq: str) -> Iterable[int]:
        """Generate binary-encoded k-mer indices using sliding window."""
        if len(seq) < self.k:
            return

        # Initialize first k-mer
        kmer_idx = 0
        for i in range(self.k):
            kmer_idx |= ArrayKmerIndexer.NUCL2BIN[seq[i]] << ((self.k - i - 1) * 2)
        yield kmer_idx

        # Slide window: shift left 2 bits, add new nucleotide
        for i in range(self.k, len(seq)):
            kmer_idx <<= 2
            kmer_idx &= self.mask
            kmer_idx |= ArrayKmerIndexer.NUCL2BIN[seq[i]]
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
        """Find indexed sequences with shared k-mers."""
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
