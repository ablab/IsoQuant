############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################
import math
from collections import defaultdict
from .common import bit_to_str, str_to_2bit


class KmerIndexer:
    # @params:
    # known_strings: collection of strings (barcodes or UMI)
    # kmer_size: K to use for indexing
    def __init__(self, known_strings, kmer_size=6):
        self.seq_list = list(known_strings)
        self.k = kmer_size
        self.index = defaultdict(list)
        self._index()

    def _get_kmers(self, seq):
        if len(seq) < self.k:
            return
        kmer = seq[:self.k]
        yield kmer
        for i in range(self.k, len(seq)):
            kmer = kmer[1:] + seq[i]
            yield kmer

    def _index(self):
        for i, barcode in enumerate(self.seq_list):
            for kmer in self._get_kmers(barcode):
                self.index[kmer].append(i)

    def append(self, barcode):
        self.seq_list.append(barcode)
        index = len(self.seq_list) - 1
        for kmer in self._get_kmers(barcode):
            self.index[kmer].append(index)

    def empty(self):
        return len(self.seq_list) == 0

    # @params:
    # sequence: a string to be searched against known strings
    # max_hits: return at most max_hits candidates
    # min_kmers: minimal number of matching k-mers
    # @return
    # a list of pairs (string, numer of common k-mers) sorted descending by the number of shared k-mers
    def get_occurrences(self, sequence, max_hits=0, min_kmers=1, hits_delta=1, ignore_equal=False):
        barcode_counts = defaultdict(int)
        barcode_positions = defaultdict(list)
        for pos, kmer in enumerate(self._get_kmers(sequence)):
            for i in self.index[kmer]:
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
            return {}

        top_hits = max(result, key=lambda x: x[1])[1]
        result = filter(lambda x: x[1] >= top_hits - hits_delta, result)
        result = sorted(result, reverse=True, key=lambda x: x[1])

        if max_hits == 0:
            return {x[0]:x for x in result}
        return {x[0]:x for x in list(result)[:max_hits]}


class ArrayKmerIndexer:
    NUCL2BIN = {'A':0, 'C':1, 'G':2, 'T':3, 'a':0, 'c':1, 'g':2, 't':3}

    # @params:
    # known_strings: collection of strings (barcodes or UMI)
    # kmer_size: K to use for indexing
    def __init__(self, known_strings, kmer_size=6):
        self.seq_list = list(known_strings)
        self.k = kmer_size
        total_kmers = int(math.pow(4, kmer_size))
        self.index = []
        for i in range(total_kmers):
            self.index.append([])
        self.mask = (1 << (2 * self.k)) - 1
        self._index()

    def _get_kmer_indexes(self, seq):
        if len(seq) < self.k:
            return
        kmer_idx = 0
        for i in range(self.k):
            # TODO: NUCL2BIN[seq[i]] -> (seq[i] & 6) >> 1
            kmer_idx |= ArrayKmerIndexer.NUCL2BIN[seq[i]] << ((self.k - i - 1) * 2)
        yield kmer_idx
        for i in range(self.k, len(seq)):
            kmer_idx <<= 2
            kmer_idx &= self.mask
            kmer_idx |= ArrayKmerIndexer.NUCL2BIN[seq[i]]
            yield kmer_idx

    def _index(self):
        for i, barcode in enumerate(self.seq_list):
            for kmer_idx in self._get_kmer_indexes(barcode):
                self.index[kmer_idx].append(i)

    def append(self, barcode):
        self.seq_list.append(barcode)
        index = len(self.seq_list) - 1
        for kmer_idx in self._get_kmer_indexes(barcode):
            self.index[kmer_idx].append(index)

    def empty(self):
        return len(self.seq_list) == 0

    # @params:
    # sequence: a string to be searched against known strings
    # max_hits: return at most max_hits candidates
    # min_kmers: minimal number of matching k-mers
    # @return
    # a list of pairs (string, numer of common k-mers) sorted descending by the number of shared k-mers
    def get_occurrences(self, sequence, max_hits=0, min_kmers=1, hits_delta=1, ignore_equal=False):
        barcode_counts = defaultdict(int)
        barcode_positions = defaultdict(list)

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
            return {}

        top_hits = max(result, key=lambda x: x[1])[1]
        result = filter(lambda x: x[1] >= top_hits - hits_delta, result)
        result = sorted(result, reverse=True, key=lambda x: x[1])

        if max_hits == 0:
            return {x[0]:x for x in result}
        return {x[0]:x for x in list(result)[:max_hits]}


class Array2BitKmerIndexer:
    # @params:
    # known_strings: collection of strings (barcodes or UMI)
    # kmer_size: K to use for indexing
    def __init__(self, known_bin_seq, kmer_size=12, seq_len=25):
        self.k = kmer_size
        total_kmers = int(math.pow(4, kmer_size))
        self.index = []
        for i in range(total_kmers):
            self.index.append([])
        self.mask = (1 << (2 * self.k)) - 1
        self.seq_len = seq_len
        self.seq_mask = (1 << (2 * self.seq_len)) - 1
        self.total_sequences = 0
        self._index(known_bin_seq)

    def _get_kmer_bin_indexes(self, bin_seq):
        kmer_idx = 0
        for i in range(self.seq_len - self.k + 1):
            yield (bin_seq >> ((self.seq_len - self.k - i) * 2)) & self.mask

    def _index(self, known_bin_seq):
        for bin_seq in known_bin_seq:
            self.total_sequences += 1
            for kmer_idx in self._get_kmer_bin_indexes(bin_seq):
                self.index[kmer_idx].append(bin_seq)

    def append(self, bin_seq):
        self.total_sequences += 1
        for kmer_idx in self._get_kmer_bin_indexes(bin_seq):
            self.index[kmer_idx].append(bin_seq)

    def empty(self):
        return self.total_sequences == 0

    # @params:
    # sequence: a string to be searched against known strings
    # max_hits: return at most max_hits candidates
    # min_kmers: minimal number of matching k-mers
    # @return
    # a list of pairs (string, numer of common k-mers) sorted descending by the number of shared k-mers
    def get_occurrences(self, sequence, max_hits=0, min_kmers=1, hits_delta=1, ignore_equal=False):
        barcode_counts = defaultdict(int)
        barcode_positions = defaultdict(list)

        seq = str_to_2bit(sequence, len(sequence))
        for pos, kmer_idx in enumerate(self._get_kmer_bin_indexes(seq)):
            for barcode in self.index[kmer_idx]:
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
            return {}

        top_hits = max(result, key=lambda x: x[1])[1]
        result = filter(lambda x: x[1] >= top_hits - hits_delta, result)
        result = sorted(result, reverse=True, key=lambda x: x[1])

        if max_hits == 0:
            return {bit_to_str(x[0], self.seq_len):x for x in result}
        return {bit_to_str(x[0], self.seq_len):x for x in list(result)[:max_hits]}
