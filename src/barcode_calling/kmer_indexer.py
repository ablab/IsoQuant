############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################
import math
from collections import defaultdict
import ctypes
import os


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


class HashKmerIndexerC:
    lib = ctypes.cdll.LoadLibrary(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'cseqlib.so'))
    # @params:
    # known_strings: collection of strings (barcodes or UMI)
    # kmer_size: K to use for indexing
    def __init__(self, barcodes, kmer_size=7):
        self._kmer_index = HashKmerIndexerC.lib.KmerIndexer_new(ctypes.c_int(kmer_size))
        for b in barcodes:
            HashKmerIndexerC.lib.KmerIndexer_add(self._kmer_index, b.encode(), ctypes.c_int(len(b)))
        self.barcodes = barcodes

    def append(self, barcode):
        HashKmerIndexerC.lib.KmerIndexer_add(self._kmer_index, barcode.encode(), ctypes.c_int(len(barcode)))

    def empty(self):
        return ctypes.cast(HashKmerIndexerC.lib.KmerIndexer_empty(self._kmer_index), bool)

    # @params:
    # sequence: a string to be searched against known strings
    # max_hits: return at most max_hits candidates
    # min_kmers: minimal number of matching k-mers
    # @return
    # a list of pairs (string, numer of common k-mers) sorted descending by the number of shared k-mers
    def get_occurrences(self, sequence, max_hits=0, min_kmers=1, hits_delta=1):
        ArrayType = ctypes.c_int * 10000
        result = ArrayType()
        HashKmerIndexerC.lib.KmerIndexer_get_occurrences(self._kmer_index, sequence.encode(),
                                                         ctypes.c_int(len(sequence)),
                                                         ctypes.c_int(max_hits), ctypes.c_int(min_kmers),
                                                         ctypes.c_int(hits_delta), result, ctypes.c_int(100))

        total_candidates = result[0]
        current_idx = 1
        result_dict = {}
        for i in range(total_candidates):
            seq = self.barcodes[result[current_idx]]
            current_idx += 1
            hit_count = result[current_idx]
            current_idx += 1
            pos_hits = result[current_idx : current_idx + hit_count]
            current_idx += hit_count
            result_dict[seq] = (seq, hit_count, pos_hits)

        return result_dict


class ArrayKmerIndexerC:
    lib = ctypes.cdll.LoadLibrary(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'cseqlib.so'))
    # @params:
    # known_strings: collection of strings (barcodes or UMI)
    # kmer_size: K to use for indexing
    def __init__(self, barcodes, kmer_size=7):
        self._kmer_index = ArrayKmerIndexerC.lib.ArrayKmerIndexer_new(ctypes.c_int(kmer_size))
        for b in barcodes:
            ArrayKmerIndexerC.lib.ArrayKmerIndexer_add(self._kmer_index, b.encode(), ctypes.c_int(len(b)))
        self.barcodes = barcodes

    def append(self, barcode):
        ArrayKmerIndexerC.lib.ArrayKmerIndexer_add(self._kmer_index, barcode.encode(), ctypes.c_int(len(barcode)))

    def empty(self):
        return ctypes.cast(ArrayKmerIndexerC.lib.ArrayKmerIndexer_empty(self._kmer_index), bool)

    # @params:
    # sequence: a string to be searched against known strings
    # max_hits: return at most max_hits candidates
    # min_kmers: minimal number of matching k-mers
    # @return
    # a list of pairs (string, numer of common k-mers) sorted descending by the number of shared k-mers
    def get_occurrences(self, sequence, max_hits=0, min_kmers=1, hits_delta=1):
        ArrayType = ctypes.c_int * 10000
        result = ArrayType()
        ArrayKmerIndexerC.lib.ArrayKmerIndexer_get_occurrences(self._kmer_index, sequence.encode(),
                                                               ctypes.c_int(len(sequence)),
                                                               ctypes.c_int(max_hits), ctypes.c_int(min_kmers),
                                                               ctypes.c_int(hits_delta), result, ctypes.c_int(100))

        total_candidates = result[0]
        current_idx = 1
        result_dict = {}
        for i in range(total_candidates):
            seq = self.barcodes[result[current_idx]]
            current_idx += 1
            hit_count = result[current_idx]
            current_idx += 1
            pos_hits = result[current_idx : current_idx + hit_count]
            current_idx += hit_count
            result_dict[seq] = (seq, hit_count, pos_hits)

        return result_dict