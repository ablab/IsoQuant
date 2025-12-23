############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################
import math
import numpy
from collections import defaultdict
from multiprocessing import shared_memory
from .common import bit_to_str, str_to_2bit


class SharedMemoryIndexInfo:
    def __init__(self, barcode_count, kmer_size, seq_len, index_size, barcodes_sm_name, index_sm_name, index_range_sm_name):
        self.barcode_count = barcode_count
        self.kmer_size = kmer_size
        self.seq_len = seq_len
        self.index_size = index_size
        self.barcodes_sm_name = barcodes_sm_name
        self.index_sm_name = index_sm_name
        self.index_range_sm_name = index_range_sm_name

    def __getstate__(self):
        return (self.barcode_count,
                self.kmer_size,
                self.seq_len,
                self.index_size,
                self.barcodes_sm_name,
                self.index_sm_name,
                self.index_range_sm_name)

    def __setstate__(self, state):
        self.barcode_count = state[0]
        self.kmer_size = state[1]
        self.seq_len = state[2]
        self.index_size = state[3]
        self.barcodes_sm_name = state[4]
        self.index_sm_name = state[5]
        self.index_range_sm_name = state[6]


class SharedMemoryArray2BitKmerIndexer:
    # @params:
    # known_bin_seq: collection of strings in binary or string format (barcodes or UMI)
    # kmer_size: K to use for indexing
    SEQ_DTYPE = numpy.uint64
    KMER_DTYPE = numpy.uint32
    INDEX_DTYPE = numpy.uint32

    def __init__(self, known_bin_seq: list, kmer_size=12, seq_len=25):
        self.main_instance = True
        self.k = kmer_size
        total_kmers = int(math.pow(4, self.k))
        tmp_index = []
        for i in range(total_kmers):
            tmp_index.append([])
        self.mask = (1 << (2 * self.k)) - 1
        self.seq_len = seq_len
        self.seq_mask = (1 << (2 * self.seq_len)) - 1
        self._index(known_bin_seq, tmp_index)
        self.index_size = sum(len(x) for x in tmp_index)
        self.total_sequences = len(known_bin_seq)
        self.barcodes_shared_memory = shared_memory.SharedMemory(create=True, size=self.total_sequences*self.SEQ_DTYPE().nbytes)
        self.known_bin_seq = numpy.ndarray(shape=(self.total_sequences, ), dtype=self.SEQ_DTYPE, buffer=self.barcodes_shared_memory.buf)
        self.known_bin_seq[:] = known_bin_seq[:]
        self.index_shared_memory = shared_memory.SharedMemory(create=True, size=self.index_size*self.KMER_DTYPE().nbytes)
        self.index = numpy.ndarray(shape=(self.index_size, ), dtype=self.KMER_DTYPE, buffer=self.index_shared_memory.buf)
        self.index_ranges_shared_memory = shared_memory.SharedMemory(create=True, size=(total_kmers+1)*self.INDEX_DTYPE().nbytes)
        self.index_ranges = numpy.ndarray(shape=(total_kmers + 1, ), dtype=self.INDEX_DTYPE, buffer=self.index_ranges_shared_memory.buf)
        self.index_ranges[0] = 0
        i = 0
        index_i = 1
        for l in tmp_index:
            for e in l:
                self.index[i] = e
                i += 1
            self.index_ranges[index_i] = i
            index_i += 1

    def __del__(self):
        self.barcodes_shared_memory.close()
        self.index_shared_memory.close()
        self.index_ranges_shared_memory.close()
        if self.main_instance:
            self.barcodes_shared_memory.unlink()
            self.index_shared_memory.unlink()
            self.index_ranges_shared_memory.unlink()

    @classmethod
    def from_sharable_info(cls, shared_mem_index_info: SharedMemoryIndexInfo):
        kmer_index = cls.__new__(cls)
        kmer_index.main_instance = False
        kmer_index.k = shared_mem_index_info.kmer_size
        kmer_index.mask = (1 << (2 * kmer_index.k)) - 1
        kmer_index.seq_len = shared_mem_index_info.seq_len
        kmer_index.seq_mask = (1 << (2 * kmer_index.seq_len)) - 1
        kmer_index.index_size = shared_mem_index_info.index_size
        kmer_index.total_sequences = shared_mem_index_info.barcode_count
        total_kmers = int(math.pow(4, kmer_index.k))
        kmer_index.barcodes_shared_memory = shared_memory.SharedMemory(create=False, name=shared_mem_index_info.barcodes_sm_name)
        kmer_index.known_bin_seq = numpy.ndarray(shape=(kmer_index.total_sequences, ), dtype=SharedMemoryArray2BitKmerIndexer.SEQ_DTYPE, buffer=kmer_index.barcodes_shared_memory.buf)
        kmer_index.index_shared_memory = shared_memory.SharedMemory(create=False, name=shared_mem_index_info.index_sm_name)
        kmer_index.index = numpy.ndarray(shape=(kmer_index.index_size, ), dtype=kmer_index.KMER_DTYPE, buffer=kmer_index.index_shared_memory.buf)
        kmer_index.index_ranges_shared_memory = shared_memory.SharedMemory(create=False, name=shared_mem_index_info.index_range_sm_name)
        kmer_index.index_ranges = numpy.ndarray(shape=(total_kmers + 1, ), dtype=kmer_index.INDEX_DTYPE, buffer=kmer_index.index_ranges_shared_memory.buf)
        return kmer_index

    def get_sharable_info(self):
        return SharedMemoryIndexInfo(self.total_sequences, self.k, self.seq_len, self.index_size,
                                     self.barcodes_shared_memory.name, self.index_shared_memory.name,
                                     self.index_ranges_shared_memory.name)

    def _get_kmer_bin_indexes(self, bin_seq):
        for i in range(self.seq_len - self.k + 1):
            yield (bin_seq >> ((self.seq_len - self.k - i) * 2)) & self.mask

    def _index(self, known_bin_seq, tmp_index):
        for i, bin_seq in enumerate(known_bin_seq):
            for kmer_idx in self._get_kmer_bin_indexes(bin_seq):
                tmp_index[kmer_idx].append(i)

    def empty(self):
        return len(self.known_bin_seq) == 0

    # @params:
    # sequence: a string to be searched against known strings
    # max_hits: return at most max_hits candidates
    # min_kmers: minimal number of matching k-mers
    # @return
    # a list of (pattern: str, number of shared kmers: int, their positions: list)
    # sorted descending by the number of shared k-m
    def get_occurrences(self, sequence, max_hits=0, min_kmers=1, hits_delta=1, ignore_equal=False):
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
            return  [(bit_to_str(x[0], self.seq_len), x[1], x[2]) for x in result]
        return [(bit_to_str(x[0], self.seq_len), x[1], x[2]) for x in list(result)[:max_hits]]
    