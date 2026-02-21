############################################################################
# Copyright (c) 2023-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
10x Genomics barcode detectors.

Supports 10x Genomics v3 single-cell and Visium HD spatial platforms.
"""

import logging
from typing import List, Optional

from .. import ArrayKmerIndexer, Array2BitKmerIndexer
from ..indexers import KmerIndexer
from ..common import (
    find_polyt_start, reverese_complement,
    find_candidate_with_max_score_ssw, find_candidate_with_max_score_ssw_var_len,
    detect_exact_positions, str_to_2bit,
)
from .base import TenXBarcodeDetectionResult

logger = logging.getLogger('IsoQuant')


class TenXBarcodeDetector:
    """10x Genomics v3 single-cell barcode detector."""

    TSO = "CCCATGTACTCTGCGTTGATACCACTGCTT"
    R1 = "CTACACGACGCTCTTCCGATCT"  # 10x 3'
    BARCODE_LEN_10X = 16
    UMI_LEN = 12

    UMI_LEN_DELTA = 2
    TERMINAL_MATCH_DELTA = 2
    STRICT_TERMINAL_MATCH_DELTA = 1

    def __init__(self, barcode_list: List[str]):
        """
        Initialize 10x detector.

        Args:
            barcode_list: List of known barcodes (whitelist)
            umi_list: Optional list of known UMIs for validation
        """
        self.r1_indexer = KmerIndexer([TenXBarcodeDetector.R1], kmer_size=7)
        if len(barcode_list) <= 100000:
            self.barcode_indexer = ArrayKmerIndexer(barcode_list, kmer_size=6)
            self.max_barcodes_hits = 20
            self.min_matching_kmers = 1
            self.min_score = 14
            self.score_diff = 0
        else:
            barcode_bit_list = [str_to_2bit(b) for b in barcode_list]
            self.barcode_indexer = Array2BitKmerIndexer(barcode_bit_list, kmer_size=9, seq_len=self.BARCODE_LEN_10X)
            self.small_barcode_indexer = Array2BitKmerIndexer(barcode_bit_list, kmer_size=6, seq_len=self.BARCODE_LEN_10X)
            self.max_barcodes_hits = 10
            self.min_matching_kmers = 2
            self.min_score = 15
            self.score_diff = 1

        logger.debug("Min score set to %d" % self.min_score)

    def find_barcode_umi(self, read_id: str, sequence: str) -> TenXBarcodeDetectionResult:
        """
        Detect barcode and UMI in a read sequence.

        Args:
            read_id: Read identifier
            sequence: Read sequence

        Returns:
            Detection result with barcode, UMI, and feature positions
        """
        read_result = self._find_barcode_umi_fwd(read_id, sequence)
        if read_result.polyT != -1:
            read_result.set_strand("+")
        if read_result.is_valid():
            return read_result

        rev_seq = reverese_complement(sequence)
        read_rev_result = self._find_barcode_umi_fwd(read_id, rev_seq)
        if read_rev_result.polyT != -1:
            read_rev_result.set_strand("-")
        if read_rev_result.is_valid():
            return read_rev_result

        return read_result if read_result.more_informative_than(read_rev_result) else read_rev_result

    def _find_barcode_umi_fwd(self, read_id: str, sequence: str) -> TenXBarcodeDetectionResult:
        polyt_start = find_polyt_start(sequence)

        r1_start, r1_end = None, None
        if polyt_start != -1:
            # use relaxed parameters if polyA is found
            r1_occurrences = self.r1_indexer.get_occurrences(sequence[0:polyt_start + 1])
            r1_start, r1_end = detect_exact_positions(sequence, 0, polyt_start + 1,
                                                      self.r1_indexer.k, self.R1,
                                                      r1_occurrences, min_score=11,
                                                      end_delta=self.TERMINAL_MATCH_DELTA)

        if r1_start is None:
            # if polyT was not found, or linker was not found to the left of polyT,
            # look for linker in the entire read
            r1_occurrences = self.r1_indexer.get_occurrences(sequence)
            r1_start, r1_end = detect_exact_positions(sequence, 0, len(sequence),
                                                      self.r1_indexer.k, self.R1,
                                                      r1_occurrences, min_score=18,
                                                      start_delta=self.STRICT_TERMINAL_MATCH_DELTA,
                                                      end_delta=self.STRICT_TERMINAL_MATCH_DELTA)

        if r1_start is None:
            return TenXBarcodeDetectionResult(read_id, polyT=polyt_start)
        logger.debug("LINKER: %d-%d" % (r1_start, r1_end))

        if polyt_start == -1 or polyt_start - r1_end > self.BARCODE_LEN_10X + self.UMI_LEN + 10:
            # if polyT was not detected earlier, use relaxed parameters once the linker is found
            presumable_polyt_start = r1_end + self.BARCODE_LEN_10X + self.UMI_LEN
            search_start = presumable_polyt_start - 4
            search_end = min(len(sequence), presumable_polyt_start + 10)
            polyt_start = find_polyt_start(sequence[search_start:search_end], window_size=5, polya_fraction=1.0)
            if polyt_start != -1:
                polyt_start += search_start

        barcode_start = r1_end + 1
        barcode_end = r1_end + self.BARCODE_LEN_10X + 1
        potential_barcode = sequence[barcode_start:barcode_end + 1]
        logger.debug("Barcode: %s" % (potential_barcode))
        matching_barcodes = self.barcode_indexer.get_occurrences(potential_barcode,
                                                                 max_hits=self.max_barcodes_hits,
                                                                 min_kmers=self.min_matching_kmers)
        barcode, bc_score, bc_start, bc_end = \
            find_candidate_with_max_score_ssw(matching_barcodes, potential_barcode,
                                              min_score=self.min_score,
                                              score_diff=self.score_diff)

        if barcode is None:
            matching_barcodes = self.small_barcode_indexer.get_occurrences(potential_barcode,
                                                                           max_hits=self.max_barcodes_hits,
                                                                           min_kmers=self.min_matching_kmers)
            barcode, bc_score, bc_start, bc_end = \
                find_candidate_with_max_score_ssw(matching_barcodes, potential_barcode,
                                                  min_score=self.min_score,
                                                  score_diff=self.score_diff)
            if barcode is None:
                return TenXBarcodeDetectionResult(read_id, polyT=polyt_start, r1=r1_end)
        logger.debug("Found: %s %d-%d" % (barcode, bc_start, bc_end))
        # position of barcode end in the reference: end of potential barcode minus bases to the alignment end
        read_barcode_end = barcode_start + bc_end - 1
        potential_umi_start = read_barcode_end + 1
        potential_umi_end = polyt_start - 1
        if potential_umi_end - potential_umi_start <= 5:
            potential_umi_end = potential_umi_start + self.UMI_LEN - 1
        potential_umi = sequence[potential_umi_start:potential_umi_end + 1]
        logger.debug("Potential UMI: %s" % potential_umi)

        good_umi = False
        umi = potential_umi
        if self.UMI_LEN - self.UMI_LEN_DELTA <= len(umi) <= self.UMI_LEN + self.UMI_LEN_DELTA:
            good_umi = True

        if not umi:
            return TenXBarcodeDetectionResult(read_id, barcode, BC_score=bc_score, polyT=polyt_start, r1=r1_end)
        return TenXBarcodeDetectionResult(read_id, barcode, umi, bc_score, good_umi, polyT=polyt_start, r1=r1_end)

    def find_barcode_umi_no_polya(self, read_id: str, sequence: str) -> TenXBarcodeDetectionResult:
        """Find barcode without polyA requirement."""
        read_result = self._find_barcode_umi_fwd(read_id, sequence)
        if read_result.polyT != -1:
            read_result.set_strand("+")
        if read_result.is_valid():
            return read_result

        rev_seq = reverese_complement(sequence)
        read_rev_result = self._find_barcode_umi_fwd(read_id, rev_seq)
        if read_rev_result.polyT != -1:
            read_rev_result.set_strand("-")
        if read_rev_result.is_valid():
            return read_rev_result

        return read_result if read_result.more_informative_than(read_rev_result) else read_rev_result

    @staticmethod
    def result_type():
        return TenXBarcodeDetectionResult

    @classmethod
    def header(cls):
        return cls.result_type().header()


class TenXv2BarcodeDetector(TenXBarcodeDetector):
    """TenX v2 barcode detector."""
    UMI_LEN = 10

    def __init__(self, barcode_list: List[str]):
        """Initialize TenX v2 detector."""
        super().__init__(barcode_list)


class VisiumHDBarcodeDetector:
    """Visium HD spatial transcriptomics barcode detector."""

    R1 = "CTACACGACGCTCTTCCGATCT"  # 10x 3'
    BARCODE1_LEN_VIS = 16
    BARCODE2_LEN_VIS = 15
    SEPARATOR_BASES = 2
    TOTAL_BARCODE_LEN_VIS = BARCODE1_LEN_VIS + BARCODE2_LEN_VIS
    UMI_LEN_VIS = 9

    UMI_LEN_DELTA = 2
    TERMINAL_MATCH_DELTA = 2
    STRICT_TERMINAL_MATCH_DELTA = 1

    def __init__(self, barcode_pair_list: List[List[str]]):
        """
        Initialize Visium HD detector.

        Args:
            barcode_pair_list: List of two lists [part1_barcodes, part2_barcodes]
        """
        assert len(barcode_pair_list) == 2
        self.r1_indexer = KmerIndexer([VisiumHDBarcodeDetector.R1], kmer_size=7)
        self.part1_list = barcode_pair_list[0]
        self.part2_list = barcode_pair_list[1]

        self.part1_barcode_indexer = ArrayKmerIndexer(self.part1_list, kmer_size=7)
        self.part2_barcode_indexer = ArrayKmerIndexer(self.part2_list, kmer_size=7)
        self.max_barcodes_hits = 20
        self.min_matching_kmers = 2
        self.min_score = 13

        logger.debug("Min score set to %d" % self.min_score)

    def find_barcode_umi(self, read_id: str, sequence: str) -> TenXBarcodeDetectionResult:
        """Detect barcode and UMI in a read sequence."""
        read_result = self._find_barcode_umi_fwd(read_id, sequence)
        if read_result.polyT != -1:
            read_result.set_strand("+")
        if read_result.is_valid():
            return read_result

        rev_seq = reverese_complement(sequence)
        read_rev_result = self._find_barcode_umi_fwd(read_id, rev_seq)
        if read_rev_result.polyT != -1:
            read_rev_result.set_strand("-")
        if read_rev_result.is_valid():
            return read_rev_result

        return read_result if read_result.more_informative_than(read_rev_result) else read_rev_result

    def _find_barcode_umi_fwd(self, read_id: str, sequence: str) -> TenXBarcodeDetectionResult:
        logger.debug("===== " + read_id)
        polyt_start = find_polyt_start(sequence)

        r1_start, r1_end = None, None
        if polyt_start != -1:
            # use relaxed parameters if polyA is found
            r1_occurrences = self.r1_indexer.get_occurrences(sequence[0:polyt_start + 1])
            r1_start, r1_end = detect_exact_positions(sequence, 0, polyt_start + 1,
                                                      self.r1_indexer.k, self.R1,
                                                      r1_occurrences, min_score=10,
                                                      end_delta=self.TERMINAL_MATCH_DELTA)

        if r1_start is None:
            # if polyT was not found, or linker was not found to the left of polyT,
            # look for linker in the entire read
            r1_occurrences = self.r1_indexer.get_occurrences(sequence)
            r1_start, r1_end = detect_exact_positions(sequence, 0, len(sequence),
                                                      self.r1_indexer.k, self.R1,
                                                      r1_occurrences, min_score=17,
                                                      start_delta=self.STRICT_TERMINAL_MATCH_DELTA,
                                                      end_delta=self.STRICT_TERMINAL_MATCH_DELTA)

        if r1_start is not None:
            logger.debug("PRIMER: %d-%d" % (r1_start, r1_end))

        if r1_end is not None and (polyt_start == -1 or polyt_start - r1_end > self.TOTAL_BARCODE_LEN_VIS + self.UMI_LEN_VIS + 10):
            # if polyT was not detected earlier, use relaxed parameters once the linker is found
            presumable_polyt_start = r1_end + self.TOTAL_BARCODE_LEN_VIS + self.UMI_LEN_VIS
            search_start = presumable_polyt_start - 4
            search_end = min(len(sequence), presumable_polyt_start + 10)
            polyt_start = find_polyt_start(sequence[search_start:search_end], window_size=5, polya_fraction=1.0)
            if polyt_start != -1:
                polyt_start += search_start

        if polyt_start == -1:
            if r1_start is None:
                return TenXBarcodeDetectionResult(read_id, polyT=polyt_start)
            # no polyT, start from the left
            potential_umi_start = r1_end + 1
            potential_umi_end = potential_umi_start + self.UMI_LEN_VIS - 1
            potential_umi = sequence[potential_umi_start:potential_umi_end + 1]
            logger.debug("Potential UMI: %s" % potential_umi)

            barcode1_start = r1_end + self.UMI_LEN_VIS + 1
            barcode1_end = barcode1_start + self.BARCODE1_LEN_VIS - 1
            potential_barcode1 = sequence[barcode1_start:barcode1_end + 1]
            matching_barcodes1 = self.part1_barcode_indexer.get_occurrences(potential_barcode1, max_hits=self.max_barcodes_hits)
            barcode1, bc1_score, bc1_start, bc1_end = \
                find_candidate_with_max_score_ssw_var_len(matching_barcodes1, potential_barcode1, min_score=self.min_score)
            logger.debug("Barcode 1: %s, %s" % (potential_barcode1, barcode1))
            real_bc1_end = barcode1_start + bc1_end

            barcode2_start = real_bc1_end + 1
            barcode2_end = barcode2_start + self.BARCODE2_LEN_VIS - 1
            potential_barcode2 = sequence[barcode2_start:barcode2_end + 1]
            matching_barcodes2 = self.part2_barcode_indexer.get_occurrences(potential_barcode2, max_hits=self.max_barcodes_hits)
            barcode2, bc2_score, bc2_start, bc2_end = \
                find_candidate_with_max_score_ssw_var_len(matching_barcodes2, potential_barcode2, min_score=self.min_score)
            logger.debug("Barcode 2: %s, %s" % (potential_barcode2, barcode2))

            if barcode1 is None or barcode2 is None:
                return TenXBarcodeDetectionResult(read_id, polyT=polyt_start, r1=r1_end)

            return TenXBarcodeDetectionResult(read_id, barcode1 + "|" + barcode2, potential_umi, bc1_score + bc2_score,
                                              UMI_good=True, polyT=polyt_start, r1=r1_end)

        barcode2_end = polyt_start - 1 - self.SEPARATOR_BASES
        barcode2_start = barcode2_end - self.BARCODE2_LEN_VIS + 1
        potential_barcode2 = sequence[barcode2_start:barcode2_end + 1]
        matching_barcodes2 = self.part2_barcode_indexer.get_occurrences(potential_barcode2, max_hits=self.max_barcodes_hits)
        barcode2, bc2_score, bc2_start, bc2_end = \
            find_candidate_with_max_score_ssw_var_len(matching_barcodes2, potential_barcode2, min_score=self.min_score)
        logger.debug("Barcode 2: %s, %s" % (potential_barcode2, barcode2))

        real_bc2_start = barcode2_start + bc2_start
        barcode1_end = real_bc2_start - 1
        barcode1_start = barcode1_end - self.BARCODE1_LEN_VIS + 1
        potential_barcode1 = sequence[barcode1_start:barcode1_end + 1]
        matching_barcodes1 = self.part1_barcode_indexer.get_occurrences(potential_barcode1, max_hits=self.max_barcodes_hits)
        barcode1, bc1_score, bc1_start, bc1_end = \
            find_candidate_with_max_score_ssw_var_len(matching_barcodes1, potential_barcode1, min_score=self.min_score)
        logger.debug("Barcode 1: %s, %s" % (potential_barcode1, barcode1))
        real_bc1_start = barcode1_start + bc1_start

        potential_umi_end = real_bc1_start - 1
        if r1_end is not None:
            potential_umi_start = r1_end + 1
        else:
            potential_umi_start = max(0, potential_umi_end - self.UMI_LEN_VIS)
        umi_good = abs(potential_umi_end - potential_umi_start + 1 - self.UMI_LEN_VIS) <= self.UMI_LEN_DELTA
        potential_umi = sequence[potential_umi_start:potential_umi_end + 1]
        logger.debug("Potential UMI: %s" % potential_umi)

        if barcode1 is None or barcode2 is None:
            return TenXBarcodeDetectionResult(read_id, polyT=polyt_start, r1=r1_end if r1_end is not None else -1)

        return TenXBarcodeDetectionResult(read_id, barcode1 + "|" + barcode2, potential_umi, bc1_score + bc2_score,
                                          UMI_good=umi_good, polyT=polyt_start, r1=r1_end if r1_end is not None else -1)

    def find_barcode_umi_no_polya(self, read_id: str, sequence: str) -> TenXBarcodeDetectionResult:
        """Find barcode without polyA requirement."""
        read_result = self._find_barcode_umi_fwd(read_id, sequence)
        if read_result.polyT != -1:
            read_result.set_strand("+")
        if read_result.is_valid():
            return read_result

        rev_seq = reverese_complement(sequence)
        read_rev_result = self._find_barcode_umi_fwd(read_id, rev_seq)
        if read_rev_result.polyT != -1:
            read_rev_result.set_strand("-")
        if read_rev_result.is_valid():
            return read_rev_result

        return read_result if read_result.more_informative_than(read_rev_result) else read_rev_result

    @staticmethod
    def result_type():
        return TenXBarcodeDetectionResult

    @classmethod
    def header(cls):
        return cls.result_type().header()
