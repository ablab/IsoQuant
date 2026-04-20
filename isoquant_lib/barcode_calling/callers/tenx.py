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
    detect_exact_positions, detect_first_exact_positions, str_to_2bit,
    find_optimal_kmer_size,
)
from .base import TenXBarcodeDetectionResult, TenXSplitBarcodeDetectionResult, SplittingBarcodeDetectionResult

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
        bc_length = len(barcode_list[0])
        bc_count = len(barcode_list)
        if bc_count > 1000000:
            logger.warning("The number of barcodes is large: barcode calling may take substantial amount of time and RAM, and affect the accuracy")
            logger.warning("We suggest to use a sub-list of barcodes derived from short-read analysis whenever possible")

        self.k = find_optimal_kmer_size(bc_length, bc_count)
        if len(barcode_list) <= 100000:
            self.barcode_indexer = ArrayKmerIndexer(barcode_list, kmer_size=self.k)
            self.max_barcodes_hits = 20
            self.min_matching_kmers = 1
            self.min_score = 14
            self.score_diff = 0
        else:
            barcode_bit_list = [str_to_2bit(b) for b in barcode_list]
            self.barcode_indexer = Array2BitKmerIndexer(barcode_bit_list, kmer_size=self.k, seq_len=self.BARCODE_LEN_10X)
            self.max_barcodes_hits = 10
            self.min_matching_kmers = 2
            self.min_score = 16
            self.score_diff = 1

        logger.info("Indexed %d barcodes of length %d with k-mer size %d" % (len(barcode_list), bc_length, self.k))
        logger.info("Minimal alignment score set to %d" % self.min_score)

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

    # Maximum search margin before polyT for direct barcode search (beyond BC+UMI region)
    BARCODE_SEARCH_MARGIN = 30

    def _find_barcode_near_polyt(self, read_id: str, sequence: str,
                                 polyt_start: int) -> TenXBarcodeDetectionResult:
        """Fallback: search for barcode directly in the region before polyT.

        Used when R1 linker is not detected (partial/error-corrupted R1).
        The polyT anchor constrains the barcode position: BC(16)...UMI(12)...polyT.
        """
        search_start = max(0, polyt_start - self.BARCODE_LEN_10X - self.UMI_LEN - self.BARCODE_SEARCH_MARGIN)
        search_region = sequence[search_start:polyt_start]
        if len(search_region) < self.BARCODE_LEN_10X:
            return TenXBarcodeDetectionResult(read_id, polyT=polyt_start)

        matching_barcodes = self.barcode_indexer.get_occurrences(search_region,
                                                                 max_hits=self.max_barcodes_hits,
                                                                 min_kmers=self.min_matching_kmers)
        barcode, bc_score, bc_start, bc_end = \
            find_candidate_with_max_score_ssw(matching_barcodes, search_region,
                                              min_score=self.min_score,
                                              score_diff=self.score_diff)
        if barcode is None:
            return TenXBarcodeDetectionResult(read_id, polyT=polyt_start)

        logger.debug("BARCODE_NEAR_POLYT: %s score=%d at %d-%d" % (barcode, bc_score, search_start + bc_start, search_start + bc_end))
        abs_bc_end = search_start + bc_end
        potential_umi_start = abs_bc_end + 1
        potential_umi_end = polyt_start - 1
        if potential_umi_end - potential_umi_start <= 5:
            potential_umi_end = potential_umi_start + self.UMI_LEN - 1
        potential_umi = sequence[potential_umi_start:potential_umi_end + 1]

        good_umi = self.UMI_LEN - self.UMI_LEN_DELTA <= len(potential_umi) <= self.UMI_LEN + self.UMI_LEN_DELTA
        if not potential_umi:
            return TenXBarcodeDetectionResult(read_id, barcode, BC_score=bc_score, polyT=polyt_start)
        return TenXBarcodeDetectionResult(read_id, barcode, potential_umi, bc_score, good_umi,
                                          polyT=polyt_start)

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
                r1_start, r1_end = detect_exact_positions(sequence, 0, polyt_start + 1,
                                                          self.r1_indexer.k, self.R1,
                                                          r1_occurrences, min_score=8,
                                                          end_delta=self.STRICT_TERMINAL_MATCH_DELTA)

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
            if polyt_start != -1:
                return self._find_barcode_near_polyt(read_id, sequence, polyt_start)
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
                r1_start, r1_end = detect_exact_positions(sequence, 0, polyt_start + 1,
                                                          self.r1_indexer.k, self.R1,
                                                          r1_occurrences, min_score=8,
                                                          end_delta=self.STRICT_TERMINAL_MATCH_DELTA)

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


class TenXSplittingBarcodeDetector(TenXBarcodeDetector):
    """10x Genomics v3 splitting detector — finds multiple barcodes in concatenated reads."""

    MIN_REMAINING_SEQ = 50
    DEFAULT_POLYT_STEP = 100
    MIN_SPLIT_STEP = 150

    def __init__(self, barcode_list: List[str]):
        super().__init__(barcode_list)
        from ..indexers import KmerIndexer as _KmerIndexer
        self.tso_indexer = _KmerIndexer([self.TSO], kmer_size=9)

    # Maximum allowed distance between R1 end and polyT for a valid split detection.
    # Expected: R1...BC(16)...UMI(12)...polyT = ~28bp, allow generous margin.
    MAX_R1_POLYT_DISTANCE = 50

    def _find_barcode_umi_fwd_local(self, read_id: str, sequence: str) -> TenXBarcodeDetectionResult:
        """Find barcode near polyT only — no fallback to full-read R1 search.

        In split mode, the fallback R1 search is wasteful because
        the consistency check rejects far-away R1 matches anyway.
        """
        polyt_start = find_polyt_start(sequence)
        if polyt_start == -1:
            return TenXBarcodeDetectionResult(read_id, polyT=-1)

        r1_occurrences = self.r1_indexer.get_occurrences(sequence[0:polyt_start + 1])
        r1_start, r1_end = detect_exact_positions(sequence, 0, polyt_start + 1,
                                                   self.r1_indexer.k, self.R1,
                                                   r1_occurrences, min_score=11,
                                                   end_delta=self.TERMINAL_MATCH_DELTA)

        if r1_start is None:
            r1_start, r1_end = detect_exact_positions(sequence, 0, polyt_start + 1,
                                                      self.r1_indexer.k, self.R1,
                                                      r1_occurrences, min_score=8,
                                                      end_delta=self.STRICT_TERMINAL_MATCH_DELTA)

        if r1_start is None:
            return TenXBarcodeDetectionResult(read_id, polyT=polyt_start)

        barcode_start = r1_end + 1
        barcode_end = r1_end + self.BARCODE_LEN_10X + 1
        potential_barcode = sequence[barcode_start:barcode_end + 1]
        matching_barcodes = self.barcode_indexer.get_occurrences(potential_barcode,
                                                                 max_hits=self.max_barcodes_hits,
                                                                 min_kmers=self.min_matching_kmers)
        barcode, bc_score, bc_start, bc_end = \
            find_candidate_with_max_score_ssw(matching_barcodes, potential_barcode,
                                              min_score=self.min_score,
                                              score_diff=self.score_diff)

        if barcode is None:
            return TenXBarcodeDetectionResult(read_id, polyT=polyt_start, r1=r1_end)

        read_barcode_end = barcode_start + bc_end - 1
        potential_umi_start = read_barcode_end + 1
        potential_umi_end = polyt_start - 1
        if potential_umi_end - potential_umi_start <= 5:
            potential_umi_end = potential_umi_start + self.UMI_LEN - 1
        potential_umi = sequence[potential_umi_start:potential_umi_end + 1]

        good_umi = self.UMI_LEN - self.UMI_LEN_DELTA <= len(potential_umi) <= self.UMI_LEN + self.UMI_LEN_DELTA
        if not potential_umi:
            return TenXBarcodeDetectionResult(read_id, barcode, BC_score=bc_score, polyT=polyt_start, r1=r1_end)
        return TenXBarcodeDetectionResult(read_id, barcode, potential_umi, bc_score, good_umi,
                                          polyT=polyt_start, r1=r1_end)

    def _find_barcode_umi_split_fwd(self, read_id: str, sequence: str,
                                     offset: int = 0) -> TenXSplitBarcodeDetectionResult:
        """Detect one barcode+UMI+TSO pattern in the forward-oriented subsequence."""
        base_result = self._find_barcode_umi_fwd_local(read_id, sequence)

        tso_start = -1
        if base_result.polyT != -1 and base_result.is_valid():
            tso_search_start = base_result.polyT
            tso_search_end = len(sequence) - 1
            tso_occurrences = self.tso_indexer.get_occurrences_substr(
                sequence, tso_search_start, tso_search_end
            )
            tso_s, tso_e = detect_first_exact_positions(
                sequence, tso_search_start, tso_search_end + 1,
                self.tso_indexer.k, self.TSO,
                tso_occurrences,
                min_score=18, start_delta=3, end_delta=3
            )
            if tso_s is not None:
                tso_start = tso_s

        result = TenXSplitBarcodeDetectionResult(
            read_id,
            barcode=base_result.get_barcode(),
            UMI=base_result.get_umi(),
            BC_score=base_result.BC_score,
            UMI_good=base_result.UMI_good,
            polyT=base_result.polyT,
            r1=base_result.r1,
            tso=tso_start,
        )
        if offset > 0:
            result.update_coordinates(offset)
        return result

    def _is_consistent_detection(self, r: 'TenXSplitBarcodeDetectionResult') -> bool:
        """Check if polyT and R1 positions are consistent (belong to the same molecule).

        The fallback R1 search in _find_barcode_umi_fwd can match R1 from a
        different molecule far away when a spurious polyT is found in cDNA.
        Reject these by checking that R1 and polyT are within expected distance.
        """
        if r.r1 == -1 or r.polyT == -1 or not r.is_valid():
            return True  # no barcode found, nothing to reject
        return r.r1 < r.polyT and r.polyT - r.r1 <= self.MAX_R1_POLYT_DISTANCE

    def _scan_strand(self, read_id: str, sequence: str, strand: str,
                     read_result: SplittingBarcodeDetectionResult) -> None:
        """Scan one strand for barcode patterns, appending results."""
        current_start = 0
        prev_start = -1
        while True:
            seq = sequence[current_start:]
            r = self._find_barcode_umi_split_fwd(read_id, seq, offset=current_start)
            if r.polyT == -1:
                break
            r.set_strand(strand)

            if self._is_consistent_detection(r):
                read_result.append(r)
                if r.tso != -1:
                    next_start = r.tso + len(self.TSO)
                else:
                    next_start = r.polyT + self.DEFAULT_POLYT_STEP
            else:
                # polyT and R1 are from different molecules — skip to before R1
                # so the next iteration can detect the molecule with proper R1/polyT proximity
                next_start = r.r1 - len(self.R1) - 10

            next_start = max(prev_start + self.MIN_SPLIT_STEP, next_start)
            prev_start = next_start
            current_start = next_start
            if len(sequence) - current_start < self.MIN_REMAINING_SEQ:
                break

    def find_barcode_umi(self, read_id: str, sequence: str) -> SplittingBarcodeDetectionResult:
        """Find all barcode+UMI patterns in a (potentially concatenated) read.

        Scans both forward and reverse strands, collecting molecules from both
        into a single result (molecules in concatenated reads can alternate orientation).
        """
        read_result = SplittingBarcodeDetectionResult(read_id)

        # Forward strand scan
        self._scan_strand(read_id, sequence, "+", read_result)

        # Reverse strand scan — collect into the same result
        rev_seq = reverese_complement(sequence)
        self._scan_strand(read_id, rev_seq, "-", read_result)

        # If no valid barcodes found, try polyT-anchored barcode search as last resort
        # (handles reads with partial/corrupted R1 where k-mer seeding fails)
        has_valid = any(r.is_valid() for r in read_result.detected_patterns)
        if not has_valid:
            for try_seq, strand in [(sequence, "+"), (rev_seq, "-")]:
                polyt = find_polyt_start(try_seq)
                if polyt != -1:
                    base_result = self._find_barcode_near_polyt(read_id, try_seq, polyt)
                    if base_result.is_valid():
                        split_result = TenXSplitBarcodeDetectionResult(
                            read_id, base_result.get_barcode(), base_result.get_umi(),
                            base_result.BC_score, base_result.UMI_good,
                            polyT=base_result.polyT, r1=-1, tso=-1)
                        split_result.set_strand(strand)
                        read_result.detected_patterns = [split_result]
                        break

        if read_result.empty():
            r = self._find_barcode_umi_split_fwd(read_id, sequence)
            read_result.append(r)

        read_result.filter()
        return read_result

    @staticmethod
    def result_type():
        return TenXSplitBarcodeDetectionResult

    @classmethod
    def header(cls):
        return cls.result_type().header()


class TenXv2SplittingBarcodeDetector(TenXSplittingBarcodeDetector):
    """10x Genomics v2 splitting detector."""

    UMI_LEN = 10
    # RC of v2 TSO oligo (AAGCAGTGGTATCAACGCAGAGTAC) as it appears in the read after polyT
    TSO = "GTACTCTGCGTTGATACCACTGCTT"

    def __init__(self, barcode_list: List[str]):
        super().__init__(barcode_list)
