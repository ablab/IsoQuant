############################################################################
# Copyright (c) 2023-2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Curio barcode detectors (formerly DoubleBarcodeDetector).

Curio platform uses double barcodes split by a linker sequence.
"""

import logging
from typing import List, Optional

from ..indexers import KmerIndexer, ArrayKmerIndexer
from ..common import (
    find_polyt_start, reverese_complement,
    find_candidate_with_max_score_ssw, detect_exact_positions,
)
from .base import LinkerBarcodeDetectionResult

logger = logging.getLogger('IsoQuant')


class CurioBarcodeDetector:
    """
    Curio platform barcode detector.

    Detects barcodes split by a linker sequence, with primer and polyT detection.
    """

    LINKER = "TCTTCAGCGTTCCCGAGA"
    PCR_PRIMER = "TACACGACGCTCTTCCGATCT"
    LEFT_BC_LENGTH = 8
    RIGHT_BC_LENGTH = 6
    BC_LENGTH = LEFT_BC_LENGTH + RIGHT_BC_LENGTH
    UMI_LEN = 9
    NON_T_UMI_BASES = 2
    UMI_LEN_DELTA = 2
    TERMINAL_MATCH_DELTA = 2
    STRICT_TERMINAL_MATCH_DELTA = 1

    def __init__(self, barcode_list: List[str], umi_list: Optional[List[str]] = None,
                 min_score: int = 13):
        """
        Initialize Curio barcode detector.

        Args:
            barcode_list: List of known barcodes, or tuple of (left, right) barcode lists
            umi_list: Optional list of known UMIs for validation
            min_score: Minimum alignment score for barcode matching
        """
        self.pcr_primer_indexer = ArrayKmerIndexer([CurioBarcodeDetector.PCR_PRIMER], kmer_size=6)
        self.linker_indexer = ArrayKmerIndexer([CurioBarcodeDetector.LINKER], kmer_size=5)
        joint_barcode_list = []
        if isinstance(barcode_list, tuple):
            # barcodes provided separately
            for b1 in barcode_list[0]:
                for b2 in barcode_list[1]:
                    joint_barcode_list.append(b1 + b2)
        else:
            joint_barcode_list = barcode_list
        self.barcode_indexer = ArrayKmerIndexer(joint_barcode_list, kmer_size=6)
        self.umi_set = None
        if umi_list:
            self.umi_set = set(umi_list)
            logger.debug("Loaded %d UMIs" % len(umi_list))
            self.umi_indexer = KmerIndexer(umi_list, kmer_size=5)
        self.min_score = min_score

    def find_barcode_umi(self, read_id: str, sequence: str) -> LinkerBarcodeDetectionResult:
        """
        Detect barcode and UMI in a read sequence.

        Tries forward strand first, then reverse complement.

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

    def _find_barcode_umi_fwd(self, read_id: str, sequence: str) -> LinkerBarcodeDetectionResult:
        """Detect barcode/UMI in forward orientation."""
        polyt_start = find_polyt_start(sequence)

        linker_start, linker_end = None, None
        if polyt_start != -1:
            # use relaxed parameters if polyA is found
            linker_occurrences = self.linker_indexer.get_occurrences(sequence[0:polyt_start + 1])
            linker_start, linker_end = detect_exact_positions(sequence, 0, polyt_start + 1,
                                                              self.linker_indexer.k, self.LINKER,
                                                              linker_occurrences, min_score=11,
                                                              start_delta=self.TERMINAL_MATCH_DELTA,
                                                              end_delta=self.TERMINAL_MATCH_DELTA)

        if linker_start is None:
            # if polyT was not found, or linker was not found to the left of polyT,
            # look for linker in the entire read
            linker_occurrences = self.linker_indexer.get_occurrences(sequence)
            linker_start, linker_end = detect_exact_positions(sequence, 0, len(sequence),
                                                              self.linker_indexer.k, self.LINKER,
                                                              linker_occurrences, min_score=14,
                                                              start_delta=self.STRICT_TERMINAL_MATCH_DELTA,
                                                              end_delta=self.STRICT_TERMINAL_MATCH_DELTA)

        if linker_start is None:
            return LinkerBarcodeDetectionResult(read_id, polyT=polyt_start)
        logger.debug("LINKER: %d-%d" % (linker_start, linker_end))

        if polyt_start == -1 or polyt_start - linker_end > self.RIGHT_BC_LENGTH + self.UMI_LEN + 10:
            # if polyT was not detected earlier, use relaxed parameters once the linker is found
            presumable_polyt_start = linker_end + self.RIGHT_BC_LENGTH + self.UMI_LEN
            search_start = presumable_polyt_start - 4
            search_end = min(len(sequence), presumable_polyt_start + 10)
            polyt_start = find_polyt_start(sequence[search_start:search_end], window_size=5, polya_fraction=1.0)
            if polyt_start != -1:
                polyt_start += search_start

        primer_occurrences = self.pcr_primer_indexer.get_occurrences(sequence[:linker_start])
        primer_start, primer_end = detect_exact_positions(sequence, 0, linker_start,
                                                          self.pcr_primer_indexer.k, self.PCR_PRIMER,
                                                          primer_occurrences, min_score=5,
                                                          end_delta=self.TERMINAL_MATCH_DELTA)
        if primer_start is not None:
            logger.debug("PRIMER: %d-%d" % (primer_start, primer_end))
        else:
            primer_start = -1
            primer_end = -1

        barcode_start = primer_end + 1 if primer_start != -1 else linker_start - self.LEFT_BC_LENGTH
        barcode_end = linker_end + self.RIGHT_BC_LENGTH + 1
        potential_barcode = sequence[barcode_start:linker_start] + sequence[linker_end + 1:barcode_end + 1]
        logger.debug("Barcode: %s" % (potential_barcode))
        matching_barcodes = self.barcode_indexer.get_occurrences(potential_barcode)
        barcode, bc_score, bc_start, bc_end = \
            find_candidate_with_max_score_ssw(matching_barcodes, potential_barcode, min_score=self.min_score)

        if barcode is None:
            return LinkerBarcodeDetectionResult(read_id, polyT=polyt_start, primer=primer_end,
                                                linker_start=linker_start, linker_end=linker_end)
        logger.debug("Found: %s %d-%d" % (barcode, bc_start, bc_end))
        # position of barcode end in the reference: end of potential barcode minus bases to the alignment end
        read_barcode_end = barcode_start + bc_end - 1 + (linker_end - linker_start + 1)

        potential_umi_start = read_barcode_end + 1
        potential_umi_end = polyt_start - 1
        if potential_umi_end - potential_umi_start <= 5:
            potential_umi_end = potential_umi_start + self.UMI_LEN - 1
        potential_umi = sequence[potential_umi_start:potential_umi_end + 1]
        logger.debug("Potential UMI: %s" % potential_umi)

        umi = None
        good_umi = False
        if self.umi_set:
            matching_umis = self.umi_indexer.get_occurrences(potential_umi)
            umi, umi_score, umi_start, umi_end = \
                find_candidate_with_max_score_ssw(matching_umis, potential_umi, min_score=7)
            logger.debug("Found UMI %s %d-%d" % (umi, umi_start, umi_end))

        if not umi:
            umi = potential_umi
            if self.UMI_LEN - self.UMI_LEN_DELTA <= len(umi) <= self.UMI_LEN + self.UMI_LEN_DELTA and \
                    all(x != "T" for x in umi[-self.NON_T_UMI_BASES:]):
                good_umi = True

        if not umi:
            return LinkerBarcodeDetectionResult(read_id, barcode, BC_score=bc_score,
                                                polyT=polyt_start, primer=primer_end,
                                                linker_start=linker_start, linker_end=linker_end)
        return LinkerBarcodeDetectionResult(read_id, barcode, umi, bc_score, good_umi,
                                            polyT=polyt_start, primer=primer_end,
                                            linker_start=linker_start, linker_end=linker_end)

    @staticmethod
    def result_type():
        return LinkerBarcodeDetectionResult


class CurioIlluminaDetector:
    """
    Curio barcode detector optimized for Illumina short reads.

    Uses different scoring thresholds suitable for higher-quality Illumina data.
    """

    LINKER = "TCTTCAGCGTTCCCGAGA"
    PCR_PRIMER = "TACACGACGCTCTTCCGATCT"
    LEFT_BC_LENGTH = 8
    RIGHT_BC_LENGTH = 6
    BC_LENGTH = LEFT_BC_LENGTH + RIGHT_BC_LENGTH
    MIN_BC_LEN = BC_LENGTH - 4
    UMI_LEN = 9
    NON_T_UMI_BASES = 2
    UMI_LEN_DELTA = 2
    SCORE_DIFF = 1

    TERMINAL_MATCH_DELTA = 1
    STRICT_TERMINAL_MATCH_DELTA = 0

    def __init__(self, joint_barcode_list: List[str], umi_list: Optional[List[str]] = None,
                 min_score: int = 14):
        """
        Initialize Illumina Curio detector.

        Args:
            joint_barcode_list: List of known barcodes
            umi_list: Optional list of known UMIs
            min_score: Minimum alignment score
        """
        self.pcr_primer_indexer = KmerIndexer([CurioBarcodeDetector.PCR_PRIMER], kmer_size=6)
        self.linker_indexer = KmerIndexer([CurioBarcodeDetector.LINKER], kmer_size=5)
        self.barcode_indexer = KmerIndexer(joint_barcode_list, kmer_size=5)
        self.umi_set = None
        if umi_list:
            self.umi_set = set(umi_list)
            logger.debug("Loaded %d UMIs" % len(umi_list))
            self.umi_indexer = KmerIndexer(umi_list, kmer_size=5)
        self.min_score = min_score

    def find_barcode_umi(self, read_id: str, sequence: str) -> LinkerBarcodeDetectionResult:
        """Detect barcode and UMI."""
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

    def _find_barcode_umi_fwd(self, read_id: str, sequence: str) -> LinkerBarcodeDetectionResult:
        # look for linker in the entire read
        linker_occurrences = self.linker_indexer.get_occurrences(sequence)
        linker_start, linker_end = detect_exact_positions(sequence, 0, len(sequence),
                                                          self.linker_indexer.k, self.LINKER,
                                                          linker_occurrences, min_score=14,
                                                          start_delta=self.TERMINAL_MATCH_DELTA,
                                                          end_delta=self.TERMINAL_MATCH_DELTA)

        if linker_start is None:
            return LinkerBarcodeDetectionResult(read_id)
        logger.debug("LINKER: %d-%d" % (linker_start, linker_end))
        primer_end = -1  # forget about primer

        # use relaxed parameters once the linker is found
        presumable_polyt_start = linker_end + self.RIGHT_BC_LENGTH + self.UMI_LEN
        search_start = presumable_polyt_start - 4
        search_end = min(len(sequence), presumable_polyt_start + 10)
        polyt_start = find_polyt_start(sequence[search_start:search_end], window_size=5, polya_fraction=1.0)
        if polyt_start != -1:
            polyt_start += search_start

        barcode_start = linker_start - self.LEFT_BC_LENGTH
        if barcode_start < 0:
            barcode_start = 0
        barcode_end = linker_end + self.RIGHT_BC_LENGTH + 1
        if barcode_end > len(sequence):
            barcode_end = len(sequence)

        potential_barcode = sequence[barcode_start:linker_start] + sequence[linker_end + 1:barcode_end]
        logger.debug("Barcode: %s" % (potential_barcode))
        if len(potential_barcode) < self.MIN_BC_LEN:
            return LinkerBarcodeDetectionResult(read_id, linker_start=linker_start, linker_end=linker_end)
        matching_barcodes = self.barcode_indexer.get_occurrences(potential_barcode)
        barcode, bc_score, bc_start, bc_end = \
            find_candidate_with_max_score_ssw(matching_barcodes, potential_barcode,
                                              min_score=len(potential_barcode) - 1, score_diff=self.SCORE_DIFF)

        if barcode is None:
            return LinkerBarcodeDetectionResult(read_id, polyT=polyt_start, primer=-1,
                                                linker_start=linker_start, linker_end=linker_end)
        logger.debug("Found: %s %d-%d" % (barcode, bc_start, bc_end))
        # position of barcode end in the reference: end of potential barcode minus bases to the alignment end
        read_barcode_end = barcode_start + bc_end - 1 + (linker_end - linker_start + 1)

        potential_umi_start = read_barcode_end + 1
        potential_umi_end = polyt_start - 1
        if polyt_start != -1 or potential_umi_end - potential_umi_start <= 5:
            potential_umi_end = potential_umi_start + self.UMI_LEN - 1
        potential_umi = sequence[potential_umi_start:min(potential_umi_end + 1, len(sequence))]
        logger.debug("Potential UMI: %s" % potential_umi)

        umi = None
        good_umi = False
        if self.umi_set:
            matching_umis = self.umi_indexer.get_occurrences(potential_umi)
            umi, umi_score, umi_start, umi_end = \
                find_candidate_with_max_score_ssw(matching_umis, potential_umi, min_score=7)
            logger.debug("Found UMI %s %d-%d" % (umi, umi_start, umi_end))

        if not umi:
            umi = potential_umi
            if self.UMI_LEN - self.UMI_LEN_DELTA <= len(umi) <= self.UMI_LEN + self.UMI_LEN_DELTA and \
                    all(x != "T" for x in umi[-self.NON_T_UMI_BASES:]):
                good_umi = True

        if not umi:
            return LinkerBarcodeDetectionResult(read_id, barcode, BC_score=bc_score,
                                                polyT=polyt_start, primer=primer_end,
                                                linker_start=linker_start, linker_end=linker_end)
        return LinkerBarcodeDetectionResult(read_id, barcode, umi, bc_score, good_umi,
                                            polyT=polyt_start, primer=primer_end,
                                            linker_start=linker_start, linker_end=linker_end)

    @staticmethod
    def result_type():
        return LinkerBarcodeDetectionResult
