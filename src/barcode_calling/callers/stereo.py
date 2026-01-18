############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Stereo-seq barcode detectors.

Stereo-seq uses spatial barcodes with linker sequences and optional TSO detection.
"""

import logging
from typing import List, Optional

from ..indexers import ArrayKmerIndexer, Dict2BitKmerIndexer, SharedMemoryArray2BitKmerIndexer
from ..common import (
    find_polyt_start, reverese_complement, batch_str_to_2bit_chunked,
    find_candidate_with_max_score_ssw, detect_exact_positions, detect_first_exact_positions,
)
from .base import LinkerBarcodeDetectionResult, TSOBarcodeDetectionResult, SplittingBarcodeDetectionResult

logger = logging.getLogger('IsoQuant')


class StereoBarcodeDetector:
    """Stereo-seq barcode detector for standard (non-splitting) mode."""

    LINKER = "TTGTCTTCCTAAGAC"
    TSO_PRIMER = "ACTGAGAGGCATGGCGACCTTATCAG"
    PC1_PRIMER = "CTTCCGATCTATGGCGACCTTATCAG"
    BC_LENGTH = 25
    UMI_LEN = 10
    NON_T_UMI_BASES = 0
    UMI_LEN_DELTA = 4
    TERMINAL_MATCH_DELTA = 3
    STRICT_TERMINAL_MATCH_DELTA = 1

    def __init__(self, barcodes: List[str], min_score: int = 21):
        """
        Initialize Stereo-seq detector.

        Args:
            barcodes: List of known barcode sequences
            min_score: Minimum alignment score for barcode matching
        """
        self.main_primer = StereoBarcodeDetector.PC1_PRIMER
        self.pcr_primer_indexer = ArrayKmerIndexer([self.main_primer], kmer_size=6)
        self.linker_indexer = ArrayKmerIndexer([StereoBarcodeDetector.LINKER], kmer_size=5)
        self.strict_linker_indexer = ArrayKmerIndexer([StereoBarcodeDetector.LINKER], kmer_size=6)

        self.barcode_indexer = None
        if barcodes:
            bit_barcodes = batch_str_to_2bit_chunked(iter(barcodes), seq_len=self.BC_LENGTH)
            self.barcode_indexer = Dict2BitKmerIndexer(bit_barcodes, kmer_size=14, seq_len=self.BC_LENGTH)
            logger.info("Indexed %d barcodes" % len(bit_barcodes))

        self.umi_set = None
        self.min_score = min_score

    def find_barcode_umi_multiple(self, read_id: str, sequence: str) -> List[LinkerBarcodeDetectionResult]:
        """Find multiple barcodes in a single read."""
        read_result = []
        r = self._find_barcode_umi_fwd(read_id, sequence)
        current_start = 0
        while r.polyT != -1:
            r.set_strand("+")
            read_result.append(r)
            new_start = r.polyT + 50
            current_start += new_start
            if len(sequence) - current_start < 50:
                break
            seq = sequence[current_start:]
            new_id = read_id + "_%d" % current_start
            r = self._find_barcode_umi_fwd(new_id, seq)

        rev_seq = reverese_complement(sequence)
        read_id += "_R"
        rr = self._find_barcode_umi_fwd(read_id, rev_seq)
        current_start = 0
        while rr.polyT != -1:
            rr.set_strand("-")
            read_result.append(rr)
            new_start = rr.polyT + 50
            current_start += new_start
            if len(rev_seq) - current_start < 50:
                break
            seq = rev_seq[current_start:]
            new_id = read_id + "_%d" % current_start
            rr = self._find_barcode_umi_fwd(new_id, seq)

        if not read_result:
            read_result.append(r)
        return read_result

    def find_barcode_umi(self, read_id: str, sequence: str) -> LinkerBarcodeDetectionResult:
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

    def _find_barcode_umi_fwd(self, read_id: str, sequence: str) -> LinkerBarcodeDetectionResult:
        polyt_start = find_polyt_start(sequence)

        linker_start, linker_end = None, None
        if polyt_start != -1:
            # use relaxed parameters if polyA is found
            linker_occurrences = self.linker_indexer.get_occurrences(sequence[0:polyt_start + 1])
            linker_start, linker_end = detect_exact_positions(sequence, 0, polyt_start + 1,
                                                              self.linker_indexer.k, StereoBarcodeDetector.LINKER,
                                                              linker_occurrences, min_score=12,
                                                              start_delta=self.TERMINAL_MATCH_DELTA,
                                                              end_delta=self.TERMINAL_MATCH_DELTA)

        if linker_start is None:
            # if polyT was not found, or linker was not found to the left of polyT,
            # look for linker in the entire read
            linker_occurrences = self.strict_linker_indexer.get_occurrences(sequence)
            linker_start, linker_end = detect_exact_positions(sequence, 0, len(sequence),
                                                              self.linker_indexer.k, StereoBarcodeDetector.LINKER,
                                                              linker_occurrences, min_score=12,
                                                              start_delta=self.STRICT_TERMINAL_MATCH_DELTA,
                                                              end_delta=self.STRICT_TERMINAL_MATCH_DELTA)

        if linker_start is None:
            return LinkerBarcodeDetectionResult(read_id, polyT=polyt_start)
        logger.debug("LINKER: %d-%d" % (linker_start, linker_end))

        if polyt_start == -1:
            # if polyT was not detected earlier, use relaxed parameters once the linker is found
            presumable_polyt_start = linker_end + self.UMI_LEN
            search_start = presumable_polyt_start - 4
            search_end = min(len(sequence), presumable_polyt_start + 10)
            polyt_start = find_polyt_start(sequence[search_start:search_end], window_size=5, polya_fraction=1.0)
            if polyt_start != -1:
                polyt_start += search_start

        primer_occurrences = self.pcr_primer_indexer.get_occurrences(sequence[:linker_start])
        primer_start, primer_end = detect_exact_positions(sequence, 0, linker_start,
                                                          self.pcr_primer_indexer.k, self.main_primer,
                                                          primer_occurrences, min_score=12,
                                                          end_delta=self.TERMINAL_MATCH_DELTA)
        if primer_start is not None:
            logger.debug("PRIMER: %d-%d" % (primer_start, primer_end))
        else:
            primer_start = -1
            primer_end = linker_start - self.BC_LENGTH - 1

        if primer_end < 0:
            return LinkerBarcodeDetectionResult(read_id, polyT=polyt_start, primer=-1,
                                                linker_start=linker_start, linker_end=linker_end)

        barcode_start = primer_end + 1
        barcode_end = linker_start - 1
        bc_len = barcode_end - barcode_start
        if abs(bc_len - self.BC_LENGTH) > 10:
            return LinkerBarcodeDetectionResult(read_id, polyT=polyt_start, primer=primer_end,
                                                linker_start=linker_start, linker_end=linker_end)

        potential_barcode = sequence[barcode_start:barcode_end + 1]
        logger.debug("Barcode: %s" % (potential_barcode))
        matching_barcodes = self.barcode_indexer.get_occurrences(potential_barcode, max_hits=10, min_kmers=2)
        barcode, bc_score, bc_start, bc_end = \
            find_candidate_with_max_score_ssw(matching_barcodes, potential_barcode,
                                              min_score=self.min_score, sufficient_score=self.BC_LENGTH - 1)

        if barcode is None:
            return LinkerBarcodeDetectionResult(read_id, polyT=polyt_start, primer=primer_end,
                                                linker_start=linker_start, linker_end=linker_end)
        logger.debug("Found: %s %d-%d" % (barcode, bc_start, bc_end))

        potential_umi_start = linker_end + 1
        potential_umi_end = polyt_start - 1
        if potential_umi_end - potential_umi_start <= 5:
            potential_umi_end = potential_umi_start + self.UMI_LEN - 1
        umi = None
        good_umi = False
        if potential_umi_start + 2 * self.UMI_LEN > potential_umi_end > potential_umi_start:
            umi = sequence[potential_umi_start:potential_umi_end + 1]
            logger.debug("Potential UMI: %s" % umi)
            good_umi = abs(len(umi) - self.UMI_LEN) <= self.UMI_LEN_DELTA

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


class SharedMemoryStereoBarcodeDetector(StereoBarcodeDetector):
    """Stereo-seq detector with shared memory for large barcode sets."""

    MIN_BARCODES_FOR_SHARED_MEM = 1000000

    def __init__(self, barcodes: List[str], min_score: int = 21):
        super().__init__([], min_score=min_score)
        # Convert to 2-bit encoding first (memory efficient)
        bit_barcodes = batch_str_to_2bit_chunked(iter(barcodes), seq_len=self.BC_LENGTH)
        self.barcode_count = len(bit_barcodes)
        self.bit_barcodes = None
        if self.barcode_count < self.MIN_BARCODES_FOR_SHARED_MEM:
            # For small/medium sets, use Dict2BitKmerIndexer
            self.bit_barcodes = bit_barcodes
            self.barcode_indexer = Dict2BitKmerIndexer(bit_barcodes, kmer_size=14, seq_len=self.BC_LENGTH)
        else:
            # For large sets, use shared memory indexer
            self.barcode_indexer = SharedMemoryArray2BitKmerIndexer(bit_barcodes, kmer_size=14,
                                                                    seq_len=self.BC_LENGTH)

        logger.info("Indexed %d barcodes" % self.barcode_count)

    def __getstate__(self):
        if self.barcode_count < self.MIN_BARCODES_FOR_SHARED_MEM:
            return (self.min_score,
                    self.barcode_count,
                    self.bit_barcodes)
        else:
            return (self.min_score,
                    self.barcode_count,
                    self.barcode_indexer.get_sharable_info())

    def __setstate__(self, state):
        self.min_score = state[0]
        super().__init__([], min_score=self.min_score)
        self.bit_barcodes = None
        self.barcode_count = state[1]
        if self.barcode_count < self.MIN_BARCODES_FOR_SHARED_MEM:
            self.bit_barcodes = state[2]
            self.barcode_indexer = Dict2BitKmerIndexer(self.bit_barcodes, kmer_size=14, seq_len=self.BC_LENGTH)
        else:
            self.barcode_indexer = SharedMemoryArray2BitKmerIndexer.from_sharable_info(state[2])


class StereoSplittingBarcodeDetector:
    """Stereo-seq detector with read splitting for concatenated reads."""

    TSO5 = "CCCGCCTCTCAGTACGTCAGCAG"
    LINKER = "TTGTCTTCCTAAGAC"
    TSO_PRIMER = "ACTGAGAGGCATGGCGACCTTATCAG"
    PC1_PRIMER = "CTTCCGATCTATGGCGACCTTATCAG"
    BC_LENGTH = 25
    UMI_LEN = 10
    NON_T_UMI_BASES = 0
    UMI_LEN_DELTA = 3
    TERMINAL_MATCH_DELTA = 3
    STRICT_TERMINAL_MATCH_DELTA = 1

    def __init__(self, barcodes: List[str], min_score: int = 21):
        """
        Initialize splitting detector.

        Args:
            barcodes: List of known barcode sequences
            min_score: Minimum alignment score
        """
        self.main_primer = self.PC1_PRIMER
        self.tso5_indexer = ArrayKmerIndexer([self.TSO5], kmer_size=8)
        self.pcr_primer_indexer = ArrayKmerIndexer([self.main_primer], kmer_size=6)
        self.linker_indexer = ArrayKmerIndexer([self.LINKER], kmer_size=5)
        self.strict_linker_indexer = ArrayKmerIndexer([StereoBarcodeDetector.LINKER], kmer_size=7)

        self.barcode_indexer = None
        if barcodes:
            bit_barcodes = batch_str_to_2bit_chunked(iter(barcodes), seq_len=self.BC_LENGTH)
            self.barcode_indexer = Dict2BitKmerIndexer(bit_barcodes, kmer_size=14, seq_len=self.BC_LENGTH)
            logger.info("Indexed %d barcodes" % len(bit_barcodes))
        self.umi_set = None
        self.min_score = min_score

    def find_barcode_umi(self, read_id: str, sequence: str) -> SplittingBarcodeDetectionResult:
        """Find multiple barcodes in a concatenated read."""
        read_result = SplittingBarcodeDetectionResult(read_id)
        logger.debug("Looking in forward direction")
        r = self._find_barcode_umi_fwd(read_id, sequence)
        prev_start = 0
        while r.polyT != -1:
            r.set_strand("+")
            read_result.append(r)
            if r.tso5 != -1:
                current_start = r.tso5 + 15
            else:
                current_start = r.polyT + 100
            # always make a step
            current_start = max(prev_start + 150, current_start)
            prev_start = current_start
            if len(sequence) - current_start < 50:
                break

            logger.debug("Looking further from %d" % current_start)
            seq = sequence[current_start:]
            r = self._find_barcode_umi_fwd(read_id, seq)
            r.update_coordinates(current_start)

        logger.debug("Looking in reverse direction")
        rev_seq = reverese_complement(sequence)
        r = self._find_barcode_umi_fwd(read_id, rev_seq)
        prev_start = 0
        while r.polyT != -1:
            r.set_strand("-")
            read_result.append(r)
            if r.tso5 != -1:
                current_start = r.tso5 + 15
            else:
                current_start = r.polyT + 100
            # always make a step
            current_start = max(prev_start + 150, current_start)
            prev_start = current_start
            if len(rev_seq) - current_start < 50:
                break

            logger.debug("Looking further from %d" % current_start)
            seq = rev_seq[current_start:]
            r = self._find_barcode_umi_fwd(read_id, seq)
            r.update_coordinates(current_start)

        if read_result.empty():
            # add empty result anyway
            read_result.append(r)

        read_result.filter()
        logger.debug("Total barcodes detected %d" % len(read_result.detected_patterns))
        return read_result

    def find_barcode_umi_single(self, read_id: str, sequence: str) -> TSOBarcodeDetectionResult:
        """Find single barcode (non-splitting mode)."""
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

    def _find_barcode_umi_fwd(self, read_id: str, sequence: str) -> TSOBarcodeDetectionResult:
        polyt_start = find_polyt_start(sequence)
        logger.debug("PolyT found right away %d" % polyt_start)

        linker_start, linker_end = None, None
        tso5_start = None
        if polyt_start != -1:
            # use relaxed parameters if polyA is found
            logger.debug("Looking for linker in %d" % len(sequence[0:polyt_start + 1]))
            linker_occurrences = self.linker_indexer.get_occurrences(sequence[0:polyt_start + 1])
            linker_start, linker_end = detect_exact_positions(sequence, 0, polyt_start + 1,
                                                              self.linker_indexer.k, self.LINKER,
                                                              linker_occurrences, min_score=10,
                                                              start_delta=self.TERMINAL_MATCH_DELTA,
                                                              end_delta=self.TERMINAL_MATCH_DELTA)

            tso5_occurrences = self.tso5_indexer.get_occurrences(sequence[polyt_start + 1:])
            tso5_start, tso5_end = detect_first_exact_positions(sequence, polyt_start + 1, len(sequence),
                                                                self.tso5_indexer.k, self.TSO5,
                                                                tso5_occurrences, min_score=15,
                                                                start_delta=self.TERMINAL_MATCH_DELTA,
                                                                end_delta=self.TERMINAL_MATCH_DELTA)

        if linker_start is None:
            # if polyT was not found, or linker was not found to the left of polyT,
            # look for linker in the entire read
            linker_occurrences = self.strict_linker_indexer.get_occurrences(sequence)
            linker_start, linker_end = detect_first_exact_positions(sequence, 0, len(sequence),
                                                                    self.linker_indexer.k, StereoBarcodeDetector.LINKER,
                                                                    linker_occurrences, min_score=12,
                                                                    start_delta=self.STRICT_TERMINAL_MATCH_DELTA,
                                                                    end_delta=self.STRICT_TERMINAL_MATCH_DELTA)

        if linker_start is None:
            return TSOBarcodeDetectionResult(read_id, polyT=polyt_start)
        logger.debug("LINKER: %d-%d" % (linker_start, linker_end))

        if polyt_start == -1 or polyt_start < linker_start:
            # if polyT was not detected earlier, use relaxed parameters once the linker is found
            presumable_polyt_start = linker_end + self.UMI_LEN
            search_start = presumable_polyt_start - 4
            search_end = min(len(sequence), presumable_polyt_start + 10)
            polyt_start = find_polyt_start(sequence[search_start:search_end], window_size=5, polya_fraction=1.0)
            if polyt_start != -1:
                polyt_start += search_start
                logger.debug("PolyT found later %d" % polyt_start)
            else:
                logger.debug("PolyT was not found %d" % polyt_start)

            tso5_occurrences = self.tso5_indexer.get_occurrences(sequence[polyt_start + 1:])
            tso5_start, tso5_end = detect_first_exact_positions(sequence, polyt_start + 1, len(sequence),
                                                                self.tso5_indexer.k, self.TSO5,
                                                                tso5_occurrences, min_score=15,
                                                                start_delta=self.TERMINAL_MATCH_DELTA,
                                                                end_delta=self.TERMINAL_MATCH_DELTA)

        if tso5_start is not None:
            logger.debug("TSO found %d" % tso5_start)
            # check that no another linker is found inbetween polyA and TSO 5'
            linker_occurrences = self.strict_linker_indexer.get_occurrences(sequence[polyt_start + 1: tso5_start])
            new_linker_start, new_linker_end = detect_exact_positions(sequence, polyt_start + 1, tso5_start,
                                                                      self.linker_indexer.k, self.LINKER,
                                                                      linker_occurrences, min_score=12,
                                                                      start_delta=self.STRICT_TERMINAL_MATCH_DELTA,
                                                                      end_delta=self.STRICT_TERMINAL_MATCH_DELTA)

            if new_linker_start is not None and new_linker_start != -1 and new_linker_start - polyt_start > 100:
                # another linker found inbetween polyT and TSO
                logger.debug("Another linker was found before TSO: %d" % new_linker_start)
                tso5_start = new_linker_start - self.BC_LENGTH - len(self.main_primer) - len(self.TSO5)
                logger.debug("TSO updated %d" % tso5_start)
        else:
            tso5_start = -1

        primer_occurrences = self.pcr_primer_indexer.get_occurrences(sequence[:linker_start])
        primer_start, primer_end = detect_exact_positions(sequence, 0, linker_start,
                                                          self.pcr_primer_indexer.k, self.main_primer,
                                                          primer_occurrences, min_score=12,
                                                          end_delta=self.TERMINAL_MATCH_DELTA)
        if primer_start is not None:
            logger.debug("PRIMER: %d-%d" % (primer_start, primer_end))
        else:
            primer_end = linker_start - self.BC_LENGTH - 1

        if primer_end < 0:
            return TSOBarcodeDetectionResult(read_id, polyT=polyt_start, primer=-1,
                                             linker_start=linker_start, linker_end=linker_end, tso=tso5_start)

        barcode_start = primer_end + 1
        barcode_end = linker_start - 1
        bc_len = barcode_end - barcode_start
        if abs(bc_len - self.BC_LENGTH) > 10:
            return TSOBarcodeDetectionResult(read_id, polyT=polyt_start, primer=primer_end,
                                             linker_start=linker_start, linker_end=linker_end, tso=tso5_start)

        potential_barcode = sequence[barcode_start:barcode_end + 1]
        logger.debug("Barcode: %s" % (potential_barcode))
        if not self.barcode_indexer:
            return TSOBarcodeDetectionResult(read_id, potential_barcode, BC_score=0,
                                             polyT=polyt_start, primer=primer_end,
                                             linker_start=linker_start, linker_end=linker_end, tso=tso5_start)

        matching_barcodes = self.barcode_indexer.get_occurrences(potential_barcode, max_hits=10, min_kmers=2)
        barcode, bc_score, bc_start, bc_end = \
            find_candidate_with_max_score_ssw(matching_barcodes, potential_barcode,
                                              min_score=self.min_score, sufficient_score=self.BC_LENGTH - 1)

        if barcode is None:
            return TSOBarcodeDetectionResult(read_id, polyT=polyt_start, primer=primer_end,
                                             linker_start=linker_start, linker_end=linker_end, tso=tso5_start)
        logger.debug("Found: %s %d-%d" % (barcode, bc_start, bc_end))

        potential_umi_start = linker_end + 1
        potential_umi_end = polyt_start - 1
        if potential_umi_end - potential_umi_start <= 5:
            potential_umi_end = potential_umi_start + self.UMI_LEN - 1
        umi = None
        good_umi = False
        if potential_umi_start + 2 * self.UMI_LEN > potential_umi_end > potential_umi_start:
            umi = sequence[potential_umi_start:potential_umi_end + 1]
            logger.debug("Potential UMI: %s" % umi)
            good_umi = abs(len(umi) - self.UMI_LEN) <= self.UMI_LEN_DELTA

        if not umi:
            return TSOBarcodeDetectionResult(read_id, barcode, BC_score=bc_score,
                                             polyT=polyt_start, primer=primer_end,
                                             linker_start=linker_start, linker_end=linker_end, tso=tso5_start)
        return TSOBarcodeDetectionResult(read_id, barcode, umi, bc_score, good_umi,
                                         polyT=polyt_start, primer=primer_end,
                                         linker_start=linker_start, linker_end=linker_end, tso=tso5_start)

    @staticmethod
    def result_type():
        return SplittingBarcodeDetectionResult


class SharedMemoryStereoSplittingBarcodeDetector(StereoSplittingBarcodeDetector):
    """Stereo-seq splitting detector with shared memory for large barcode sets."""

    MIN_BARCODES_FOR_SHARED_MEM = 1000000

    def __init__(self, barcodes: List[str], min_score: int = 21):
        super().__init__([], min_score=min_score)
        # Convert to 2-bit encoding first (memory efficient)
        bit_barcodes = batch_str_to_2bit_chunked(iter(barcodes), seq_len=self.BC_LENGTH)
        self.barcode_count = len(bit_barcodes)
        self.bit_barcodes = None
        if self.barcode_count < self.MIN_BARCODES_FOR_SHARED_MEM:
            # For small/medium sets, use Dict2BitKmerIndexer
            self.bit_barcodes = bit_barcodes
            self.barcode_indexer = Dict2BitKmerIndexer(bit_barcodes, kmer_size=14, seq_len=self.BC_LENGTH)
        else:
            # For large sets, use shared memory indexer
            self.barcode_indexer = SharedMemoryArray2BitKmerIndexer(bit_barcodes, kmer_size=14,
                                                                    seq_len=self.BC_LENGTH)

        logger.info("Indexed %d barcodes" % self.barcode_count)

    def __getstate__(self):
        if self.barcode_count < self.MIN_BARCODES_FOR_SHARED_MEM:
            return (self.min_score,
                    self.barcode_count,
                    self.bit_barcodes)
        else:
            return (self.min_score,
                    self.barcode_count,
                    self.barcode_indexer.get_sharable_info())

    def __setstate__(self, state):
        self.min_score = state[0]
        super().__init__([], min_score=self.min_score)
        self.bit_barcodes = None
        self.barcode_count = state[1]
        if self.barcode_count < self.MIN_BARCODES_FOR_SHARED_MEM:
            self.bit_barcodes = state[2]
            self.barcode_indexer = Dict2BitKmerIndexer(self.bit_barcodes, kmer_size=14, seq_len=self.BC_LENGTH)
        else:
            self.barcode_indexer = SharedMemoryArray2BitKmerIndexer.from_sharable_info(state[2])


class SharedMemoryWrapper:
    """Generic wrapper to add shared memory support to any barcode detector."""

    def __init__(self, barcode_detector_class, barcodes: List[str], min_score: int = 21):
        """
        Initialize wrapper.

        Args:
            barcode_detector_class: Detector class to wrap
            barcodes: List of known barcodes
            min_score: Minimum alignment score
        """
        self.barcode_detector_class = barcode_detector_class
        self.min_score = min_score
        self.barcode_detector = self.barcode_detector_class([], self.min_score)
        # Convert to 2-bit and use shared memory indexer
        barcode_iter = iter(barcodes) if isinstance(barcodes, list) else barcodes
        bit_barcodes = batch_str_to_2bit_chunked(barcode_iter, seq_len=self.barcode_detector_class.BC_LENGTH)
        self.barcode_detector.barcode_indexer = SharedMemoryArray2BitKmerIndexer(bit_barcodes, kmer_size=14,
                                                                                 seq_len=self.barcode_detector_class.BC_LENGTH)

        logger.info("Indexed %d barcodes" % self.barcode_detector.barcode_indexer.total_sequences)

    def __getstate__(self):
        return (self.barcode_detector_class,
                self.min_score,
                self.barcode_detector.barcode_indexer.get_sharable_info())

    def __setstate__(self, state):
        self.barcode_detector_class = state[0]
        self.min_score = state[1]
        self.barcode_detector = self.barcode_detector_class([], self.min_score)
        self.barcode_detector.barcode_indexer = SharedMemoryArray2BitKmerIndexer.from_sharable_info(state[2])
