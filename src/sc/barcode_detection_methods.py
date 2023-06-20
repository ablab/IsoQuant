############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging

from .kmer_index import KmerIndexer
from .sequence_common import (
    reverse_complement,
    detect_exact_positions,
    find_candidate_with_max_score_ssw,
    find_polyt_start,
)
from reports import *


logger = logging.getLogger('IsoQuant')


class TenXBarcodeDetector:
    def __init__(self):
        pass


class DoubleBarcodeDetector:
    LINKER = "TCTTCAGCGTTCCCGAGA"
    PCR_PRIMER = "TACACGACGCTCTTCCGATCT"
    LEFT_BC_LENGTH = 8
    RIGHT_BC_LENGTH = 6
    BC_LENGTH = LEFT_BC_LENGTH + RIGHT_BC_LENGTH
    MIN_BC_LEN = BC_LENGTH - 4
    UMI_LENGTH = 9
    NON_T_UMI_BASES = 2
    UMI_LEN_DELTA = 2
    SCORE_DIFF = 1

    TERMINAL_MATCH_DELTA = 1
    STRICT_TERMINAL_MATCH_DELTA = 0

    def __init__(self, joint_barcode_list, umi_list=None, min_score=14):
        self.pcr_primer_indexer = KmerIndexer([DoubleBarcodeDetector.PCR_PRIMER], kmer_size=6)
        self.linker_indexer = KmerIndexer([DoubleBarcodeDetector.LINKER], kmer_size=5)
        self.barcode_indexer = KmerIndexer(joint_barcode_list, kmer_size=5)
        self.umi_set = None
        if umi_list:
            self.umi_set =  set(umi_list)
            logger.debug("Loaded %d UMIs" % len(umi_list))
            self.umi_indexer = KmerIndexer(umi_list, kmer_size=5)
        self.min_score = min_score

    def find_barcode_umi(self, read_id, sequence):
        read_result = self._find_barcode_umi_fwd(read_id, sequence)
        if read_result.polyT != -1:
            read_result.set_strand("+")
        if read_result.is_valid():
            return read_result

        rev_seq = reverse_complement(sequence)
        read_rev_result = self._find_barcode_umi_fwd(read_id, rev_seq)
        if read_rev_result.polyT != -1:
            read_rev_result.set_strand("-")
        if read_rev_result.is_valid():
            return read_rev_result

        return read_result if read_result.more_informative_than(read_rev_result) else read_rev_result

    def _find_barcode_umi_fwd(self, read_id, sequence):
        # look for linker in the entire read
        linker_occurrences = self.linker_indexer.get_occurrences(sequence)
        linker_start, linker_end = detect_exact_positions(sequence, 0, len(sequence),
                                                          self.linker_indexer.k, self.LINKER,
                                                          linker_occurrences, min_score=14,
                                                          start_delta=self.TERMINAL_MATCH_DELTA,
                                                          end_delta=self.TERMINAL_MATCH_DELTA)

        if linker_start is None:
            return BarcodeDetectionResult(read_id)
        logger.debug("LINKER: %d-%d" % (linker_start, linker_end))
        primer_end = -1 # forget about primer

        # use relaxed parameters once the linker is found
        presumable_polyt_start = linker_end + self.RIGHT_BC_LENGTH + self.UMI_LENGTH
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
            return BarcodeDetectionResult(read_id, linker_start=linker_start, linker_end=linker_end)
        matching_barcodes = self.barcode_indexer.get_occurrences(potential_barcode)
        barcode, bc_score, bc_start, bc_end = \
            find_candidate_with_max_score_ssw(matching_barcodes, potential_barcode,
                                              min_score=len(potential_barcode) - 1, score_diff=self.SCORE_DIFF)

        if barcode is None:
            return BarcodeDetectionResult(read_id, polyt_start, -1, linker_start, linker_end)
        logger.debug("Found: %s %d-%d" % (barcode, bc_start, bc_end))
        # position of barcode end in the reference: end of potential barcode minus bases to the alignment end
        read_barcode_end = linker_end + self.RIGHT_BC_LENGTH + 1 - (len(potential_barcode) - bc_end - 1)

        potential_umi_start = read_barcode_end + 1
        potential_umi_end = polyt_start - 1
        if polyt_start != -1 or potential_umi_end - potential_umi_start <= 5:
            potential_umi_end = potential_umi_start + self.UMI_LENGTH - 1
        potential_umi = sequence[potential_umi_start:min(potential_umi_end + 1, len(sequence))]
        logger.debug("Potential UMI: %s" % potential_umi)

        umi = None
        good_umi = False
        if self.umi_set:
            matching_umis = self.umi_indexer.get_occurrences(potential_umi)
            umi, umi_score, umi_start, umi_end = \
                find_candidate_with_max_score_ssw(matching_umis, potential_umi, min_score=7)
            logger.debug("Found UMI %s %d-%d" % (umi, umi_start, umi_end))

        if not umi :
            umi = potential_umi
            if self.UMI_LENGTH - self.UMI_LEN_DELTA <= len(umi) <= self.UMI_LENGTH + self.UMI_LEN_DELTA and \
                    all(x != "T" for x in umi[-self.NON_T_UMI_BASES:]):
                good_umi = True

        if not umi:
            return BarcodeDetectionResult(read_id, polyt_start, primer_end, linker_start, linker_end, barcode, bc_score)
        return BarcodeDetectionResult(read_id, polyt_start, primer_end, linker_start, linker_end,
                                      barcode, bc_score, umi, good_umi)