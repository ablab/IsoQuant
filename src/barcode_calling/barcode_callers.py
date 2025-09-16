###########################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict

from .kmer_indexer import KmerIndexer, ArrayKmerIndexer
from .common import find_polyt_start, reverese_complement, find_candidate_with_max_score_ssw, detect_exact_positions

logger = logging.getLogger('IsoQuant')


class BarcodeDetectionResult:
    NOSEQ = "*"

    def __init__(self, read_id, barcode=NOSEQ, UMI=NOSEQ, BC_score=-1, UMI_good=False, strand=".", additional_info=None):
        self.read_id = read_id
        self.barcode = barcode
        self.UMI = UMI
        self.BC_score = BC_score
        self.UMI_good = UMI_good
        self.strand = strand

    def is_valid(self):
        return self.barcode != BarcodeDetectionResult.NOSEQ

    def more_informative_than(self, that):
        raise NotImplemented()

    def get_additional_attributes(self):
        raise NotImplemented()

    def set_strand(self, strand):
        self.strand = strand

    def __str__(self):
        return "%s\t%s\t%s\t%d\t%s\t%s" % (self.read_id, self.barcode, self.UMI,
                                           self.BC_score, self.UMI_good, self.strand)

    @staticmethod
    def header():
        return "#read_id\tbarcode\tUMI\tBC_score\tvalid_UMI\tstrand"


class DoubleBarcodeDetectionResult(BarcodeDetectionResult):
    def __init__(self, read_id, barcode=BarcodeDetectionResult.NOSEQ, UMI=BarcodeDetectionResult.NOSEQ,
                 BC_score=-1, UMI_good=False, strand=".",
                 polyT=-1, primer=-1, linker_start=-1, linker_end=-1):
        BarcodeDetectionResult.__init__(self, read_id, barcode, UMI, BC_score, UMI_good, strand)
        self.primer = primer
        self.linker_start = linker_start
        self.linker_end = linker_end
        self.polyT = polyT

    def is_valid(self):
        return self.barcode != BarcodeDetectionResult.NOSEQ

    def more_informative_than(self, that):
        if self.polyT != that.polyT:
            return self.polyT > that.polyT
        if self.primer != that.primer:
            return self.primer > that.primer
        if self.linker_start != that.linker_start:
            return self.linker_start > that.linker_start
        return self.BC_score > that.BC_score

    def get_additional_attributes(self):
        attr = []
        if self.polyT != -1:
            attr.append("PolyT detected")
        if self.primer != -1:
            attr.append("Primer detected")
        if self.linker_start != -1:
            attr.append("Linker detected")
        return attr

    def set_strand(self, strand):
        self.strand = strand

    def __str__(self):
        return (BarcodeDetectionResult.__str__(self) +
                "\t%d\t%d\t%d\t%d" % (self.polyT, self.primer, self.linker_start, self.linker_end))

    @staticmethod
    def header():
        return BarcodeDetectionResult.header() + "\tpolyT_start\tprimer_end\tlinker_start\tlinker_end"


class TenXBarcodeDetectionResult(BarcodeDetectionResult):
    def __init__(self, read_id, barcode=BarcodeDetectionResult.NOSEQ, UMI=BarcodeDetectionResult.NOSEQ,
                 BC_score=-1, UMI_good=False, strand=".",
                 polyT=-1, r1=-1):
        BarcodeDetectionResult.__init__(self, read_id, barcode, UMI, BC_score, UMI_good, strand)
        self.r1 = r1
        self.polyT = polyT

    def is_valid(self):
        return self.barcode != BarcodeDetectionResult.NOSEQ

    def more_informative_than(self, that):
        if self.polyT != that.polyT:
            return self.polyT > that.polyT
        if self.r1 != that.r1:
            return self.r1 > that.r1
        return self.BC_score > that.BC_score

    def get_additional_attributes(self):
        attr = []
        if self.polyT != -1:
            attr.append("PolyT detected")
        if self.r1 != -1:
            attr.append("R1 detected")
        return attr

    def set_strand(self, strand):
        self.strand = strand

    def __str__(self):
        return (BarcodeDetectionResult.__str__(self) +
                "\t%d\t%d" % (self.polyT, self.r1))

    @staticmethod
    def header():
        return BarcodeDetectionResult.header() + "\tpolyT_start\tR1_end"


class ReadStats:
    def __init__(self):
        self.read_count = 0
        self.bc_count = 0
        self.umi_count = 0
        self.additional_attributes_counts = defaultdict(int)

    def add_read(self, barcode_detection_result):
        self.read_count += 1
        for a in barcode_detection_result.get_additional_attributes():
            self.additional_attributes_counts[a] += 1
        if barcode_detection_result.barcode != BarcodeDetectionResult.NOSEQ:
            self.bc_count += 1
        if barcode_detection_result.UMI_good:
            self.umi_count += 1

    def __str__(self):
        human_readable_str =  ("Total reads:\t%d\nBarcode detected:\t%d\nReliable UMI:\t%d\n" %
                      (self.read_count, self.bc_count, self.umi_count))
        for a in self.additional_attributes_counts:
            human_readable_str += "%s\t%d\n" % (a, self.additional_attributes_counts[a])
        return human_readable_str


class DoubleBarcodeDetector:
    LINKER = "TCTTCAGCGTTCCCGAGA"
    PCR_PRIMER = "TACACGACGCTCTTCCGATCT"
    LEFT_BC_LENGTH = 8
    RIGHT_BC_LENGTH = 6
    BC_LENGTH = LEFT_BC_LENGTH + RIGHT_BC_LENGTH
    UMI_LENGTH = 9
    NON_T_UMI_BASES = 2
    UMI_LEN_DELTA = 2
    TERMINAL_MATCH_DELTA = 2
    STRICT_TERMINAL_MATCH_DELTA = 1

    def __init__(self, joint_barcode_list, umi_list=None, min_score=13):
        self.pcr_primer_indexer = ArrayKmerIndexer([DoubleBarcodeDetector.PCR_PRIMER], kmer_size=6)
        self.linker_indexer = ArrayKmerIndexer([DoubleBarcodeDetector.LINKER], kmer_size=5)
        self.barcode_indexer = ArrayKmerIndexer(joint_barcode_list, kmer_size=6)
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

        rev_seq = reverese_complement(sequence)
        read_rev_result = self._find_barcode_umi_fwd(read_id, rev_seq)
        if read_rev_result.polyT != -1:
            read_rev_result.set_strand("-")
        if read_rev_result.is_valid():
            return read_rev_result

        return read_result if read_result.more_informative_than(read_rev_result) else read_rev_result

    def _find_barcode_umi_fwd(self, read_id, sequence):
        polyt_start = find_polyt_start(sequence)

        linker_start, linker_end = None, None
        if polyt_start != -1:
            # use relaxed parameters is polyA is found
            linker_occurrences = self.linker_indexer.get_occurrences(sequence[0:polyt_start + 1])
            linker_start, linker_end = detect_exact_positions(sequence, 0, polyt_start + 1,
                                                              self.linker_indexer.k, self.LINKER,
                                                              linker_occurrences, min_score=11,
                                                              start_delta=self.TERMINAL_MATCH_DELTA,
                                                              end_delta=self.TERMINAL_MATCH_DELTA)

        if linker_start is None:
            # if polyT was not found, or linker was not found to the left of polyT, look for linker in the entire read
            linker_occurrences = self.linker_indexer.get_occurrences(sequence)
            linker_start, linker_end = detect_exact_positions(sequence, 0, len(sequence),
                                                              self.linker_indexer.k, self.LINKER,
                                                              linker_occurrences, min_score=14,
                                                              start_delta=self.STRICT_TERMINAL_MATCH_DELTA,
                                                              end_delta=self.STRICT_TERMINAL_MATCH_DELTA)

        if linker_start is None:
            return DoubleBarcodeDetectionResult(read_id, polyT=polyt_start)
        logger.debug("LINKER: %d-%d" % (linker_start, linker_end))

        if polyt_start == -1 or polyt_start - linker_end > self.RIGHT_BC_LENGTH + self.UMI_LENGTH + 10:
            # if polyT was not detected earlier, use relaxed parameters once the linker is found
            presumable_polyt_start = linker_end + self.RIGHT_BC_LENGTH + self.UMI_LENGTH
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
            return DoubleBarcodeDetectionResult(read_id, polyT=polyt_start, primer=primer_end,
                                                linker_start=linker_start, linker_end=linker_end)
        logger.debug("Found: %s %d-%d" % (barcode, bc_start, bc_end))
        # position of barcode end in the reference: end of potential barcode minus bases to the alignment end
        read_barcode_end = linker_end + self.RIGHT_BC_LENGTH + 1 - (len(potential_barcode) - bc_end - 1)

        potential_umi_start = read_barcode_end + 1
        potential_umi_end = polyt_start - 1
        if potential_umi_end - potential_umi_start <= 5:
            potential_umi_end = potential_umi_start + self.UMI_LENGTH - 1
        potential_umi = sequence[potential_umi_start:potential_umi_end + 1]
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
            return DoubleBarcodeDetectionResult(read_id, barcode, BC_score=bc_score,
                                            polyT=polyt_start, primer=primer_end,
                                            linker_start=linker_start, linker_end=linker_end)
        return DoubleBarcodeDetectionResult(read_id, barcode, umi, bc_score, good_umi,
                                            polyT=polyt_start, primer=primer_end,
                                            linker_start=linker_start, linker_end=linker_end)

    @staticmethod
    def result_type():
        return DoubleBarcodeDetectionResult


class BruteForceDoubleBarcodeDetector:
    LINKER = "TCTTCAGCGTTCCCGAGA"
    LEFT_BC_LENGTH = 2
    RIGHT_BC_LENGTH = 12
    BC_LENGTH = LEFT_BC_LENGTH + RIGHT_BC_LENGTH

    def __init__(self, joint_barcode_list):
        self.barcode_set = set(joint_barcode_list)

    def find_barcode_umi(self, read_id, sequence):
        linker_found, barcode = self._find_barcode_umi_fwd(read_id, sequence)
        if linker_found:
            return linker_found, barcode

        rev_seq = reverese_complement(sequence)
        return self._find_barcode_umi_fwd(read_id, rev_seq)

    def _find_barcode_umi_fwd(self, read_id, sequence):
        pos = sequence.find(self.LINKER)
        if pos == -1:
            return DoubleBarcodeDetectionResult(read_id)

        bc_start = max(0, pos - self.LEFT_BC_LENGTH)
        barcode = sequence[bc_start:pos]
        linker_end = pos + len(self.LINKER)
        bc_end = min(len(sequence), linker_end + 1 + self.RIGHT_BC_LENGTH)
        barcode += sequence[linker_end + 1:bc_end]
        if len(barcode) != self.BC_LENGTH or barcode not in self.barcode_set:
            return DoubleBarcodeDetectionResult(read_id, linker_start=pos)
        return DoubleBarcodeDetectionResult(read_id, barcode, BC_score=len(barcode), linker_start=pos)

    @staticmethod
    def result_type():
        return DoubleBarcodeDetectionResult


class IlluminaDoubleBarcodeDetector:
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

        rev_seq = reverese_complement(sequence)
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
            return DoubleBarcodeDetectionResult(read_id)
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
            return DoubleBarcodeDetectionResult(read_id, linker_start=linker_start, linker_end=linker_end)
        matching_barcodes = self.barcode_indexer.get_occurrences(potential_barcode)
        barcode, bc_score, bc_start, bc_end = \
            find_candidate_with_max_score_ssw(matching_barcodes, potential_barcode,
                                              min_score=len(potential_barcode) - 1, score_diff=self.SCORE_DIFF)

        if barcode is None:
            return DoubleBarcodeDetectionResult(read_id, polyT=polyt_start, primer=-1,
                                                linker_start=linker_start, linker_end=linker_end)
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
            return DoubleBarcodeDetectionResult(read_id, barcode, BC_score=bc_score,
                                                polyT=polyt_start, primer=primer_end,
                                                linker_start=linker_start, linker_end=linker_end)
        return DoubleBarcodeDetectionResult(read_id, barcode, umi, bc_score, good_umi,
                                            polyT=polyt_start, primer=primer_end,
                                            linker_start=linker_start, linker_end=linker_end)

    @staticmethod
    def result_type():
        return DoubleBarcodeDetectionResult


class TenXBarcodeDetector:
    TSO = "CCCATGTACTCTGCGTTGATACCACTGCTT"
    # R1 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"  #
    R1 = "CTACACGACGCTCTTCCGATCT" # 10x 3'
    BARCODE_LEN_10X = 16
    UMI_LEN_10X = 12

    UMI_LEN_DELTA = 2
    TERMINAL_MATCH_DELTA = 2
    STRICT_TERMINAL_MATCH_DELTA = 1

    def __init__(self, barcode_list, umi_list=None):
        self.r1_indexer = KmerIndexer([TenXBarcodeDetector.R1], kmer_size=7)
        self.barcode_indexer = KmerIndexer(barcode_list, kmer_size=6)
        self.umi_set = None
        if umi_list:
            self.umi_set =  set(umi_list)
            logger.debug("Loaded %d UMIs" % len(umi_list))
            self.umi_indexer = KmerIndexer(umi_list, kmer_size=5)
        self.min_score = 14
        if len(barcode_list) > 100000:
            self.min_score = 16
        logger.debug("Min score set to %d" % self.min_score)

    def find_barcode_umi(self, read_id, sequence):
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

    def _find_barcode_umi_fwd(self, read_id, sequence):
        polyt_start = find_polyt_start(sequence)

        r1_start, r1_end = None, None
        if polyt_start != -1:
            # use relaxed parameters is polyA is found
            r1_occurrences = self.r1_indexer.get_occurrences(sequence[0:polyt_start + 1])
            r1_start, r1_end = detect_exact_positions(sequence, 0, polyt_start + 1,
                                                      self.r1_indexer.k, self.R1,
                                                      r1_occurrences, min_score=11,
                                                      end_delta=self.TERMINAL_MATCH_DELTA)

        if r1_start is None:
            # if polyT was not found, or linker was not found to the left of polyT, look for linker in the entire read
            r1_occurrences = self.r1_indexer.get_occurrences(sequence)
            r1_start, r1_end = detect_exact_positions(sequence, 0, len(sequence),
                                                      self.r1_indexer.k, self.R1,
                                                      r1_occurrences, min_score=18,
                                                      start_delta=self.STRICT_TERMINAL_MATCH_DELTA,
                                                      end_delta=self.STRICT_TERMINAL_MATCH_DELTA)

        if r1_start is None:
            return TenXBarcodeDetectionResult(read_id, polyT=polyt_start)
        logger.debug("LINKER: %d-%d" % (r1_start, r1_end))

        if polyt_start == -1 or polyt_start - r1_end > self.BARCODE_LEN_10X + self.UMI_LEN_10X + 10:
            # if polyT was not detected earlier, use relaxed parameters once the linker is found
            presumable_polyt_start = r1_end + self.BARCODE_LEN_10X + self.UMI_LEN_10X
            search_start = presumable_polyt_start - 4
            search_end = min(len(sequence), presumable_polyt_start + 10)
            polyt_start = find_polyt_start(sequence[search_start:search_end], window_size=5, polya_fraction=1.0)
            if polyt_start != -1:
                polyt_start += search_start

        barcode_start = r1_end + 1
        barcode_end = r1_end + self.BARCODE_LEN_10X + 1
        potential_barcode = sequence[barcode_start:barcode_end + 1]
        logger.debug("Barcode: %s" % (potential_barcode))
        matching_barcodes = self.barcode_indexer.get_occurrences(potential_barcode)
        barcode, bc_score, bc_start, bc_end = \
            find_candidate_with_max_score_ssw(matching_barcodes, potential_barcode, min_score=self.min_score)

        if barcode is None:
            return TenXBarcodeDetectionResult(read_id, polyT=polyt_start, r1=r1_end)
        logger.debug("Found: %s %d-%d" % (barcode, bc_start, bc_end))
        # position of barcode end in the reference: end of potential barcode minus bases to the alignment end
        read_barcode_end = r1_end + self.BARCODE_LEN_10X + 1 - (len(potential_barcode) - bc_end - 1)
        potential_umi_start = read_barcode_end + 1
        potential_umi_end = polyt_start - 1
        if potential_umi_end - potential_umi_start <= 5:
            potential_umi_end = potential_umi_start + self.UMI_LEN_10X - 1
        potential_umi = sequence[potential_umi_start:potential_umi_end + 1]
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
            if self.UMI_LEN_10X - self.UMI_LEN_DELTA <= len(umi) <= self.UMI_LEN_10X + self.UMI_LEN_DELTA:
                good_umi = True

        if not umi:
            return TenXBarcodeDetectionResult(read_id, barcode, BC_score=bc_score, polyT=polyt_start, r1=r1_end)
        return TenXBarcodeDetectionResult(read_id, barcode, umi, bc_score, good_umi, polyT=polyt_start, r1=r1_end)

    def find_barcode_umi_no_polya(self, read_id, sequence):
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
