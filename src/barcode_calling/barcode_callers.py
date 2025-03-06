###########################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict

from .kmer_indexer import KmerIndexer, ArrayKmerIndexer, Array2BitKmerIndexer
from .common import find_polyt_start, reverese_complement, find_candidate_with_max_score_ssw, detect_exact_positions, \
    detect_first_exact_positions, str_to_2bit

logger = logging.getLogger('IsoQuant')


def increase_if_valid(val, delta):
    if val and val != -1:
        return val + delta
    return val


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

    def update_coordinates(self, delta):
        pass

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

    def update_coordinates(self, delta):
        self.primer = increase_if_valid(self.primer, delta)
        self.linker_start = increase_if_valid(self.linker_start, delta)
        self.linker_end = increase_if_valid(self.linker_end, delta)
        self.polyT = increase_if_valid(self.polyT, delta)

    def more_informative_than(self, that):
        if self.BC_score != that.BC_score:
            return self.BC_score > that.BC_score
        if self.linker_start != that.linker_start:
            return self.linker_start > that.linker_start
        if self.primer != that.primer:
            return self.primer > that.primer
        return self.polyT > that.polyT

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


class StereoBarcodeDetectionResult(DoubleBarcodeDetectionResult):
    def __init__(self, read_id, barcode=BarcodeDetectionResult.NOSEQ, UMI=BarcodeDetectionResult.NOSEQ,
                 BC_score=-1, UMI_good=False, strand=".",
                 polyT=-1, primer=-1, linker_start=-1, linker_end=-1, tso=-1):
        DoubleBarcodeDetectionResult.__init__(self, read_id, barcode, UMI, BC_score, UMI_good, strand,
                                              polyT, primer, linker_start, linker_end)
        self.tso5 = tso

    def update_coordinates(self, delta):
        self.tso5 = increase_if_valid(self.tso5, delta)
        DoubleBarcodeDetectionResult.update_coordinates(self, delta)

    def __str__(self):
        return (DoubleBarcodeDetectionResult.__str__(self) +
                "\t%d" % self.tso5)

    @staticmethod
    def header():
        return DoubleBarcodeDetectionResult.header() + "\tTSO5"


class TenXBarcodeDetectionResult(BarcodeDetectionResult):
    def __init__(self, read_id, barcode=BarcodeDetectionResult.NOSEQ, UMI=BarcodeDetectionResult.NOSEQ,
                 BC_score=-1, UMI_good=False, strand=".",
                 polyT=-1, r1=-1):
        BarcodeDetectionResult.__init__(self, read_id, barcode, UMI, BC_score, UMI_good, strand)
        self.r1 = r1
        self.polyT = polyT

    def is_valid(self):
        return self.barcode != BarcodeDetectionResult.NOSEQ

    def update_coordinates(self, delta):
        self.r1 = increase_if_valid(self.r1, delta)
        self.polyT = increase_if_valid(self.polyT, delta)

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


class SplittingBarcodeDetectionResult:
    def __init__(self, read_id):
        self.read_id = read_id
        self.detected_patterns = []

    def append(self, barcode_detection_result):
        self.detected_patterns.append(barcode_detection_result)

    def filter(self):
        if not self.detected_patterns: return
        barcoded_results = []
        for r in self.detected_patterns:
            if r.barcode != BarcodeDetectionResult.NOSEQ:
                barcoded_results.append(r)

        if not barcoded_results:
            self.detected_patterns = [self.detected_patterns[0]]
        else:
            self.detected_patterns = barcoded_results

    @staticmethod
    def header():
        return StereoBarcodeDetectionResult.header()


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

    def add_custom_stats(self, stat_name, val: int):
        self.additional_attributes_counts[stat_name] += val

    def __str__(self):
        human_readable_str =  ("Total reads\t%d\nBarcode detected\t%d\nReliable UMI\t%d\n" %
                      (self.read_count, self.bc_count, self.umi_count))
        for a in self.additional_attributes_counts:
            human_readable_str += "%s\t%d\n" % (a, self.additional_attributes_counts[a])
        return human_readable_str


class StereoBarcodeDetector:
    LINKER = "TTGTCTTCCTAAGAC"
    TSO_PRIMER = "ACTGAGAGGCATGGCGACCTTATCAG"
    PC1_PRIMER = "CTTCCGATCTATGGCGACCTTATCAG"
    BC_LENGTH = 25
    UMI_LENGTH = 10
    NON_T_UMI_BASES = 0
    UMI_LEN_DELTA = 4
    TERMINAL_MATCH_DELTA = 3
    STRICT_TERMINAL_MATCH_DELTA = 1

    def __init__(self, barcodes, min_score=21, primer=1):
        if primer == 1:
            self.MAIN_PRIMER = StereoBarcodeDetector.TSO_PRIMER
        else:
            self.MAIN_PRIMER = StereoBarcodeDetector.PC1_PRIMER
        self.pcr_primer_indexer = ArrayKmerIndexer([self.MAIN_PRIMER], kmer_size=6)
        self.linker_indexer = ArrayKmerIndexer([StereoBarcodeDetector.LINKER], kmer_size=5)
        self.strict_linker_indexer = ArrayKmerIndexer([StereoBarcodeDetector.LINKER], kmer_size=6)
        bit_barcodes = map(str_to_2bit, barcodes)
        self.barcode_indexer = Array2BitKmerIndexer(bit_barcodes, kmer_size=14, seq_len=self.BC_LENGTH)
        logger.info("Indexed %d barcodes" % self.barcode_indexer.total_sequences)
        self.umi_set = None
        self.min_score = min_score

    def find_barcode_umi_multiple(self, read_id, sequence):
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
                                                              self.linker_indexer.k, StereoBarcodeDetector.LINKER,
                                                              linker_occurrences, min_score=12,
                                                              start_delta=self.TERMINAL_MATCH_DELTA,
                                                              end_delta=self.TERMINAL_MATCH_DELTA)

        if linker_start is None:
            # if polyT was not found, or linker was not found to the left of polyT, look for linker in the entire read
            linker_occurrences = self.strict_linker_indexer.get_occurrences(sequence)
            linker_start, linker_end = detect_exact_positions(sequence, 0, len(sequence),
                                                              self.linker_indexer.k, StereoBarcodeDetector.LINKER,
                                                              linker_occurrences, min_score=12,
                                                              start_delta=self.STRICT_TERMINAL_MATCH_DELTA,
                                                              end_delta=self.STRICT_TERMINAL_MATCH_DELTA)

        if linker_start is None:
            return DoubleBarcodeDetectionResult(read_id, polyT=polyt_start)
        logger.debug("LINKER: %d-%d" % (linker_start, linker_end))

        if polyt_start == -1:
            # if polyT was not detected earlier, use relaxed parameters once the linker is found
            presumable_polyt_start = linker_end + self.UMI_LENGTH
            search_start = presumable_polyt_start - 4
            search_end = min(len(sequence), presumable_polyt_start + 10)
            polyt_start = find_polyt_start(sequence[search_start:search_end], window_size=5, polya_fraction=1.0)
            if polyt_start != -1:
                polyt_start += search_start

        primer_occurrences = self.pcr_primer_indexer.get_occurrences(sequence[:linker_start])
        primer_start, primer_end = detect_exact_positions(sequence, 0, linker_start,
                                                          self.pcr_primer_indexer.k, self.MAIN_PRIMER,
                                                          primer_occurrences, min_score=12,
                                                          end_delta=self.TERMINAL_MATCH_DELTA)
        if primer_start is not None:
            logger.debug("PRIMER: %d-%d" % (primer_start, primer_end))
        else:
            primer_start = -1
            primer_end = linker_start - self.BC_LENGTH - 1

        if primer_end < 0:
            return DoubleBarcodeDetectionResult(read_id, polyT=polyt_start, primer=-1,
                                                linker_start=linker_start, linker_end=linker_end)

        barcode_start = primer_end + 1
        barcode_end = linker_start - 1
        bc_len = barcode_end - barcode_start
        if abs(bc_len - self.BC_LENGTH) > 10:
            return DoubleBarcodeDetectionResult(read_id, polyT=polyt_start, primer=primer_end,
                                                linker_start=linker_start, linker_end=linker_end)

        potential_barcode = sequence[barcode_start:barcode_end + 1]
        logger.debug("Barcode: %s" % (potential_barcode))
        matching_barcodes = self.barcode_indexer.get_occurrences(potential_barcode, max_hits=10, min_kmers=2)
        barcode, bc_score, bc_start, bc_end = \
            find_candidate_with_max_score_ssw(matching_barcodes, potential_barcode,
                                              min_score=self.min_score, sufficient_score=self.BC_LENGTH - 1)

        if barcode is None:
            return DoubleBarcodeDetectionResult(read_id, polyT=polyt_start, primer=primer_end,
                                                linker_start=linker_start, linker_end=linker_end)
        logger.debug("Found: %s %d-%d" % (barcode, bc_start, bc_end))

        potential_umi_start = linker_end + 1
        potential_umi_end = polyt_start - 1
        umi = None
        good_umi = False
        if potential_umi_start + 2 * self.UMI_LENGTH > potential_umi_end > potential_umi_start:
            umi = sequence[potential_umi_start:potential_umi_end + 1]
            logger.debug("Potential UMI: %s" % umi)
            good_umi = abs(len(umi) - self.UMI_LENGTH) <= self.UMI_LEN_DELTA

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


class StereoBarcodeDetectorTSO(StereoBarcodeDetector):
    def __init__(self, barcode_list, min_score=21):
        StereoBarcodeDetector.__init__(self, barcode_list, min_score, primer=1)


class StereoBarcodeDetectorPC(StereoBarcodeDetector):
    def __init__(self, barcode_list, min_score=21):
        StereoBarcodeDetector.__init__(self, barcode_list, min_score, primer=2)


class StereoSplttingBarcodeDetector:
    TSO5 = "CCCGCCTCTCAGTACGTCAGCAG"
    LINKER = "TTGTCTTCCTAAGAC"
    TSO_PRIMER = "ACTGAGAGGCATGGCGACCTTATCAG"
    PC1_PRIMER = "CTTCCGATCTATGGCGACCTTATCAG"
    BC_LENGTH = 25
    UMI_LENGTH = 10
    NON_T_UMI_BASES = 0
    UMI_LEN_DELTA = 3
    TERMINAL_MATCH_DELTA = 3
    STRICT_TERMINAL_MATCH_DELTA = 1

    def __init__(self, barcodes, min_score=21, primer=1):
        if primer == 1:
            self.MAIN_PRIMER = self.TSO_PRIMER
        else:
            self.MAIN_PRIMER = self.PC1_PRIMER
        self.tso5_indexer = ArrayKmerIndexer([self.TSO5], kmer_size=6)
        self.pcr_primer_indexer = ArrayKmerIndexer([self.MAIN_PRIMER], kmer_size=6)
        self.linker_indexer = ArrayKmerIndexer([self.LINKER], kmer_size=5)
        self.strict_linker_indexer = ArrayKmerIndexer([StereoBarcodeDetector.LINKER], kmer_size=6)
        #self.barcode_indexer = KmerIndexer(barcodes, kmer_size=14)
        #logger.info("Indexed %d barcodes" % len(self.barcode_indexer.seq_list))
        bit_barcodes = map(str_to_2bit, barcodes)
        self.barcode_indexer = Array2BitKmerIndexer(bit_barcodes, kmer_size=14, seq_len=self.BC_LENGTH)
        logger.info("Indexed %d barcodes" % self.barcode_indexer.total_sequences)
        self.umi_set = None
        self.min_score = min_score

    def find_barcode_umi(self, read_id, sequence):
        read_result = SplittingBarcodeDetectionResult(read_id)
        logger.debug("Looking in forward direction")
        r = self._find_barcode_umi_fwd(read_id, sequence)
        while r.polyT != -1:
            r.set_strand("+")
            read_result.append(r)
            if r.tso5 != -1:
                current_start = r.tso5 + 15
            else:
                current_start = r.polyT + 100
            if len(sequence) - current_start < 50:
                break

            logger.debug("Looking further from %d" % current_start)
            seq = sequence[current_start:]
            r = self._find_barcode_umi_fwd(read_id, seq)
            r.update_coordinates(current_start)

        logger.debug("Looking in reverse direction")
        rev_seq = reverese_complement(sequence)
        r = self._find_barcode_umi_fwd(read_id, rev_seq)
        while r.polyT != -1:
            r.set_strand("-")
            read_result.append(r)
            if r.tso5 != -1:
                current_start = r.tso5 + 15
            else:
                current_start = r.polyT + 100
            if len(rev_seq) - current_start < 50:
                break

            logger.debug("Looking further from %d" % current_start)
            seq = rev_seq[current_start:]
            r = self._find_barcode_umi_fwd(read_id, seq)
            r.update_coordinates(current_start)

        read_result.filter()
        logger.debug("Total barcodes detected %d" % len(read_result.detected_patterns))
        return read_result

    def find_barcode_umi_single(self, read_id, sequence):
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
        logger.debug("PolyT found right away %d" % polyt_start)

        linker_start, linker_end = None, None
        tso5_start = None
        if polyt_start != -1:
            # use relaxed parameters is polyA is found
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
            # if polyT was not found, or linker was not found to the left of polyT, look for linker in the entire read
            linker_occurrences = self.strict_linker_indexer.get_occurrences(sequence)
            linker_start, linker_end = detect_first_exact_positions(sequence, 0, len(sequence),
                                                                    self.linker_indexer.k, StereoBarcodeDetector.LINKER,
                                                                    linker_occurrences, min_score=12,
                                                                    start_delta=self.STRICT_TERMINAL_MATCH_DELTA,
                                                                    end_delta=self.STRICT_TERMINAL_MATCH_DELTA)

        if linker_start is None:
            return StereoBarcodeDetectionResult(read_id, polyT=polyt_start)
        logger.debug("LINKER: %d-%d" % (linker_start, linker_end))

        if polyt_start == -1 or polyt_start < linker_start:
            # if polyT was not detected earlier, use relaxed parameters once the linker is found
            presumable_polyt_start = linker_end + self.UMI_LENGTH
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

        if tso5_start:
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
                tso5_start = new_linker_start - self.BC_LENGTH - len(self.MAIN_PRIMER) - len(self.TSO5)
                logger.debug("TSO updated %d" % tso5_start)
        else:
            tso5_start = -1

        primer_occurrences = self.pcr_primer_indexer.get_occurrences(sequence[:linker_start])
        primer_start, primer_end = detect_exact_positions(sequence, 0, linker_start,
                                                          self.pcr_primer_indexer.k, self.MAIN_PRIMER,
                                                          primer_occurrences, min_score=12,
                                                          end_delta=self.TERMINAL_MATCH_DELTA)
        if primer_start is not None:
            logger.debug("PRIMER: %d-%d" % (primer_start, primer_end))
        else:
            primer_end = linker_start - self.BC_LENGTH - 1

        if primer_end < 0:
            return StereoBarcodeDetectionResult(read_id, polyT=polyt_start, primer=-1,
                                                linker_start=linker_start, linker_end=linker_end, tso=tso5_start)

        barcode_start = primer_end + 1
        barcode_end = linker_start - 1
        bc_len = barcode_end - barcode_start
        if abs(bc_len - self.BC_LENGTH) > 10:
            return StereoBarcodeDetectionResult(read_id, polyT=polyt_start, primer=primer_end,
                                                linker_start=linker_start, linker_end=linker_end, tso=tso5_start)

        potential_barcode = sequence[barcode_start:barcode_end + 1]
        logger.debug("Barcode: %s" % (potential_barcode))
        matching_barcodes = self.barcode_indexer.get_occurrences(potential_barcode, max_hits=10, min_kmers=2)
        barcode, bc_score, bc_start, bc_end = \
            find_candidate_with_max_score_ssw(matching_barcodes, potential_barcode,
                                              min_score=self.min_score, sufficient_score=self.BC_LENGTH - 1)

        if barcode is None:
            return StereoBarcodeDetectionResult(read_id, polyT=polyt_start, primer=primer_end,
                                                linker_start=linker_start, linker_end=linker_end, tso=tso5_start)
        logger.debug("Found: %s %d-%d" % (barcode, bc_start, bc_end))

        potential_umi_start = linker_end + 1
        potential_umi_end = polyt_start - 1
        umi = None
        good_umi = False
        if potential_umi_start + 2 * self.UMI_LENGTH > potential_umi_end > potential_umi_start:
            umi = sequence[potential_umi_start:potential_umi_end + 1]
            logger.debug("Potential UMI: %s" % umi)
            good_umi = abs(len(umi) - self.UMI_LENGTH) <= self.UMI_LEN_DELTA

        if not umi:
            return StereoBarcodeDetectionResult(read_id, barcode, BC_score=bc_score,
                                            polyT=polyt_start, primer=primer_end,
                                            linker_start=linker_start, linker_end=linker_end, tso=tso5_start)
        return StereoBarcodeDetectionResult(read_id, barcode, umi, bc_score, good_umi,
                                            polyT=polyt_start, primer=primer_end,
                                            linker_start=linker_start, linker_end=linker_end, tso=tso5_start)

    @staticmethod
    def result_type():
        return SplittingBarcodeDetectionResult


class StereoSplitBarcodeDetectorTSO(StereoSplttingBarcodeDetector):
    def __init__(self, barcode_list, min_score=21):
        StereoSplttingBarcodeDetector.__init__(self, barcode_list, min_score, primer=1)


class StereoSplitBarcodeDetectorPC(StereoSplttingBarcodeDetector):
    def __init__(self, barcode_list, min_score=21):
        StereoSplttingBarcodeDetector.__init__(self, barcode_list, min_score, primer=2)


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
