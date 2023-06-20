############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


NOSEQ = "*"


class BarcodeDetectionResult:
    def __init__(self, read_id, polyT=-1, primer=-1, linker_start=-1, linker_end=-1,
                 barcode=NOSEQ, BC_score=-1, UMI=NOSEQ, UMI_good=False, strand="."):
        self.read_id = read_id
        self.polyT = polyT
        self.primer = primer
        self.linker_start = linker_start
        self.linker_end = linker_end
        self.barcode = barcode
        self.UMI = UMI
        self.BC_score = BC_score
        self.UMI_good = UMI_good
        self.strand = strand

    def is_valid(self):
        return self.barcode != NOSEQ

    def more_informative_than(self, that):
        if self.polyT != that.polyT:
            return self.polyT > that.polyT
        if self.primer != that.primer:
            return self.primer > that.primer
        if self.linker_start != that.linker_start:
            return self.linker_start > that.linker_start
        return self.BC_score > that.BC_score

    def set_strand(self, strand):
        self.strand = strand

    def __str__(self):
        return "%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%s" % (self.read_id, self.polyT, self.primer,
                                                           self.linker_start, self.linker_end, self.strand,
                                                           self.barcode, self.UMI, self.BC_score, self.UMI_good)


class ReadStats:
    def __init__(self):
        self.read_count = 0
        self.polyT_count = 0
        self.primer_count = 0
        self.linker_count = 0
        self.bc_count = 0
        self.umi_count = 0

    def add_read(self, barcode_detection_result):
        self.read_count += 1
        if barcode_detection_result.polyT != -1:
            self.polyT_count += 1
        if barcode_detection_result.primer != -1:
            self.primer_count += 1
        if barcode_detection_result.linker_start != -1:
            self.linker_count += 1
        if barcode_detection_result.barcode != NOSEQ:
            self.bc_count += 1
        if barcode_detection_result.UMI_good:
            self.umi_count += 1

    def __str__(self):
        return "Total reads:\t%d\nPolyT found:\t%d\nPrimer found:\t%d\n" \
               "Linker found:\t%d\nBarcode detected:\t%d\nReliable UMI:\t%d\n" % \
            (self.read_count, self.polyT_count, self.primer_count, self.linker_count, self.bc_count, self.umi_count)