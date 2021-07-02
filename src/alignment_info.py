############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import pysam
from Bio import SeqIO

from src.polya_finder import *
from src.polya_verification import *

logger = logging.getLogger('IsoQuant')


class AlignmentInfo:
    def __init__(self, alignment):
        self.alignment = alignment
        # concat indels
        self.read_exons, self.read_blocks, self.cigar_blocks = get_read_blocks(alignment.reference_start,
                                                                               alignment.cigartuples)
        self.aligned_pairs = None
        self.aligned_pairs_index = None
        if not self.read_exons:
            return
        self.read_start = self.read_exons[0][0]
        self.read_end = self.read_exons[-1][1]
        self.polya_info = None
        self.exons_changed = False
        self.cage_hits = []
        self.combined_profile = None

    def construct_profiles(self, profile_constructor):
        self.combined_profile = profile_constructor.construct_profiles(self.read_exons, self.polya_info, self.cage_hits)

    def set_aligned_pairs(self):
        self.aligned_pairs = self.alignment.get_aligned_pairs()
        self.aligned_pairs_index = []
        exon_index = 0
        for i, read_pos, ref_pos in enumerate(self.aligned_pairs):
            if ref_pos == self.read_exons[exon_index][0]:
                self.aligned_pairs_index.append(i)
                exon_index += 1
        self.aligned_pairs_index.append(len(self.aligned_pairs))

    def get_seq(self, ref_start, ref_end, intron_region=None):
        if not self.aligned_pairs:
            self.set_aligned_pairs()

        if intron_region is None:
            selected_pairs = self.aligned_pairs
        else:
            selected_pairs = self.aligned_pairs[self.aligned_pairs_index[intron_region[0]]:
                                                self.aligned_pairs_index[intron_region[1]+2]]
        seq = ''
        last_ref_pos = 0
        for read_pos, ref_pos in selected_pairs:
            if ref_start <= last_ref_pos <= ref_end:
                if read_pos is not None:
                    seq += self.alignment.seq[read_pos]
            last_ref_pos = ref_pos or last_ref_pos
        return seq

    def add_cage_info(self, cage_finder):
        self.cage_hits = cage_finder.find_cage_peak(self.alignment)

    def add_polya_info(self, polya_finder, polya_fixer):
        self.polya_info = polya_finder.detect_polya(self.alignment)
        polya_exon_count, polyt_exon_count = polya_fixer.correct_read_info(self.read_exons, self.polya_info)

        if polya_exon_count > 0:
            self.polya_info.internal_polya_pos = shift_polya(self.read_exons, polya_exon_count,
                                                             self.polya_info.internal_polya_pos)
            self.polya_info.external_polya_pos = shift_polya(self.read_exons, polya_exon_count,
                                                             self.polya_info.external_polya_pos)
            self.read_exons = self.read_exons[:-polya_exon_count]
            self.read_blocks = self.read_blocks[:-polya_exon_count]
            self.cigar_blocks = self.cigar_blocks[:-polya_exon_count]
            logger.debug("Trimming polyA exons %d: %s" % (polya_exon_count, str(self.read_exons)))
            self.exons_changed = True

        if polyt_exon_count > 0:
            self.polya_info.internal_polyt_pos = shift_polyt(self.read_exons, polyt_exon_count,
                                                             self.polya_info.internal_polyt_pos)
            self.polya_info.external_polyt_pos = shift_polyt(self.read_exons, polyt_exon_count,
                                                             self.polya_info.external_polyt_pos)
            self.read_exons = self.read_exons[polyt_exon_count:]
            self.read_blocks = self.read_blocks[polyt_exon_count:]
            self.cigar_blocks = self.cigar_blocks[polyt_exon_count:]
            logger.debug("Trimming polyT exons %d: %s" % (polyt_exon_count, str(self.read_exons)))
            self.exons_changed = True

        if self.exons_changed:
            self.read_start = self.read_exons[0][0]
            self.read_end = self.read_exons[-1][1]
