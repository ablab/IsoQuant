############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import namedtuple

import pysam
from Bio import pairwise2

from src.isoform_assignment import *
from src.gene_info import *
from src.long_read_profiles import *

logger = logging.getLogger('IsoQuant')
pairwise2.MAX_ALIGNMENTS = 1


class AlignmentRefiner:
    def __init__(self, gene_info, params):
        self.gene_info = gene_info
        self.params = params
        self.delta = self.params.delta
        self.scores = (2, -1, -1, -.2)

    def classify_skipped_exon(self, isoform_junctions, isoform_cregion, read_junctions, read_cregion,
                              total_intron_len_diff, read_introns_known, total_exon_len, alignment):

        assert len(set(read_cregion)) == 1

        read_cpos = read_cregion[0]
        read_left, read_right = read_junctions[read_cpos]
        l, r = isoform_cregion
        iso_left, iso_right = isoform_junctions[l][0], isoform_junctions[r][1]
        exon_left, exon_right = isoform_junctions[l][1], isoform_junctions[r][0]
        if total_intron_len_diff < 2 * self.delta and total_exon_len <= self.max_missed_exon_len:
            if alignment is None or self.reference is None:
                return MatchEventSubtype.exon_misallignment

            if iso_left - self.delta <= read_left and iso_right + self.delta >= read_right:
                seq = self.get_read_sequence(alignment, (iso_left, read_left)) \
                      + self.get_read_sequence(alignment, (read_right, iso_right))
                ref_seq = self.reference.fetch(alignment.reference_name, exon_left, exon_right)
                if self.sequences_match(seq, ref_seq):
                    return MatchEventSubtype.exon_misallignment
                else:
                    return MatchEventSubtype.exon_skipping_novel_intron

        if read_introns_known:
            return MatchEventSubtype.exon_skipping_known_intron
        return MatchEventSubtype.exon_skipping_novel_intron

    def classify_intron_misalignment(self, read_junctions, isoform_junctions, read_cpos, isoform_cpos, alignment=None):

        read_left, read_right = read_junctions[read_cpos]
        iso_left, iso_right = isoform_junctions[isoform_cpos]

        if abs(read_left - iso_left) + abs(read_right - iso_right) <= 2 * self.delta:
            return MatchEventSubtype.intron_shift

        if alignment is None or self.reference is None:
            if abs(read_left - iso_left) <= self.max_intron_shift:
                return MatchEventSubtype.intron_shift
            return MatchEventSubtype.intron_alternation_novel

        read_region, ref_region = self.get_aligned_regions_intron(read_left, read_right, iso_left, iso_right)
        seq = self.get_read_sequence(alignment, read_region)
        ref_seq = self.reference.fetch(alignment.reference_name, *ref_region)

        if self.sequences_match(seq, ref_seq):
            return MatchEventSubtype.intron_shift
        return MatchEventSubtype.intron_alternation_novel

    def sequences_match(self, seq, ref_seq):
        score, size = 0, 1
        for a in pairwise2.align.globalms(seq, ref_seq, *self.scores):
            score, size = a[2], a[4]
        if score > 0.7 * size:
            return True
        return False

    @staticmethod
    def get_aligned_regions_intron(read_left, read_right, iso_left, iso_right):
        if iso_left <= read_left:
            return (iso_left, read_left), (iso_right, read_right)
        return (read_right, iso_right), (read_left, iso_left)

