
############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import pytest

from src.alignment_info import AlignmentInfo
from src.polya_finder import PolyAFinder
from src.polya_verification import PolyAFixer

SEQ1 = 'CTCAAGACCAAGAAGGACGACATGACCATGGCTTAAAAGAGTCTGCTCCCCACAGCCCCCTGCGAT' \
       'GGATGGACGGAGGAACCAGGGTCGGACGACCTCCGATGCTAAGAGCACTCCAACTGCTGCAAACCG' \
       'AAGAAGCAGGCATCGGACDCCTTTTTTTTTTAGGCGCCAAAAAAAAAAAACCAAAAAAAAAAAAA'
TUPLES1 = [(0, 167), (4, 30)]
TUPLES2 = [(0, 167), (4, 30), (5, 5), (4, 16), (0, 30), (5, 70), (0, 201)]
PAIRS = [(0, 75897), (1, 76063), (2, 75995)]


class Alignment:
    def __init__(self, query_name, cigartuples, seq, reference_start, reference_end, reference_name, reference_id,
                 aligned_pairs):
        self.query_name = query_name
        self.cigartuples = cigartuples
        self.seq = seq
        self.reference_start = reference_start
        self.reference_end = reference_end
        self.reference_name = reference_name
        self.reference_id = reference_id
        self.aligned_pairs = aligned_pairs

    def get_aligned_pairs(self):
        return self.aligned_pairs


class TestAlignmentInfo:

    @pytest.mark.parametrize("read_exons, read_blocks, cigar_blocks",
                             [
                                 ([(75897, 76063), (76064, 76294)], [(0, 166), (213, 443)], [(0, 0), (4, 6)])
                             ])
    def test_basic(self, read_exons, read_blocks, cigar_blocks):
        alignment = Alignment('query', TUPLES2, SEQ1, 75896, 113577, 'reference', 789, PAIRS)
        alignment_info = AlignmentInfo(alignment)

        assert alignment_info.alignment == alignment
        assert alignment_info.read_start == read_exons[0][0]
        assert alignment_info.read_end == read_exons[-1][-1]
        assert alignment_info.read_exons == read_exons
        assert alignment_info.read_blocks == read_blocks
        assert alignment_info.cigar_blocks == cigar_blocks

        assert alignment_info.aligned_pairs is None
        assert alignment_info.aligned_pairs_start_index is None
        assert alignment_info.aligned_pairs_end_index is None
        assert alignment_info.polya_info.external_polyt_pos == -1
        assert alignment_info.polya_info.external_polya_pos == -1
        assert alignment_info.polya_info.internal_polya_pos == -1
        assert alignment_info.polya_info.internal_polyt_pos == -1
        assert alignment_info.combined_profile is None
        assert alignment_info.exons_changed is False
        assert len(alignment_info.cage_hits) == 0

    def test_alignment_pairs(self):
        alignment = Alignment('query', TUPLES1, SEQ1, 75896, 113577, 'reference', 789, PAIRS)
        alignment_info = AlignmentInfo(alignment)
        alignment_info.set_aligned_pairs()

        assert len(PAIRS) in alignment_info.aligned_pairs_start_index
        assert len(alignment_info.aligned_pairs_start_index) == 3

    @pytest.mark.parametrize("ref_start, ref_end, intron_index, left_site, chr_record, expected",
                             [
                                 (75999, 81000, None, None, None, (0, 0)),
                                 (75999, 81000, 0, False, None, (0, 0)),
                                 (75999, 81000, 0, True, None, (0, 0))
                             ])
    def test_get_error_count(self, ref_start, ref_end, intron_index, left_site, chr_record, expected):
        alignment = Alignment('query', TUPLES2, SEQ1, 75896, 113577, 'reference', 789, PAIRS)
        alignment_info = AlignmentInfo(alignment)

        indel_count, mismatch_count = alignment_info.get_error_count(ref_start, ref_end, intron_index, left_site,
                                                                     chr_record)
        assert (indel_count, mismatch_count) == expected

    @pytest.mark.parametrize("read_start, read_end",
                             [
                                 (51055, 51221)
                             ])
    def test_add_polya_info(self, read_start, read_end):
        alignment = Alignment('query', TUPLES2, SEQ1, 51054, 35000, 'reference', 789, PAIRS)
        alignment_info = AlignmentInfo(alignment)
        polya_finder = PolyAFinder()
        polya_fixer = PolyAFixer(1)

        alignment_info.add_polya_info(polya_finder, polya_fixer)

        assert alignment_info.exons_changed is True
        assert alignment_info.polya_info.internal_polya_pos == read_end
        assert alignment_info.polya_info.internal_polyt_pos == -1
        assert alignment_info.polya_info.external_polya_pos == -1
        assert alignment_info.polya_info.external_polyt_pos == -1
        assert alignment_info.read_exons == [(read_start, read_end)]
        assert alignment_info.read_blocks == [(0, (read_end - read_start))]
        assert alignment_info.cigar_blocks == [(0, 0)]
        assert alignment_info.read_start == read_start
        assert alignment_info.read_end == read_end

    def test_get_seq(self):
        alignment = Alignment('query', TUPLES1, SEQ1, 75896, 113577, 'reference', 789, PAIRS)
        alignment_info = AlignmentInfo(alignment)

        with pytest.raises(AssertionError):
            alignment_info.get_seq(76000, 81000)
