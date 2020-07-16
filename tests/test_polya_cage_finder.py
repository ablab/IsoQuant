import os
from collections import namedtuple, Counter

import pytest

from src.polya_finder import PolyAFinder, CagePeakFinder


PolyAAlignment = namedtuple('PolyAAlignment', ('query_name', 'cigartuples', 'seq', 'reference_start', 'reference_end'))
CageAlignment = namedtuple('CageAlignment', ('query_name', 'reference_name', 'query_alignment_start',
                                             'query_alignment_end', 'is_reverse'))


TUPLES1 = [(4, 2), (0, 15), (1, 1), (0, 5), (2, 1), (0, 16), (3, 4642), (0, 5), (1, 1), (0, 1), (1, 1),
           (0, 15), (2, 1), (0, 4), (1, 1), (0, 6), (2, 1), (0, 19), (1, 1), (0, 32), (1, 1), (0, 3),
           (3, 7649), (0, 25), (1, 1), (0, 13), (1, 1), (0, 29)]

TUPLES2 = [(0, 167), (4, 30)]

TUPLES3 = [(0, 167), (4, 30), (5, 30)]

SEQ1 = 'CTCAAGACCAAGAAGGACGACATGACCATGGCTTAAAAGAGTCTGCTCCCCACAGCCCCCTGCGAT' \
       'GGATGGACGGAGGAACCAGGGTCGGACGACCTCCGATGCTAAGAGCACTCCAACTGCTGCAAACCG' \
       'AAGAAGCAGGCATCGGAGACACTCCCAAACCAGGAGCGACCAAGCACTGGCGCATGTGACTCAAG'

SEQ2 = 'CTCAAGACCAAGAAGGACGACATGACCATGGCTTAAAAGAGTCTGCTCCCCACAGCCCCCTGCGAT' \
       'GGATGGACGGAGGAACCAGGGTCGGACGACCTCCGATGCTAAGAGCACTCCAACTGCTGCAAACCG' \
       'AAGAAGCAGGCATCGGACDCCTTTTTTTTTTAGGCGCCAAAAAAAAAAAACCAAAAAAAAAAAAA'

SEQ3 = 'TTTTTTTTTTTTTTTTTTTTCTCAAGACCAAGAAGGACGACATGACCATGGCTTAAAAGAGTCTGC' \
       'TCCCCACAGCCCCCTGCGATGGATGGACGGAGGAACCAGGGTCGGACGACCTCCGATGCTAAGAGC' \
       'ACTCCAACTGCTGCAAACCGAAGAAGCAGGCATCGGACDCCAAAAAAAAACCAAAAAAAAAAAAA'


class TestPolyAFinder:
    polya_finder = PolyAFinder()

    @pytest.mark.parametrize("alignment",
                             [PolyAAlignment('aligned_segment1', TUPLES1, SEQ1, 51054, 63535)])
    def test_no_t_head(self, alignment):
        assert self.polya_finder.find_polyt_head(alignment) == -1

    @pytest.mark.parametrize("alignment",
                             [PolyAAlignment('aligned_segment1', TUPLES1, SEQ1, 51054, 63535)])
    def test_no_a_tail(self, alignment):
        assert self.polya_finder.find_polya_tail(alignment) == -1

    @pytest.mark.parametrize("alignment, expected",
                             [(PolyAAlignment('aligned_segment1', TUPLES2, SEQ2, 51054, 63535), 63535),
                              (PolyAAlignment('aligned_segment1', TUPLES3, SEQ2, 51054, 63535), 63535)])
    def test_find_a_tail(self, alignment, expected):
        assert self.polya_finder.find_polya_tail(alignment) == expected

    @pytest.mark.parametrize("alignment",
                             [PolyAAlignment('aligned_segment1', list(reversed(TUPLES2)), SEQ3, 51054, 63535)])
    def test_find_t_head(self, alignment):
        assert self.polya_finder.find_polyt_head(alignment) == 51054

    @pytest.mark.parametrize("seq, expected",
                             [('TTTATATTTTTTTTTTGCC', -1),
                              ('AAAAAAAAAAACGAAAAAAAAAAAAAAAA', 0)])
    def test_find_polya_pos(self, seq, expected):
        assert self.polya_finder.find_polya(seq) == expected


class TestCagePeakFinder:
    source_dir = os.path.dirname(os.path.realpath(__file__))
    cage_peak_finder = CagePeakFinder(os.path.join(source_dir, 'toy_data/MAPT.Mouse.CAGE.bed'), 50)

    @pytest.mark.parametrize("alignment",
                             [CageAlignment('aligned_segment1', 'chr11', 24, 2165, False),
                              CageAlignment('aligned_segment2', 'chr11', 0, 1296, False),
                              CageAlignment('aligned_segment3', 'chr11', 80, 372, True)])
    def test_no_intersection(self, alignment):
        assert self.cage_peak_finder.find_cage_peak(alignment) == []

    @pytest.mark.parametrize("alignment",
                             [CageAlignment('aligned_segment1', 'chr11', 57, 522, False),
                              CageAlignment('aligned_segment2', 'chr11', 69, 1499, False)])
    def test_one_intersection(self, alignment):
        assert len(self.cage_peak_finder.find_cage_peak(alignment)) == 1

    @pytest.mark.parametrize("alignment",
                             [CageAlignment('aligned_segment1', 'chr11', 134, 1322, False),
                              CageAlignment('aligned_segment2', 'chr11', 102, 2191, False),
                              CageAlignment('aligned_segment2', 'chr11', 11, 126, True)])
    def test_several_hits(self, alignment):
        assert len(self.cage_peak_finder.find_cage_peak(alignment)) > 1

