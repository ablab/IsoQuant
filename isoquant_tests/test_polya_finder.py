
############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from collections import namedtuple

import pytest

from isoquant_lib.terminal_prediction.polya_finder import PolyAFinder


PolyAAlignment = namedtuple('PolyAAlignment', ('query_name', 'cigartuples', 'seq', 'reference_start', 'reference_end'))


TUPLES1 = [(4, 2), (0, 15), (1, 1), (0, 5), (2, 1), (0, 16), (3, 4642), (0, 5), (1, 1), (0, 1), (1, 1),
           (0, 15), (2, 1), (0, 4), (1, 1), (0, 6), (2, 1), (0, 19), (1, 1), (0, 32), (1, 1), (0, 3),
           (3, 7649), (0, 25), (1, 1), (0, 13), (1, 1), (0, 29)]

TUPLES2 = [(0, 167), (4, 30)]

TUPLES3 = [(0, 167), (4, 30), (5, 5)]

TUPLES4 = [(0, 150), (3, 100), (0, 10), (3, 100), (0, 11), (4, 27)]

TUPLES5 = [(0, 147), (3, 100), (0, 5), (2, 2), (0, 4), (3, 100), (0, 5), (1, 2), (0, 5), (4, 30)]

TUPLES6 = [(0, 147), (3, 100), (0, 5), (2, 2), (0, 4), (3, 100), (0, 5), (1, 2), (0, 5), (4, 30), (5, 5)]

SEQ1 = 'CTCAAGACCAAGAAGGACGACATGACCATGGCTTAAAAGAGTCTGCTCCCCACAGCCCCCTGCGAT' \
       'GGATGGACGGAGGAACCAGGGTCGGACGACCTCCGATGCTAAGAGCACTCCAACTGCTGCAAACCG' \
       'AAGAAGCAGGCATCGGAGACACTCCCAAACCAGGAGCGACCAAGCACTGGCGCATGTGACTCAAG'

SEQ2 = 'CTCAAGACCAAGAAGGACGACATGACCATGGCTTAAAAGAGTCTGCTCCCCACAGCCCCCTGCGAT' \
       'GGATGGACGGAGGAACCAGGGTCGGACGACCTCCGATGCTAAGAGCACTCCAACTGCTGCAAACCG' \
       'AAGAAGCAGGCATCGGACDCCTTTTTTTTTTAGGCGCCAAAAAAAAAAAACCAAAAAAAAAAAAA'

SEQ3 = 'TTTTTTTTTTTTTTTTTTTTCTCAAGACCAAGAAGGACGACATGACCATGGCTTAAAAGAGTCTGC' \
       'TCCCCACAGCCCCCTGCGATGGATGGACGGAGGAACCAGGGTCGGACGACCTCCGATGCTAAGAGC' \
       'ACTCCAACTGCTGCAAACCGAAGAAGCAGGCATCGGACGCCAAAAAAAAACCAAAAAAAAAAAAA'

SEQ4 = 'TTTTTTTATTTGTTTTTTTTTTTTTTTCTCAAGACCAAGAAGGACGACATGACCATGGCTTAAAAGAGTCTGC' \
       'TCCCCACAGCCCCCTGCGATGGATGGACGGAGGAACCAGGGTCGGACGACCTCCGATGCTAAGAGC' \
       'ACTCCAACTGCTGCAAACCGAAGAAGCAGGCATCGGACGCCAAAAAAAAACCAAAAAAAAAAAAA'


class TestPolyAFinder:
    polya_finder = PolyAFinder()

    @pytest.mark.parametrize("alignment",
                             [PolyAAlignment('aligned_segment1', TUPLES1, SEQ1, 51054, 63535)])
    def test_no_t_head(self, alignment):
        assert self.polya_finder.find_polyt_external(alignment) == -1

    @pytest.mark.parametrize("alignment",
                             [PolyAAlignment('aligned_segment1', TUPLES1, SEQ1, 51054, 63535)])
    def test_no_a_tail(self, alignment):
        assert self.polya_finder.find_polya_external(alignment) == -1

    @pytest.mark.parametrize("alignment, expected",
                             [(PolyAAlignment('aligned_segment1', TUPLES2, SEQ2, 51054, 63535), 63538),
                              (PolyAAlignment('aligned_segment1', TUPLES3, SEQ2, 51054, 63535), 63538)])
    def test_find_a_tail(self, alignment, expected):
        assert self.polya_finder.find_polya_external(alignment) == expected

    @pytest.mark.parametrize("alignment, expected",
                             [(PolyAAlignment('aligned_segment1', TUPLES4, SEQ2, 1000, 1300), 1300),
                              (PolyAAlignment('aligned_segment2', TUPLES5, SEQ2, 1000, 1300), 1303)])
    def test_find_fake_a_tail(self, alignment, expected):
        assert self.polya_finder.find_polya_external(alignment) == expected

    @pytest.mark.parametrize("alignment",
                             [PolyAAlignment('aligned_segment1', list(reversed(TUPLES2)), SEQ3, 51054, 63535),
                              PolyAAlignment('aligned_segment1', list(reversed(TUPLES3)), SEQ3, 51054, 63535)])
    def test_find_t_head(self, alignment):
        assert self.polya_finder.find_polyt_external(alignment) == 51043

    @pytest.mark.parametrize("alignment, expected",
                             [(PolyAAlignment('aligned_segment1', list(reversed(TUPLES5)), SEQ4, 1000, 1300), 996),
                              (PolyAAlignment('aligned_segment2', list(reversed(TUPLES6)), SEQ4, 1000, 1300), 996)])
    def test_find_fake_t_head(self, alignment, expected):
        assert self.polya_finder.find_polyt_external(alignment) == expected

    @pytest.mark.parametrize("seq, expected",
                             [('TTTATATTGGTATTTTGCC', -1),
                              ('AAAAAAAAAAACGAAAAAAAAAAAAAAAA', 0),
                              ('TGGCGGTACAAAAAAAAAAACGAAAAAAAAAAAAAAAA', 9)])
    def test_find_polya_pos(self, seq, expected):
        assert self.polya_finder.find_polya(seq) == expected
