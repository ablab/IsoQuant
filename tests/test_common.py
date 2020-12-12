import pytest

from src.common import *


class TestMiscFunction:
    def test_get_first(self):
        assert get_first_best_from_sorted([(1,1), (2,2), (3,3)]) == [1]
        assert get_first_best_from_sorted([(1,1), (2,1), (3,3)]) == [1, 2]

    @pytest.mark.parametrize("collection, elem, expected",
                             [([0, 1, 2], 2, 2), ([0, 1, 2], 0, 0), ([0, 1, 2], 1, 1), ([1, 1, 0], 1, 1),
                              ([1, 1, 1], 1, 2)],
                             ids=("last", "first", "middle", "several", "all"))
    def test_correct(self, collection, elem, expected):
        assert expected == rindex(collection, elem)

    @pytest.mark.parametrize("collection", [(), ("b", "c")])
    def test_no_such_elem(self, collection):
        with pytest.raises(ValueError, match="a is not in list"):
            rindex(collection, "a")


class TestRanges:
    @pytest.mark.parametrize("range1, range2, delta, expected",
                             [((1, 20), (1, 20), 0, True), ((1, 20), (2, 21), 1, True), ((1, 10), (1, 20), 5, False)],
                             ids=("equals", "less_than_delta", "bigger_than_delta"))
    def test_equal(self, range1, range2, delta, expected):
        assert expected == equal_ranges(range1, range2, delta)

    @pytest.mark.parametrize("range1, range2, expected",
                             [((1, 20), (1, 20), True), ((1, 20), (20, 121), True), ((1, 10), (11, 20), False),
                              ((111, 210), (11, 110), False), ((120, 210), (101, 121), True)])
    def test_overlaps(self, range1, range2, expected):
        assert expected == overlaps(range1, range2)

    @pytest.mark.parametrize("range1, range2, delta, expected",
                             [((1, 20), (20, 25), 0, True), ((18, 120), (2, 20), 3, True), ((1, 10), (6, 20), 6, False)])
    def test_overlaps_at_least(self, range1, range2, delta, expected):
        assert expected == overlaps_at_least(range1, range2, delta)

    @pytest.mark.parametrize("range1, range2, delta, expected",
                             [((1, 25), (20, 25), 0, True), ((18, 40), (16, 41), 3, True), ((10, 25), (5, 30), 7, False)])
    def test_contains_approx(self, range1, range2, delta, expected):
        assert expected == contains_approx(range1, range2, delta)

    @pytest.mark.parametrize("range1, range2, expected",
                             [((1, 20), (1, 20), True), ((1, 120), (20, 101), True), ((1, 10), (11, 20), False),
                              ((111, 210), (111, 211), False), ((120, 210), (119, 121), False)])
    def test_contains(self, range1, range2, expected):
        assert expected == contains(range1, range2)


    @pytest.mark.parametrize("range1, range2, expected",
                             [((1, 20), (1, 20), True), ((1, 120), (20, 201), False), ((1, 10), (11, 20), False),
                              ((111, 210), (111, 211), False), ((120, 210), (119, 201), True)])
    def test_covers_start(self, range1, range2, expected):
        assert expected == covers_start(range1, range2)

    @pytest.mark.parametrize("range1, range2, expected",
                             [((1, 20), (1, 20), True), ((1, 120), (20, 201), True), ((1, 10), (11, 20), False),
                              ((111, 210), (111, 211), True), ((120, 210), (121, 211), True)])
    def test_covers_end(self, range1, range2, expected):
        assert expected == covers_end(range1, range2)

    @pytest.mark.parametrize("range, expected",
                             [((1, 20), 20), ((1, 1), 1)])
    def test_interval_len(self, range, expected):
        assert expected == interval_len(range)

    @pytest.mark.parametrize("ranges, expected",
                             [([(1, 20), (40, 50)], 31), ([(1, 1), (5, 5)], 2)])
    def test_intervals_len(self, ranges, expected):
        assert expected == intervals_total_length(ranges)

    @pytest.mark.parametrize("ranges, point, expected",
                             [([(1, 20), (40, 50)], 60, 31), ([(1, 1), (5, 5)], 10, 2),
                              ([(1, 20), (40, 50)], 39, 20), ([(1, 20), (40, 50)], 41, 21)])
    def test_intervals_len_to_point(self, ranges, point, expected):
        assert expected == sum_intervals_to_point(ranges, point)

    @pytest.mark.parametrize("ranges, point, expected",
                             [([(1, 20), (40, 50)], 1, 30), ([(1, 1), (5, 5)], 0, 2),
                              ([(1, 20), (40, 50)], 39, 11), ([(1, 20), (40, 50)], 11, 20)])
    def test_intervals_len_from_point(self, ranges, point, expected):
        assert expected == sum_intervals_from_point(ranges, point)


class TestSimilarityFunctions:
    @pytest.mark.parametrize("ranges1, ranges2, expected",
                             [([(1, 70)], [(31, 100)], 0.4),
                              ([(1, 20), (41, 50), (61, 70)], [(1, 20), (41, 50), (61, 70)], 1.0),
                              ([(1, 20), (41, 50), (61, 70)], [(41, 50), (61, 70), (101, 160)], 0.2),
                              ([(1, 20), (41, 50), (61, 70)], [(46, 55), (60, 62), (70, 74), (101, 150)], 0.08),
                              ([(1, 20), (41, 50), (61, 70)], [(146, 155), (160, 172)], 0.0),
                              ([(146, 155), (160, 172)], [(1, 20), (41, 50), (61, 70)], 0.0),
                              ([(11, 20), (41, 50), (61, 70)], [(1, 4), (10, 12), (14, 16), (19, 30), (30, 34), (60, 70), (101, 149)], 0.17)])
    def test_jaccard_index(self, ranges1, ranges2, expected):
        assert expected == jaccard_similarity(ranges1, ranges2)

    @pytest.mark.parametrize("read_exons, isoform_exons, expected",
                             [([(1, 100)], [(31, 110)], 0.7),
                              ([(1, 20), (41, 50), (61, 70)], [(1, 20), (41, 50), (61, 70)], 1.0),
                              ([(1, 20), (41, 50), (61, 80)], [(41, 50), (61, 70), (101, 110)], 0.4),
                              ([(1, 20), (41, 50), (61, 80)], [(0, 0), (46, 55), (60, 62), (70, 74), (80,82), (101, 110)], 0.26),
                              ([(1, 20), (41, 50), (61, 80)], [(146, 155), (160, 172)], 0.0),
                              ([(146, 155), (160, 172)], [(1, 20), (41, 50), (61, 70)], 0.0),
                              ([(11, 20), (41, 50), (61, 70), (201, 220)], [(1, 4), (10, 12), (14, 16), (19, 30), (30, 34), (60, 70), (101, 109)], 0.34)])
    def test_read_coverage_fraction(self, read_exons, isoform_exons, expected):
        assert expected == read_coverage_fraction(read_exons, isoform_exons)

    @pytest.mark.parametrize("read_exons, isoform_region, expected",
                             [([(1, 100)], (31, 110), 0.3),
                              ([(1, 100)], (31, 70), 0.6),
                              ([(1, 20), (41, 50), (61, 70)], (1, 70), 0.0),
                              ([(1, 20), (41, 50), (61, 70)], (0, 71), 0.0),
                              ([(1, 20), (41, 50), (61, 80)], (35, 55), 0.8),
                              ([(1, 20), (41, 50), (61, 80)], (146, 172), 1.0),
                              ([(11, 20), (41, 50), (61, 90)], (16, 70), 0.5)])
    def test_extra_flanking(self, read_exons, isoform_region, expected):
        assert expected == extra_exon_percentage(isoform_region, read_exons)


class TestExonFunctions:
    @pytest.mark.parametrize("exons, expected",
                             [([(1, 100)], []),
                              ([(1, 20), (41, 50), (61, 70)], [(21, 40), (51, 60)])])
    def test_junctions_from_blocks(self, exons, expected):
        assert expected == junctions_from_blocks(exons)

    @pytest.mark.parametrize("region, introns, intron_position, expected",
                             [((0, 100), [(20, 80)], 0, (81, 100)),
                              ((0, 100), [(20, 40), (60, 80)], 0, (41, 59)),
                              ((0, 100), [(20, 40), (60, 80)], 1, (81, 100))])
    def test_following_exon(self, region, introns, intron_position, expected):
        assert expected == get_following_exon_from_junctions(region, introns, intron_position)

    @pytest.mark.parametrize("region, introns, intron_position, expected",
                             [((0, 100), [(20, 80)], 0, (0, 19)),
                              ((0, 100), [(20, 40), (60, 80)], 0, (0, 19)),
                              ((0, 100), [(20, 40), (60, 80)], 1, (41, 59))])
    def test_preceding_exon(self, region, introns, intron_position, expected):
        assert expected == get_preceding_exon_from_junctions(region, introns, intron_position)

    @pytest.mark.parametrize("region, introns, exon_position, expected",
                             [((0, 100), [(20, 80)], 0, (0, 19)),
                              ((0, 100), [(20, 80)], 1, (81, 100)),
                              ((0, 100), [(20, 40), (60, 80)], 0, (0, 19)),
                              ((0, 100), [(20, 40), (60, 80)], 1, (41, 59)),
                              ((0, 100), [(20, 40), (60, 80)], 2, (81, 100))])
    def test_get_exon(self, region, introns, exon_position, expected):
        assert expected == get_exon(region, introns, exon_position)

    @pytest.mark.parametrize("exons, expected",
                             [([(1, 100)], [(2, 100)]),
                              ([(1, 20), (41, 50), (61, 70)], [(2, 20), (42, 50), (62, 70)])])
    def test_correct_bam_coords(self, exons, expected):
        assert expected == correct_bam_coords(exons)

    @pytest.mark.parametrize("exons, cigar_tuples, expected",
                             [([(1, 100)], [(0, 99)], [(1, 100)]),
                              ([(10, 20), (20, 30), (41, 50), (51, 60), (61, 70), (72, 80), (100, 110), (110, 114), (117, 120)],
                               [(4, 10), (0, 10), (1, 2), (0, 10), (3, 11), (0, 9), (2, 1), (0, 9), (2, 1), (0, 9), (2, 2), (0, 8), (3, 20), (0, 10), (1, 3), (0, 4),  (2, 3), (0, 3), (5, 10)],
                               [(10, 30), (41, 80), (100, 120)])])
    def test_concat(self, exons, cigar_tuples, expected):
        assert expected == concat_gapless_blocks(exons, cigar_tuples)


class TestProfiles:
    @pytest.mark.parametrize("profile1, profile2, expected", [([-2, 1], [1, 1], True), ([-2, -1, -2], [-1, -1, -1], True),
                                                              ([1], [1], True),
                                                              ([-2, 1, -1, -2], [-1, 1, 1, 1], False),
                                                              ([-2, 1, -1, -2], [-1, 1, -1, 1], True)])
    def test_is_subprofile(self, profile1, profile2, expected):
        assert expected == is_subprofile(profile1, profile2)

    @pytest.mark.parametrize("profile1, profile2, expected", [([], [], 0), ([0], [1], 0), ([1], [1], 0),
                                                              ([1, -1], [1, 1], 1),
                                                              ([0], [0], 0), ([0, 1], [1, 1], 0),
                                                              ([-1, 1], [1, 0], 1), ([1, 0, -1], [0, 1, 1], 1)])
    def test_difference_in_present_features(self, profile1, profile2, expected):
        assert expected == difference_in_present_features(profile1, profile2)

    @pytest.mark.parametrize("profile1, profile2, expected", [([0], [1], 0), ([0, 1], [1, 1], 1),
                                                              ([0, 1], [1, -2], 0), ([0, 1, -1], [-2, 1, 1], 1)])
    def test_count_both_present_features(self, profile1, profile2, expected):
        assert expected == count_both_present_features(profile1, profile2)

    @pytest.mark.parametrize("profile1, profile2, expected", [([1], [1], True), ([1, 1], [0, 1], False),
                                                              ([-1, 1, 1, -2], [0, 1, 1, -2], True),
                                                              ([-2, 1, 1], [0, 1, -1], False)])
    def test_all_features_present(self, profile1, profile2, expected):
        assert expected == all_features_present(profile1, profile2)

    @pytest.mark.parametrize("profile1, profile2, expected", [([1], [1], [1]), ([0, 1], [1, 1], [0, 1]),
                                                              ([0, 1, 1, 0], [-1, 1, 1, -2], [0, 1, 1, 0]),
                                                              ([0, 1, -1], [-2, 1, 1], [0, 1, 0])])
    def test_find_matching_positions(self, profile1, profile2, expected):
        assert expected == find_matching_positions(profile1, profile2)

    @pytest.mark.parametrize("profile1, profile2, expected", [([1], [1], True), ([0, 1], [1, 1], True),
                                                              ([0, 1, 1, 0], [-1, 1, 1, -2], True),
                                                              ([1, -1, 0, 0], [-1, 1, 1, -2], False),
                                                              ([0, 1, -1], [-2, 1, 1], True)])
    def test_has_overlapping_features(self, profile1, profile2, expected):
        assert expected == has_overlapping_features(profile1, profile2)

    @pytest.mark.parametrize("profile1, profile2, expected", [([1], [1], False), ([0, 1], [1, 1], False),
                                                              ([0, 1, 1, 0], [-1, 1, 1, -2], False),
                                                              ([1, -1, 0, 0], [-1, 1, 1, -2], True),
                                                              ([0, 1, -1], [-2, 1, 1], True)])
    def test_has_inconsistent_features(self, profile1, profile2, expected):
        assert expected == has_inconsistent_features(profile1, profile2)

    @pytest.mark.parametrize("profile1, mask, expected", [([1], [1], [1]), ([0, 1], [1, 1], [0, 1]),
                                                              ([0, 1, 1, 0], [1, 1, 0, 0], [0, 1, 0, 0]),
                                                              ([-2, -1, 1, 0], [1, 1, 0, 1], [-2, -1, 0, 0])])
    def test_mask_profile(self, profile1, mask, expected):
        assert expected == mask_profile(profile1, mask)

    @pytest.mark.parametrize("profile1, features, expected", [([1], [(10, 20)], [(10, 20)]), ([0, 1], [(10, 20), (30, 40)], [(30, 40)]),
                                                              ([0, 1, 0, 1], [(10, 20), (30, 40), (110, 120), (130, 140)], [(30, 40), (130, 140)])])
    def test_get_blocks_from_profile(self, profile1, features, expected):
        assert expected == get_blocks_from_profile(features, profile1)

    @pytest.mark.parametrize("profile1, profile2, expected", [([1], [1], False), ([0, 1], [1, 1], True),
                                                              ([0, 1, 1, 0], [-1, 1, 1, -2], False),
                                                              ([-1, -1, 1, 0], [-1, 1, 1, -2], True),
                                                              ([0, 1, -1], [-2, 1, 1], False)])
    def test_left_truncated(self, profile1, profile2, expected):
        assert expected == left_truncated(profile1, profile2)

    @pytest.mark.parametrize("profile1, profile2, expected", [([1], [1], False), ([1, 0], [1, 1], True),
                                                              ([0, 1, 1, 0], [-1, 1, 1, -2], False),
                                                              ([-1, 1, -1, 0], [-1, 1, 1, -2], True),
                                                              ([0, 1, -1], [-2, 1, 1], True)])
    def test_right_truncated(self, profile1, profile2, expected):
        assert expected == right_truncated(profile1, profile2)


class TestListToStr:
    def test_empty(self):
        assert "." == list_to_str([])

    @pytest.mark.parametrize("element_list, sep, expected", [([1], ",", "1"), ([1, 0, 2], ":", "1:0:2")])
    def test_not_empty(self, element_list, sep, expected):
        assert expected == list_to_str(element_list, sep)

