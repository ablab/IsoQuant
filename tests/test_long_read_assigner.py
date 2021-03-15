from collections import namedtuple, Counter

import pytest
import gffutils
import os

from src.long_read_assigner import *
from src.gene_info import *
from src.polya_finder import PolyAInfo


class Params:
    def __init__(self, delta):
        self.delta = delta
        self.minor_exon_extension = 15
        self.major_exon_extension = 100
        self.min_abs_exon_overlap = 10
        self.min_rel_exon_overlap = 0.2
        self.max_suspicious_intron_abs_len = 0
        self.max_suspicious_intron_rel_len = 0
        self.max_fake_terminal_exon_len = 0
        self.micro_intron_length = 0
        self.max_intron_abs_diff = 5
        self.max_intron_rel_diff = 0.1
        self.apa_delta = 50
        self.minimal_exon_overlap = 5
        self.minimal_intron_absence_overlap = 5
        self.max_intron_shift = 10
        self.max_missed_exon_len = 10
        self.resolve_ambiguous = AmbiguityResolvingMethod.monoexon_only
        self.correct_minor_errors = True

IntronProfiles = namedtuple("IntronProfiles", ("features", ))
GeneInfoTuple = namedtuple("GeneInfo", ("intron_profiles", "start", "end"))


class TestMatchProfileAndFindMatchingIsoforms:
    gene_info = GeneInfoTuple(IntronProfiles([(50, 60), (80, 100), (80, 110), (200, 210)]), 0, 300)
    assigner = LongReadAssigner(gene_info, Params(0))  # we need no params here

    @pytest.mark.parametrize("read_gene_profile, isoform_profiles, hint, expected",
                             [([], dict(id1=[1], id2=[-1]), None, [IsoformDiff("id2", 0), IsoformDiff("id1", 0)])])
    def test_empty(self, read_gene_profile, isoform_profiles, hint, expected):
        with pytest.raises(AssertionError):
            self.check(read_gene_profile, isoform_profiles, hint, expected)

    @pytest.mark.parametrize("read_gene_profile, isoform_profiles, hint, expected",
                             [([-1, 1, -1, 0], dict(id1=[-1, 1, 1, -1], id2=[-1, 1 ,-2], id3=[-1, 1, 1, 1]),
                               {"id2", "id3"}, [IsoformDiff("id2", -1), IsoformDiff("id3", 1)])])
    def test_different_length(self, read_gene_profile, isoform_profiles, hint, expected):
        with pytest.raises(AssertionError):
            self.check(read_gene_profile, isoform_profiles, hint, expected)

    def check(self, read_gene_profile, isoform_profiles, hint, expected):
        assert expected == self.assigner.match_profile(read_gene_profile, isoform_profiles, hint)
        expected = sorted([x[0] for x in expected if x[1] == 0])
        assert expected == self.assigner.find_matching_isoforms(read_gene_profile, isoform_profiles, hint)

    @pytest.mark.parametrize("read_gene_profile, isoform_profiles, hint, expected",
                             [([-1, 1, -1, 0], dict(id1=[-1, 1, -1, -1], id2=[-1, 1, -1, -2]),
                               None, [IsoformDiff("id1", 0), IsoformDiff("id2", 0)]),
                              ([0, 1, -1, 1], dict(id1=[-1, 1, -1, 1], id2=[-2, 1, -1, 1]),
                               None, [IsoformDiff("id1", 0), IsoformDiff("id2", 0)])])
    def test_all_equals(self, read_gene_profile, isoform_profiles, hint, expected):
        self.check(read_gene_profile, isoform_profiles, hint, expected)

    @pytest.mark.parametrize("read_gene_profile, isoform_profiles, hint, expected",
                             [([-1, 1, -1, 0], dict(id1=[-1, 1, -1, -1], id2=[1, 1, 1, -1]),
                               None, [IsoformDiff("id1", 0), IsoformDiff("id2", 2)]),
                              ([0, 1, 1, -1, 0], dict(id1=[-1, 1, 1,  1, -2], id2=[-2, 1, 1, -1, -2],
                                                       id3=[ 1, -1, 1, -2, -2]),
                               None, [IsoformDiff("id2", 0), IsoformDiff("id1", 1), IsoformDiff("id3", 2)])])
    def test_some_equals(self, read_gene_profile, isoform_profiles, hint, expected):
        self.check(read_gene_profile, isoform_profiles, hint, expected)

    @pytest.mark.parametrize("read_gene_profile, isoform_profiles, hint, expected",
                             [([-1, 1, -1, 0], dict(id1=[-2, 1, -1, -1], id2=[1, 1, 1, 1]),
                               None, [IsoformDiff("id1", 1), IsoformDiff("id2", 2)]),
                              ([-1, 1, 1, -1, 0], dict(id1=[-1, 1, 1, 1, 1], id2=[1, 1, -1, -1, -1],
                                                       id3=[-2, -1, -1, 1, -2]),
                               None, [IsoformDiff("id1", 1), IsoformDiff("id2", 2), IsoformDiff("id3", 4)])])
    def test_no_equals(self, read_gene_profile, isoform_profiles, hint, expected):
        self.check(read_gene_profile, isoform_profiles, hint, expected)

    @pytest.mark.parametrize("read_gene_profile, isoform_profiles, hint, expected",
                             [([-1, 1, 1, 0], dict(id1=[-1, 1, 1, -1], id2=[1, 1, 1, -1], id3=[-1, -1, -1, 1]),
                               {"id2", "id3"}, [IsoformDiff("id2", 1), IsoformDiff("id3", 2)]),  # skip matched return another one
                              ([-1, 1, -1, 0], dict(id1=[-1, 1, 1, -2], id2=[-1, -1, 1, -2], id3=[-1, 1, -1, -2]),
                               {"id3"}, [IsoformDiff("id3", 0)])])  # matched
    def test_hint(self, read_gene_profile, isoform_profiles, hint, expected):
        self.check(read_gene_profile, isoform_profiles, hint, expected)


class TestCompareJunctions:
    gene_info = GeneInfoTuple(IntronProfiles([(50, 60), (80, 100), (80, 110), (80, 150), (200, 210)]), 0, 300)

    @pytest.mark.parametrize("junctions, region, delta, expexted",
                             [([(50, 60), (80, 100), (200, 210)], (1, 1), 3, [0, 1, -1, 0]),
                              ([(50, 60), (80, 100), (200, 210)], (0, 1), 3, [1, 1, -1, 0])])
    def test_profile_for_junctions_introns(self, junctions, region, delta, expexted):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        profile = assigner.intron_comparator.profile_for_junctions_introns(junctions, region)
        assert expexted == profile.gene_profile

    @pytest.mark.parametrize("junctions, region, delta, expexted",
                             [([(50, 60), (82, 100), (200, 210)], (1, 1), 3, True),
                              ([(48, 61), (80, 98), (200, 210)], (0, 1), 3, True),
                              ([(48, 61), (80, 98), (160, 220)], (0, 1), 3, True),
                              ([(48, 66), (80, 99), (200, 210)], (0, 1), 3, False)])
    def test_profile_for_junctions_introns(self, junctions, region, delta, expexted):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        assert expexted == assigner.intron_comparator.are_known_introns(junctions, region)

    @pytest.mark.parametrize("read_intron_read_profile, read_region, read_introns, delta, expexted",
                             [([0, 1, 1], (0, 300), [(30, 60), (82, 100), (200, 210)], 3, MatchEventSubtype.extra_intron_flanking_left),
                              ([0, 1, 1], (0, 300), [(5, 60), (82, 100), (200, 210)], 3, MatchEventSubtype.fake_terminal_exon_left),
                              ([1, 1, 0], (0, 300), [(50, 60), (82, 100), (200, 210)], 3, MatchEventSubtype.extra_intron_flanking_right),
                              ([1, 1, 0], (0, 300), [(50, 60), (82, 100), (200, 292)], 3, MatchEventSubtype.fake_terminal_exon_right)])
    def test_add_extra_out_exon_events(self, read_intron_read_profile, read_region, read_introns, delta, expexted):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        assigner.params.max_fake_terminal_exon_len = 10
        events = []
        assigner.intron_comparator.add_extra_out_exon_events(events, read_intron_read_profile, read_region, read_introns, 100)
        assert expexted == events[0].event_type

    @pytest.mark.parametrize("read_region, junctions, region, delta, expexted",
                             [((0, 300), [(50, 60), (82, 100), (200, 210)], (0, 0), 3, True),
                              ((0, 300), [(48, 61), (80, 98), (200, 210)], (0, 1), 3, True),
                              ((0, 300), [(38, 61), (80, 98), (160, 220)], (0, 0), 3, False),
                              ((50, 300), [(58, 76), (80, 97), (130, 210)], (0, 1), 3, False)])
    def test_profile_for_junctions_introns(self, read_region, junctions, region, delta, expexted):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        assigner.params.max_suspicious_intron_abs_len = 20
        assigner.params.max_suspicious_intron_rel_len = 0.2
        assert expexted == assigner.intron_comparator.are_suspicious_introns(read_region, junctions, region)

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([], (20, 200), [], (20, 290), 1)])
    def test_monoexon_match(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == 1
        assert match_events[0].event_type == MatchEventSubtype.mono_exon_match

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([], (20, 200), [(150, 170)], (20, 290), 1)])
    def test_unspliced_intron_retention(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == 1
        assert match_events[0].event_type == MatchEventSubtype.unspliced_intron_retention

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([], (20, 100), [(50, 170)], (20, 290), 1)])
    def test_incomplete_intron_retention_right(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == 1
        assert match_events[0].event_type == MatchEventSubtype.incomplete_intron_retention_right

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([], (150, 320), [(50, 170)], (20, 290), 1)])
    def test_incomplete_intron_retention_left(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == 1
        assert match_events[0].event_type == MatchEventSubtype.incomplete_intron_retention_left

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([], (100, 150), [(50, 170)], (20, 290), 1),
                              ([], (20, 55), [(50, 170)], (20, 290), 1)])
    def test_mono_exonic(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == 1
        assert match_events[0].event_type == MatchEventSubtype.mono_exonic

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([(1, 10), (15,  20)], (0, 30), [(2, 10), (15,  19)], (0, 40), 3),
                              ([(1, 10), (15, 20)], (0, 30), [(1, 10), (15, 20)], (0, 50), 0),
                              ([(1, 10), (15, 20)], (0, 30), [(1, 10), (15, 21)], (0, 30), 1),
                              ([(15, 20), (25, 35)], (10, 40), [(1, 10), (15, 21), (25, 34)], (0, 40), 2)])
    def test_no_contradiction(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == 1
        assert match_events[0].event_type == MatchEventSubtype.none

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta, expected_len",
                             [([(1, 100), (150,  200)], (0, 300), [(2, 101)], (0, 120), 1, 1),
                              ([(1, 100), (150, 200), (250, 360), (380, 390)], (0, 400), [(3, 100), (150, 201)], (0, 249), 3, 2)])
    def test_extra_intron_out_right(self, read_junctions, read_region, isoform_junctions, isoform_region, delta, expected_len):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == expected_len
        assert match_events[0].event_type == MatchEventSubtype.extra_intron_flanking_right

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta, expected_len",
                             [([(1, 100), (150, 200)], (0, 300), [(150, 201)], (110, 220), 1, 1),
                              ([(1, 100), (150, 200), (250, 360)], (0, 400), [(251, 361)], (201, 405), 3, 2)])
    def test_extra_intron_out_left(self, read_junctions, read_region, isoform_junctions, isoform_region, delta, expected_len):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == expected_len
        assert match_events[0].event_type == MatchEventSubtype.extra_intron_flanking_left

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta, expected",
                             [([(20, 50), (60, 100), (150,  200)], (0, 300), [(20, 51), (150, 201)], (0, 290), 2,
                               MatchEventSubtype.extra_intron_novel),
                              ([(20, 40), (50, 60), (150, 200)], (0, 300), [(20, 41), (150, 201)], (0, 290), 1,
                               MatchEventSubtype.extra_intron_known),
                              ([(20, 40), (78, 112), (150, 200)], (0, 300), [(20, 41), (150, 201)], (0, 290), 3,
                               MatchEventSubtype.extra_intron_known),
                              ([(20, 40), (98, 102), (150, 200)], (0, 300), [(20, 41), (150, 201)], (0, 290), 3,
                               MatchEventSubtype.none),
                              ([(20, 50), (60, 100), (150, 200)], (15, 300), [(60, 100), (150, 201)], (45, 290), 2,
                               MatchEventSubtype.fake_terminal_exon_left),
                              ([(20, 50), (60, 100), (150, 200)], (0, 205), [(20, 50), (60, 100)], (0, 170), 2,
                               MatchEventSubtype.fake_terminal_exon_right)
                              ])
    def test_extra_intron(self, read_junctions, read_region, isoform_junctions, isoform_region, delta, expected):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        assigner.params.max_suspicious_intron_abs_len = 10
        assigner.params.max_suspicious_intron_rel_len = 0.1
        assigner.params.max_fake_terminal_exon_len = 10
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == 1
        assert match_events[0].event_type == expected

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([(10, 50), (150,  200)], (0, 300), [(10, 51), (80, 100), (150, 200), (225, 240)], (0, 310), 3)])
    def test_missed_intron_in(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        event_types = [match_events[i].event_type for i in range(len(match_events))]
        assert len(match_events) == 2
        assert set(event_types) == {MatchEventSubtype.intron_retention}

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta, expected",
                             [([(10, 50)], (0, 100), [(10, 25), (40, 49)], (0, 99), 3, MatchEventSubtype.exon_skipping_novel),
                              ([(80, 110)], (50, 150), [(80, 90), (105, 110)], (50, 149), 3, MatchEventSubtype.exon_skipping_known),
                              ([(80, 110)], (20, 150), [(65, 90), (105, 110)], (20, 149), 3, MatchEventSubtype.exon_merge_known),
                              ([(70, 110)], (20, 150), [(55, 90), (105, 110)], (20, 149), 3, MatchEventSubtype.exon_merge_novel),
                              ([(81, 109)], (50, 150), [(80, 98), (102, 110)], (50, 149), 3, MatchEventSubtype.exon_misalignment)])
    def test_exon_skipping(self, read_junctions, read_region, isoform_junctions, isoform_region, delta, expected):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == 1
        assert match_events[0].event_type == expected

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta, expected",
                             [([(10, 30), (50, 100)], (0, 120), [(10, 60), (80, 100)], (0, 119), 3,
                               MatchEventSubtype.mutually_exclusive_exons_novel),
                              ([(50, 60), (81, 149)], (0, 200), [(50, 100), (119, 150)], (0, 219), 3,
                               MatchEventSubtype.mutually_exclusive_exons_known),
                              ([(50, 60), (81, 130), (140, 181)], (0, 200), [(50, 100), (119, 200)], (0, 219), 3,
                               MatchEventSubtype.alternative_structure_novel),
                              ([(50, 60), (81, 150), (200, 210)], (0, 250), [(50, 99), (115, 210)], (0, 279), 3,
                               MatchEventSubtype.alternative_structure_known)
                              ])
    def test_complex_alteration(self, read_junctions, read_region, isoform_junctions, isoform_region, delta, expected):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == 1
        assert match_events[0].event_type == expected

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta, expected",
                             [([(10, 30), (55, 100)], (0, 150), [(10, 100)], (9, 149), 2,
                               MatchEventSubtype.exon_gain_novel),
                              ([(50, 60), (80, 100)], (0, 150), [(50, 100)], (9, 149), 2,
                               MatchEventSubtype.exon_gain_known),
                              ([(10, 30), (55, 100)], (0, 150), [(10, 75)], (9, 149), 2,
                               MatchEventSubtype.exon_detatch_novel),
                              ([(50, 60), (80, 110)], (0, 150), [(50, 90)], (9, 149), 2,
                               MatchEventSubtype.exon_detatch_known),
                              ([(50, 60), (80, 110)], (20, 150), [(20, 60), (80, 110)], (0, 150), 2,
                               MatchEventSubtype.terminal_exon_shift_known),
                              ([(45, 60), (80, 110)], (22, 150), [(10, 60), (80, 110)], (0, 150), 2,
                               MatchEventSubtype.terminal_exon_shift_novel),
                              ([(50, 60), (80, 110)], (30, 150), [(20, 60), (80, 110)], (0, 150), 2,
                               MatchEventSubtype.terminal_exon_misalignment_left),
                              ([(10, 13), (95, 100)], (0, 150), [(10, 100)], (9, 149), 2,
                               MatchEventSubtype.intron_retention)])
    def test_exon_gain(self, read_junctions, read_region, isoform_junctions, isoform_region, delta, expected):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        assigner.params.max_suspicious_intron_abs_len = 10
        assigner.params.max_suspicious_intron_rel_len = 0.1
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == 1
        assert match_events[0].event_type == expected

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta, expected",
                             [([(50, 70), (80, 100)], (0, 150), [(50, 60), (80, 100)], (9, 149), 2,
                               MatchEventSubtype.alt_right_site_novel),
                              ([(50, 60), (80, 110)], (0, 150), [(50, 60), (80, 100)], (9, 149), 2,
                               MatchEventSubtype.alt_right_site_known),
                              ([(50, 60), (90, 100)], (0, 150), [(50, 60), (80, 100)], (9, 149), 2,
                               MatchEventSubtype.alt_left_site_novel),
                              ([(50, 60), (80, 100)], (0, 150), [(40, 60), (80, 100)], (9, 149), 2,
                               MatchEventSubtype.alt_left_site_known),
                              ([(50, 60), (80, 100)], (0, 150), [(50, 60), (95, 115)], (9, 149), 2,
                               MatchEventSubtype.intron_migration),
                              ([(50, 60), (85, 105)], (0, 150), [(50, 60), (80, 100)], (9, 149), 2,
                               MatchEventSubtype.intron_shift),
                              ([(50, 60), (68, 85)], (0, 150), [(50, 60), (80, 100)], (9, 149), 2,
                               MatchEventSubtype.intron_alternation_novel),
                              ([(50, 60), (75, 115)], (0, 150), [(50, 60), (80, 100)], (9, 149), 2,
                               MatchEventSubtype.intron_alternation_novel),
                              ([(50, 60), (80, 110)], (0, 150), [(50, 60), (65, 100)], (9, 149), 2,
                               MatchEventSubtype.intron_alternation_known),
                              ([(50, 60), (98, 100)], (0, 150), [(50, 60), (80, 100)], (9, 149), 2,
                               MatchEventSubtype.intron_retention),
                              ([(50, 60), (93, 96)], (0, 150), [(50, 60), (80, 100)], (9, 149), 2,
                               MatchEventSubtype.intron_retention)
                              ])
    def test_single_introns(self, read_junctions, read_region, isoform_junctions, isoform_region, delta, expected):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        assigner.params.max_suspicious_intron_abs_len = 5
        assigner.params.max_suspicious_intron_rel_len = 0.1
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == 1
        assert match_events[0].event_type == expected

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta, expected",
                             [([(80, 100)], (0, 150), [(50, 60), (80, 100)], (9, 149), 2,
                               MatchEventSubtype.intron_retention),
                              ([(50, 60)], (0, 95), [(50, 60), (80, 100)], (9, 149), 2,
                               MatchEventSubtype.incomplete_intron_retention_right),
                              ([(80, 100)], (40, 150), [(20, 60), (80, 100)], (9, 149), 2,
                               MatchEventSubtype.incomplete_intron_retention_left),
                              ])
    def test_ir(self, read_junctions, read_region, isoform_junctions, isoform_region, delta, expected):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        assigner.params.max_suspicious_intron_abs_len = 5
        assigner.params.max_suspicious_intron_rel_len = 0.1
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == 1
        assert match_events[0].event_type == expected

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([(50, 60)], (0, 85), [(50, 60), (80, 100)], (9, 149), 2),
                              ([(80, 100)], (57, 150), [(20, 60), (80, 100)], (9, 149), 2),
                              ])
    def test_no_ir(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        assigner.params.max_suspicious_intron_abs_len = 5
        assigner.params.max_suspicious_intron_rel_len = 0.1
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert match_events[0].event_type == MatchEventSubtype.none


class TestDetectContradictionType:
    gene_info = GeneInfoTuple(IntronProfiles([(1, 1)]), 1, 100)

    def test(self):
        assigner = LongReadAssigner(self.gene_info, Params(3))
        match_events = assigner.intron_comparator.detect_contradiction_type((0, 200), [(50, 75)], (0, 200), [(45, 70)], [((0, 0), (0, 0))])
        assert len(match_events) == 1
        assert match_events[0].event_type == MatchEventSubtype.intron_shift


class TestAssignIsoform:
    params = Params(3)
    source_dir = os.path.dirname(os.path.realpath(__file__))
    gffutils_db = gffutils.FeatureDB(os.path.join(source_dir, 'toy_data/synth.db'), keep_order=True)
    gene_db = gffutils_db['ENSMUSG00000020196.10']
    gene_info = GeneInfo([gene_db], gffutils_db)
    assigner = LongReadAssigner(gene_info, params)
    gene_region = (gene_info.start, gene_info.end)

    intron_profile_constructor = \
        OverlappingFeaturesProfileConstructor(gene_info.intron_profiles.features, gene_region,
                                              comparator=partial(equal_ranges, delta=params.delta),
                                              absence_condition=partial(overlaps_at_least,
                                                                        delta=params.minimal_intron_absence_overlap),
                                              delta=params.delta)
    exon_profile_constructor = \
        OverlappingFeaturesProfileConstructor(gene_info.exon_profiles.features, gene_region,
                                              comparator=partial(equal_ranges, delta=params.delta),
                                              delta=params.delta)
    split_exon_profile_constructor = \
        NonOverlappingFeaturesProfileConstructor(gene_info.split_exon_profiles.features,
                                                 comparator=partial(overlaps_at_least,
                                                                    delta=params.minimal_exon_overlap),
                                                 delta=params.delta)

    @pytest.mark.parametrize("sorted_blocks, polya_pos, polyt_pos, isoform_id, expected_event",
                             [([(1000, 1100), (2000, 2100), (2300, 2400), (3000, 3300), (9500, 10000)], -1, -1,
                               "ENSMUST00000001712.7", MatchEventSubtype.fsm),
                              ([(1000, 1100), (2000, 2100), (2300, 2400), (3000, 3300)], -1, -1,
                               "ENSMUST00000001712.7", MatchEventSubtype.ism_right),
                              ([(2000, 2098), (2301, 2400), (3001, 3300), (9500, 10003)], -1, -1,
                               "ENSMUST00000001712.7", MatchEventSubtype.ism_left),
                              ([(2000, 2100), (2300, 2400), (3000, 3300)], -1, -1,
                               "ENSMUST00000001712.7", MatchEventSubtype.ism_internal),
                              ([(1000, 1100), (2000, 2200), (2500, 2600), (3000, 3300), (6000, 6010), (9500, 10000)],
                               -1, -1, "ENSMUST00000001713.7", MatchEventSubtype.fsm),
                              ([(1040, 1100), (2000, 2201), (2500, 2602), (2998, 3300), (6000, 6011), (9500, 9600)],
                               -1, -1, "ENSMUST00000001713.7", MatchEventSubtype.fsm),
                              ([(7998, 8201), (8500, 8800)], -1, -1,
                               "ENSMUST00000001715.7", MatchEventSubtype.fsm),
                              ([(7000, 7500)], -1, -1,
                               "ENSMUST00000001714.7", MatchEventSubtype.mono_exon_match),
                              ([(7100, 7300)], -1, -1,
                               "ENSMUST00000001714.7", MatchEventSubtype.mono_exon_match),
                              ([(6999, 7503)], -1, -1,
                               "ENSMUST00000001714.7", MatchEventSubtype.mono_exon_match),
                              ])
    def test_assign_unique(self, sorted_blocks, polya_pos, polyt_pos, isoform_id, expected_event):
        intron_profile = self.intron_profile_constructor.construct_intron_profile(sorted_blocks, polya_pos, polyt_pos)
        exon_profile = self.exon_profile_constructor.construct_exon_profile(sorted_blocks, polya_pos, polyt_pos)
        split_exon_profile = self.split_exon_profile_constructor.construct_profile(sorted_blocks, polya_pos, polyt_pos)
        combined_profile = CombinedReadProfiles(intron_profile, exon_profile, split_exon_profile,
                                                PolyAInfo(polya_pos, polyt_pos, -1, -1))

        read_assignment = self.assigner.assign_to_isoform("read_id", combined_profile)
        assert read_assignment.assignment_type == ReadAssignmentType.unique
        assert len(read_assignment.isoform_matches) == 1
        assert read_assignment.isoform_matches[0].assigned_transcript == isoform_id
        assert expected_event in {e.event_type for e in read_assignment.isoform_matches[0].match_subclassifications}

    @pytest.mark.parametrize("sorted_blocks, polya_pos, polyt_pos, isoform_id, expected_event",
                             [([(1000, 1108), (2007, 2100), (2300, 2400), (3000, 3300), (9500, 10000)], -1, -1,
                               "ENSMUST00000001712.7", MatchEventSubtype.intron_shift),
                              ([(1000, 1100), (2000, 2092), (2295, 2400), (3000, 3300)], -1, -1,
                               "ENSMUST00000001712.7", MatchEventSubtype.intron_shift),
                              ([(2000, 2098), (2301, 2400), (3001, 3300), (9500, 10012)], -1, -1,
                               "ENSMUST00000001712.7", MatchEventSubtype.exon_elongation_right),
                              ([(1000, 1100), (2000, 2200), (2500, 2600), (3000, 3303), (9496, 10000)],
                               -1, -1, "ENSMUST00000001713.7", MatchEventSubtype.exon_misalignment),
                              ([(988, 1100), (2000, 2200), (2500, 2602), (2998, 3300), (6000, 6011), (9500, 9600)],
                               -1, -1, "ENSMUST00000001713.7", MatchEventSubtype.exon_elongation_left),
                              ([(7998, 8201), (8500, 8800), (9900, 9907)], -1, -1,
                               "ENSMUST00000001715.7", MatchEventSubtype.fake_terminal_exon_right),
                              ([(4000, 4010), (7998, 8201), (8500, 8800)], -1, -1,
                               "ENSMUST00000001715.7", MatchEventSubtype.fake_terminal_exon_left),
                              ])
    def test_assign_unique_minor_error(self, sorted_blocks, polya_pos, polyt_pos, isoform_id, expected_event):
        self.assigner.params.max_missed_exon_len = 20
        self.assigner.params.max_fake_terminal_exon_len = 20
        intron_profile = self.intron_profile_constructor.construct_intron_profile(sorted_blocks, polya_pos, polyt_pos)
        exon_profile = self.exon_profile_constructor.construct_exon_profile(sorted_blocks, polya_pos, polyt_pos)
        split_exon_profile = self.split_exon_profile_constructor.construct_profile(sorted_blocks, polya_pos, polyt_pos)
        combined_profile = CombinedReadProfiles(intron_profile, exon_profile, split_exon_profile,
                                                PolyAInfo(polya_pos, polyt_pos, -1, -1))

        read_assignment = self.assigner.assign_to_isoform("read_id", combined_profile)
        assert read_assignment.assignment_type == ReadAssignmentType.unique_minor_difference
        assert len(read_assignment.isoform_matches) == 1
        assert read_assignment.isoform_matches[0].assigned_transcript == isoform_id
        assert expected_event in {e.event_type for e in read_assignment.isoform_matches[0].match_subclassifications}

    @pytest.mark.parametrize("sorted_blocks, polya_pos, polyt_pos, isoform_id",
                             [([(1000, 1108), (2007, 2100)], -1, -1,
                               "ENSMUST00000001712.7"),
                              ([(1050, 1100), (2001, 2092)], -1, -1,
                               "ENSMUST00000001712.7"),
                              ])
    def test_assign_unique_ambiguous(self, sorted_blocks, polya_pos, polyt_pos, isoform_id):
        self.assigner.params.max_missed_exon_len = 20
        self.assigner.params.max_fake_terminal_exon_len = 20
        intron_profile = self.intron_profile_constructor.construct_intron_profile(sorted_blocks, polya_pos, polyt_pos)
        exon_profile = self.exon_profile_constructor.construct_exon_profile(sorted_blocks, polya_pos, polyt_pos)
        split_exon_profile = self.split_exon_profile_constructor.construct_profile(sorted_blocks, polya_pos, polyt_pos)
        combined_profile = CombinedReadProfiles(intron_profile, exon_profile, split_exon_profile,
                                                PolyAInfo(polya_pos, polyt_pos, -1, -1))

        read_assignment = self.assigner.assign_to_isoform("read_id", combined_profile)
        assert read_assignment.assignment_type == ReadAssignmentType.ambiguous
        assert len(read_assignment.isoform_matches) > 1
        assert isoform_id in {im.assigned_transcript for im in read_assignment.isoform_matches}

    @pytest.mark.parametrize("sorted_blocks, polya_pos, polyt_pos, isoform_id, expected_event",
                             [([(1000, 1200), (2000, 2100), (2300, 2400), (3000, 3300), (9500, 10000)], -1, -1,
                               "ENSMUST00000001712.7", MatchEventSubtype.alt_left_site_novel),
                              ([(1000, 1100), (2000, 2400), (3000, 3300), (9500, 10000)], -1, -1,
                               "ENSMUST00000001712.7", MatchEventSubtype.intron_retention),
                              ([(500, 600), (1000, 1100), (2000, 2100), (2300, 2400), (3000, 3300), (9500, 10000)], -1, -1,
                               "ENSMUST00000001712.7", MatchEventSubtype.extra_intron_flanking_left),
                              ([(1000, 1100), (2000, 2100), (2300, 2400), (2700, 2750), (3000, 3300)], -1, -1,
                               "ENSMUST00000001712.7", MatchEventSubtype.exon_gain_novel),
                              ([(2000, 2098), (2301, 2400), (3001, 3300), (9500, 11003)], -1, -1,
                               "ENSMUST00000001712.7", MatchEventSubtype.major_exon_elongation_right),
                              ([(1000, 1100), (2000, 2200), (2500, 2600), (3000, 3300), (5000, 5010), (9500, 10000)],
                               -1, -1, "ENSMUST00000001713.7", MatchEventSubtype.mutually_exclusive_exons_novel),
                              ([(1040, 1100), (2000, 2201), (2500, 2602), (2998, 3300), (6000, 6011), (9500, 9600), (9700, 9999)],
                               -1, -1, "ENSMUST00000001713.7", MatchEventSubtype.extra_intron_novel),
                              ([(7998, 8201), (8300, 8333), (8500, 8800)], -1, -1,
                               "ENSMUST00000001715.7", MatchEventSubtype.exon_gain_novel),
                              ([(7100, 7800)], -1, -1,
                               "ENSMUST00000001714.7", MatchEventSubtype.major_exon_elongation_right),
                              ([(1040, 1100), (2000, 2201), (2500, 2602), (2998, 3900)], -1, -1,
                               "ENSMUST00000001713.7", MatchEventSubtype.incomplete_intron_retention_right),
                              ([(2500, 6010)], -1, -1,
                               "ENSMUST00000001713.7", MatchEventSubtype.unspliced_intron_retention),
                              ([(900, 1100), (2000, 2100), (2300, 2400), (3000, 3300), (9500, 10000)], -1, 900,
                               "ENSMUST00000001712.7", MatchEventSubtype.alternative_polya_site_left),
                              ([(1000, 1100), (1900, 2100)], -1, -1,
                               "ENSMUST00000001713.7", MatchEventSubtype.alt_right_site_novel),
                              ([(1000, 2100)], -1, -1,
                               "ENSMUST00000001713.7", MatchEventSubtype.unspliced_intron_retention),
                              ([(1000, 1100), (2000, 3300)], -1, -1,
                               "ENSMUST00000001713.7", MatchEventSubtype.intron_retention),
                              ])
    def test_assign_inconsistent(self, sorted_blocks, polya_pos, polyt_pos, isoform_id, expected_event):
        self.assigner.params.max_missed_exon_len = 20
        self.assigner.params.max_fake_terminal_exon_len = 20
        intron_profile = self.intron_profile_constructor.construct_intron_profile(sorted_blocks, polya_pos, polyt_pos)
        exon_profile = self.exon_profile_constructor.construct_exon_profile(sorted_blocks, polya_pos, polyt_pos)
        split_exon_profile = self.split_exon_profile_constructor.construct_profile(sorted_blocks, polya_pos, polyt_pos)
        combined_profile = CombinedReadProfiles(intron_profile, exon_profile, split_exon_profile,
                                                PolyAInfo(polya_pos, polyt_pos, -1, -1))

        read_assignment = self.assigner.assign_to_isoform("read_id", combined_profile)
        assert read_assignment.assignment_type == ReadAssignmentType.inconsistent
        if isoform_id:
            assert isoform_id in {im.assigned_transcript for im in read_assignment.isoform_matches}
        if expected_event:
            assert expected_event in {e.event_type for e in read_assignment.isoform_matches[0].match_subclassifications}
