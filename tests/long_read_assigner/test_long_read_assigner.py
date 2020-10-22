from collections import namedtuple, Counter

import pytest

from src.long_read_assigner import LongReadAssigner, MatchEventSubtype


Params = namedtuple("Params", ("delta", ))
IntronProfiles = namedtuple("IntronProfiles", ("features", ))
GeneInfo = namedtuple("GeneInfo", ("intron_profiles", "start", "end"))


class TestMatchProfileAndFindMatchingIsoforms:
    gene_info = GeneInfo(IntronProfiles([(1, 1)]), 1, 100)
    assigner = LongReadAssigner(gene_info, Params(0))  # we need no params here

    @pytest.mark.parametrize("read_gene_profile, isoform_profiles, hint, expected",
                             [([], dict(id1=[1], id2=[0]), None, [("id2", -1), ("id1", 0)])])
    def test_empty(self, read_gene_profile, isoform_profiles, hint, expected):
        with pytest.raises(AssertionError):
            self.check(read_gene_profile, isoform_profiles, hint, expected)

    @pytest.mark.parametrize("read_gene_profile, isoform_profiles, hint, expected",
                             [([-1, 1, 0, -1], dict(id1=[-1, 1, 0, -1], id2=[-1, 1, 0], id3=[-1, 1, 0, 1]),
                               {"id2", "id3"}, [("id2", -1), ("id3", 1)])])
    def test_different_length(self, read_gene_profile, isoform_profiles, hint, expected):
        with pytest.raises(AssertionError):
            self.check(read_gene_profile, isoform_profiles, hint, expected)

    def check(self, read_gene_profile, isoform_profiles, hint, expected):
        assert expected == self.assigner.match_profile(read_gene_profile, isoform_profiles, hint)
        expected = {x[0] for x in expected if x[1] == 0}
        assert expected == self.assigner.find_matching_isoforms(read_gene_profile, isoform_profiles, hint)

    @pytest.mark.parametrize("read_gene_profile, isoform_profiles, hint, expected",
                             [([-1, 1, 0, -1], dict(id1=[-1, 1, 0, -1], id2=[-1, 1, 0, -1]),
                               None, [("id1", 0), ("id2", 0)]),
                              ([-1, 1, 0, -1], dict(id1=[-1, 1, 1, -1], id2=[-1, 1, 0, -1]),
                               None, [("id1", 0), ("id2", 0)])])
    def test_all_equals(self, read_gene_profile, isoform_profiles, hint, expected):
        self.check(read_gene_profile, isoform_profiles, hint, expected)

    @pytest.mark.parametrize("read_gene_profile, isoform_profiles, hint, expected",
                             [([-1, 1, 0, -1], dict(id1=[-1, 1, -1, -1], id2=[1, 1, 1, -1]),
                               None, [("id1", 0), ("id2", 1)]),
                              ([-1, 1, 1, 0, -1], dict(id1=[-1, 1, 1, 1, -1], id2=[-1, -1, 0, 1, -1],
                                                       id3=[-1, 0, 1, -1, -1]),
                               None, [("id1", 0), ("id3", 0), ("id2", 1)])])
    def test_some_equals(self, read_gene_profile, isoform_profiles, hint, expected):
        self.check(read_gene_profile, isoform_profiles, hint, expected)

    @pytest.mark.parametrize("read_gene_profile, isoform_profiles, hint, expected",
                             [([-1, 1, 0, -1], dict(id1=[-1, -1, 1, -1], id2=[1, 0, 1, -1]),
                               None, [("id1", 1), ("id2", 1)]),
                              ([-1, 1, 1, 1, -1], dict(id1=[-1, 0, 1, 1, 1], id2=[1, 0, 0, -1, -1],
                                                       id3=[1, 0, 0, 0, 1]),
                               None, [("id1", 1), ("id2", 2), ("id3", 2)])])
    def test_no_equals(self, read_gene_profile, isoform_profiles, hint, expected):
        self.check(read_gene_profile, isoform_profiles, hint, expected)

    @pytest.mark.parametrize("read_gene_profile, isoform_profiles, hint, expected",
                             [([-1, 1, 0, -1], dict(id1=[-1, 1, 0, -1], id2=[0, 1, 0, -1], id3=[-1, 1, 0, 1]),
                               {"id2", "id3"}, [("id2", 0), ("id3", 1)]),  # skip matched return another one
                              ([-1, 1, 0, -1], dict(id1=[-1, 0, 0, -1], id2=[-1, 1, -1], id3=[-1, 1, 0, -1]),
                               {"id3"}, [("id3", 0)])])  # matched
    def test_hint(self, read_gene_profile, isoform_profiles, hint, expected):
        self.check(read_gene_profile, isoform_profiles, hint, expected)


class TestCompareJunctions:
    gene_info = GeneInfo(IntronProfiles([(1, 1)]), 1, 100)

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([(1, 10), (15,  20)], (1, 10), [(2, 10), (15,  19)], (1, 10), 1),
                              ([(1, 10), (15, 20)], (1, 10), [(1, 10), (15, 20)], (15, 20), 0),
                              ([(1, 10), (15, 20)], (1, 10), [(1, 10), (15, 21), (25, 30)], (1, 10), 1),
                              ([(15, 20), (25, 35)], (15, 20), [(1, 10), (15, 21), (25, 34)], (1, 10), 1)])
    def test_no_contradiction(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == 1
        assert match_events[0].event_type == MatchEventSubtype.none

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([(1, 10), (15,  20)], (1, 10), [(2, 10)], (1, 10), 1),
                              ([(1, 10), (15, 20), (25, 36)], (25, 30), [(1, 10), (15, 21)], (1, 10), 1)])
    def test_extra_intron_out_right(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == 1
        assert match_events[0].event_type == MatchEventSubtype.extra_intron_flanking_right

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([(1, 10), (15,  20)], (1, 10), [(15, 19)], (15, 19), 1),
                              ([(1, 5), (7, 10), (15, 20)], (1, 10), [(15, 19)], (15, 19), 1)],
                             ids=("extra_intron", "two_extra_intron"))
    def test_extra_intron_out_left(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == 1
        assert match_events[0].event_type == MatchEventSubtype.extra_intron_flanking_left

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([(1, 5), (8, 10), (15,  20)], (15, 20), [(1, 5), (15, 19)], (15, 19), 1)])
    def test_extra_intron(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == 1
        assert match_events[0].event_type == MatchEventSubtype.extra_intron

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([(1, 5), (25,  30)], (15, 20), [(1, 5), (8, 10), (15, 19), (25, 30)], (15, 19), 1)])
    def test_missed_intron_in(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        event_types = [match_events[i].event_type for i in range(len(match_events))]
        assert len(match_events) == 2
        assert set(event_types) == {MatchEventSubtype.intron_retention}

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([(1, 50)], (1, 50), [(1, 5), (45, 49)], (45, 49), 1)])
    def test_exon_skipping_novel_intron(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == 1
        assert match_events[0].event_type == MatchEventSubtype.exon_skipping_novel_intron

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([(1, 10), (15, 50)], (15, 50), [(1, 49)], (1, 49), 1)])
    def test_exon_gain_novel(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(self.gene_info, Params(delta))
        match_events = assigner.intron_comparator.compare_junctions(read_junctions, read_region,
                                                                    isoform_junctions, isoform_region)
        assert len(match_events) == 1
        assert match_events[0].event_type == MatchEventSubtype.exon_gain_novel


class TestDetectContradictionType:
    Params = namedtuple("Params", ("delta", "max_intron_shift"))
    gene_info = GeneInfo(IntronProfiles([(1, 1)]), 1, 100)

    def test(self):
        assigner = LongReadAssigner(self.gene_info, self.Params(1, 1))
        match_events = assigner.intron_comparator.detect_contradiction_type((0, 200), [(50, 75)], (0, 200), [(45, 70)], [((0, 0), (0, 0))])
        assert len(match_events) == 1
        assert match_events[0].event_type == MatchEventSubtype.intron_shift
