from collections import namedtuple

import pytest

from src.long_read_assigner import LongReadAssigner, MatchingEvent


class TestMatchProfileAndFindMatchingIsoforms:
    assigner = LongReadAssigner(None, None)  # we need no params here

    @pytest.mark.parametrize("read_gene_profile, isoform_profiles, hint, expected",
                             [([], dict(id1=[], id2=[0]), None, [("id2", -1), ("id1", 0)])])
    def test_empty(self, read_gene_profile, isoform_profiles, hint, expected):
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
                             [([-1, 1, 0, -1], dict(id1=[-1, 1, 0, -1], id2=[-1, 1, 0], id3=[-1, 1, 0, 1]),
                               None, [("id2", -1), ("id1", 0), ("id3", 1)]),
                              ([-1, 1, 0, -1], dict(id1=[-1, 0, 0, -1], id2=[-1, 1, -1], id3=[-1, 1, 1, 1]),
                               None, [("id2", -1), ("id1", 0), ("id3", 1)])])
    def test_has_diff_length(self, read_gene_profile, isoform_profiles, hint, expected):
        self.check(read_gene_profile, isoform_profiles, hint, expected)

    @pytest.mark.parametrize("read_gene_profile, isoform_profiles, hint, expected",
                             [([-1, 1, 0, -1], dict(id1=[-1, 1, 0, -1], id2=[-1, 1, 0], id3=[-1, 1, 0, 1]),
                               {"id2", "id3"}, [("id2", -1), ("id3", 1)]),  # skip matched return nothing
                              ([-1, 1, 0, -1], dict(id1=[-1, 1, 0, -1], id2=[0, 1, 0, -1], id3=[-1, 1, 0, 1]),
                               {"id2", "id3"}, [("id2", 0), ("id3", 1)]),  # skip matched return another one
                              ([-1, 1, 0, -1], dict(id1=[-1, 0, 0, -1], id2=[-1, 1, -1], id3=[-1, 1, 0, -1]),
                               {"id3"}, [("id3", 0)])])  # matched
    def test_hint(self, read_gene_profile, isoform_profiles, hint, expected):
        self.check(read_gene_profile, isoform_profiles, hint, expected)


class TestCompareJunctions:
    Params = namedtuple("Params", "delta")

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([(1, 10), (15,  20)], (1, 10), [(2, 10), (15,  19)], (1, 10), 1),
                              ([(1, 10), (15, 20)], (1, 10), [(1, 10), (15, 20)], (15, 20), 0),
                              ([(1, 10), (15, 20)], (1, 10), [(1, 10), (15, 21), (25, 30)], (1, 10), 1),
                              ([(15, 20), (25, 35)], (15, 20), [(1, 10), (15, 21), (25, 34)], (1, 10), 1)])
    def test_no_contradiction(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(None, self.Params(delta))
        assert (MatchingEvent.no_contradiction
                == assigner.compare_junctions(read_junctions, read_region, isoform_junctions, isoform_region))

    @pytest.mark.parametrize("isoform_junctions, isoform_region, delta",
                             [([(2, 10)], (1, 10), 1),
                              ([(1, 10), (15, 21)], (1, 10), 1)])
    def test_unspliced(self, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(None, self.Params(delta))
        assert (MatchingEvent.unspliced == assigner.compare_junctions([], [], isoform_junctions, isoform_region))

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([(1, 10), (15,  20)], (1, 10), [(2, 10)], (1, 10), 1),
                              ([(1, 10), (15, 20), (25, 36)], (25, 30), [(1, 10), (15, 21)], (1, 10), 1)])
    def test_extra_intron_out(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(None, self.Params(delta))
        assert (MatchingEvent.extra_intron_out
                == assigner.compare_junctions(read_junctions, read_region, isoform_junctions, isoform_region))

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([(1, 10), (15,  20)], (1, 10), [(15, 19)], (15, 19), 1),
                              ([(1, 5), (7, 10), (15, 20)], (1, 10), [(15, 19)], (15, 19), 1)],
                              ids=("extra_intron", "two_extra_intron"))
    def test_extra_introns_before(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(None, self.Params(delta))
        assert (MatchingEvent.extra_intron_out
                == assigner.compare_junctions(read_junctions, read_region, isoform_junctions, isoform_region))

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([(1, 5), (8, 10), (15,  20)], (15, 20), [(1, 5), (15, 19)], (15, 19), 1)])
    def test_extra_intron_in(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(None, self.Params(delta))
        assert (MatchingEvent.extra_intron_out
                == assigner.compare_junctions(read_junctions, read_region, isoform_junctions, isoform_region))

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([(1, 5), (25,  30)], (15, 20), [(1, 5), (8, 10), (15, 19)], (15, 19), 1),
                              ([(1, 5), (25,  30)], (15, 20), [(1, 5), (8, 10), (15, 19), (25, 30)], (15, 19), 1)])
    def test_missed_intron_in(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(None, self.Params(delta))
        assert (MatchingEvent.intron_retention
                == assigner.compare_junctions(read_junctions, read_region, isoform_junctions, isoform_region))

    @pytest.mark.parametrize("read_junctions, read_region, isoform_junctions, isoform_region, delta",
                             [([(1, 50)], (1, 50), [(1, 5), (45, 49)], (45, 49), 1),
                              ([(1, 10), (15, 50)], (15, 50), [(1, 49)], (1, 49), 1)],
                             ids=("one_in_read_two_in_isoforms", "two_in_read_one_in_isoforms"))
    def test_overlap(self, read_junctions, read_region, isoform_junctions, isoform_region, delta):
        assigner = LongReadAssigner(None, self.Params(delta))
        assert (MatchingEvent.intron_retention
                == assigner.compare_junctions(read_junctions, read_region, isoform_junctions, isoform_region))


class TestDetectContradictionType:
    Params = namedtuple("Params", ("delta", "max_intron_shift"))
    IntronProfiles = namedtuple("IntronProfiles", ("features"))
    GeneInfo = namedtuple("GeneInfo", ("intron_profiles", "start", "end"))

    def test(self):
        assigner = LongReadAssigner(self.GeneInfo(self.IntronProfiles([(1, 1)]), 1, 100), self.Params(1, 1))
        assert "intron_shift" is assigner.detect_contradiction_type([(50, 75)], [(50, 75)], [((0, 0), (0, 0))])
