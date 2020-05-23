import pytest

from src.long_read_assigner import LongReadAssigner


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
