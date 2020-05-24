import pytest

from src.common import difference_in_present_features, list_to_str, rindex


class TestDifferenceInPresentFeatures:
    @pytest.mark.parametrize("profile1, profile2, expected", [([], [], 0), ([], [1], -1), ([1], [1], 0),
                                                              ([1, -1], [1, 1], 1)])
    def test_ones(self, profile1, profile2, expected):
        assert expected == difference_in_present_features(profile1, profile2)

    @pytest.mark.parametrize("profile1, profile2, expected", [([0], [0], 0), ([0, 1], [1, 1], 0),
                                                              ([-1, 1], [1, 0], 1), ([1, 0, -1], [0, 1, 1], 1)])
    def test_zeroes_doesnot_affect(self, profile1, profile2, expected):
        assert expected == difference_in_present_features(profile1, profile2)


class TestListToStr:
    def test_empty(self):
        assert "." == list_to_str([])

    @pytest.mark.parametrize("element_list, sep, expected", [([1], ",", "1"), ([1, 0, 2], ":", "1:0:2")])
    def test_not_empty(self, element_list, sep, expected):
        assert expected == list_to_str(element_list, sep)


class TestRindex:
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
