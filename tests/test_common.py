import pytest

from src.common import difference_in_present_features


class TestDifferenceInPresentFeatures:
    @pytest.mark.parametrize("profile1, profile2, expected", [([], [], 0), ([], [1], -1), ([1], [1], 0),
                                                              ([1, -1], [1, 1], 1)])
    def test_ones(self, profile1, profile2, expected):
        assert expected == difference_in_present_features(profile1, profile2)

    @pytest.mark.parametrize("profile1, profile2, expected", [([0], [0], 0), ([0, 1], [1, 1], 0),
                                                              ([-1, 1], [1, 0], 1), ([1, 0, -1], [0, 1, 1], 1)])
    def test_zeroes_doesnot_affect(self, profile1, profile2, expected):
        assert expected == difference_in_present_features(profile1, profile2)
