import pytest

from IlluminaExonCorrector import IlluminaExonCorrector, VoidExonCorrector

class TestExonCorrection:
    
    @pytest.mark.parametrize('short_introns', 'long_exons',
                             [({(5154787, 5159231),(5159335, 5163585),(5163676, 5165837),(5163717, 5165837),(5165952, 5168251),(5168357, 5171292),(5171347, 5187612),(5187711, 5194499),(5194693, 5203306),(5203486, 5206034),(5206161, 5213972),(5214075, 5220170),
                             (5219146, 5324823),(5220285, 5232327)},
                             [(5154647, 5154786), (5159232, 5159334), (5163586, 5163675),(5165838, 5165951), (5168252, 5168356), (5171293, 5171346), (5187613, 5187710),(5194500, 5194692), (5203307, 5203485), (5206035, 5206160), (5213973, 5214074),
                             (5220171, 5220284), (5232328, 5232505)])
                             ])
    def test_no_correction:
        # not using short exons as those are never calculated
        # from the bam file I get introns directly
        corrector = IlluminaExonCorrector.from_data(short_introns)
        corrected_exons = corrector.correct_exons(long_exons)
        assert corrected_exons == long_exons
    
