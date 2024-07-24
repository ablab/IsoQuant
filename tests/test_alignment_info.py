from src.alignment_info import AlignmentInfo
import pytest

from src.polya_finder import PolyAFinder
from src.polya_verification import PolyAFixer

SEQ1 = 'CTCAAGACCAAGAAGGACGACATGACCATGGCTTAAAAGAGTCTGCTCCCCACAGCCCCCTGCGAT'
TUPLES1 = [(0, 99)]
TUPLES2 = [(1, 475), (75, 899)]
PAIRS = [(0, 75897), (1, None), (2, 75995)]


class Alignment:
    def __init__(self, query_name, cigartuples, seq, reference_start, reference_end, reference_name, reference_id,
                 aligned_pairs):
        self.query_name = query_name
        self.cigartuples = cigartuples
        self.seq = seq
        self.reference_start = reference_start
        self.reference_end = reference_end
        self.reference_name = reference_name
        self.reference_id = reference_id
        self.aligned_pairs = aligned_pairs

    def get_aligned_pairs(self):
        return self.aligned_pairs


class TestAlignmentInfo:

    def test_basic(self):
        alignment = Alignment('query', TUPLES2, SEQ1, 75896, 113577, 'reference', 789, PAIRS)
        alignment_info = AlignmentInfo(alignment)

        assert alignment_info.alignment == alignment
        # TODO: give magic number a name
        assert alignment_info.region_to_check == 12
        # TODO: add test for common.get_read_blocks()
        # TODO: change the testdata in order to get read_start and _end

    def test_alignment_pairs(self):
        alignment = Alignment('query', TUPLES1, SEQ1, 75896, 113577, 'reference', 789, PAIRS)
        alignment_info = AlignmentInfo(alignment)
        alignment_info.set_aligned_pairs()

        assert len(PAIRS) in alignment_info.aligned_pairs_start_index
        assert len(alignment_info.aligned_pairs_start_index) == 3

    # TODO: parametrize to cover all conditions branches
    def test_get_error_count(self):
        alignment = Alignment('query', TUPLES1, SEQ1, 75896, 113577, 'reference', 789, PAIRS)
        alignment_info = AlignmentInfo(alignment)

        indel_count, mismatch_count = alignment_info.get_error_count(75999, 81000)
        assert indel_count == 0
        assert mismatch_count == 0

    # TODO: parametrize to cover all conditions branches
    def test_add_polya_info(self):
        alignment = Alignment('query', TUPLES1, SEQ1, 51054, 63535, 'reference', 789, PAIRS)
        alignment_info = AlignmentInfo(alignment)
        polya_finder = PolyAFinder()
        polya_fixer = PolyAFixer(35)

        alignment_info.add_polya_info(polya_finder, polya_fixer)
        # TODO: assert internal fields: internal_polya_pos, external_polya_pos, read_exons, read_blocks, cigar_blocks, exons_changed, read_start, read_end

    def test_get_seq(self):
        alignment = Alignment('query', TUPLES1, SEQ1, 75896, 113577, 'reference', 789, PAIRS)
        alignment_info = AlignmentInfo(alignment)

        with pytest.raises(AssertionError):
            alignment_info.get_seq(76000, 81000)
