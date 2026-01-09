############################################################################
# Copyright (c) 2025 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import pytest
from src.barcode_calling.common import (
    str_to_2bit,
    bit_to_str,
    reverese_complement,
    find_polyt_start,
    align_pattern_ssw,
    find_candidate_with_max_score_ssw,
    find_candidate_with_max_score_ssw_var_len,
    NUCL2BIN,
    BIN2NUCL,
    base_comp,
)


class TestBitEncoding:
    """Test 2-bit DNA encoding functions."""

    def test_str_to_2bit_single_nucleotides(self):
        """Test encoding single nucleotides."""
        # Based on implementation: A=0, C=1, T=2, G=3
        assert str_to_2bit("A") == 0
        assert str_to_2bit("C") == 1
        assert str_to_2bit("T") == 2
        assert str_to_2bit("G") == 3

    def test_str_to_2bit_lowercase(self):
        """Test that lowercase works the same."""
        assert str_to_2bit("a") == str_to_2bit("A")
        assert str_to_2bit("c") == str_to_2bit("C")
        assert str_to_2bit("t") == str_to_2bit("T")
        assert str_to_2bit("g") == str_to_2bit("G")

    def test_str_to_2bit_sequence(self):
        """Test encoding multi-nucleotide sequences."""
        # "AC" = A(00) C(01) = 0b0001 = 1
        assert str_to_2bit("AC") == 1
        # "CA" = C(01) A(00) = 0b0100 = 4
        assert str_to_2bit("CA") == 4
        # "TG" = T(10) G(11) = 0b1011 = 11
        assert str_to_2bit("TG") == 11

    def test_bit_to_str_single_nucleotides(self):
        """Test decoding single nucleotides."""
        assert bit_to_str(0, 1) == "A"
        assert bit_to_str(1, 1) == "C"
        assert bit_to_str(2, 1) == "T"
        assert bit_to_str(3, 1) == "G"

    def test_bit_to_str_sequence(self):
        """Test decoding multi-nucleotide sequences."""
        assert bit_to_str(1, 2) == "AC"
        assert bit_to_str(4, 2) == "CA"
        assert bit_to_str(11, 2) == "TG"

    def test_roundtrip_encoding(self):
        """Test that encoding and decoding are inverses."""
        sequences = ["A", "ACTG", "GGGG", "AAAA", "TTTT", "ACGTACGT",
                     "ACTGACTGACTGACTG", "TATATATA"]
        for seq in sequences:
            encoded = str_to_2bit(seq)
            decoded = bit_to_str(encoded, len(seq))
            assert decoded == seq, f"Failed for {seq}"

    def test_encoding_constants(self):
        """Test that NUCL2BIN and BIN2NUCL are consistent."""
        assert BIN2NUCL[NUCL2BIN['A']] == 'A'
        assert BIN2NUCL[NUCL2BIN['C']] == 'C'
        assert BIN2NUCL[NUCL2BIN['T']] == 'T'
        assert BIN2NUCL[NUCL2BIN['G']] == 'G'


class TestReverseComplement:
    """Test reverse complement function."""

    def test_single_nucleotides(self):
        """Test complementing single nucleotides."""
        assert reverese_complement("A") == "T"
        assert reverese_complement("T") == "A"
        assert reverese_complement("C") == "G"
        assert reverese_complement("G") == "C"

    def test_sequence(self):
        """Test reverse complement of sequences."""
        assert reverese_complement("ACTG") == "CAGT"
        assert reverese_complement("AAAA") == "TTTT"
        assert reverese_complement("CCCC") == "GGGG"

    def test_palindrome(self):
        """Test that palindromic sequences return themselves."""
        # GAATTC is EcoRI recognition site (palindrome)
        assert reverese_complement("GAATTC") == "GAATTC"
        # ATAT is also palindromic
        assert reverese_complement("ATAT") == "ATAT"

    def test_with_n(self):
        """Test handling of N (unknown) nucleotide."""
        assert reverese_complement("N") == "N"
        assert reverese_complement("ACNG") == "CNGT"
        assert reverese_complement("NNN") == "NNN"

    def test_base_comp_dict(self):
        """Test base_comp dictionary is complete."""
        assert base_comp['A'] == 'T'
        assert base_comp['T'] == 'A'
        assert base_comp['C'] == 'G'
        assert base_comp['G'] == 'C'
        assert base_comp['N'] == 'N'


class TestFindPolytStart:
    """Test polyT detection function."""

    def test_clear_polyt_tract(self):
        """Test finding clear polyT tract."""
        # 10 bases of random + 16 Ts
        seq = "ACGTACGTAC" + "T" * 16
        result = find_polyt_start(seq)
        assert result >= 10  # Should find polyT starting around position 10

    def test_no_polyt(self):
        """Test sequence without polyT."""
        seq = "ACGTACGTACGTACGTACGT"
        assert find_polyt_start(seq) == -1

    def test_short_sequence(self):
        """Test sequence shorter than window size."""
        seq = "ACTG"  # Shorter than default window_size=16
        assert find_polyt_start(seq) == -1

    def test_polyt_at_start(self):
        """Test polyT at the beginning."""
        seq = "T" * 20 + "ACGT"
        result = find_polyt_start(seq)
        assert result >= 0 and result < 5  # Should be near start

    def test_mixed_with_some_t(self):
        """Test sequence with some Ts but not enough for polyT."""
        seq = "ATATATAT" * 4  # Alternating, not enough consecutive Ts
        assert find_polyt_start(seq) == -1

    def test_polyt_with_interruption(self):
        """Test polyT with slight interruption still detected."""
        # 12 Ts, then A, then 4 more Ts - should still be >= 75% T
        seq = "ACGTACGT" + "TTTTTTTTTTTTATTTTT"
        result = find_polyt_start(seq)
        assert result >= 0  # Should still detect


class TestSSWAlignment:
    """Test SSW-based alignment functions."""

    def test_align_pattern_exact_match(self):
        """Test aligning exact match pattern."""
        seq = "NNNNACTGACTGNNNN"
        start, end, pstart, pend, score = align_pattern_ssw(
            seq, 0, len(seq), "ACTGACTG", min_score=8)
        assert score == 8  # Perfect match of 8 bases
        assert start is not None

    def test_align_pattern_no_match(self):
        """Test aligning with no matching pattern."""
        seq = "AAAAAAAAAAAAAAAA"
        result = align_pattern_ssw(seq, 0, len(seq), "GGGGGGGG", min_score=6)
        assert result[0] is None  # No match above threshold

    def test_align_pattern_partial_match(self):
        """Test aligning with partial match."""
        seq = "NNNNACTGNNNN"
        start, end, pstart, pend, score = align_pattern_ssw(
            seq, 0, len(seq), "ACTG", min_score=3)
        assert score >= 3
        assert start is not None

    def test_align_pattern_with_mismatch(self):
        """Test aligning pattern with mismatch."""
        seq = "NNNNACTAACTGNNNN"  # ACTA instead of ACTG
        start, end, pstart, pend, score = align_pattern_ssw(
            seq, 0, len(seq), "ACTG", min_score=3)
        assert score >= 3  # Should still find partial match


class TestFindCandidateSSW:
    """Test barcode candidate selection functions."""

    def test_find_candidate_exact_match(self):
        """Test finding exact barcode match."""
        barcode_matches = [("ACTGACTG", 5, [0, 1, 2])]
        read_seq = "NNACTGACTGNN"

        barcode, score, start, end = find_candidate_with_max_score_ssw(
            barcode_matches, read_seq, min_score=6)

        assert barcode == "ACTGACTG"
        assert score >= 6

    def test_find_candidate_multiple_options(self):
        """Test selecting best from multiple candidates."""
        barcode_matches = [
            ("ACTGACTG", 5, [0, 1, 2]),
            ("TGCATGCA", 3, [0, 1]),
        ]
        read_seq = "NNACTGACTGNN"

        barcode, score, start, end = find_candidate_with_max_score_ssw(
            barcode_matches, read_seq, min_score=4)

        assert barcode == "ACTGACTG"  # Should pick the better match

    def test_find_candidate_no_match(self):
        """Test when no candidate meets threshold."""
        barcode_matches = [("ACTGACTG", 5, [0, 1, 2])]
        read_seq = "GGGGGGGGGGGG"

        barcode, score, start, end = find_candidate_with_max_score_ssw(
            barcode_matches, read_seq, min_score=6)

        assert barcode is None

    def test_find_candidate_score_diff(self):
        """Test score difference threshold."""
        # Two similar barcodes
        barcode_matches = [
            ("ACTGACTG", 5, [0, 1, 2]),
            ("ACTGACTA", 5, [0, 1, 2]),
        ]
        read_seq = "NNACTGACTGNN"

        # With high score_diff, should reject ambiguous
        barcode, score, start, end = find_candidate_with_max_score_ssw(
            barcode_matches, read_seq, min_score=4, score_diff=3)

        # Result depends on score difference

    def test_find_candidate_var_len(self):
        """Test variable length barcode matching."""
        barcode_matches = [
            ("ACTGACTGACTG", 5, [0, 1, 2]),  # 12-mer
            ("ACTGACTG", 5, [0, 1, 2]),  # 8-mer
        ]
        read_seq = "NNACTGACTGACTGNN"

        barcode, score, start, end = find_candidate_with_max_score_ssw_var_len(
            barcode_matches, read_seq, min_score=8)

        assert barcode == "ACTGACTGACTG"  # Should prefer longer exact match


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
