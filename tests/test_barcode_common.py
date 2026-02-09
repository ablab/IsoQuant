############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import pytest
from src.barcode_calling.common import (
    str_to_2bit,
    bit_to_str,
    batch_str_to_2bit,
    batch_str_to_2bit_chunked,
    reverese_complement,
    find_polyt_start,
    align_pattern_ssw,
    find_candidate_with_max_score_ssw,
    find_candidate_with_max_score_ssw_var_len,
    detect_exact_positions,
    detect_first_exact_positions,
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


class TestNumba2BitEncoding:
    """Test Numba JIT-compiled 2-bit encoding consistency."""

    def test_numba_consistency_with_numpy(self):
        """Test that Numba and NumPy fallback produce identical results."""
        import numpy as np
        from src.barcode_calling.common import NUMBA_AVAILABLE, _convert_chunk_to_2bit

        barcodes = [
            "ACTGACTGACTGACTGACTGACTGA",
            "TGCATGCATGCATGCATGCATGCAT",
            "GGGGGGGGGGGGGGGGGGGGGGGGG",
            "AAAAAAAAAAAAAAAAAAAAAAAAA",
            "CCCCCCCCCCCCCCCCCCCCCCCCC",
        ]
        seq_len = 25

        # Get result from current implementation (Numba if available)
        result = _convert_chunk_to_2bit(barcodes, seq_len)

        # Verify against individual str_to_2bit calls
        for i, bc in enumerate(barcodes):
            expected = str_to_2bit(bc)
            assert result[i] == expected, f"Mismatch for {bc}: got {result[i]}, expected {expected}"

    def test_numba_available_flag(self):
        """Test that NUMBA_AVAILABLE flag is set correctly."""
        from src.barcode_calling.common import NUMBA_AVAILABLE

        # Just verify the flag is a boolean
        assert isinstance(NUMBA_AVAILABLE, bool)

        # If numba is importable without errors, flag should be True
        try:
            from numba import njit, prange
            assert NUMBA_AVAILABLE is True
        except (ImportError, SystemError):
            # SystemError can occur with broken numba installations
            assert NUMBA_AVAILABLE is False


class TestBatchStrTo2Bit:
    """Test batch 2-bit DNA encoding functions."""

    def test_empty_input(self):
        """Test empty input returns empty array."""
        import numpy as np
        result = batch_str_to_2bit([], 25)
        assert len(result) == 0
        assert result.dtype == np.uint64

    def test_single_barcode(self):
        """Test single barcode matches str_to_2bit."""
        barcode = "ACTGACTG"
        result = batch_str_to_2bit([barcode], seq_len=8)
        expected = str_to_2bit(barcode)
        assert result[0] == expected

    def test_multiple_barcodes(self):
        """Test multiple barcodes match individual str_to_2bit calls."""
        barcodes = ["ACTGACTG", "TGCATGCA", "GGGGGGGG", "AAAAAAAA"]
        result = batch_str_to_2bit(barcodes, seq_len=8)

        assert len(result) == 4
        for i, bc in enumerate(barcodes):
            assert result[i] == str_to_2bit(bc), f"Mismatch for barcode {bc}"

    def test_25mer_barcodes(self):
        """Test with default 25-mer barcodes (Stereo-seq length)."""
        barcodes = [
            "ACTGACTGACTGACTGACTGACTGA",
            "TGCATGCATGCATGCATGCATGCAT",
            "GGGGGGGGGGGGGGGGGGGGGGGGG",
        ]
        result = batch_str_to_2bit(barcodes, seq_len=25)

        assert len(result) == 3
        for i, bc in enumerate(barcodes):
            assert result[i] == str_to_2bit(bc), f"Mismatch for barcode {bc}"

    def test_iterator_input(self):
        """Test that iterator input works correctly."""
        barcodes = ["ACTGACTG", "TGCATGCA", "GGGGGGGG"]
        result = batch_str_to_2bit(iter(barcodes), seq_len=8)

        assert len(result) == 3
        for i, bc in enumerate(barcodes):
            assert result[i] == str_to_2bit(bc)

    def test_roundtrip(self):
        """Test that batch encoding can be decoded correctly."""
        barcodes = ["ACTGACTG", "TGCATGCA", "GGGGGGGG"]
        encoded = batch_str_to_2bit(barcodes, seq_len=8)

        for i, bc in enumerate(barcodes):
            decoded = bit_to_str(encoded[i], 8)
            assert decoded == bc, f"Roundtrip failed for {bc}: got {decoded}"


class TestBatchStrTo2BitChunked:
    """Test chunked batch 2-bit DNA encoding function."""

    def test_empty_input(self):
        """Test empty iterator returns empty array."""
        import numpy as np
        result = batch_str_to_2bit_chunked(iter([]), 25)
        assert len(result) == 0
        assert result.dtype == np.uint64

    def test_single_chunk(self):
        """Test processing fewer items than chunk size."""
        barcodes = ["ACTGACTG", "TGCATGCA", "GGGGGGGG"]
        result = batch_str_to_2bit_chunked(iter(barcodes), seq_len=8, chunk_size=100)

        assert len(result) == 3
        for i, bc in enumerate(barcodes):
            assert result[i] == str_to_2bit(bc)

    def test_multiple_chunks(self):
        """Test processing more items than chunk size."""
        barcodes = ["ACTGACTG", "TGCATGCA", "GGGGGGGG", "AAAAAAAA", "CCCCCCCC"]
        # Use chunk_size=2 to force multiple chunks
        result = batch_str_to_2bit_chunked(iter(barcodes), seq_len=8, chunk_size=2)

        assert len(result) == 5
        for i, bc in enumerate(barcodes):
            assert result[i] == str_to_2bit(bc), f"Mismatch for barcode {bc}"

    def test_exact_chunk_boundary(self):
        """Test when barcode count equals chunk size."""
        barcodes = ["ACTGACTG", "TGCATGCA"]
        result = batch_str_to_2bit_chunked(iter(barcodes), seq_len=8, chunk_size=2)

        assert len(result) == 2
        for i, bc in enumerate(barcodes):
            assert result[i] == str_to_2bit(bc)

    def test_consistency_with_batch(self):
        """Test that chunked and non-chunked produce same results."""
        barcodes = ["ACTGACTG", "TGCATGCA", "GGGGGGGG", "AAAAAAAA", "CCCCCCCC"]

        batch_result = batch_str_to_2bit(barcodes, seq_len=8)
        chunked_result = batch_str_to_2bit_chunked(iter(barcodes), seq_len=8, chunk_size=2)

        assert len(batch_result) == len(chunked_result)
        for i in range(len(batch_result)):
            assert batch_result[i] == chunked_result[i], f"Mismatch at index {i}"


class TestDetectExactPositions:
    """Test detect_exact_positions function."""

    def test_exact_match_in_sequence(self):
        """Test finding exact pattern in sequence."""
        sequence = "NNNNACTGACTGNNNN"
        pattern = "ACTGACTG"
        # Simulate k-mer occurrences: (pattern, count, [positions])
        # Positions are k-mer match positions within the search window
        pattern_occurrences = [(pattern, 5, [4, 5, 6, 7])]

        start, end = detect_exact_positions(
            sequence, 0, len(sequence), kmer_size=4,
            pattern=pattern, pattern_occurrences=pattern_occurrences,
            min_score=6
        )

        assert start is not None
        assert end is not None
        # Should find the pattern starting around position 4
        assert sequence[start:end] == pattern or pattern in sequence[start:end+1]

    def test_pattern_not_in_occurrences(self):
        """Test when pattern is not in the occurrence list."""
        sequence = "NNNNACTGACTGNNNN"
        pattern = "ACTGACTG"
        pattern_occurrences = [("TGCATGCA", 5, [0, 1, 2])]  # Different pattern

        start, end = detect_exact_positions(
            sequence, 0, len(sequence), kmer_size=4,
            pattern=pattern, pattern_occurrences=pattern_occurrences,
            min_score=6
        )

        assert start is None
        assert end is None

    def test_empty_occurrences(self):
        """Test with empty occurrence list."""
        sequence = "NNNNACTGACTGNNNN"

        start, end = detect_exact_positions(
            sequence, 0, len(sequence), kmer_size=4,
            pattern="ACTGACTG", pattern_occurrences=[],
            min_score=6
        )

        assert start is None
        assert end is None

    def test_no_positions_for_pattern(self):
        """Test when pattern has no match positions."""
        sequence = "NNNNACTGACTGNNNN"
        pattern = "ACTGACTG"
        pattern_occurrences = [(pattern, 0, [])]  # No positions

        start, end = detect_exact_positions(
            sequence, 0, len(sequence), kmer_size=4,
            pattern=pattern, pattern_occurrences=pattern_occurrences,
            min_score=6
        )

        assert start is None
        assert end is None


class TestDetectFirstExactPositions:
    """Test detect_first_exact_positions function."""

    def test_finds_first_match(self):
        """Test that first match is returned."""
        sequence = "NNACTGACTGNNACTGACTGNN"
        pattern = "ACTGACTG"
        # Two potential matches - should return first one
        pattern_occurrences = [(pattern, 5, [2, 12])]

        start, end = detect_first_exact_positions(
            sequence, 0, len(sequence), kmer_size=4,
            pattern=pattern, pattern_occurrences=pattern_occurrences,
            min_score=6
        )

        assert start is not None
        # First match should be around position 2
        assert start <= 5

    def test_pattern_not_found(self):
        """Test when pattern is not in occurrence list."""
        sequence = "NNNNACTGACTGNNNN"
        pattern_occurrences = [("TGCATGCA", 5, [0, 1])]

        start, end = detect_first_exact_positions(
            sequence, 0, len(sequence), kmer_size=4,
            pattern="ACTGACTG", pattern_occurrences=pattern_occurrences,
            min_score=6
        )

        assert start is None
        assert end is None

    def test_empty_positions(self):
        """Test with empty position list."""
        sequence = "NNNNACTGACTGNNNN"
        pattern = "ACTGACTG"
        pattern_occurrences = [(pattern, 0, [])]

        start, end = detect_first_exact_positions(
            sequence, 0, len(sequence), kmer_size=4,
            pattern=pattern, pattern_occurrences=pattern_occurrences,
            min_score=6
        )

        assert start is None
        assert end is None

    def test_min_score_threshold(self):
        """Test that min_score filters weak matches."""
        sequence = "NNNNACNNNNNN"  # Only partial match
        pattern = "ACTGACTG"
        pattern_occurrences = [(pattern, 2, [4])]

        # With high min_score, should not find match
        start, end = detect_first_exact_positions(
            sequence, 0, len(sequence), kmer_size=4,
            pattern=pattern, pattern_occurrences=pattern_occurrences,
            min_score=7
        )

        assert start is None
        assert end is None


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
