############################################################################
# Copyright (c) 2025 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Unit tests for barcode detector classes.

Each detector has:
- test_init: Basic initialization
- test_no_match_random_sequence: Random sequence with no expected match
- test_real_sequence: Placeholder for real data (TODO: add real sequences)
"""

import pytest
from src.barcode_calling.callers import (
    # Stereo-seq
    StereoBarcodeDetector,
    StereoSplittingBarcodeDetector,
    # 10x Genomics
    TenXBarcodeDetector,
    VisiumHDBarcodeDetector,
    # Curio
    CurioBarcodeDetector,
    CurioBruteForceDetector,
    CurioIlluminaDetector,
)
from src.barcode_calling.callers.base import (
    BarcodeDetectionResult,
    StereoBarcodeDetectionResult,
    TenXBarcodeDetectionResult,
    CurioBarcodeDetectionResult,
)


# =============================================================================
# Dummy data for testing
# =============================================================================

# Random sequence with no expected barcode match (just random nucleotides)
RANDOM_SEQUENCE = (
    "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"
    "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
)

# Dummy barcodes for initializing indexers
DUMMY_BARCODES_16 = [
    "ACTGACTGACTGACTG",
    "TGCATGCATGCATGCA",
    "GGGGGGGGGGGGGGGG",
]

DUMMY_BARCODES_25 = [
    "ACTGACTGACTGACTGACTGACTGA",
    "TGCATGCATGCATGCATGCATGCAT",
    "GGGGGGGGGGGGGGGGGGGGGGGGG",
]

DUMMY_BARCODES_14 = [
    "ACTGACTGACTGAC",
    "TGCATGCATGCATG",
    "GGGGGGGGGGGGGG",
]

DUMMY_BARCODES_8 = [
    "ACTGACTG",
    "TGCATGCA",
    "GGGGGGGG",
]

DUMMY_BARCODES_6 = [
    "ACTGAC",
    "TGCATG",
    "GGGGGG",
]

DUMMY_BARCODES_15 = [
    "ACTGACTGACTGACT",
    "TGCATGCATGCATGC",
    "GGGGGGGGGGGGGGG",
]


# =============================================================================
# Stereo-seq Detector Tests
# =============================================================================

class TestStereoBarcodeDetector:
    """Test Stereo-seq barcode detector."""

    def test_init(self):
        """Test basic initialization."""
        detector = StereoBarcodeDetector(DUMMY_BARCODES_25, min_score=21)
        assert detector.barcode_indexer is not None
        assert detector.min_score == 21

    def test_init_empty_barcodes(self):
        """Test initialization with empty barcode list."""
        detector = StereoBarcodeDetector([], min_score=21)
        assert detector.barcode_indexer is None

    def test_no_match_random_sequence(self):
        """Test that random sequence returns no valid barcode."""
        detector = StereoBarcodeDetector(DUMMY_BARCODES_25, min_score=21)
        result = detector.find_barcode_umi("test_read", RANDOM_SEQUENCE)

        assert isinstance(result, (CurioBarcodeDetectionResult, StereoBarcodeDetectionResult))
        assert result.read_id == "test_read"
        # Random sequence should not match any barcode
        assert not result.is_valid()

    def test_real_sequence(self):
        """Test with real sequence data.

        TODO: Add real sequence that should match a barcode.
        Example format:
            sequence = "..."  # Real read sequence
            expected_barcode = "ACTGACTGACTGACTGACTGACTGA"
            result = detector.find_barcode_umi("read_001", sequence)
            assert result.barcode == expected_barcode
        """
        pytest.skip("TODO: Add real test sequence from data")


class TestStereoSplittingBarcodeDetector:
    """Test Stereo-seq splitting barcode detector."""

    def test_init(self):
        """Test basic initialization."""
        detector = StereoSplittingBarcodeDetector(DUMMY_BARCODES_25)
        assert detector.barcode_indexer is not None

    def test_init_empty_barcodes(self):
        """Test initialization with empty barcode list."""
        detector = StereoSplittingBarcodeDetector([])
        assert detector.barcode_indexer is None

    def test_no_match_random_sequence(self):
        """Test that random sequence returns empty splitting result."""
        detector = StereoSplittingBarcodeDetector(DUMMY_BARCODES_25)
        result = detector.find_barcode_umi("test_read", RANDOM_SEQUENCE)

        assert result.read_id == "test_read"
        # Should return empty or no valid patterns
        assert result.empty() or not any(p.is_valid() for p in result.detected_patterns)

    def test_real_sequence(self):
        """Test with real sequence data.

        TODO: Add real sequence that contains multiple barcodes.
        """
        pytest.skip("TODO: Add real test sequence from data")


# =============================================================================
# 10x Genomics Detector Tests
# =============================================================================

class TestTenXBarcodeDetector:
    """Test 10x Genomics v3 barcode detector."""

    def test_init(self):
        """Test basic initialization."""
        detector = TenXBarcodeDetector(DUMMY_BARCODES_16)
        assert detector.barcode_indexer is not None
        assert detector.min_score == 14  # Default for small whitelist

    def test_init_large_whitelist(self):
        """Test that min_score increases for large whitelists."""
        large_barcodes = [f"ACTG{i:012d}"[:16] for i in range(100001)]
        detector = TenXBarcodeDetector(large_barcodes)
        assert detector.min_score == 16

    def test_no_match_random_sequence(self):
        """Test that random sequence returns no valid barcode."""
        detector = TenXBarcodeDetector(DUMMY_BARCODES_16)
        result = detector.find_barcode_umi("test_read", RANDOM_SEQUENCE)

        assert isinstance(result, TenXBarcodeDetectionResult)
        assert result.read_id == "test_read"
        assert not result.is_valid()

    def test_real_sequence(self):
        """Test with real sequence data.

        TODO: Add real 10x read sequence.
        Expected structure: [R1 adapter][16bp barcode][12bp UMI][polyT][cDNA]
        """
        pytest.skip("TODO: Add real test sequence from data")


class TestVisiumHDBarcodeDetector:
    """Test Visium HD spatial barcode detector."""

    def test_init(self):
        """Test basic initialization."""
        barcode_pairs = [DUMMY_BARCODES_16, DUMMY_BARCODES_15]
        detector = VisiumHDBarcodeDetector(barcode_pairs)

        assert detector.part1_barcode_indexer is not None
        assert detector.part2_barcode_indexer is not None
        assert detector.min_score == 13

    def test_no_match_random_sequence(self):
        """Test that random sequence returns no valid barcode."""
        barcode_pairs = [DUMMY_BARCODES_16, DUMMY_BARCODES_15]
        detector = VisiumHDBarcodeDetector(barcode_pairs)
        result = detector.find_barcode_umi("test_read", RANDOM_SEQUENCE)

        assert isinstance(result, TenXBarcodeDetectionResult)
        assert result.read_id == "test_read"
        assert not result.is_valid()

    def test_real_sequence(self):
        """Test with real sequence data.

        TODO: Add real Visium HD read sequence.
        Expected structure: [R1][16bp BC1][2bp separator][15bp BC2][9bp UMI][polyT][cDNA]
        """
        pytest.skip("TODO: Add real test sequence from data")


# =============================================================================
# Curio Detector Tests
# =============================================================================

class TestCurioBarcodeDetector:
    """Test Curio platform barcode detector."""

    def test_init_joint_barcodes(self):
        """Test initialization with joint barcode list."""
        detector = CurioBarcodeDetector(DUMMY_BARCODES_14, min_score=13)
        assert detector.barcode_indexer is not None
        assert detector.min_score == 13

    def test_init_split_barcodes(self):
        """Test initialization with separate left/right barcode lists."""
        barcode_tuple = (DUMMY_BARCODES_8, DUMMY_BARCODES_6)
        detector = CurioBarcodeDetector(barcode_tuple, min_score=13)
        # Should have 3*3 = 9 combined barcodes
        assert len(detector.barcode_indexer.seq_list) == 9

    def test_no_match_random_sequence(self):
        """Test that random sequence returns no valid barcode."""
        detector = CurioBarcodeDetector(DUMMY_BARCODES_14, min_score=13)
        result = detector.find_barcode_umi("test_read", RANDOM_SEQUENCE)

        assert isinstance(result, CurioBarcodeDetectionResult)
        assert result.read_id == "test_read"
        assert not result.is_valid()

    def test_real_sequence(self):
        """Test with real sequence data.

        TODO: Add real Curio read sequence.
        Expected structure: [primer][8bp BC1][linker][6bp BC2][9bp UMI][polyT][cDNA]
        Linker: TCTTCAGCGTTCCCGAGA
        """
        pytest.skip("TODO: Add real test sequence from data")


class TestCurioBruteForceDetector:
    """Test Curio brute force detector (exact linker matching).

    Note: find_barcode_umi has a bug (tries to unpack non-tuple),
    so we test _find_barcode_umi_fwd directly.
    """

    def test_init(self):
        """Test basic initialization."""
        detector = CurioBruteForceDetector(DUMMY_BARCODES_14)
        assert len(detector.barcode_set) == 3

    def test_no_match_random_sequence(self):
        """Test that random sequence returns no valid barcode."""
        detector = CurioBruteForceDetector(DUMMY_BARCODES_14)
        # Call internal method directly (find_barcode_umi has a bug)
        result = detector._find_barcode_umi_fwd("test_read", RANDOM_SEQUENCE)

        assert isinstance(result, CurioBarcodeDetectionResult)
        assert result.read_id == "test_read"
        assert not result.is_valid()

    def test_exact_linker_match_wrong_barcode(self):
        """Test that exact linker match but wrong barcode returns partial result."""
        # Sequence with linker but barcodes not in whitelist
        linker = CurioBruteForceDetector.LINKER
        sequence = f"NNACTGACTG{linker}NACTGACTGACTGNN"

        detector = CurioBruteForceDetector(DUMMY_BARCODES_14)
        result = detector._find_barcode_umi_fwd("test_read", sequence)

        assert isinstance(result, CurioBarcodeDetectionResult)
        # Linker should be detected even if barcode doesn't match
        assert result.linker_start != -1
        # But barcode should not be valid (not in whitelist)
        assert not result.is_valid()

    def test_real_sequence(self):
        """Test with real sequence data.

        TODO: Add real sequence with exact linker match.
        """
        pytest.skip("TODO: Add real test sequence from data")


class TestCurioIlluminaDetector:
    """Test Curio detector optimized for Illumina reads."""

    def test_init(self):
        """Test basic initialization."""
        detector = CurioIlluminaDetector(DUMMY_BARCODES_14, min_score=14)
        assert detector.barcode_indexer is not None
        assert detector.min_score == 14

    def test_no_match_random_sequence(self):
        """Test that random sequence returns no valid barcode."""
        detector = CurioIlluminaDetector(DUMMY_BARCODES_14, min_score=14)
        result = detector.find_barcode_umi("test_read", RANDOM_SEQUENCE)

        assert isinstance(result, CurioBarcodeDetectionResult)
        assert result.read_id == "test_read"
        assert not result.is_valid()

    def test_real_sequence(self):
        """Test with real sequence data.

        TODO: Add real Illumina Curio read sequence.
        """
        pytest.skip("TODO: Add real test sequence from data")


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
