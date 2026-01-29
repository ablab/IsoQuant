############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import pytest
from src.barcode_calling.callers import (
    BarcodeDetectionResult,
    LinkerBarcodeDetectionResult,
    TSOBarcodeDetectionResult,
    TenXBarcodeDetectionResult,
    SplittingBarcodeDetectionResult,
    ReadStats,
    increase_if_valid
)


class TestIncreaseIfValid:
    """Test coordinate increment utility function."""

    def test_increment_valid_positive(self):
        """Test incrementing valid positive value."""
        assert increase_if_valid(10, 5) == 15

    def test_increment_valid_zero(self):
        """Test incrementing zero (treated as invalid)."""
        assert increase_if_valid(0, 5) == 0

    def test_increment_invalid_negative(self):
        """Test incrementing -1 (invalid marker)."""
        assert increase_if_valid(-1, 5) == -1

    def test_increment_none(self):
        """Test incrementing None."""
        assert increase_if_valid(None, 5) is None


class TestBarcodeDetectionResult:
    """Test base barcode detection result class."""

    def test_init_default(self):
        """Test initialization with defaults."""
        result = BarcodeDetectionResult("read_001")

        assert result.read_id == "read_001"
        assert result.barcode == BarcodeDetectionResult.NOSEQ
        assert result.UMI == BarcodeDetectionResult.NOSEQ
        assert result.BC_score == -1
        assert result.UMI_good is False
        assert result.strand == "."

    def test_init_with_values(self):
        """Test initialization with all values."""
        result = BarcodeDetectionResult(
            read_id="read_001",
            barcode="ACTGACTG",
            UMI="GGGGGGGG",
            BC_score=16,
            UMI_good=True,
            strand="+"
        )

        assert result.read_id == "read_001"
        assert result.barcode == "ACTGACTG"
        assert result.UMI == "GGGGGGGG"
        assert result.BC_score == 16
        assert result.UMI_good is True
        assert result.strand == "+"

    def test_is_valid_with_barcode(self):
        """Test validity check with detected barcode."""
        result = BarcodeDetectionResult("read_001", barcode="ACTG")
        assert result.is_valid() is True

    def test_is_valid_without_barcode(self):
        """Test validity check without barcode."""
        result = BarcodeDetectionResult("read_001")
        assert result.is_valid() is False

    def test_set_strand(self):
        """Test setting strand."""
        result = BarcodeDetectionResult("read_001")
        result.set_strand("+")
        assert result.strand == "+"

    def test_str_format(self):
        """Test string formatting."""
        result = BarcodeDetectionResult(
            "read_001", "ACTG", "GGGG", 15, True, "+"
        )
        output = str(result)

        assert "read_001" in output
        assert "ACTG" in output
        assert "GGGG" in output
        assert "15" in output
        assert "True" in output
        assert "+" in output

    def test_header(self):
        """Test TSV header."""
        header = BarcodeDetectionResult.header()
        assert "#read_id" in header
        assert "barcode" in header
        assert "UMI" in header
        assert "BC_score" in header


class TestCurioBarcodeDetectionResult:
    """Test double barcode detection result class."""

    def test_init_default(self):
        """Test initialization with defaults."""
        result = LinkerBarcodeDetectionResult("read_001")

        assert result.read_id == "read_001"
        assert result.polyT == -1
        assert result.primer == -1
        assert result.linker_start == -1
        assert result.linker_end == -1

    def test_init_with_positions(self):
        """Test initialization with all positions."""
        result = LinkerBarcodeDetectionResult(
            read_id="read_001",
            barcode="ACTGACTG",
            UMI="GGGG",
            BC_score=14,
            UMI_good=True,
            strand="+",
            polyT=100,
            primer=50,
            linker_start=60,
            linker_end=75
        )

        assert result.polyT == 100
        assert result.primer == 50
        assert result.linker_start == 60
        assert result.linker_end == 75

    def test_update_coordinates(self):
        """Test coordinate shifting."""
        result = LinkerBarcodeDetectionResult(
            "read_001",
            polyT=100,
            primer=50,
            linker_start=60,
            linker_end=75
        )

        result.update_coordinates(10)

        assert result.polyT == 110
        assert result.primer == 60
        assert result.linker_start == 70
        assert result.linker_end == 85

    def test_update_coordinates_invalid(self):
        """Test coordinate shifting with invalid values."""
        result = LinkerBarcodeDetectionResult("read_001")
        result.update_coordinates(10)

        # Invalid values (-1) should remain unchanged
        assert result.polyT == -1
        assert result.primer == -1

    def test_more_informative_than_by_score(self):
        """Test comparison by barcode score."""
        result1 = LinkerBarcodeDetectionResult("read_001", BC_score=16)
        result2 = LinkerBarcodeDetectionResult("read_001", BC_score=14)

        assert result1.more_informative_than(result2) is True
        assert result2.more_informative_than(result1) is False

    def test_more_informative_than_by_linker(self):
        """Test comparison by linker position (tie-breaker)."""
        result1 = LinkerBarcodeDetectionResult(
            "read_001", BC_score=16, linker_start=80
        )
        result2 = LinkerBarcodeDetectionResult(
            "read_001", BC_score=16, linker_start=60
        )

        # Higher linker position is more informative
        assert result1.more_informative_than(result2) is True

    def test_get_additional_attributes_all(self):
        """Test attribute detection with all features."""
        result = LinkerBarcodeDetectionResult(
            "read_001",
            polyT=100,
            primer=50,
            linker_start=60
        )

        attrs = result.get_additional_attributes()

        assert "PolyT detected" in attrs
        assert "Primer detected" in attrs
        assert "Linker detected" in attrs

    def test_get_additional_attributes_partial(self):
        """Test attribute detection with some features."""
        result = LinkerBarcodeDetectionResult(
            "read_001",
            polyT=100,
            # No primer or linker
        )

        attrs = result.get_additional_attributes()

        assert "PolyT detected" in attrs
        assert "Primer detected" not in attrs
        assert "Linker detected" not in attrs

    def test_str_format(self):
        """Test string formatting includes positions."""
        result = LinkerBarcodeDetectionResult(
            "read_001", "ACTG", "GGGG", 16, True, "+",
            polyT=100, primer=50, linker_start=60, linker_end=75
        )
        output = str(result)

        assert "100" in output  # polyT
        assert "50" in output   # primer
        assert "60" in output   # linker_start
        assert "75" in output   # linker_end


class TestStereoBarcodeDetectionResult:
    """Test Stereo-seq detection result class."""

    def test_init_with_tso(self):
        """Test initialization with TSO position."""
        result = TSOBarcodeDetectionResult(
            "read_001",
            polyT=100,
            tso=150
        )

        assert result.tso5 == 150

    def test_update_coordinates_includes_tso(self):
        """Test coordinate shifting includes TSO."""
        result = TSOBarcodeDetectionResult(
            "read_001",
            polyT=100,
            tso=150
        )

        result.update_coordinates(10)

        assert result.tso5 == 160
        assert result.polyT == 110

    def test_get_additional_attributes_with_tso(self):
        """Test attribute detection includes TSO."""
        result = TSOBarcodeDetectionResult(
            "read_001",
            polyT=100,
            tso=150
        )

        attrs = result.get_additional_attributes()

        assert "PolyT detected" in attrs
        assert "TSO detected" in attrs


class TestTenXBarcodeDetectionResult:
    """Test 10x Genomics detection result class."""

    def test_init_with_r1(self):
        """Test initialization with R1 position."""
        result = TenXBarcodeDetectionResult(
            "read_001",
            polyT=100,
            r1=20
        )

        assert result.r1 == 20
        assert result.polyT == 100

    def test_update_coordinates(self):
        """Test coordinate shifting."""
        result = TenXBarcodeDetectionResult(
            "read_001",
            polyT=100,
            r1=20
        )

        result.update_coordinates(10)

        assert result.r1 == 30
        assert result.polyT == 110

    def test_more_informative_than_by_polyt(self):
        """Test comparison prioritizes polyT."""
        result1 = TenXBarcodeDetectionResult("read_001", polyT=100)
        result2 = TenXBarcodeDetectionResult("read_001", polyT=80)

        assert result1.more_informative_than(result2) is True

    def test_more_informative_than_by_r1(self):
        """Test comparison uses R1 as tie-breaker."""
        result1 = TenXBarcodeDetectionResult("read_001", polyT=100, r1=30)
        result2 = TenXBarcodeDetectionResult("read_001", polyT=100, r1=20)

        assert result1.more_informative_than(result2) is True

    def test_get_additional_attributes(self):
        """Test attribute detection."""
        result = TenXBarcodeDetectionResult(
            "read_001",
            polyT=100,
            r1=20
        )

        attrs = result.get_additional_attributes()

        assert "PolyT detected" in attrs
        assert "R1 detected" in attrs


class TestSplittingBarcodeDetectionResult:
    """Test splitting barcode detection result class."""

    def test_init(self):
        """Test initialization."""
        result = SplittingBarcodeDetectionResult("read_001")

        assert result.read_id == "read_001"
        assert result.detected_patterns == []

    def test_append(self):
        """Test appending detection patterns."""
        result = SplittingBarcodeDetectionResult("read_001")

        pattern1 = TSOBarcodeDetectionResult("read_001", barcode="ACTG")
        pattern2 = TSOBarcodeDetectionResult("read_001", barcode="TGCA")

        result.append(pattern1)
        result.append(pattern2)

        assert len(result.detected_patterns) == 2

    def test_empty_true(self):
        """Test empty detection."""
        result = SplittingBarcodeDetectionResult("read_001")
        assert result.empty() is True

    def test_empty_false(self):
        """Test non-empty detection."""
        result = SplittingBarcodeDetectionResult("read_001")
        pattern = TSOBarcodeDetectionResult("read_001", barcode="ACTG")
        result.append(pattern)

        assert result.empty() is False

    def test_filter_keeps_barcoded(self):
        """Test filter keeps results with barcodes."""
        result = SplittingBarcodeDetectionResult("read_001")

        barcoded = TSOBarcodeDetectionResult("read_001", barcode="ACTG")
        unbarcoded = TSOBarcodeDetectionResult("read_001")  # No barcode

        result.append(barcoded)
        result.append(unbarcoded)
        result.filter()

        # Should keep only barcoded result
        assert len(result.detected_patterns) == 1
        assert result.detected_patterns[0].barcode == "ACTG"

    def test_filter_keeps_first_if_none_barcoded(self):
        """Test filter keeps first result if none have barcodes."""
        result = SplittingBarcodeDetectionResult("read_001")

        unbarcoded1 = TSOBarcodeDetectionResult("read_001")
        unbarcoded2 = TSOBarcodeDetectionResult("read_001")

        result.append(unbarcoded1)
        result.append(unbarcoded2)
        result.filter()

        # Should keep only first result
        assert len(result.detected_patterns) == 1


class TestReadStats:
    """Test statistics tracker class."""

    def test_init(self):
        """Test initialization."""
        stats = ReadStats()

        assert stats.read_count == 0
        assert stats.bc_count == 0
        assert stats.umi_count == 0
        assert len(stats.additional_attributes_counts) == 0

    def test_add_read_with_barcode(self):
        """Test adding read with valid barcode."""
        stats = ReadStats()
        result = LinkerBarcodeDetectionResult(
            "read_001",
            barcode="ACTG",
            UMI_good=True,
            polyT=100
        )

        stats.add_read(result)

        assert stats.read_count == 1
        assert stats.bc_count == 1
        assert stats.umi_count == 1
        assert stats.additional_attributes_counts["PolyT detected"] == 1

    def test_add_read_without_barcode(self):
        """Test adding read without barcode."""
        stats = ReadStats()
        result = LinkerBarcodeDetectionResult("read_001")  # No barcode

        stats.add_read(result)

        assert stats.read_count == 1
        assert stats.bc_count == 0
        assert stats.umi_count == 0

    def test_add_multiple_reads(self):
        """Test adding multiple reads."""
        stats = ReadStats()

        result1 = LinkerBarcodeDetectionResult("read_001", barcode="ACTG", polyT=100)
        result2 = LinkerBarcodeDetectionResult("read_002", barcode="TGCA", primer=50)
        result3 = LinkerBarcodeDetectionResult("read_003")  # No barcode

        stats.add_read(result1)
        stats.add_read(result2)
        stats.add_read(result3)

        assert stats.read_count == 3
        assert stats.bc_count == 2
        assert stats.additional_attributes_counts["PolyT detected"] == 1
        assert stats.additional_attributes_counts["Primer detected"] == 1

    def test_add_custom_stats(self):
        """Test adding custom statistics."""
        stats = ReadStats()

        stats.add_custom_stats("Custom feature", 5)
        stats.add_custom_stats("Custom feature", 3)

        assert stats.additional_attributes_counts["Custom feature"] == 8

    def test_str_format(self):
        """Test string formatting."""
        stats = ReadStats()
        result = LinkerBarcodeDetectionResult("read_001", barcode="ACTG", UMI_good=True)
        stats.add_read(result)

        output = str(stats)

        assert "Total reads\t1" in output
        assert "Barcode detected\t1" in output
        assert "Reliable UMI\t1" in output

    def test_iter(self):
        """Test iteration over statistics."""
        stats = ReadStats()
        result = LinkerBarcodeDetectionResult("read_001", barcode="ACTG", polyT=100)
        stats.add_read(result)

        lines = list(stats)

        assert "Total reads: 1" in lines
        assert "Barcode detected: 1" in lines
        assert "PolyT detected: 1" in lines


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
