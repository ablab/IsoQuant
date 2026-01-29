############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Tests for universal barcode extraction.

Tests the universal barcode calling infrastructure:
- DetectedElement container class
- ExtractionResult dict-based detection result
- UniversalSingleMoleculeExtractor for pattern detection
- BarcodeResult protocol compliance
"""

import os
import pytest

from src.barcode_calling.callers.extraction_result import (
    DetectedElement,
    ExtractionResult,
    ReadStats
)
from src.barcode_calling.callers.molecule_structure import (
    ElementType,
    MoleculeElement,
    MoleculeStructure
)
from src.barcode_calling.callers.universal_extraction import (
    UniversalSingleMoleculeExtractor
)
from src.barcode_calling.callers.protocol import BarcodeResult

# Test data directory
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "universal_data")


# =============================================================================
# Test Fixtures
# =============================================================================

def create_simple_structure():
    """Create simple molecule structure: Barcode:UMI:PolyT:cDNA"""
    mdf_content = """Barcode:UMI:PolyT:cDNA
Barcode\tVAR_LIST\tATCCTTAGTGTTTGTC,CTTAGGAGTGGTTAGC,TGAGGGCCAACGTGCT
UMI\tVAR_ANY\t12
"""
    return MoleculeStructure(iter(mdf_content.strip().split('\n')))


def create_10x_structure():
    """Create 10x-like molecule structure with R1 adapter and TSO.

    Uses inline barcodes instead of file reference to avoid path issues.
    """
    # Use inline barcodes to avoid relative path issues with VAR_FILE
    mdf_content = """R1:Barcode:UMI:PolyT:cDNA:TSO
R1\tCONST\tCTACACGACGCTCTTCCGATCT
Barcode\tVAR_LIST\tATCCTTAGTGTTTGTC,CTTAGGAGTGGTTAGC,TGAGGGCCAACGTGCT,AAACCCGGGTTTAAAC
UMI\tVAR_ANY\t12
TSO\tCONST\tCCCATGTACTCTGCGTTGATACCACTGCTT
"""
    return MoleculeStructure(iter(mdf_content.strip().split('\n')))


def create_detected_elements():
    """Create sample detected elements for testing."""
    return {
        "Barcode": DetectedElement(22, 37, 16, "ATCCTTAGTGTTTGTC"),
        "UMI": DetectedElement(38, 49, 0, "AACCTTAGACCG"),
        "PolyT": DetectedElement(50, 79, 0)
    }


# =============================================================================
# Test DetectedElement
# =============================================================================

class TestDetectedElement:
    """Test DetectedElement container class."""

    def test_default_init(self):
        """Test initialization with default values."""
        element = DetectedElement()

        assert element.start == -1
        assert element.end == -1
        assert element.score == -1
        assert element.seq is None

    def test_custom_init(self):
        """Test initialization with custom values."""
        element = DetectedElement(start=10, end=25, score=15, seq="ACTGACTG")

        assert element.start == 10
        assert element.end == 25
        assert element.score == 15
        assert element.seq == "ACTGACTG"

    def test_partial_init(self):
        """Test initialization with some values."""
        element = DetectedElement(start=100, end=150)

        assert element.start == 100
        assert element.end == 150
        assert element.score == -1
        assert element.seq is None

    def test_coordinates_only(self):
        """Test element with only coordinates (no sequence)."""
        element = DetectedElement(50, 70, 0)

        assert element.start == 50
        assert element.end == 70
        assert element.seq is None


# =============================================================================
# Test ExtractionResult
# =============================================================================

class TestExtractionResult:
    """Test ExtractionResult dict-based detection result."""

    def test_init(self):
        """Test initialization."""
        structure = create_simple_structure()
        detected = create_detected_elements()

        result = ExtractionResult(structure, "read_001", "+", detected)

        assert result.read_id == "read_001"
        assert result.strand == "+"
        assert result.molecule_structure is structure
        assert "Barcode" in result.detected_results

    def test_barcode_identification_by_prefix(self):
        """Test barcode elements are identified by 'barcode' prefix."""
        structure = create_simple_structure()
        detected = create_detected_elements()

        result = ExtractionResult(structure, "read_001", "+", detected)

        # Should identify 'Barcode' element as barcode (case-insensitive prefix)
        assert "Barcode" in result._barcode_elements

    def test_umi_identification_by_prefix(self):
        """Test UMI elements are identified by 'umi' prefix."""
        structure = create_simple_structure()
        detected = create_detected_elements()

        result = ExtractionResult(structure, "read_001", "+", detected)

        # Should identify 'UMI' element as UMI (case-insensitive prefix)
        assert "UMI" in result._umi_elements

    def test_get_barcode_single(self):
        """Test getting single barcode."""
        structure = create_simple_structure()
        detected = create_detected_elements()

        result = ExtractionResult(structure, "read_001", "+", detected)

        assert result.get_barcode() == "ATCCTTAGTGTTTGTC"

    def test_get_barcode_missing(self):
        """Test getting barcode when not detected."""
        structure = create_simple_structure()
        detected = {
            "UMI": DetectedElement(38, 49, 0, "AACCTTAGACCG"),
            "PolyT": DetectedElement(50, 79, 0)
        }

        result = ExtractionResult(structure, "read_001", "+", detected)

        assert result.get_barcode() == ExtractionResult.NOSEQ

    def test_get_umi(self):
        """Test getting UMI."""
        structure = create_simple_structure()
        detected = create_detected_elements()

        result = ExtractionResult(structure, "read_001", "+", detected)

        assert result.get_umi() == "AACCTTAGACCG"

    def test_get_umi_missing(self):
        """Test getting UMI when not detected."""
        structure = create_simple_structure()
        detected = {
            "Barcode": DetectedElement(22, 37, 16, "ATCCTTAGTGTTTGTC"),
            "PolyT": DetectedElement(50, 79, 0)
        }

        result = ExtractionResult(structure, "read_001", "+", detected)

        assert result.get_umi() == ExtractionResult.NOSEQ

    def test_is_valid_with_barcode(self):
        """Test is_valid returns True when barcode detected."""
        structure = create_simple_structure()
        detected = create_detected_elements()

        result = ExtractionResult(structure, "read_001", "+", detected)

        assert result.is_valid() is True

    def test_is_valid_without_barcode(self):
        """Test is_valid returns False when no barcode."""
        structure = create_simple_structure()
        detected = {"PolyT": DetectedElement(50, 79, 0)}

        result = ExtractionResult(structure, "read_001", "+", detected)

        assert result.is_valid() is False

    def test_has_barcode(self):
        """Test has_barcode method."""
        structure = create_simple_structure()
        detected = create_detected_elements()

        result = ExtractionResult(structure, "read_001", "+", detected)

        assert result.has_barcode() is True

    def test_has_umi(self):
        """Test has_umi method."""
        structure = create_simple_structure()
        detected = create_detected_elements()

        result = ExtractionResult(structure, "read_001", "+", detected)

        assert result.has_umi() is True

    def test_has_umi_false(self):
        """Test has_umi returns False when no UMI."""
        structure = create_simple_structure()
        detected = {
            "Barcode": DetectedElement(22, 37, 16, "ATCCTTAGTGTTTGTC"),
            "PolyT": DetectedElement(50, 79, 0)
        }

        result = ExtractionResult(structure, "read_001", "+", detected)

        assert result.has_umi() is False

    def test_get_barcode_score(self):
        """Test getting barcode score."""
        structure = create_simple_structure()
        detected = create_detected_elements()

        result = ExtractionResult(structure, "read_001", "+", detected)

        assert result.get_barcode_score() == 16

    def test_get_barcode_score_no_barcode(self):
        """Test barcode score when no barcode detected."""
        structure = create_simple_structure()
        detected = {"PolyT": DetectedElement(50, 79, 0)}

        result = ExtractionResult(structure, "read_001", "+", detected)

        assert result.get_barcode_score() == -1

    def test_set_strand(self):
        """Test setting strand."""
        structure = create_simple_structure()
        detected = create_detected_elements()

        result = ExtractionResult(structure, "read_001", ".", detected)
        result.set_strand("-")

        assert result.strand == "-"

    def test_update_coordinates(self):
        """Test shifting all coordinates."""
        structure = create_simple_structure()
        detected = {
            "Barcode": DetectedElement(22, 37, 16, "ATCCTTAGTGTTTGTC"),
            "UMI": DetectedElement(38, 49, 0, "AACCTTAGACCG"),
        }

        result = ExtractionResult(structure, "read_001", "+", detected)
        result.update_coordinates(100)

        assert result.detected_results["Barcode"].start == 122
        assert result.detected_results["Barcode"].end == 137
        assert result.detected_results["UMI"].start == 138
        assert result.detected_results["UMI"].end == 149

    def test_update_coordinates_invalid(self):
        """Test coordinate shift preserves -1 values."""
        structure = create_simple_structure()
        detected = {
            "Barcode": DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ),
        }

        result = ExtractionResult(structure, "read_001", "+", detected)
        result.update_coordinates(100)

        # -1 values (invalid/not detected) should remain -1
        assert result.detected_results["Barcode"].start == -1
        assert result.detected_results["Barcode"].end == -1

    def test_more_informative_than_by_count(self):
        """Test comparison by detected element count."""
        structure = create_simple_structure()

        detected1 = {
            "Barcode": DetectedElement(22, 37, 16, "ATCCTTAGTGTTTGTC"),
            "UMI": DetectedElement(38, 49, 0, "AACCTTAGACCG"),
            "PolyT": DetectedElement(50, 79, 0),
        }
        detected2 = {
            "Barcode": DetectedElement(22, 37, 16, "ATCCTTAGTGTTTGTC"),
        }

        result1 = ExtractionResult(structure, "read_001", "+", detected1)
        result2 = ExtractionResult(structure, "read_001", "+", detected2)

        assert result1.more_informative_than(result2) is True
        assert result2.more_informative_than(result1) is False

    def test_more_informative_than_by_score(self):
        """Test comparison by total score when counts equal."""
        structure = create_simple_structure()

        detected1 = {
            "Barcode": DetectedElement(22, 37, 20, "ATCCTTAGTGTTTGTC"),
        }
        detected2 = {
            "Barcode": DetectedElement(22, 37, 10, "ATCCTTAGTGTTTGTC"),
        }

        result1 = ExtractionResult(structure, "read_001", "+", detected1)
        result2 = ExtractionResult(structure, "read_001", "+", detected2)

        assert result1.more_informative_than(result2) is True

    def test_str_format(self):
        """Test string formatting for TSV output."""
        structure = create_simple_structure()
        detected = create_detected_elements()

        result = ExtractionResult(structure, "read_001", "+", detected)
        output = str(result)

        # Standard format: read_id, barcode, UMI, BC_score, valid_UMI, strand
        assert "read_001" in output
        assert "ATCCTTAGTGTTTGTC" in output
        assert "AACCTTAGACCG" in output
        assert "16" in output
        assert "True" in output
        assert "+" in output

    def test_header_format(self):
        """Test header generation."""
        structure = create_simple_structure()
        detected = create_detected_elements()

        result = ExtractionResult(structure, "read_001", "+", detected)
        header = result.header()

        assert "#read_id" in header
        assert "barcode" in header
        assert "UMI" in header
        assert "BC_score" in header
        assert "valid_UMI" in header
        assert "strand" in header

    def test_get_additional_attributes(self):
        """Test getting additional attributes."""
        structure = create_simple_structure()
        detected = create_detected_elements()

        result = ExtractionResult(structure, "read_001", "+", detected)
        attrs = result.get_additional_attributes()

        assert "Barcode detected" in attrs
        assert "UMI detected" in attrs
        assert "PolyT detected" in attrs

    def test_noseq_constant(self):
        """Test NOSEQ constant."""
        assert ExtractionResult.NOSEQ == "*"


class TestExtractionResultProtocolCompliance:
    """Test that ExtractionResult conforms to BarcodeResult protocol."""

    def test_isinstance_check(self):
        """Test that ExtractionResult is recognized as BarcodeResult."""
        structure = create_simple_structure()
        detected = create_detected_elements()

        result = ExtractionResult(structure, "read_001", "+", detected)

        # Protocol is runtime_checkable
        assert isinstance(result, BarcodeResult)

    def test_has_required_attributes(self):
        """Test that all protocol attributes exist."""
        structure = create_simple_structure()
        detected = create_detected_elements()

        result = ExtractionResult(structure, "read_001", "+", detected)

        # Required attributes
        assert hasattr(result, 'NOSEQ')
        assert hasattr(result, 'read_id')
        assert hasattr(result, 'strand')

    def test_has_required_methods(self):
        """Test that all protocol methods exist."""
        structure = create_simple_structure()
        detected = create_detected_elements()

        result = ExtractionResult(structure, "read_001", "+", detected)

        # Required methods
        assert callable(getattr(result, 'get_barcode', None))
        assert callable(getattr(result, 'get_umi', None))
        assert callable(getattr(result, 'is_valid', None))
        assert callable(getattr(result, 'has_barcode', None))
        assert callable(getattr(result, 'has_umi', None))
        assert callable(getattr(result, 'set_strand', None))
        assert callable(getattr(result, 'update_coordinates', None))
        assert callable(getattr(result, 'more_informative_than', None))
        assert callable(getattr(result, 'get_additional_attributes', None))
        assert callable(getattr(result, 'header', None))


# =============================================================================
# Test ExtractionResult ReadStats
# =============================================================================

class TestExtractionReadStats:
    """Test ReadStats class for ExtractionResult."""

    def test_init(self):
        """Test initialization."""
        stats = ReadStats()

        assert stats.read_count == 0
        assert stats.bc_count == 0
        assert stats.umi_count == 0
        assert len(stats.pattern_counts) == 0

    def test_add_read_with_barcode(self):
        """Test adding read with barcode detected."""
        stats = ReadStats()
        structure = create_simple_structure()
        detected = create_detected_elements()
        result = ExtractionResult(structure, "read_001", "+", detected)

        stats.add_read(result)

        assert stats.read_count == 1
        assert stats.bc_count == 1
        assert stats.umi_count == 1
        assert stats.pattern_counts["Barcode"] == 1
        assert stats.pattern_counts["UMI"] == 1
        assert stats.pattern_counts["PolyT"] == 1

    def test_add_read_without_barcode(self):
        """Test adding read without barcode."""
        stats = ReadStats()
        structure = create_simple_structure()
        detected = {"PolyT": DetectedElement(50, 79, 0)}
        result = ExtractionResult(structure, "read_001", "+", detected)

        stats.add_read(result)

        assert stats.read_count == 1
        assert stats.bc_count == 0
        assert stats.umi_count == 0
        assert stats.pattern_counts["PolyT"] == 1

    def test_add_multiple_reads(self):
        """Test adding multiple reads."""
        stats = ReadStats()
        structure = create_simple_structure()

        # Read 1: has barcode, UMI, polyT
        detected1 = create_detected_elements()
        result1 = ExtractionResult(structure, "read_001", "+", detected1)

        # Read 2: has only barcode
        detected2 = {"Barcode": DetectedElement(22, 37, 16, "ATCCTTAGTGTTTGTC")}
        result2 = ExtractionResult(structure, "read_002", "+", detected2)

        # Read 3: no barcode
        detected3 = {"PolyT": DetectedElement(50, 79, 0)}
        result3 = ExtractionResult(structure, "read_003", "+", detected3)

        stats.add_read(result1)
        stats.add_read(result2)
        stats.add_read(result3)

        assert stats.read_count == 3
        assert stats.bc_count == 2  # reads 1 and 2
        assert stats.umi_count == 1  # only read 1
        assert stats.pattern_counts["Barcode"] == 2
        assert stats.pattern_counts["PolyT"] == 2

    def test_str_format(self):
        """Test string formatting."""
        stats = ReadStats()
        structure = create_simple_structure()
        detected = create_detected_elements()
        result = ExtractionResult(structure, "read_001", "+", detected)

        stats.add_read(result)
        output = str(stats)

        assert "Total reads:\t1" in output
        assert "Barcode detected:\t1" in output
        assert "UMI detected:\t1" in output


# =============================================================================
# Test UniversalSingleMoleculeExtractor
# =============================================================================

class TestUniversalSingleMoleculeExtractor:
    """Test UniversalSingleMoleculeExtractor class."""

    def test_init_simple_structure(self):
        """Test initialization with simple structure."""
        structure = create_simple_structure()
        extractor = UniversalSingleMoleculeExtractor(structure)

        assert extractor.has_polyt is True
        assert extractor.has_cdna is True
        assert "UMI" in extractor.elements_to_extract

    def test_init_10x_structure(self):
        """Test initialization with 10x structure."""
        structure = create_10x_structure()
        extractor = UniversalSingleMoleculeExtractor(structure)

        assert extractor.has_polyt is True
        assert extractor.has_cdna is True
        assert "R1" in extractor.constant_elements_to_detect
        assert "TSO" in extractor.constant_elements_to_detect

    def test_result_type(self):
        """Test result_type returns ExtractionResult."""
        structure = create_simple_structure()
        extractor = UniversalSingleMoleculeExtractor(structure)

        assert extractor.result_type() is ExtractionResult

    def test_header(self):
        """Test header generation."""
        structure = create_simple_structure()
        extractor = UniversalSingleMoleculeExtractor(structure)

        header = extractor.header()

        assert "#read_id" in header
        assert "strand" in header


class TestPatternDetection:
    """Test actual barcode/UMI detection with synthetic sequences."""

    def test_find_barcode_umi_forward(self):
        """Test detection in forward strand read."""
        structure = create_10x_structure()
        extractor = UniversalSingleMoleculeExtractor(structure)

        # Test read with known structure (from test_reads.fq READ_1)
        # Structure: R1(22bp) + Barcode(16bp) + UMI(12bp) + PolyT + cDNA
        read_id = "READ_1"
        sequence = "TACACGACGCTCTTCCGATCTATCCTTAGTGTTTGTCAACCTTAGACCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGGGCCTCAGCGTAATCCTATTGTGCATG"

        result = extractor.find_barcode_umi(read_id, sequence)

        assert result.read_id == read_id
        assert result.strand == "+"
        # Should detect polyT
        assert "PolyT" in result.detected_results
        polyt = result.detected_results["PolyT"]
        assert polyt.start != -1

    def test_find_barcode_umi_reverse(self):
        """Test detection in reverse strand read."""
        structure = create_10x_structure()
        extractor = UniversalSingleMoleculeExtractor(structure)

        # Test read with TSO at beginning (from test_reads.fq READ_2)
        # This indicates reverse strand - TSO is at 3' end
        read_id = "READ_2"
        # Starts with TSO adapter, ends with R1
        sequence = "AAGCAGTGGTATCAACGCAGAGTACATGGGAAGCAAGGAGGCAAGGCCCGCGCCAAGGCCAAGTCGCGGTCTTCCCGGGCCCGGGCTACCTATTCCCGGTGGGGCGTGTGCACCGGCTGCTGCGCAAGGGCAACTACGCGGAGCGTGTGGGCGCCGGCGCGCCGGTATACATGGCGGCGGTGCTGGAGTACCTAACGGCCGAGATCCTGGAGCTGGCGGGCAACGCGGCCCGCGACAACAAGAAGACTCGCATCATCCCGCGCCACCTGCAGCTGGCCATCCGCAACGACGAGGAGCTCAACAAGCTGCTGGGCAAAGTGACGATCGCGCAGGGCGGCGTCCTGCCCAACATCCAGGCCGTGCTGCTGCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATCAAGGAGGCTCGCTAACCACTCCTAAGAGATCGGAAGAGCGTCGTGTA"

        result = extractor.find_barcode_umi(read_id, sequence)

        assert result.read_id == read_id
        # Should detect in reverse orientation (TSO at start means polyA/polyT would be internal)

    def test_no_match_random_sequence(self):
        """Test detection with random sequence (no patterns)."""
        structure = create_simple_structure()
        extractor = UniversalSingleMoleculeExtractor(structure)

        read_id = "RANDOM"
        sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" * 5

        result = extractor.find_barcode_umi(read_id, sequence)

        assert result.read_id == read_id
        # May or may not detect anything, but should not crash


class TestMultipleBarcodeNaming:
    """Test barcode/UMI identification with various naming conventions."""

    def test_barcode_prefix_lowercase(self):
        """Test 'barcode' prefix detection (lowercase)."""
        mdf_content = """barcode1:cDNA
barcode1\tVAR_LIST\tAAAA,CCCC,GGGG
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))
        detected = {"barcode1": DetectedElement(0, 3, 10, "AAAA")}

        result = ExtractionResult(structure, "read_001", "+", detected)

        assert result.get_barcode() == "AAAA"
        assert "barcode1" in result._barcode_elements

    def test_barcode_prefix_mixedcase(self):
        """Test 'Barcode' prefix detection (mixed case)."""
        mdf_content = """Barcode:cDNA
Barcode\tVAR_LIST\tAAAA,CCCC,GGGG
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))
        detected = {"Barcode": DetectedElement(0, 3, 10, "CCCC")}

        result = ExtractionResult(structure, "read_001", "+", detected)

        assert result.get_barcode() == "CCCC"

    def test_umi_prefix_lowercase(self):
        """Test 'umi' prefix detection (lowercase)."""
        mdf_content = """umi:cDNA
umi\tVAR_ANY\t10
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))
        detected = {"umi": DetectedElement(0, 9, 0, "ACGTACGTAC")}

        result = ExtractionResult(structure, "read_001", "+", detected)

        assert result.get_umi() == "ACGTACGTAC"
        assert "umi" in result._umi_elements

    def test_multiple_barcodes_concatenated(self):
        """Test multiple barcode elements are concatenated with '|'."""
        mdf_content = """Barcode1:Barcode2:cDNA
Barcode1\tVAR_LIST\tAAAA,CCCC
Barcode2\tVAR_LIST\tGGGG,TTTT
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))
        detected = {
            "Barcode1": DetectedElement(0, 3, 10, "AAAA"),
            "Barcode2": DetectedElement(4, 7, 10, "GGGG")
        }

        result = ExtractionResult(structure, "read_001", "+", detected)

        assert result.get_barcode() == "AAAA|GGGG"

    def test_multiple_umis_concatenated(self):
        """Test multiple UMI elements are concatenated with '|'."""
        mdf_content = """UMI1:UMI2:cDNA
UMI1\tVAR_ANY\t6
UMI2\tVAR_ANY\t6
"""
        structure = MoleculeStructure(iter(mdf_content.strip().split('\n')))
        detected = {
            "UMI1": DetectedElement(0, 5, 0, "AAAAAA"),
            "UMI2": DetectedElement(6, 11, 0, "CCCCCC")
        }

        result = ExtractionResult(structure, "read_001", "+", detected)

        assert result.get_umi() == "AAAAAA|CCCCCC"


class TestIntegration:
    """Integration tests with test data files."""

    def test_process_test_reads_file(self):
        """Test processing reads from test_reads.fq."""
        structure = create_10x_structure()
        extractor = UniversalSingleMoleculeExtractor(structure)

        # Load test reads
        test_file = os.path.join(TEST_DATA_DIR, "test_reads.fq")
        reads = []
        with open(test_file) as f:
            lines = f.readlines()
            for i in range(0, len(lines), 4):
                if i + 1 < len(lines):
                    read_id = lines[i].strip()[1:].split()[0]  # Remove @ and take first part
                    sequence = lines[i + 1].strip()
                    reads.append((read_id, sequence))

        # Process each read
        results = []
        for read_id, sequence in reads:
            result = extractor.find_barcode_umi(read_id, sequence)
            results.append(result)

        # Should process all reads without error
        assert len(results) == len(reads)

        # Check results are ExtractionResult instances
        for result in results:
            assert isinstance(result, ExtractionResult)
            assert result.read_id in [r[0] for r in reads]


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
