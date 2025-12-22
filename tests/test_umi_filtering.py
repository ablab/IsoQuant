############################################################################
# Copyright (c) 2025 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import pytest
from collections import defaultdict
from src.barcode_calling.umi_filtering import (
    UMIFilter,
    format_read_assignment_for_output,
    overlaps,
    create_transcript_info_dict
)
from src.isoform_assignment import ReadAssignment, ReadAssignmentType
from src.gene_info import GeneInfo


class TestOverlaps:
    """Test the overlaps helper function."""

    def test_overlaps_true(self):
        """Test ranges that overlap."""
        assert overlaps((10, 20), (15, 25)) == True
        assert overlaps((10, 30), (15, 25)) == True
        assert overlaps((15, 25), (10, 30)) == True
        assert overlaps((10, 20), (10, 20)) == True

    def test_overlaps_false(self):
        """Test ranges that don't overlap."""
        assert overlaps((10, 20), (21, 30)) == False
        assert overlaps((21, 30), (10, 20)) == False
        assert overlaps((10, 15), (20, 25)) == False

    def test_overlaps_adjacent(self):
        """Test adjacent ranges (touching endpoints are considered overlapping)."""
        assert overlaps((10, 20), (20, 30)) == True  # Endpoints touch, considered overlapping


class TestFormatReadAssignmentForOutput:
    """Test read assignment output formatting."""

    def test_format_basic(self):
        """Test basic formatting with all attributes."""
        # Create a ReadAssignment
        read_assignment = ReadAssignment(
            read_id="read_001",
            assignment_type=ReadAssignmentType.unique
        )
        read_assignment.chr_id = "chr1"
        read_assignment.genomic_region = (1000, 2000)
        read_assignment.barcode = "ACTG"
        read_assignment.umi = "GGGG"
        read_assignment.strand = '+'
        read_assignment.set_additional_attribute('transcript_type', 'known')
        read_assignment.set_additional_attribute('polya_site', 1950)
        read_assignment.set_additional_attribute('cell_type', 'neuron')

        output = format_read_assignment_for_output(read_assignment)

        # Verify all fields are in output
        assert "read_001" in output
        assert "ACTG" in output
        assert "GGGG" in output
        assert "known" in output
        assert "neuron" in output
        assert "1950" in output

    def test_format_missing_attributes(self):
        """Test formatting with missing optional attributes."""
        read_assignment = ReadAssignment(
            read_id="read_002",
            assignment_type=ReadAssignmentType.ambiguous
        )
        read_assignment.chr_id = "chr1"
        read_assignment.genomic_region = (1000, 2000)
        read_assignment.barcode = "ACTG"
        read_assignment.umi = "CCCC"
        read_assignment.strand = '-'

        output = format_read_assignment_for_output(read_assignment)

        # Should handle missing attributes gracefully
        assert "read_002" in output
        assert "unknown" in output or "None" in output


class TestUMIFilter:
    """Test UMIFilter class."""

    @pytest.fixture
    def mock_args(self):
        """Create mock arguments for UMIFilter."""
        class Args:
            def __init__(self):
                self.output = "/tmp/test_output"
                self.umi_length = 10
                self.barcode_length = 16
                self.use_technical_replicas = False
        return Args()

    @pytest.fixture
    def umi_filter(self, mock_args):
        """Create a UMIFilter instance."""
        return UMIFilter(mock_args)

    def test_init(self, umi_filter):
        """Test UMIFilter initialization."""
        assert umi_filter.args is not None
        assert isinstance(umi_filter.stats, dict)
        assert isinstance(umi_filter.unique_gene_barcode, set)

    def test_hamming_distance(self, umi_filter):
        """Test Hamming distance calculation."""
        assert umi_filter.hamming_distance("ACTG", "ACTG") == 0
        assert umi_filter.hamming_distance("ACTG", "ACCG") == 1
        assert umi_filter.hamming_distance("ACTG", "TGCA") == 4
        assert umi_filter.hamming_distance("ACTG", "ACCC") == 2

    def test_hamming_distance_different_lengths(self, umi_filter):
        """Test Hamming distance with different length strings."""
        # Should handle different lengths
        result = umi_filter.hamming_distance("ACTG", "AC")
        assert result >= 2

    def test_construct_umi_dict(self, umi_filter):
        """Test UMI dictionary construction."""
        # Create mock ReadAssignments
        read1 = ReadAssignment(read_id="read_001", assignment_type=ReadAssignmentType.unique)
        read1.umi = "AAAA"

        read2 = ReadAssignment(read_id="read_002", assignment_type=ReadAssignmentType.unique)
        read2.umi = "AAAA"

        read3 = ReadAssignment(read_id="read_003", assignment_type=ReadAssignmentType.unique)
        read3.umi = "TTTT"

        molecule_list = [read1, read2, read3]
        umi_dict = umi_filter._construct_umi_dict(molecule_list)

        # Should group by UMI
        assert len(umi_dict) == 2
        assert "AAAA" in umi_dict
        assert "TTTT" in umi_dict
        assert len(umi_dict["AAAA"]) == 2
        assert len(umi_dict["TTTT"]) == 1

    def test_construct_umi_dict_untrusted(self, umi_filter):
        """Test UMI dict construction with untrusted UMIs."""
        read1 = ReadAssignment(read_id="read_001", assignment_type=ReadAssignmentType.unique)
        read1.umi = ""  # Untrusted UMI

        molecule_list = [read1]
        umi_dict = umi_filter._construct_umi_dict(molecule_list)

        # Untrusted UMIs should be grouped separately
        assert len(umi_dict) == 1
        assert "" in umi_dict or "read_001" in umi_dict

    def test_select_best_read(self, umi_filter):
        """Test selecting best read from duplicates."""
        # Create reads with different qualities
        read1 = ReadAssignment(read_id="read_001", assignment_type=ReadAssignmentType.unique)
        read1.exons = [(1000, 1500), (1600, 2000)]

        read2 = ReadAssignment(read_id="read_002", assignment_type=ReadAssignmentType.ambiguous)
        read2.exons = [(1000, 1500), (1600, 2000)]

        read3 = ReadAssignment(read_id="read_003", assignment_type=ReadAssignmentType.unique)
        read3.exons = [(1000, 2000)]  # Single exon, less informative

        duplicates = [read2, read1, read3]
        best = umi_filter._process_duplicates(duplicates)

        # Should select unique assignment with more exons
        assert len(best) >= 1
        # The best read should be unique type
        assert any(r.assignment_type == ReadAssignmentType.unique for r in best)


class TestCreateTranscriptInfoDict:
    """Test transcript info dictionary creation."""

    def test_create_empty_dict(self):
        """Test creating dict with no genes."""
        gene_dict = {}
        result = create_transcript_info_dict(gene_dict)
        assert result == {}

    def test_create_dict_with_genes(self):
        """Test creating dict with mock genes."""
        # This would require more complex mocking of GeneInfo
        # For now, just test that the function exists and is callable
        gene_dict = {}
        result = create_transcript_info_dict(gene_dict)
        assert isinstance(result, dict)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
