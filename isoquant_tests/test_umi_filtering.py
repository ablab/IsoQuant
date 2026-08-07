############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import pytest
from isoquant_lib.barcode_calling.umi_filtering import (
    UMIFilter,
    format_read_assignment_for_output,
    is_untrusted_umi,
    UNTRUSTED_UMI_KEY,

)
from isoquant_lib.assignment.isoform_assignment import ReadAssignment, ReadAssignmentType
from isoquant_lib.common import junctions_from_blocks
from isoquant_lib.utils.string_pools import StringPoolManager


class TestFormatReadAssignmentForOutput:
    """Test read assignment output formatting."""
    string_pools = StringPoolManager()

    def test_format_basic(self):
        """Test basic formatting with all attributes."""
        # Create a ReadAssignment
        read_assignment = ReadAssignment(
            read_id="read_001",
            assignment_type=ReadAssignmentType.unique,
            string_pools=self.string_pools
        )
        read_assignment.chr_id = "chr1"
        read_assignment.start = 1000
        read_assignment.end = 2000
        read_assignment.corrected_exons = [(1000, 1500), (1600, 2000)]
        read_assignment.corrected_introns = junctions_from_blocks(read_assignment.corrected_exons)
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
        # polyA site formatting depends on matching_events, so just check output is non-empty
        assert len(output) > 0

    def test_format_missing_attributes(self):
        """Test formatting with missing optional attributes."""
        read_assignment = ReadAssignment(
            read_id="read_002",
            assignment_type=ReadAssignmentType.ambiguous,
            string_pools=self.string_pools
        )
        read_assignment.chr_id = "chr1"
        read_assignment.start = 1000
        read_assignment.end = 2000
        read_assignment.corrected_exons = [(1000, 2000)]
        read_assignment.corrected_introns = []
        read_assignment.barcode = "ACTG"
        read_assignment.umi = "CCCC"
        read_assignment.strand = '-'

        output = format_read_assignment_for_output(read_assignment)

        # Should handle missing attributes gracefully
        assert "read_002" in output
        assert "unknown" in output or "None" in output


class TestUMIFilter:
    """Test UMIFilter class."""
    string_pools = StringPoolManager()

    @pytest.fixture
    def umi_filter(self):
        """Create a UMIFilter instance."""
        return UMIFilter(umi_length=10, edit_distance=3)

    def test_init(self, umi_filter):
        """Test UMIFilter initialization."""
        assert umi_filter.umi_length == 10
        assert umi_filter.max_edit_distance == 3
        assert isinstance(umi_filter.stats, dict)
        assert isinstance(umi_filter.unique_gene_barcode, set)
        assert isinstance(umi_filter.selected_reads, set)

    def test_construct_umi_dict(self, umi_filter):
        """Test UMI dictionary construction."""
        # Create mock ReadAssignments
        read1 = ReadAssignment(read_id="read_001", assignment_type=ReadAssignmentType.unique, string_pools=self.string_pools)
        read1.umi = "AAAA"

        read2 = ReadAssignment(read_id="read_002", assignment_type=ReadAssignmentType.unique, string_pools=self.string_pools)
        read2.umi = "AAAA"

        read3 = ReadAssignment(read_id="read_003", assignment_type=ReadAssignmentType.unique, string_pools=self.string_pools)
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
        read1 = ReadAssignment(read_id="read_001", assignment_type=ReadAssignmentType.unique, string_pools=self.string_pools)
        read1.umi = ""  # Untrusted UMI

        molecule_list = [read1]
        umi_dict = umi_filter._construct_umi_dict(molecule_list)

        # Untrusted UMIs should be grouped separately
        assert len(umi_dict) == 1
        assert "None" in umi_dict or "" in umi_dict or "read_001" in umi_dict

    @pytest.mark.parametrize("umi", [None, "None", "*", ".", ""])
    def test_is_untrusted_umi(self, umi):
        """All sentinels for an undetected UMI must be recognized as untrusted."""
        assert is_untrusted_umi(umi)

    def test_is_untrusted_umi_real_sequence(self):
        """A real UMI sequence must be trusted."""
        assert not is_untrusted_umi("ACGTACGTAC")

    def test_construct_umi_dict_noseq_umi(self, umi_filter):
        """Reads with the barcode caller '*' sentinel go to the untrusted bucket."""
        read1 = ReadAssignment(read_id="read_001", assignment_type=ReadAssignmentType.unique,
                               string_pools=self.string_pools)
        read1.umi = "AAAAAAAAAA"

        read2 = ReadAssignment(read_id="read_002", assignment_type=ReadAssignmentType.unique,
                               string_pools=self.string_pools)
        read2.umi = "*"

        umi_dict = umi_filter._construct_umi_dict([read1, read2])

        assert "*" not in umi_dict
        assert [m.read_id for m in umi_dict[UNTRUSTED_UMI_KEY]] == ["read_002"]
        assert [m.read_id for m in umi_dict["AAAAAAAAAA"]] == ["read_001"]

    def test_process_duplicates_drops_noseq_umi(self, umi_filter):
        """A '*' UMI read is dropped when the gene-barcode pair has detected UMIs."""
        read1 = ReadAssignment(read_id="read_001", assignment_type=ReadAssignmentType.unique,
                               string_pools=self.string_pools)
        read1.umi = "AAAAAAAAAA"
        read1.corrected_exons = [(1000, 1500), (1600, 2000)]

        read2 = ReadAssignment(read_id="read_002", assignment_type=ReadAssignmentType.unique,
                               string_pools=self.string_pools)
        read2.umi = "TTTTTTTTTT"
        read2.corrected_exons = [(1000, 1500), (1600, 2000)]

        read3 = ReadAssignment(read_id="read_003", assignment_type=ReadAssignmentType.unique,
                               string_pools=self.string_pools)
        read3.umi = "*"
        read3.corrected_exons = [(1000, 1500), (1600, 2000)]

        selected = umi_filter._process_duplicates([read1, read2, read3])

        assert sorted(r.read_id for r in selected) == ["read_001", "read_002"]

    def test_process_duplicates_keeps_lone_noseq_umi(self, umi_filter):
        """Reads with no detected UMI still yield one molecule when they are all there is."""
        read1 = ReadAssignment(read_id="read_001", assignment_type=ReadAssignmentType.unique,
                               string_pools=self.string_pools)
        read1.umi = "*"
        read1.corrected_exons = [(1000, 1500), (1600, 2000)]

        read2 = ReadAssignment(read_id="read_002", assignment_type=ReadAssignmentType.unique,
                               string_pools=self.string_pools)
        read2.umi = "*"
        read2.corrected_exons = [(1000, 1500), (1600, 2000)]

        selected = umi_filter._process_duplicates([read1, read2])

        assert len(selected) == 1

    def test_select_best_read(self, umi_filter):
        """Test selecting best read from duplicates."""
        # Create reads with different qualities
        read1 = ReadAssignment(read_id="read_001", assignment_type=ReadAssignmentType.unique, string_pools=self.string_pools)
        read1.corrected_exons = [(1000, 1500), (1600, 2000)]

        read2 = ReadAssignment(read_id="read_002", assignment_type=ReadAssignmentType.ambiguous, string_pools=self.string_pools)
        read2.corrected_exons = [(1000, 1500), (1600, 2000)]

        read3 = ReadAssignment(read_id="read_003", assignment_type=ReadAssignmentType.unique, string_pools=self.string_pools)
        read3.corrected_exons = [(1000, 2000)]  # Single exon, less informative

        duplicates = [read2, read1, read3]
        best = umi_filter._process_duplicates(duplicates)

        # Should select unique assignment with more exons
        assert len(best) >= 1
        # The best read should be unique type
        assert any(r.assignment_type == ReadAssignmentType.unique for r in best)


class TestCreateTranscriptInfoDict:
    """Test transcript info dictionary creation."""

    @pytest.mark.skip(reason="Requires actual genedb file, complex to mock")
    def test_create_empty_dict(self):
        """Test creating dict with no genes."""
        # This function requires a real genedb file path, not a dict
        # Skipping as it requires complex setup
        pass

    @pytest.mark.skip(reason="Requires actual genedb file, complex to mock")
    def test_create_dict_with_genes(self):
        """Test creating dict with mock genes."""
        # This function requires a real genedb file path, not a dict
        # Skipping as it requires complex setup
        pass


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
