############################################################################
# Copyright (c) 2025 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import pytest
import tempfile
import os
from src.read_groups import (
    AlignmentTagReadGrouper,
    ReadIdSplitReadGrouper,
    ReadTableGrouper,
)


# Mock alignment object
class MockAlignment:
    def __init__(self, read_id="", tags=None):
        self.query_name = read_id
        if tags is None:
            tags = {}
        self._tags = tags

    def get_tag(self, tag_name):
        return self._tags.get(tag_name)

    def has_tag(self, tag_name):
        return tag_name in self._tags


class TestAlignmentTagReadGrouper:
    """Test AlignmentTagReadGrouper class."""

    def test_init_default_tag(self):
        """Test AlignmentTagReadGrouper with default RG tag."""
        grouper = AlignmentTagReadGrouper()
        assert grouper.tag == "RG"

    def test_init_custom_tag(self):
        """Test AlignmentTagReadGrouper with custom tag."""
        grouper = AlignmentTagReadGrouper("CB")
        assert grouper.tag == "CB"

    def test_get_group_with_tag(self):
        """Test getting group when read has tag."""
        grouper = AlignmentTagReadGrouper("CB")
        alignment = MockAlignment(tags={"CB": "ACTGACTG"})
        assert grouper.get_group_id(alignment, None) == "ACTGACTG"

    def test_get_group_missing_tag(self):
        """Test getting group when read lacks tag."""
        grouper = AlignmentTagReadGrouper("CB")
        alignment = MockAlignment()
        result = grouper.get_group_id(alignment, None)
        # Should return None
        assert result is None


class TestReadIdSplitReadGrouper:
    """Test ReadIdSplitReadGrouper class."""

    def test_init(self):
        """Test ReadIdSplitReadGrouper initialization."""
        grouper = ReadIdSplitReadGrouper("_")
        assert grouper.delim == "_"

    def test_get_group_with_delimiter(self):
        """Test getting group from read ID with delimiter."""
        grouper = ReadIdSplitReadGrouper("_")

        alignment1 = MockAlignment(read_id="read_001_groupA")
        alignment2 = MockAlignment(read_id="read_002_groupB")

        assert grouper.get_group_id(alignment1) == "groupA"
        assert grouper.get_group_id(alignment2) == "groupB"

    def test_get_group_no_delimiter(self):
        """Test getting group from read ID without delimiter."""
        grouper = ReadIdSplitReadGrouper("_")
        alignment = MockAlignment(read_id="read001")

        # Returns empty string if no delimiter found
        result = grouper.get_group_id(alignment)
        assert result == ""

    def test_get_group_multiple_delimiters(self):
        """Test getting group with multiple delimiters."""
        grouper = ReadIdSplitReadGrouper("_")
        alignment = MockAlignment(read_id="prefix_middle_suffix")

        # Should return last part after delimiter
        assert grouper.get_group_id(alignment) == "suffix"


class TestReadTableGrouper:
    """Test ReadTableGrouper class."""

    def test_init_and_get_group(self):
        """Test ReadTableGrouper initialization and file loading."""
        # Create temporary file with read-group mapping
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("read_001\tgroupA\n")
            f.write("read_002\tgroupB\n")
            f.write("read_003\tgroupA\n")
            temp_file = f.name

        try:
            from src.read_groups import SharedTableData, ReadTableGrouper

            # Create shared data with single column
            shared_data = SharedTableData(temp_file,
                                         read_id_column_index=0,
                                         group_id_column_indices=[1],
                                         delim='\t')
            grouper = ReadTableGrouper(shared_data, 0)

            alignment1 = MockAlignment(read_id="read_001")
            alignment2 = MockAlignment(read_id="read_002")
            alignment3 = MockAlignment(read_id="read_003")

            assert grouper.get_group_id(alignment1) == "groupA"
            assert grouper.get_group_id(alignment2) == "groupB"
            assert grouper.get_group_id(alignment3) == "groupA"
        finally:
            os.unlink(temp_file)

    def test_missing_read(self):
        """Test getting group for read not in file."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("read_001\tgroupA\n")
            temp_file = f.name

        try:
            from src.read_groups import SharedTableData, ReadTableGrouper

            # Create shared data with single column
            shared_data = SharedTableData(temp_file,
                                         read_id_column_index=0,
                                         group_id_column_indices=[1],
                                         delim='\t')
            grouper = ReadTableGrouper(shared_data, 0)

            alignment = MockAlignment(read_id="read_999")
            # Should return "NA" for missing reads
            result = grouper.get_group_id(alignment)
            assert result == "NA"
        finally:
            os.unlink(temp_file)


class TestSharedTableData:
    """Test SharedTableData and ReadTableGrouper classes."""

    def test_shared_data_multicolumn(self):
        """Test shared table data with multiple columns."""
        # Create temporary file with multi-column read-group mapping
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("read_001\tcell_A\tspot_1\ttissue_X\n")
            f.write("read_002\tcell_B\tspot_2\ttissue_Y\n")
            f.write("read_003\tcell_A\tspot_1\ttissue_X\n")
            temp_file = f.name

        try:
            from src.read_groups import SharedTableData, ReadTableGrouper

            # Create shared data with 3 group columns (indices 1, 2, 3)
            shared_data = SharedTableData(temp_file,
                                         read_id_column_index=0,
                                         group_id_column_indices=[1, 2, 3],
                                         delim='\t')

            # Create separate groupers for each column
            grouper_col1 = ReadTableGrouper(shared_data, 0)  # cell type
            grouper_col2 = ReadTableGrouper(shared_data, 1)  # spot
            grouper_col3 = ReadTableGrouper(shared_data, 2)  # tissue

            alignment1 = MockAlignment(read_id="read_001")
            alignment2 = MockAlignment(read_id="read_002")
            alignment3 = MockAlignment(read_id="read_003")

            # Test column 1 (cell type)
            assert grouper_col1.get_group_id(alignment1) == "cell_A"
            assert grouper_col1.get_group_id(alignment2) == "cell_B"
            assert grouper_col1.get_group_id(alignment3) == "cell_A"

            # Test column 2 (spot)
            assert grouper_col2.get_group_id(alignment1) == "spot_1"
            assert grouper_col2.get_group_id(alignment2) == "spot_2"
            assert grouper_col2.get_group_id(alignment3) == "spot_1"

            # Test column 3 (tissue)
            assert grouper_col3.get_group_id(alignment1) == "tissue_X"
            assert grouper_col3.get_group_id(alignment2) == "tissue_Y"
            assert grouper_col3.get_group_id(alignment3) == "tissue_X"

            # Verify all three groupers share the same data object
            assert grouper_col1.shared_data is shared_data
            assert grouper_col2.shared_data is shared_data
            assert grouper_col3.shared_data is shared_data

        finally:
            os.unlink(temp_file)

    def test_shared_data_missing_read(self):
        """Test shared data with missing read."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("read_001\tgroupA\tgroupX\n")
            temp_file = f.name

        try:
            from src.read_groups import SharedTableData, ReadTableGrouper

            shared_data = SharedTableData(temp_file,
                                         read_id_column_index=0,
                                         group_id_column_indices=[1, 2],
                                         delim='\t')

            grouper = ReadTableGrouper(shared_data, 0)
            alignment = MockAlignment(read_id="read_999")

            # Should return "NA" for missing reads
            result = grouper.get_group_id(alignment)
            assert result == "NA"

        finally:
            os.unlink(temp_file)


class TestParseGroupingSpec:
    """Test parse_grouping_spec function."""

    def test_parse_multicolumn_creates_separate_groupers(self):
        """Test that comma-separated columns create separate groupers sharing data."""
        # Create temporary file with multi-column data
        # Note: parse_grouping_spec appends "_spec{index}_chr_id" to the file path
        import tempfile
        temp_dir = tempfile.mkdtemp()
        base_file = os.path.join(temp_dir, "test_groups")
        chr_id = "chr1"
        spec_index = 0
        chr_file = base_file + "_spec" + str(spec_index) + "_" + chr_id

        with open(chr_file, 'w') as f:
            f.write("read_001\tcell_A\tspot_1\n")
            f.write("read_002\tcell_B\tspot_2\n")

        try:
            from src.read_groups import parse_grouping_spec, ReadTableGrouper

            # Mock args and sample objects
            class MockArgs:
                pass

            class MockSample:
                def __init__(self, file_path):
                    self.read_group_file = file_path

            args = MockArgs()
            sample = MockSample(base_file)  # Use base file, function will append spec_index and chr_id

            # Parse specification with comma-separated columns: file:path:0:1,2
            # Note: the file path in spec_string is ignored, sample.read_group_file is used
            spec_string = f"file:dummy:0:1,2"
            result = parse_grouping_spec(spec_string, args, sample, chr_id, spec_index)

            # Should return a list of two groupers
            assert isinstance(result, list)
            assert len(result) == 2

            # Both should be ReadTableGrouper instances
            assert isinstance(result[0], ReadTableGrouper)
            assert isinstance(result[1], ReadTableGrouper)

            # Both should share the same SharedTableData
            assert result[0].shared_data is result[1].shared_data

            # Column indices should be 0 and 1
            assert result[0].column_index == 0
            assert result[1].column_index == 1

            # Verify they return correct values
            alignment1 = MockAlignment(read_id="read_001")
            alignment2 = MockAlignment(read_id="read_002")

            assert result[0].get_group_id(alignment1) == "cell_A"
            assert result[0].get_group_id(alignment2) == "cell_B"
            assert result[1].get_group_id(alignment1) == "spot_1"
            assert result[1].get_group_id(alignment2) == "spot_2"

        finally:
            import shutil
            shutil.rmtree(temp_dir)

    def test_multiple_different_files(self):
        """Test that multiple different table files use separate split files."""
        import tempfile
        temp_dir = tempfile.mkdtemp()
        base_file = os.path.join(temp_dir, "test_groups")
        chr_id = "chr1"

        # Create split files for two different specs
        spec0_file = base_file + "_spec0_" + chr_id
        spec1_file = base_file + "_spec1_" + chr_id

        with open(spec0_file, 'w') as f:
            f.write("read_001\tgroupA\n")
            f.write("read_002\tgroupB\n")

        with open(spec1_file, 'w') as f:
            f.write("read_001\tgroupX\n")
            f.write("read_002\tgroupY\n")

        try:
            from src.read_groups import parse_grouping_spec, ReadTableGrouper

            class MockArgs:
                pass

            class MockSample:
                def __init__(self, file_path):
                    self.read_group_file = file_path

            args = MockArgs()
            sample = MockSample(base_file)

            # Parse two different file specs
            spec0_string = "file:groups1.tsv:0:1"
            spec1_string = "file:groups2.tsv:0:1"

            grouper0 = parse_grouping_spec(spec0_string, args, sample, chr_id, spec_index=0)
            grouper1 = parse_grouping_spec(spec1_string, args, sample, chr_id, spec_index=1)

            # Both should be ReadTableGrouper instances
            assert isinstance(grouper0, ReadTableGrouper)
            assert isinstance(grouper1, ReadTableGrouper)

            # They should use different underlying data
            assert grouper0.shared_data is not grouper1.shared_data

            # Verify they return different values for the same read
            alignment1 = MockAlignment(read_id="read_001")
            alignment2 = MockAlignment(read_id="read_002")

            assert grouper0.get_group_id(alignment1) == "groupA"
            assert grouper0.get_group_id(alignment2) == "groupB"

            assert grouper1.get_group_id(alignment1) == "groupX"
            assert grouper1.get_group_id(alignment2) == "groupY"

        finally:
            import shutil
            shutil.rmtree(temp_dir)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
