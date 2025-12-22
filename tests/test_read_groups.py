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
            grouper = ReadTableGrouper(temp_file,
                                       read_id_column_index=0,
                                       group_id_column_index=1,
                                       delim='\t')

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
            grouper = ReadTableGrouper(temp_file,
                                       read_id_column_index=0,
                                       group_id_column_index=1,
                                       delim='\t')

            alignment = MockAlignment(read_id="read_999")
            # Should return "NA" for missing reads
            result = grouper.get_group_id(alignment)
            assert result == "NA"
        finally:
            os.unlink(temp_file)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
