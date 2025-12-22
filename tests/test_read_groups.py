############################################################################
# Copyright (c) 2025 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import pytest
import tempfile
import os
from src.read_groups import (
    FileNameGrouper,
    ReadTableGrouper,
    BarcodeSpotGrouper,
    AlignmentTagReadGrouper,
    ReadIdSplitReadGrouper,
    parse_grouping_spec,
    get_grouping_strategy_names
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


class TestFileNameGrouper:
    """Test FileNameGrouper class."""

    def test_init(self):
        """Test FileNameGrouper initialization."""
        file_dict = {"file1": "label1", "file2": "label2"}
        grouper = FileNameGrouper(file_dict)
        assert grouper.file_dict == file_dict

    def test_get_group(self):
        """Test getting group from filename."""
        file_dict = {"file1.bam": "sample1", "file2.bam": "sample2"}
        grouper = FileNameGrouper(file_dict)

        assert grouper.get_group("read1", "file1.bam") == "sample1"
        assert grouper.get_group("read2", "file2.bam") == "sample2"

    def test_get_group_missing_file(self):
        """Test getting group from missing file."""
        file_dict = {"file1.bam": "sample1"}
        grouper = FileNameGrouper(file_dict)

        # Should return filename if not in dict
        assert grouper.get_group("read1", "file3.bam") == "file3.bam"


class TestAlignmentTagReadGrouper:
    """Test AlignmentTagReadGrouper class."""

    def test_init_default_tag(self):
        """Test AlignmentTagReadGrouper with default RG tag."""
        grouper = AlignmentTagReadGrouper(None)
        assert grouper.tag == "RG"

    def test_init_custom_tag(self):
        """Test AlignmentTagReadGrouper with custom tag."""
        grouper = AlignmentTagReadGrouper("CB")
        assert grouper.tag == "CB"

    def test_get_group_with_tag(self):
        """Test getting group when read has tag."""
        grouper = AlignmentTagReadGrouper("CB")
        alignment = MockAlignment(tags={"CB": "ACTGACTG"})
        assert grouper.get_group_id(alignment, None, ) == "ACTGACTG"

    def test_get_group_missing_tag(self):
        """Test getting group when read lacks tag."""
        grouper = AlignmentTagReadGrouper("CB")
        alignment = MockAlignment()
        result = grouper.get_group_id(alignment, None, )
        # Should return None or some default value
        assert result is None or result == ""


class TestReadIdSplitReadGrouper:
    """Test ReadIdSplitReadGrouper class."""

    def test_init(self):
        """Test ReadIdSplitReadGrouper initialization."""
        grouper = ReadIdSplitReadGrouper("_")
        assert grouper.delim == "_"

    def test_get_group_with_delimiter(self):
        """Test getting group from read ID with delimiter."""
        grouper = ReadIdSplitReadGrouper("_")

        assert grouper.get_group_id(MockAlignment(read_id="read_001_groupA")) == "groupA"
        assert grouper.get_group_id(MockAlignment(read_id="read_002_groupB")) == "groupB"

    def test_get_group_no_delimiter(self):
        """Test getting group from read ID without delimiter."""
        grouper = ReadIdSplitReadGrouper("_")

        # Should return read ID if no delimiter found
        assert grouper.get_group_id(MockAlignment(read_id="read001")) == "read001"

    def test_get_group_multiple_delimiters(self):
        """Test getting group with multiple delimiters."""
        grouper = ReadIdSplitReadGrouper("_")

        # Should return last part after delimiter
        assert grouper.get_group_id(MockAlignment(read_id="prefix_middle_suffix")) == "suffix"


class TestReadTableGrouper:
    """Test ReadTableGrouper class."""

    def test_init_and_load(self):
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

            assert grouper.get_group_id("read_001") == "groupA"
            assert grouper.get_group_id("read_002") == "groupB"
            assert grouper.get_group_id("read_003") == "groupA"
        finally:
            os.unlink(temp_file)

    def test_missing_read(self):
        """Test getting group for read not in file."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("read_001\tgroupA\n")
            temp_file = f.name

        try:
            grouper = ReadTableGrouper([temp_file], read_col=0, group_col=1, delimiter='\t')
            grouper.load()

            # Should return None or read ID for missing reads
            result = grouper.get_group("read_999")
            assert result is None or result == ""
        finally:
            os.unlink(temp_file)


class TestBarcodeSpotGrouper:
    """Test BarcodeSpotGrouper class."""

    def test_init(self):
        """Test BarcodeSpotGrouper initialization."""
        # Create temp files
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("read_001\tACTG\tGGGG\n")
            barcode_file = f.name

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("ACTG\tcellTypeA\n")
            spot_file = f.name

        try:
            grouper = BarcodeSpotGrouper([barcode_file], [spot_file])
            assert grouper.barcode_files == [barcode_file]
            assert grouper.spot_files == [spot_file]
        finally:
            os.unlink(barcode_file)
            os.unlink(spot_file)

    def test_load_and_get_group(self):
        """Test loading barcode-spot mappings."""
        # Create temp barcode file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("read_001\tACTG\tGGGG\n")
            f.write("read_002\tTGCA\tCCCC\n")
            barcode_file = f.name

        # Create temp spot mapping file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("ACTG\tcellTypeA\n")
            f.write("TGCA\tcellTypeB\n")
            spot_file = f.name

        try:
            grouper = BarcodeSpotGrouper([barcode_file], [spot_file])
            grouper.load()

            # Should map read -> barcode -> cell type
            assert grouper.get_group("read_001") == "cellTypeA"
            assert grouper.get_group("read_002") == "cellTypeB"
        finally:
            os.unlink(barcode_file)
            os.unlink(spot_file)

    def test_missing_barcode_in_spot_map(self):
        """Test read with barcode not in spot mapping."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("read_001\tACTG\tGGGG\n")
            barcode_file = f.name

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("TGCA\tcellTypeB\n")  # Different barcode
            spot_file = f.name

        try:
            grouper = BarcodeSpotGrouper([barcode_file], [spot_file])
            grouper.load()

            # Should return None or barcode for unmapped barcodes
            result = grouper.get_group("read_001")
            assert result is None or result == "ACTG"
        finally:
            os.unlink(barcode_file)
            os.unlink(spot_file)


class TestParseGroupingSpec:
    """Test grouping specification parsing."""

    def test_parse_file_name(self):
        """Test parsing file_name grouping."""
        file_dict = {"file1": "label1"}
        groupers = parse_grouping_spec(["file_name"], file_dict)

        assert len(groupers) == 1
        assert isinstance(groupers[0], FileNameGrouper)

    def test_parse_tag(self):
        """Test parsing tag grouping."""
        groupers = parse_grouping_spec(["tag:CB"], {})

        assert len(groupers) == 1
        assert isinstance(groupers[0], AlignmentTagReadGrouper)
        assert groupers[0].tag_name == "CB"

    def test_parse_read_id(self):
        """Test parsing read_id grouping."""
        groupers = parse_grouping_spec(["read_id:_"], {})

        assert len(groupers) == 1
        assert isinstance(groupers[0], ReadIdSplitReadGrouper)
        assert groupers[0].delimiter == "_"

    def test_parse_multiple_strategies(self):
        """Test parsing multiple grouping strategies."""
        file_dict = {"file1": "label1"}
        groupers = parse_grouping_spec(["file_name", "tag:CB"], file_dict)

        assert len(groupers) == 2
        assert isinstance(groupers[0], FileNameGrouper)
        assert isinstance(groupers[1], AlignmentTagReadGrouper)


class TestGetGroupingStrategyNames:
    """Test getting grouping strategy names."""

    def test_file_name_strategy(self):
        """Test file_name strategy name."""
        file_dict = {"file1": "label1"}
        groupers = parse_grouping_spec(["file_name"], file_dict)
        names = get_grouping_strategy_names(groupers)

        assert names == ["file_name"]

    def test_tag_strategy(self):
        """Test tag strategy name."""
        groupers = parse_grouping_spec(["tag:CB"], {})
        names = get_grouping_strategy_names(groupers)

        assert names == ["tag"]

    def test_multiple_strategies(self):
        """Test multiple strategy names."""
        file_dict = {"file1": "label1"}
        groupers = parse_grouping_spec(["file_name", "tag:RG"], file_dict)
        names = get_grouping_strategy_names(groupers)

        assert len(names) == 2
        assert "file_name" in names
        assert "tag" in names


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
