############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import pytest
import tempfile
import os
from src.string_pools import StringPool, StringPoolManager


class TestStringPool:
    """Test StringPool class."""

    def test_init(self):
        """Test initialization."""
        pool = StringPool()
        assert len(pool) == 0
        assert pool.next_id == 0

    def test_add_single(self):
        """Test adding a single string."""
        pool = StringPool()
        id1 = pool.add("gene1")
        assert id1 == 0
        assert len(pool) == 1
        assert "gene1" in pool

    def test_add_multiple(self):
        """Test adding multiple unique strings."""
        pool = StringPool()
        id1 = pool.add("gene1")
        id2 = pool.add("gene2")
        id3 = pool.add("gene3")

        assert id1 == 0
        assert id2 == 1
        assert id3 == 2
        assert len(pool) == 3

    def test_add_duplicate(self):
        """Test adding duplicate string returns same ID."""
        pool = StringPool()
        id1 = pool.add("gene1")
        id2 = pool.add("gene1")

        assert id1 == id2
        assert len(pool) == 1

    def test_get_int(self):
        """Test get_int method."""
        pool = StringPool()
        id1 = pool.get_int("gene1")
        id2 = pool.get_int("gene2")
        id3 = pool.get_int("gene1")  # Duplicate

        assert id1 == 0
        assert id2 == 1
        assert id3 == 0  # Same as first
        assert len(pool) == 2

    def test_get_str(self):
        """Test get_str method."""
        pool = StringPool()
        pool.add("gene1")
        pool.add("gene2")
        pool.add("gene3")

        assert pool.get_str(0) == "gene1"
        assert pool.get_str(1) == "gene2"
        assert pool.get_str(2) == "gene3"

    def test_get_str_invalid(self):
        """Test get_str with invalid ID."""
        pool = StringPool()
        pool.add("gene1")

        with pytest.raises(IndexError):
            pool.get_str(999)

    def test_contains(self):
        """Test __contains__ method."""
        pool = StringPool()
        pool.add("gene1")

        assert "gene1" in pool
        assert "gene2" not in pool

    def test_bidirectional_consistency(self):
        """Test that str->int->str round trip works."""
        pool = StringPool()
        strings = ["ENSG00000001", "ENSG00000002", "ENSG00000003"]

        for s in strings:
            int_id = pool.add(s)
            assert pool.get_str(int_id) == s

    def test_deterministic_order(self):
        """Test that order of addition determines IDs."""
        pool = StringPool()
        strings = ["zebra", "apple", "mango"]

        # Add in specific order
        ids = [pool.add(s) for s in strings]

        # IDs should be 0, 1, 2 regardless of alphabetical order
        assert ids == [0, 1, 2]
        assert pool.get_str(0) == "zebra"
        assert pool.get_str(1) == "apple"
        assert pool.get_str(2) == "mango"


class TestStringPoolManager:
    """Test StringPoolManager class."""

    def test_init(self):
        """Test initialization."""
        manager = StringPoolManager()

        assert len(manager.gene_pool) == 0
        assert len(manager.transcript_pool) == 0
        assert len(manager.chromosome_pool) == 0
        assert len(manager.barcode_pool) == 0
        assert len(manager.umi_pool) == 0

    def test_build_file_name_pool(self):
        """Test building file name pool."""
        manager = StringPoolManager()

        # Mock sample object without labels
        class MockSample:
            file_list = [
                ("/path/to/file1.bam", None),
                ("/path/to/file2.bam", None),
                ("/another/path/file1.bam", None),  # Same basename
            ]
            readable_names_dict = None

        sample = MockSample()
        manager.build_file_name_pool(sample)

        # Pool stores basenames without extensions (matching FileNameGrouper logic)
        assert len(manager.file_name_pool) == 2  # file1, file2
        assert "file1" in manager.file_name_pool
        assert "file2" in manager.file_name_pool

    def test_build_file_name_pool_with_labels(self):
        """Test building file name pool when -l labels are used."""
        manager = StringPoolManager()

        # Mock sample object with labels
        class MockSample:
            file_list = [
                ("/path/to/file1.bam", None),
                ("/path/to/file2.bam", None),
                ("/path/to/file3.bam", None),
            ]
            # readable_names_dict maps file paths to labels (set via -l option)
            readable_names_dict = {
                "/path/to/file1.bam": "g1",
                "/path/to/file2.bam": "g2",
                "/path/to/file3.bam": "g3",
            }

        sample = MockSample()
        manager.build_file_name_pool(sample)

        # Pool should contain labels, not basenames
        assert len(manager.file_name_pool) == 3  # g1, g2, g3
        assert "g1" in manager.file_name_pool
        assert "g2" in manager.file_name_pool
        assert "g3" in manager.file_name_pool
        # Basenames should NOT be in the pool
        assert "file1" not in manager.file_name_pool

    def test_build_barcode_spot_pool(self):
        """Test building barcode-to-spot pool."""
        manager = StringPoolManager()

        # Create temporary barcode2spot file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("ACTG\tcell_typeA\n")
            f.write("TGCA\tcell_typeB\n")
            f.write("GGGG\tcell_typeA\n")  # Duplicate spot
            temp_file = f.name

        try:
            manager.build_barcode_spot_pool(temp_file, column_index=0, file_column=1)

            assert len(manager.barcode_spot_pools[0]) == 3  # cell_typeA, cell_typeB, NA
            assert "cell_typeA" in manager.barcode_spot_pools[0]
            assert "cell_typeB" in manager.barcode_spot_pools[0]
        finally:
            os.unlink(temp_file)

    def test_load_barcode_pool(self):
        """Test loading barcode pool from file."""
        manager = StringPoolManager()

        # Create temporary barcode file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("read_001\tACTG\tGGGG\n")
            f.write("read_002\tTGCA\tAAAA\n")
            f.write("read_003\tACTG\tCCCC\n")  # Duplicate barcode
            temp_file = f.name

        try:
            manager.load_barcode_pool(temp_file)

            assert len(manager.barcode_pool) == 3  # ACTG, TGCA
            assert "ACTG" in manager.barcode_pool
            assert "TGCA" in manager.barcode_pool

            assert len(manager.umi_pool) == 3  # GGGG, AAAA, CCCC
            assert "GGGG" in manager.umi_pool
        finally:
            os.unlink(temp_file)

    def test_load_read_group_tsv_pool(self):
        """Test loading read group TSV pool."""
        manager = StringPoolManager()

        # Create temporary read group file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
            f.write("read_001\tgroupA\n")
            f.write("read_002\tgroupB\n")
            f.write("read_003\tgroupA\n")  # Duplicate group
            temp_file = f.name

        try:
            pool_key = "0:1"  # Format: spec_index:col_index
            manager.load_read_group_tsv_pool(temp_file, pool_key, col_index=1, delimiter='\t')

            assert pool_key in manager.read_group_tsv_pools
            pool = manager.read_group_tsv_pools[pool_key]
            assert len(pool) == 3  # NA (default), groupA, groupB
            assert "NA" in pool  # Default group ID always included
            assert "groupA" in pool
            assert "groupB" in pool
        finally:
            os.unlink(temp_file)

    def test_load_multiple_read_group_tsv_pools(self):
        """Test loading multiple read group TSV pools."""
        manager = StringPoolManager()

        # Create two temporary files
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f1:
            f1.write("read_001\tgroupA\n")
            temp_file1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f2:
            f2.write("read_001\tgroupX\n")
            temp_file2 = f2.name

        try:
            pool_key1 = "0:1"  # spec 0, column 1
            pool_key2 = "1:1"  # spec 1, column 1
            manager.load_read_group_tsv_pool(temp_file1, pool_key1, col_index=1, delimiter='\t')
            manager.load_read_group_tsv_pool(temp_file2, pool_key2, col_index=1, delimiter='\t')

            assert len(manager.read_group_tsv_pools) == 2
            assert "groupA" in manager.read_group_tsv_pools[pool_key1]
            assert "groupX" in manager.read_group_tsv_pools[pool_key2]
            # Both should also have NA
            assert "NA" in manager.read_group_tsv_pools[pool_key1]
            assert "NA" in manager.read_group_tsv_pools[pool_key2]
        finally:
            os.unlink(temp_file1)
            os.unlink(temp_file2)

    def test_load_multicolumn_tsv_pool(self):
        """Test loading multi-column TSV with different columns creating separate pools."""
        manager = StringPoolManager()

        # Create temporary multi-column file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
            f.write("read_001,groupA,typeX\n")
            f.write("read_002,groupB,typeY\n")
            f.write("read_003,groupA,typeX\n")  # Duplicates
            temp_file = f.name

        try:
            # Load column 1 (groups)
            pool_key1 = "0:1"
            manager.load_read_group_tsv_pool(temp_file, pool_key1, col_index=1, delimiter=',')

            # Load column 2 (types)
            pool_key2 = "0:2"
            manager.load_read_group_tsv_pool(temp_file, pool_key2, col_index=2, delimiter=',')

            assert len(manager.read_group_tsv_pools) == 2

            # Column 1 pool
            pool1 = manager.read_group_tsv_pools[pool_key1]
            assert len(pool1) == 3  # NA, groupA, groupB
            assert "NA" in pool1
            assert "groupA" in pool1
            assert "groupB" in pool1

            # Column 2 pool
            pool2 = manager.read_group_tsv_pools[pool_key2]
            assert len(pool2) == 3  # NA, typeX, typeY
            assert "NA" in pool2
            assert "typeX" in pool2
            assert "typeY" in pool2
        finally:
            os.unlink(temp_file)

    def test_load_nonexistent_file(self):
        """Test loading from nonexistent file doesn't crash."""
        manager = StringPoolManager()

        # Should not raise exception
        manager.load_barcode_pool("/nonexistent/file.tsv")
        manager.load_read_group_tsv_pool("/nonexistent/file.tsv", pool_key="0:1", col_index=1, delimiter='\t')

        assert len(manager.barcode_pool) == 0
        assert len(manager.read_group_tsv_pools) == 0

    def test_get_stats(self):
        """Test get_stats method."""
        manager = StringPoolManager()
        manager.gene_pool.add("gene1")
        manager.transcript_pool.add("tx1")
        manager.barcode_pool.add("ACTG")

        stats = manager.get_stats()

        assert "Genes: 1" in stats
        assert "Transcripts: 1" in stats
        assert "Barcodes: 1" in stats

    def test_build_chromosome_pool(self):
        """Test building chromosome pool."""
        manager = StringPoolManager()

        chr_ids = ["chr3", "chr1", "chr2", "chrX"]
        manager.build_chromosome_pool(chr_ids)

        # Pool should have all chromosomes
        assert len(manager.chromosome_pool) == 4
        assert "chr1" in manager.chromosome_pool
        assert "chr2" in manager.chromosome_pool
        assert "chr3" in manager.chromosome_pool
        assert "chrX" in manager.chromosome_pool

        # IDs should be assigned in sorted order (deterministic across workers)
        assert manager.chromosome_pool.get_int("chr1") == 0
        assert manager.chromosome_pool.get_int("chr2") == 1
        assert manager.chromosome_pool.get_int("chr3") == 2
        assert manager.chromosome_pool.get_int("chrX") == 3

    def test_set_and_get_group_spec_pool_type_file_name(self):
        """Test setting and getting pool type for file_name grouping."""
        manager = StringPoolManager()

        # Set up file_name pool
        manager.file_name_pool.add("sample1")
        manager.file_name_pool.add("sample2")

        # Set pool type for spec index 0
        manager.set_group_spec_pool_type(0, 'file_name')

        # get_read_group_pool should return file_name_pool
        pool = manager.get_read_group_pool(0)
        assert pool is manager.file_name_pool
        assert "sample1" in pool
        assert "sample2" in pool

    def test_set_and_get_group_spec_pool_type_barcode_spot(self):
        """Test setting and getting pool type for barcode_spot grouping."""
        manager = StringPoolManager()

        # Set up barcode_spot pool at column index 0
        pool = StringPool()
        pool.add("cell_typeA")
        pool.add("cell_typeB")
        manager.barcode_spot_pools[0] = pool

        # Set pool type for spec index 1 (legacy format)
        manager.set_group_spec_pool_type(1, 'barcode_spot')

        # get_read_group_pool should return barcode_spot_pools[0]
        returned_pool = manager.get_read_group_pool(1)
        assert returned_pool is manager.barcode_spot_pools[0]
        assert "cell_typeA" in returned_pool

    def test_set_and_get_group_spec_pool_type_barcode_spot_multi_column(self):
        """Test setting and getting pool type for multi-column barcode_spot grouping."""
        manager = StringPoolManager()

        # Set up barcode_spot pools for multiple columns
        pool0 = StringPool()
        pool0.add("organ_liver")
        pool0.add("organ_kidney")
        manager.barcode_spot_pools[0] = pool0

        pool1 = StringPool()
        pool1.add("cell_hepatocyte")
        pool1.add("cell_nephron")
        manager.barcode_spot_pools[1] = pool1

        # Set pool types for spec indices with explicit column format
        manager.set_group_spec_pool_type(0, 'barcode_spot:0')
        manager.set_group_spec_pool_type(1, 'barcode_spot:1')

        # get_read_group_pool should return correct pools
        assert manager.get_read_group_pool(0) is manager.barcode_spot_pools[0]
        assert manager.get_read_group_pool(1) is manager.barcode_spot_pools[1]
        assert "organ_liver" in manager.get_read_group_pool(0)
        assert "cell_hepatocyte" in manager.get_read_group_pool(1)

    def test_set_and_get_group_spec_pool_type_tsv(self):
        """Test setting and getting pool type for TSV grouping."""
        manager = StringPoolManager()

        # Create a TSV pool
        pool_key = "0:1"
        tsv_pool = StringPool()
        tsv_pool.add("groupA")
        tsv_pool.add("groupB")
        manager.read_group_tsv_pools[pool_key] = tsv_pool

        # Set pool type for spec index 0 (format: tsv:spec_idx:col_idx:delimiter)
        manager.set_group_spec_pool_type(0, 'tsv:0:1:\t')

        # get_read_group_pool should return the TSV pool
        pool = manager.get_read_group_pool(0)
        assert pool is tsv_pool
        assert "groupA" in pool

    def test_get_read_group_pool_dynamic_fallback(self):
        """Test that get_read_group_pool creates dynamic pool if no type set."""
        manager = StringPoolManager()

        # No pool type set for spec index 5
        pool = manager.get_read_group_pool(5)

        # Should create a new dynamic pool
        assert pool is not None
        assert 5 in manager.read_group_dynamic_pools
        assert pool is manager.read_group_dynamic_pools[5]

    def test_get_or_create_dynamic_pool(self):
        """Test dynamic pool creation."""
        manager = StringPoolManager()

        # First call creates the pool
        pool1 = manager._get_or_create_dynamic_pool(3)
        assert pool1 is not None
        assert 3 in manager.read_group_dynamic_pools

        # Second call returns the same pool
        pool2 = manager._get_or_create_dynamic_pool(3)
        assert pool1 is pool2

        # Different index creates different pool
        pool3 = manager._get_or_create_dynamic_pool(4)
        assert pool3 is not pool1
        assert 4 in manager.read_group_dynamic_pools

    def test_read_group_to_ids_and_from_ids(self):
        """Test converting between group strings and IDs."""
        manager = StringPoolManager()

        # Set up pools for different spec indices
        manager.set_group_spec_pool_type(0, 'file_name')
        manager.file_name_pool.add("sample1")
        manager.file_name_pool.add("sample2")

        manager.set_group_spec_pool_type(1, 'dynamic')
        dynamic_pool = manager._get_or_create_dynamic_pool(1)
        dynamic_pool.add("groupA")
        dynamic_pool.add("groupB")

        # Convert strings to IDs
        group_strings = ["sample1", "groupB"]
        ids = manager.read_group_to_ids(group_strings)

        assert ids[0] == manager.file_name_pool.get_int("sample1")
        assert ids[1] == dynamic_pool.get_int("groupB")

        # Convert IDs back to strings
        recovered = manager.read_group_from_ids(ids)
        assert recovered == group_strings

    def test_read_group_to_ids_with_none(self):
        """Test that None values are handled correctly."""
        manager = StringPoolManager()

        manager.set_group_spec_pool_type(0, 'file_name')
        manager.file_name_pool.add("sample1")

        # Convert with None value
        group_strings = ["sample1", None]
        ids = manager.read_group_to_ids(group_strings)

        assert ids[0] == manager.file_name_pool.get_int("sample1")
        assert ids[1] == -1  # None -> -1

        # Convert back
        recovered = manager.read_group_from_ids(ids)
        assert recovered[0] == "sample1"
        assert recovered[1] is None  # -1 -> None

    def test_read_group_to_ids_empty(self):
        """Test empty list handling."""
        manager = StringPoolManager()

        ids = manager.read_group_to_ids([])
        assert ids == []

        strings = manager.read_group_from_ids([])
        assert strings == []

    def test_has_dynamic_pools(self):
        """Test has_dynamic_pools method."""
        manager = StringPoolManager()

        # No dynamic pools yet
        assert manager.has_dynamic_pools() is False

        # Create empty dynamic pool
        manager._get_or_create_dynamic_pool(0)
        assert manager.has_dynamic_pools() is False  # Empty pool

        # Add data to dynamic pool
        manager.read_group_dynamic_pools[0].add("value1")
        assert manager.has_dynamic_pools() is True

    def test_serialize_and_deserialize_dynamic_pools(self):
        """Test serialization and deserialization of dynamic pools."""
        import io

        manager = StringPoolManager()

        # Create dynamic pools with data
        pool0 = manager._get_or_create_dynamic_pool(0)
        pool0.add("group_a")
        pool0.add("group_b")
        pool0.add("group_c")

        pool2 = manager._get_or_create_dynamic_pool(2)
        pool2.add("type_x")
        pool2.add("type_y")

        # Serialize
        buffer = io.BytesIO()
        manager.serialize_dynamic_pools(buffer)

        # Create new manager and deserialize
        manager2 = StringPoolManager()
        buffer.seek(0)
        manager2.deserialize_dynamic_pools(buffer)

        # Verify pools were restored
        assert 0 in manager2.read_group_dynamic_pools
        assert 2 in manager2.read_group_dynamic_pools

        # Verify pool contents
        restored_pool0 = manager2.read_group_dynamic_pools[0]
        assert len(restored_pool0) == 3
        assert "group_a" in restored_pool0
        assert "group_b" in restored_pool0
        assert "group_c" in restored_pool0

        restored_pool2 = manager2.read_group_dynamic_pools[2]
        assert len(restored_pool2) == 2
        assert "type_x" in restored_pool2
        assert "type_y" in restored_pool2

        # Verify IDs are preserved (critical for deserialization consistency)
        assert restored_pool0.get_int("group_a") == pool0.get_int("group_a")
        assert restored_pool0.get_int("group_b") == pool0.get_int("group_b")
        assert restored_pool2.get_int("type_x") == pool2.get_int("type_x")

    def test_serialize_empty_dynamic_pools(self):
        """Test serialization with no dynamic pools."""
        import io

        manager = StringPoolManager()

        # Serialize empty pools
        buffer = io.BytesIO()
        manager.serialize_dynamic_pools(buffer)

        # Deserialize into new manager
        manager2 = StringPoolManager()
        buffer.seek(0)
        manager2.deserialize_dynamic_pools(buffer)

        # Should have no dynamic pools
        assert len(manager2.read_group_dynamic_pools) == 0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
