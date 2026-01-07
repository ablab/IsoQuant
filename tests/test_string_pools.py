############################################################################
# Copyright (c) 2025 University of Helsinki
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
            manager.build_barcode_spot_pool([temp_file])

            assert len(manager.barcode_spot_pool) == 2  # cell_typeA, cell_typeB
            assert "cell_typeA" in manager.barcode_spot_pool
            assert "cell_typeB" in manager.barcode_spot_pool
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

            assert len(manager.barcode_pool) == 2  # ACTG, TGCA
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


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
