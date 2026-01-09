############################################################################
# Copyright (c) 2025 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import pytest
from src.barcode_calling.shared_mem_index import (
    SharedMemoryArray2BitKmerIndexer,
    SharedMemoryIndexInfo,
)
from src.barcode_calling.kmer_indexer import Array2BitKmerIndexer
from src.barcode_calling.common import str_to_2bit


class TestSharedMemoryArray2BitKmerIndexer:
    """Test shared memory k-mer indexer."""

    def test_init(self):
        """Test basic initialization."""
        barcodes = ["ACTGACTGACTGACTGACTGACTGA"]  # 25-mer
        bin_seqs = [str_to_2bit(b) for b in barcodes]
        indexer = SharedMemoryArray2BitKmerIndexer(bin_seqs, kmer_size=10, seq_len=25)

        assert indexer.total_sequences == 1
        assert indexer.k == 10
        assert not indexer.empty()

    def test_empty(self):
        """Test empty index."""
        indexer = SharedMemoryArray2BitKmerIndexer([], kmer_size=10, seq_len=25)
        assert indexer.empty()
        assert indexer.total_sequences == 0
        assert indexer.index_size == 0

    def test_exact_match(self):
        """Test finding exact match."""
        seq = "ACTGACTGACTGACTGACTGACTGA"
        bin_seqs = [str_to_2bit(seq)]
        indexer = SharedMemoryArray2BitKmerIndexer(bin_seqs, kmer_size=10, seq_len=25)

        results = indexer.get_occurrences(seq)
        assert len(results) > 0
        assert results[0][0] == seq

    def test_similar_match(self):
        """Test finding similar sequence (1 mismatch)."""
        barcode = "ACTGACTGACTGACTGACTGACTGA"
        query = "ACTGACTGACTGACTGACTGACTGT"  # 1 mismatch at end
        bin_seqs = [str_to_2bit(barcode)]
        indexer = SharedMemoryArray2BitKmerIndexer(bin_seqs, kmer_size=10, seq_len=25)

        results = indexer.get_occurrences(query)
        assert len(results) > 0
        assert results[0][0] == barcode

    def test_consistency_with_array2bit(self):
        """Test results match Array2BitKmerIndexer."""
        barcodes = [
            "ACTGACTGACTGACTGACTGACTGA",
            "TGCATGCATGCATGCATGCATGCAT",
            "GGGGGGGGGGGGGGGGGGGGGGGGG"
        ]
        bin_seqs = [str_to_2bit(b) for b in barcodes]

        shared = SharedMemoryArray2BitKmerIndexer(bin_seqs, kmer_size=10, seq_len=25)
        regular = Array2BitKmerIndexer(bin_seqs, kmer_size=10, seq_len=25)

        query = "ACTGACTGACTGACTGACTGACTGT"
        shared_results = shared.get_occurrences(query)
        regular_results = regular.get_occurrences(query)

        # Should find same sequences
        shared_seqs = set(seq for seq, _, _ in shared_results)
        regular_seqs = set(seq for seq, _, _ in regular_results)
        assert shared_seqs == regular_seqs

    def test_pickling_sharable_info(self):
        """Test serialization via get_sharable_info/from_sharable_info."""
        barcodes = ["ACTGACTGACTGACTGACTGACTGA"]
        bin_seqs = [str_to_2bit(b) for b in barcodes]
        indexer = SharedMemoryArray2BitKmerIndexer(bin_seqs, kmer_size=10, seq_len=25)

        # Get sharable info (simulates pickling for multiprocessing)
        info = indexer.get_sharable_info()

        assert isinstance(info, SharedMemoryIndexInfo)
        assert info.barcode_count == 1
        assert info.kmer_size == 10
        assert info.seq_len == 25

        # Reconstruct from info (simulates worker process)
        reconstructed = SharedMemoryArray2BitKmerIndexer.from_sharable_info(info)

        # Both should return same results
        query = "ACTGACTGACTGACTGACTGACTGA"
        orig_results = indexer.get_occurrences(query)
        recon_results = reconstructed.get_occurrences(query)

        assert orig_results == recon_results

    def test_multiple_barcodes(self):
        """Test with multiple similar barcodes."""
        barcodes = [
            "ACTGACTGACTGACTGACTGACTGA",
            "ACTGACTGACTGACTGACTGACTGC",
            "ACTGACTGACTGACTGACTGACTGG",
        ]
        bin_seqs = [str_to_2bit(b) for b in barcodes]
        indexer = SharedMemoryArray2BitKmerIndexer(bin_seqs, kmer_size=10, seq_len=25)

        results = indexer.get_occurrences("ACTGACTGACTGACTGACTGACTGA", max_hits=2)
        assert len(results) <= 2

    def test_index_size_calculation(self):
        """Test that index_size matches actual entries."""
        barcodes = ["ACTGACTGACTGACTGACTGACTGA", "TGCATGCATGCATGCATGCATGCAT"]
        bin_seqs = [str_to_2bit(b) for b in barcodes]
        indexer = SharedMemoryArray2BitKmerIndexer(bin_seqs, kmer_size=10, seq_len=25)

        # Each barcode generates 16 k-mers (25 - 10 + 1)
        expected_entries = 2 * 16
        assert indexer.index_size == expected_entries

    def test_max_hits(self):
        """Test maximum hits limit."""
        barcodes = [
            "ACTGACTGACTGACTGACTGACTGA",
            "ACTGACTGACTGACTGACTGACTGC",
            "ACTGACTGACTGACTGACTGACTGG",
            "ACTGACTGACTGACTGACTGACTGT",
        ]
        bin_seqs = [str_to_2bit(b) for b in barcodes]
        indexer = SharedMemoryArray2BitKmerIndexer(bin_seqs, kmer_size=10, seq_len=25)

        results = indexer.get_occurrences("ACTGACTGACTGACTGACTGACTGA", max_hits=2)
        assert len(results) <= 2

    def test_min_kmers_filter(self):
        """Test minimum k-mer threshold."""
        barcodes = ["ACTGACTGACTGACTGACTGACTGA", "TGCATGCATGCATGCATGCATGCAT"]
        bin_seqs = [str_to_2bit(b) for b in barcodes]
        indexer = SharedMemoryArray2BitKmerIndexer(bin_seqs, kmer_size=10, seq_len=25)

        # Very different sequence - should not match with high min_kmers
        results = indexer.get_occurrences("CCCCCCCCCCCCCCCCCCCCCCCCC", min_kmers=10)
        assert len(results) == 0

    def test_ignore_equal(self):
        """Test ignoring exact matches."""
        barcodes = ["ACTGACTGACTGACTGACTGACTGA", "TGCATGCATGCATGCATGCATGCAT"]
        bin_seqs = [str_to_2bit(b) for b in barcodes]
        indexer = SharedMemoryArray2BitKmerIndexer(bin_seqs, kmer_size=10, seq_len=25)

        query = "ACTGACTGACTGACTGACTGACTGA"
        results = indexer.get_occurrences(query, ignore_equal=True)

        # Should not return the exact match
        for seq, _, _ in results:
            assert seq != query


class TestSharedMemoryIndexInfo:
    """Test SharedMemoryIndexInfo serialization."""

    def test_getstate_setstate(self):
        """Test pickling of SharedMemoryIndexInfo."""
        info = SharedMemoryIndexInfo(
            barcode_count=100,
            kmer_size=10,
            seq_len=25,
            index_size=1600,
            barcodes_sm_name="test_barcodes",
            index_sm_name="test_index",
            index_range_sm_name="test_ranges"
        )

        # Simulate pickle
        state = info.__getstate__()
        assert isinstance(state, tuple)
        assert len(state) == 7

        # Simulate unpickle
        new_info = SharedMemoryIndexInfo.__new__(SharedMemoryIndexInfo)
        new_info.__setstate__(state)

        assert new_info.barcode_count == 100
        assert new_info.kmer_size == 10
        assert new_info.seq_len == 25
        assert new_info.index_size == 1600
        assert new_info.barcodes_sm_name == "test_barcodes"
        assert new_info.index_sm_name == "test_index"
        assert new_info.index_range_sm_name == "test_ranges"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
