############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import numpy
import pytest
from isoquant_lib.barcode_calling.indexers import (
    SharedMemoryArray2BitKmerIndexer,
    SharedMemoryIndexInfo,
    SharedMemorySparseAnchorIndexer,
    Array2BitKmerIndexer,
)
from isoquant_lib.barcode_calling.common import str_to_2bit


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


class TestSharedMemorySparseAnchorIndexer:
    """Test sparse 3-anchor (first/last/minimizer) shared-memory indexer."""

    def test_init(self):
        barcodes = ["ACTGACTGACTGACTGACTGACTGA"]
        bin_seqs = [str_to_2bit(b) for b in barcodes]
        indexer = SharedMemorySparseAnchorIndexer(bin_seqs, kmer_size=14, seq_len=25)

        assert indexer.total_sequences == 1
        assert indexer.k == 14
        assert not indexer.empty()
        # Exactly 3 anchors per sequence.
        assert indexer.index_size == 3

    def test_empty(self):
        indexer = SharedMemorySparseAnchorIndexer([], kmer_size=14, seq_len=25)
        assert indexer.empty()
        assert indexer.total_sequences == 0
        assert indexer.index_size == 0
        # An empty indexer should return [] for any query.
        assert indexer.get_occurrences("ACTGACTGACTGACTGACTGACTGA") == []

    def test_seq_len_too_short(self):
        # seq_len must leave at least one internal window for the minimizer.
        # For k=14 we need seq_len >= 16.
        with pytest.raises(ValueError):
            SharedMemorySparseAnchorIndexer([], kmer_size=14, seq_len=15)

    def test_exact_match(self):
        seq = "ACTGACTGACTGACTGACTGACTGA"
        bin_seqs = [str_to_2bit(seq)]
        indexer = SharedMemorySparseAnchorIndexer(bin_seqs, kmer_size=14, seq_len=25)

        results = indexer.get_occurrences(seq)
        assert len(results) == 1
        assert results[0][0] == seq
        # Exact match scores at most 3 anchors hit, Hamming 0.
        assert results[0][1] == 1000 * 3 - 0

    def test_recovery_from_front_error(self):
        # An error in the first 14-mer should still be recoverable via the
        # last anchor (positions 11..24 are untouched).
        barcode = "ACTGACTGACTGCTGACTGACTGAC"
        # Mutate position 2 (still in the front 0..10 region, outside the
        # 11..13 overlap).
        query = "ACAGACTGACTGCTGACTGACTGAC"
        assert barcode != query
        bin_seqs = [str_to_2bit(barcode)]
        indexer = SharedMemorySparseAnchorIndexer(bin_seqs, kmer_size=14, seq_len=25)

        results = indexer.get_occurrences(query, max_hits=5, min_kmers=1)
        # The last anchor covers positions 11..24 and matches exactly.
        assert any(r[0] == barcode for r in results)

    def test_recovery_from_back_error(self):
        barcode = "ACTGACTGACTGCTGACTGACTGAC"
        # Mutate position 20 (in the back region, outside the 11..13 overlap).
        query = "ACTGACTGACTGCTGACTGAATGAC"
        assert barcode != query
        bin_seqs = [str_to_2bit(barcode)]
        indexer = SharedMemorySparseAnchorIndexer(bin_seqs, kmer_size=14, seq_len=25)

        results = indexer.get_occurrences(query, max_hits=5, min_kmers=1)
        assert any(r[0] == barcode for r in results)

    def test_no_recovery_in_overlap(self):
        # Errors at positions 11..13 (the 3-base first/last overlap) break
        # every 14-mer in the barcode, so all 3 anchors fail. This is a
        # structural limit of k=14 in 25-bp barcodes (true for the dense
        # indexer too).
        barcode = "ACTGACTGACTGCTGACTGACTGAC"
        query_list = list(barcode)
        # Flip position 12.
        query_list[12] = "A" if query_list[12] != "A" else "C"
        query = "".join(query_list)
        bin_seqs = [str_to_2bit(barcode)]
        indexer = SharedMemorySparseAnchorIndexer(bin_seqs, kmer_size=14, seq_len=25)

        results = indexer.get_occurrences(query, max_hits=5, min_kmers=1)
        assert all(r[0] != barcode for r in results)

    def test_minimizer_distinct_from_terminals(self):
        # The minimizer is drawn from internal windows (p in {1..10}) so its
        # encoded value cannot equal the first or last anchor's encoded slice
        # for typical barcodes. We assert by inspecting the anchors directly.
        from isoquant_lib.barcode_calling.indexers.shared_memory import (
            _compute_anchors,
        )
        import numpy

        barcode = "ACTGACTGACTGCTGACTGACTGAC"
        seq = numpy.array([str_to_2bit(barcode)], dtype=numpy.uint64)
        first_shift = numpy.uint64((25 - 14) * 2)
        last_shift = numpy.uint64(0)
        k_mask = numpy.uint64((1 << 28) - 1)
        mid_shifts = numpy.array(
            [(25 - 14 - p) * 2 for p in range(1, 25 - 14)], dtype=numpy.uint64,
        )
        anchors = numpy.empty((1, 3), dtype=numpy.uint64)
        _compute_anchors(seq, first_shift, last_shift, k_mask, mid_shifts, anchors)
        first14, last14, min14 = (
            int(anchors[0, 0]), int(anchors[0, 1]), int(anchors[0, 2]),
        )
        # No structural collision: distinct slice content.
        assert min14 != first14
        assert min14 != last14

    def test_pickling_sharable_info(self):
        barcodes = ["ACTGACTGACTGACTGACTGACTGA"]
        bin_seqs = [str_to_2bit(b) for b in barcodes]
        indexer = SharedMemorySparseAnchorIndexer(bin_seqs, kmer_size=14, seq_len=25)

        info = indexer.get_sharable_info()
        assert isinstance(info, SharedMemoryIndexInfo)
        assert info.barcode_count == 1
        assert info.kmer_size == 14
        assert info.seq_len == 25

        reconstructed = SharedMemorySparseAnchorIndexer.from_sharable_info(info)
        query = "ACTGACTGACTGACTGACTGACTGA"
        assert (
            indexer.get_occurrences(query)
            == reconstructed.get_occurrences(query)
        )

    def test_max_hits(self):
        barcodes = [
            "ACTGACTGACTGCTGACTGACTGAC",
            "ACTGACTGACTGCTGACTGACTGAA",
            "ACTGACTGACTGCTGACTGACTGAG",
            "ACTGACTGACTGCTGACTGACTGAT",
        ]
        bin_seqs = [str_to_2bit(b) for b in barcodes]
        indexer = SharedMemorySparseAnchorIndexer(bin_seqs, kmer_size=14, seq_len=25)

        results = indexer.get_occurrences(
            "ACTGACTGACTGCTGACTGACTGAC", max_hits=2, min_kmers=1,
        )
        assert len(results) <= 2

    def test_min_kmers_filter(self):
        # A sequence with no shared anchors should not return any hits.
        barcodes = ["ACTGACTGACTGACTGACTGACTGA", "TGCATGCATGCATGCATGCATGCAT"]
        bin_seqs = [str_to_2bit(b) for b in barcodes]
        indexer = SharedMemorySparseAnchorIndexer(bin_seqs, kmer_size=14, seq_len=25)

        results = indexer.get_occurrences(
            "CCCCCCCCCCCCCCCCCCCCCCCCC", min_kmers=2,
        )
        assert results == []

    def test_ignore_equal(self):
        barcodes = ["ACTGACTGACTGACTGACTGACTGA", "TGCATGCATGCATGCATGCATGCAT"]
        bin_seqs = [str_to_2bit(b) for b in barcodes]
        indexer = SharedMemorySparseAnchorIndexer(bin_seqs, kmer_size=14, seq_len=25)

        query = "ACTGACTGACTGACTGACTGACTGA"
        results = indexer.get_occurrences(query, ignore_equal=True, min_kmers=1)
        for seq, _, _ in results:
            assert seq != query

    def test_score_ordering(self):
        # First-and-last (Tier A) candidates must outrank single-anchor
        # (Tier C) candidates regardless of Hamming.
        exact = "ACTGACTGACTGCTGACTGACTGAC"
        # A second barcode that shares only the first 14-mer with the query
        # but differs everywhere else.
        partial = "ACTGACTGACTGCTAAAAAAAAAAA"
        bin_seqs = [str_to_2bit(exact), str_to_2bit(partial)]
        indexer = SharedMemorySparseAnchorIndexer(bin_seqs, kmer_size=14, seq_len=25)

        results = indexer.get_occurrences(exact, max_hits=10, min_kmers=1)
        assert results, "expected at least one match"
        # The exact-match (Tier A) entry must be first.
        assert results[0][0] == exact

    def test_consistency_with_dense_for_exact(self):
        # For exact queries, sparse and dense must agree on the matching seq.
        barcodes = [
            "ACTGACTGACTGACTGACTGACTGA",
            "TGCATGCATGCATGCATGCATGCAT",
            "GGGGGGGGGGGGGGGGGGGGGGGGG",
        ]
        bin_seqs = [str_to_2bit(b) for b in barcodes]
        sparse = SharedMemorySparseAnchorIndexer(bin_seqs, kmer_size=14, seq_len=25)
        dense = SharedMemoryArray2BitKmerIndexer(bin_seqs, kmer_size=14, seq_len=25)

        for q in barcodes:
            sparse_seqs = {seq for seq, _, _ in sparse.get_occurrences(q, min_kmers=1)}
            dense_seqs = {seq for seq, _, _ in dense.get_occurrences(q, min_kmers=1)}
            assert q in sparse_seqs
            assert q in dense_seqs

    def test_low_mem_parity(self):
        # The materialized and fused build paths must produce identical
        # index_ranges and index byte content given the same input.
        barcodes = [
            "ACTGACTGACTGACTGACTGACTGA",
            "TGCATGCATGCATGCATGCATGCAT",
            "GGGGGGGGGGGGGGGGGGGGGGGGG",
            "AAAAAAAAAAAAAAAAAAAAAAAAA",
            "CGTACGTACGTACGTACGTACGTAC",
        ]
        bin_seqs = [str_to_2bit(b) for b in barcodes]

        regular = SharedMemorySparseAnchorIndexer(
            bin_seqs, kmer_size=14, seq_len=25, low_mem=False)
        low_mem = SharedMemorySparseAnchorIndexer(
            bin_seqs, kmer_size=14, seq_len=25, low_mem=True)

        assert regular.index_size == low_mem.index_size
        assert numpy.array_equal(regular.index_ranges, low_mem.index_ranges)
        # The inverted list contents per bucket must match as multisets.
        for k in range(len(regular.index_ranges) - 1):
            start = int(regular.index_ranges[k])
            end = int(regular.index_ranges[k + 1])
            if start == end:
                continue
            reg_bytes = bytes(regular.index[start * 5:end * 5])
            low_bytes = bytes(low_mem.index[start * 5:end * 5])
            reg_ids = sorted(
                int.from_bytes(reg_bytes[i * 5:(i + 1) * 5], 'little')
                for i in range(end - start))
            low_ids = sorted(
                int.from_bytes(low_bytes[i * 5:(i + 1) * 5], 'little')
                for i in range(end - start))
            assert reg_ids == low_ids

        # Both must give the same query results too.
        for q in barcodes:
            r_seqs = {s for s, _, _ in regular.get_occurrences(q, min_kmers=1)}
            l_seqs = {s for s, _, _ in low_mem.get_occurrences(q, min_kmers=1)}
            assert r_seqs == l_seqs


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
