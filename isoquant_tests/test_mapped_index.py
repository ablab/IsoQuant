############################################################################
# Copyright (c) 2025-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import json
import os

import numpy
import pytest

from isoquant_lib.barcode_calling.common import str_to_2bit
from isoquant_lib.barcode_calling.indexers import (
    MappedIndexInfo,
    MappedSparseAnchorIndexer,
    SharedMemorySparseAnchorIndexer,
)


BARCODES = [
    "ACTGACTGACTGACTGACTGACTGA",
    "TGCATGCATGCATGCATGCATGCAT",
    "GGGGGGGGGGGGGGGGGGGGGGGGG",
    "AAAAAAAAAAAAAAAAAAAAAAAAA",
    "CGTACGTACGTACGTACGTACGTAC",
]


def _bin_seqs(strs):
    return numpy.array([str_to_2bit(s) for s in strs], dtype=numpy.uint64)


class TestMappedSparseAnchorIndexer:

    def test_init(self, tmp_path):
        indexer = MappedSparseAnchorIndexer(
            _bin_seqs(BARCODES[:1]),
            kmer_size=14, seq_len=25, cache_dir=str(tmp_path),
        )
        assert indexer.total_sequences == 1
        assert indexer.k == 14
        assert not indexer.empty()
        # Inverted list holds 2 postings per barcode (first anchor is in the
        # permutation table, not the inverted list).
        assert indexer.index_size == 2
        # All five files exist (Step B adds first_anchor_starts.bin).
        for fname in (
            MappedSparseAnchorIndexer.BARCODES_FILENAME,
            MappedSparseAnchorIndexer.FIRST_ANCHOR_FILENAME,
            MappedSparseAnchorIndexer.RANGES_FILENAME,
            MappedSparseAnchorIndexer.INDEX_FILENAME,
            MappedSparseAnchorIndexer.META_FILENAME,
        ):
            assert os.path.exists(os.path.join(str(tmp_path), fname))

    def test_meta_contents(self, tmp_path):
        indexer = MappedSparseAnchorIndexer(
            _bin_seqs(BARCODES[:3]),
            kmer_size=14, seq_len=25, cache_dir=str(tmp_path),
        )
        with open(os.path.join(str(tmp_path), "meta.json")) as fh:
            meta = json.load(fh)
        assert meta["layout_version"] == MappedSparseAnchorIndexer.LAYOUT_VERSION
        assert meta["barcode_count"] == 3
        assert meta["kmer_size"] == 14
        assert meta["seq_len"] == 25
        # Step B: inverted list holds only last + minimizer postings → 2 * N.
        assert meta["index_size"] == 6
        assert meta.get("has_first_anchor_table") is True
        # File sizes match meta.
        assert os.path.getsize(
            os.path.join(str(tmp_path), "known_bin_seq.bin")) == 3 * 8
        assert os.path.getsize(
            os.path.join(str(tmp_path), "inverted_list.bin")) == 6 * 5
        # 4^14 + 1 cumulative offsets, uint64, for both range tables.
        assert os.path.getsize(
            os.path.join(str(tmp_path), "index_ranges.bin")) == (4**14 + 1) * 8
        assert os.path.getsize(
            os.path.join(str(tmp_path), "first_anchor_starts.bin")
        ) == (4**14 + 1) * 8
        # Persistence: directory should NOT be unlinked while indexer alive,
        # and we explicitly opt out of cleanup by keeping the indexer reachable.
        del indexer
        # Default persistent=False, so cleanup happened.
        assert not os.path.exists(str(tmp_path)) or not os.listdir(str(tmp_path))

    def test_empty(self, tmp_path):
        indexer = MappedSparseAnchorIndexer(
            [], kmer_size=14, seq_len=25, cache_dir=str(tmp_path),
        )
        assert indexer.empty()
        assert indexer.total_sequences == 0
        assert indexer.index_size == 0
        assert indexer.get_occurrences("ACTGACTGACTGACTGACTGACTGA") == []

    def test_seq_len_too_short(self, tmp_path):
        with pytest.raises(ValueError):
            MappedSparseAnchorIndexer(
                [], kmer_size=14, seq_len=15, cache_dir=str(tmp_path),
            )

    def test_exact_match(self, tmp_path):
        seq = "ACTGACTGACTGACTGACTGACTGA"
        indexer = MappedSparseAnchorIndexer(
            _bin_seqs([seq]),
            kmer_size=14, seq_len=25, cache_dir=str(tmp_path),
        )
        results = indexer.get_occurrences(seq)
        assert len(results) == 1
        assert results[0][0] == seq
        assert results[0][1] == 1000 * 3 - 0

    def test_recovery_from_front_error(self, tmp_path):
        barcode = "ACTGACTGACTGCTGACTGACTGAC"
        query = "ACAGACTGACTGCTGACTGACTGAC"
        indexer = MappedSparseAnchorIndexer(
            _bin_seqs([barcode]),
            kmer_size=14, seq_len=25, cache_dir=str(tmp_path),
        )
        results = indexer.get_occurrences(query, max_hits=5, min_kmers=1)
        assert any(r[0] == barcode for r in results)

    def test_recovery_from_back_error(self, tmp_path):
        barcode = "ACTGACTGACTGCTGACTGACTGAC"
        query = "ACTGACTGACTGCTGACTGAATGAC"
        indexer = MappedSparseAnchorIndexer(
            _bin_seqs([barcode]),
            kmer_size=14, seq_len=25, cache_dir=str(tmp_path),
        )
        results = indexer.get_occurrences(query, max_hits=5, min_kmers=1)
        assert any(r[0] == barcode for r in results)

    def test_no_recovery_in_overlap(self, tmp_path):
        barcode = "ACTGACTGACTGCTGACTGACTGAC"
        q = list(barcode)
        q[12] = "A" if q[12] != "A" else "C"
        query = "".join(q)
        indexer = MappedSparseAnchorIndexer(
            _bin_seqs([barcode]),
            kmer_size=14, seq_len=25, cache_dir=str(tmp_path),
        )
        results = indexer.get_occurrences(query, max_hits=5, min_kmers=1)
        assert all(r[0] != barcode for r in results)

    def test_pickling_sharable_info(self, tmp_path):
        barcodes = BARCODES[:1]
        indexer = MappedSparseAnchorIndexer(
            _bin_seqs(barcodes),
            kmer_size=14, seq_len=25, cache_dir=str(tmp_path),
            persistent=True,  # keep files alive across the worker simulation
        )
        info = indexer.get_sharable_info()
        assert isinstance(info, MappedIndexInfo)
        assert info.barcode_count == 1
        assert info.kmer_size == 14
        assert info.seq_len == 25
        assert info.cache_dir == str(tmp_path)
        assert info.layout_version == MappedSparseAnchorIndexer.LAYOUT_VERSION

        # Round-trip through getstate/setstate (pickling for multiprocessing).
        state = info.__getstate__()
        recovered = MappedIndexInfo.__new__(MappedIndexInfo)
        recovered.__setstate__(state)

        # Worker-side reconstruction in the same process.
        worker = MappedSparseAnchorIndexer.from_sharable_info(recovered)
        query = "ACTGACTGACTGACTGACTGACTGA"
        assert (
            indexer.get_occurrences(query)
            == worker.get_occurrences(query)
        )
        # Workers must not unlink the cache directory.
        del worker
        assert os.path.exists(str(tmp_path))

    def test_layout_version_mismatch(self, tmp_path):
        info = MappedIndexInfo(
            cache_dir=str(tmp_path),
            barcode_count=1, kmer_size=14, seq_len=25, index_size=3,
            layout_version=MappedSparseAnchorIndexer.LAYOUT_VERSION + 999,
        )
        with pytest.raises(ValueError):
            MappedSparseAnchorIndexer.from_sharable_info(info)

    def test_max_hits(self, tmp_path):
        barcodes = [
            "ACTGACTGACTGCTGACTGACTGAC",
            "ACTGACTGACTGCTGACTGACTGAA",
            "ACTGACTGACTGCTGACTGACTGAG",
            "ACTGACTGACTGCTGACTGACTGAT",
        ]
        indexer = MappedSparseAnchorIndexer(
            _bin_seqs(barcodes),
            kmer_size=14, seq_len=25, cache_dir=str(tmp_path),
        )
        results = indexer.get_occurrences(
            "ACTGACTGACTGCTGACTGACTGAC", max_hits=2, min_kmers=1,
        )
        assert len(results) <= 2

    def test_min_kmers_filter(self, tmp_path):
        indexer = MappedSparseAnchorIndexer(
            _bin_seqs(BARCODES[:2]),
            kmer_size=14, seq_len=25, cache_dir=str(tmp_path),
        )
        results = indexer.get_occurrences(
            "CCCCCCCCCCCCCCCCCCCCCCCCC", min_kmers=2,
        )
        assert results == []

    def test_ignore_equal(self, tmp_path):
        indexer = MappedSparseAnchorIndexer(
            _bin_seqs(BARCODES[:2]),
            kmer_size=14, seq_len=25, cache_dir=str(tmp_path),
        )
        query = "ACTGACTGACTGACTGACTGACTGA"
        results = indexer.get_occurrences(query, ignore_equal=True, min_kmers=1)
        for seq, _, _ in results:
            assert seq != query

    def test_score_ordering(self, tmp_path):
        exact = "ACTGACTGACTGCTGACTGACTGAC"
        partial = "ACTGACTGACTGCTAAAAAAAAAAA"
        indexer = MappedSparseAnchorIndexer(
            _bin_seqs([exact, partial]),
            kmer_size=14, seq_len=25, cache_dir=str(tmp_path),
        )
        results = indexer.get_occurrences(exact, max_hits=10, min_kmers=1)
        assert results
        assert results[0][0] == exact

    def test_low_mem_parity(self, tmp_path):
        # Same input, low_mem on vs off — index files and queries must match.
        seqs = _bin_seqs(BARCODES)
        dir_a = tmp_path / "a"
        dir_b = tmp_path / "b"
        regular = MappedSparseAnchorIndexer(
            seqs, kmer_size=14, seq_len=25,
            cache_dir=str(dir_a), low_mem=False,
        )
        low_mem = MappedSparseAnchorIndexer(
            seqs, kmer_size=14, seq_len=25,
            cache_dir=str(dir_b), low_mem=True,
        )
        assert regular.index_size == low_mem.index_size
        assert numpy.array_equal(regular.index_ranges, low_mem.index_ranges)
        for k in range(len(regular.index_ranges) - 1):
            start = int(regular.index_ranges[k])
            end = int(regular.index_ranges[k + 1])
            if start == end:
                continue
            reg_ids = sorted(
                int.from_bytes(bytes(regular.index[i * 5:(i + 1) * 5]), 'little')
                for i in range(start, end))
            low_ids = sorted(
                int.from_bytes(bytes(low_mem.index[i * 5:(i + 1) * 5]), 'little')
                for i in range(start, end))
            assert reg_ids == low_ids
        for q in BARCODES:
            r_seqs = {s for s, _, _ in regular.get_occurrences(q, min_kmers=1)}
            l_seqs = {s for s, _, _ in low_mem.get_occurrences(q, min_kmers=1)}
            assert r_seqs == l_seqs

    def test_parity_with_shared_memory(self, tmp_path):
        # Mapped and SHM variants must agree on query results.
        seqs = _bin_seqs(BARCODES)
        mapped = MappedSparseAnchorIndexer(
            seqs, kmer_size=14, seq_len=25,
            cache_dir=str(tmp_path), low_mem=True,
        )
        shm = SharedMemorySparseAnchorIndexer(
            seqs, kmer_size=14, seq_len=25, low_mem=True,
        )
        try:
            for q in BARCODES:
                mapped_results = mapped.get_occurrences(q, min_kmers=1)
                shm_results = shm.get_occurrences(q, min_kmers=1)
                m_seqs = {s for s, _, _ in mapped_results}
                s_seqs = {s for s, _, _ in shm_results}
                assert m_seqs == s_seqs
                # Scores must match too (tier ordering is deterministic).
                assert sorted(mapped_results) == sorted(shm_results)
        finally:
            del shm

    def test_persistent_keeps_files(self, tmp_path):
        seqs = _bin_seqs(BARCODES[:2])
        indexer = MappedSparseAnchorIndexer(
            seqs, kmer_size=14, seq_len=25,
            cache_dir=str(tmp_path), persistent=True,
        )
        del indexer
        # Files survived.
        assert os.path.exists(os.path.join(str(tmp_path), "meta.json"))
        assert os.path.exists(
            os.path.join(str(tmp_path), "known_bin_seq.bin"))
        assert os.path.exists(
            os.path.join(str(tmp_path), "first_anchor_starts.bin"))

    def test_layout_version_is_two(self):
        # Sanity check: Step B bumped the layout version.
        assert MappedSparseAnchorIndexer.LAYOUT_VERSION == 2

    def test_first_anchor_table_sorted_by_first_anchor(self, tmp_path):
        # After Step B, known_bin_seq is stable-sorted by the first anchor
        # k-mer, so first_anchor_starts gives the slice bounds directly.
        seqs = _bin_seqs(BARCODES)
        indexer = MappedSparseAnchorIndexer(
            seqs, kmer_size=14, seq_len=25, cache_dir=str(tmp_path),
        )
        first_shift = (25 - 14) * 2
        k_mask = (1 << 28) - 1
        stored = numpy.asarray(indexer.known_bin_seq, dtype=numpy.uint64)
        first_anchors = (stored >> numpy.uint64(first_shift)) & numpy.uint64(k_mask)
        # The first anchors of the stored barcodes must be non-decreasing.
        assert numpy.all(first_anchors[:-1] <= first_anchors[1:])
        # And first_anchor_starts must be a cumulative count of first anchors.
        for a in numpy.unique(first_anchors):
            start = int(indexer.first_anchor_starts[int(a)])
            end = int(indexer.first_anchor_starts[int(a) + 1])
            assert end - start == int((first_anchors == a).sum())
            # Every entry in the slice has first anchor == a.
            assert numpy.all(first_anchors[start:end] == a)

    def test_first_anchor_table_collision(self, tmp_path):
        # Construct multiple barcodes that share the same first anchor
        # (positions 0..13). They must all land in a contiguous slice of the
        # sorted known_bin_seq array.
        shared_prefix = "ACTGACTGACTGAC"  # first 14 nt → shared anchor
        barcodes = [
            shared_prefix + "TGACTGACTGA",
            shared_prefix + "AAACTGACTGA",
            shared_prefix + "GGGCTGACTGA",
            shared_prefix + "CCCCTGACTGA",
        ]
        # And one that has a different first anchor.
        outsider = "TTTTACTGACTGACTGACTGACTGA"
        barcodes.append(outsider)

        seqs = _bin_seqs(barcodes)
        indexer = MappedSparseAnchorIndexer(
            seqs, kmer_size=14, seq_len=25, cache_dir=str(tmp_path),
        )
        # Each barcode must be findable via an exact-match query.
        for q in barcodes:
            results = indexer.get_occurrences(q, max_hits=10, min_kmers=1)
            assert any(r[0] == q for r in results), q
            # Exact match scores 3000 (3 anchors, 0 hamming).
            assert any(r[0] == q and r[1] == 3000 for r in results), q

    def test_first_anchor_table_distinct(self, tmp_path):
        # Barcodes with totally distinct first anchors — each bucket should
        # hold exactly one entry, and queries return the matching barcode.
        barcodes = [
            "AAAAAAAAAAAAAA" + "ACTGACTGACT",
            "CCCCCCCCCCCCCC" + "ACTGACTGACT",
            "GGGGGGGGGGGGGG" + "ACTGACTGACT",
            "TTTTTTTTTTTTTT" + "ACTGACTGACT",
        ]
        seqs = _bin_seqs(barcodes)
        indexer = MappedSparseAnchorIndexer(
            seqs, kmer_size=14, seq_len=25, cache_dir=str(tmp_path),
        )
        first_shift = (25 - 14) * 2
        k_mask = (1 << 28) - 1
        for q in barcodes:
            anchor = (int(_bin_seqs([q])[0]) >> first_shift) & k_mask
            start = int(indexer.first_anchor_starts[anchor])
            end = int(indexer.first_anchor_starts[anchor + 1])
            assert end - start == 1
            # Exact-match query must return this barcode.
            results = indexer.get_occurrences(q, max_hits=5, min_kmers=1)
            assert any(r[0] == q and r[1] == 3000 for r in results)

    def test_inverted_list_shrunk_vs_three_anchors(self, tmp_path):
        # Without Step B the inverted list would be 3 * N entries; Step B
        # collapses the first-anchor third into a permutation table, so the
        # on-disk file is exactly 2/3 the size of the pre-Step-B layout.
        seqs = _bin_seqs(BARCODES)
        indexer = MappedSparseAnchorIndexer(
            seqs, kmer_size=14, seq_len=25, cache_dir=str(tmp_path),
        )
        expected = 2 * len(BARCODES) * 5  # 2 anchors * 5 bytes / entry
        actual = os.path.getsize(
            os.path.join(str(tmp_path), "inverted_list.bin"))
        assert actual == expected
        # Keep `indexer` reachable until after the size check.
        assert indexer.total_sequences == len(BARCODES)

    def test_step_b_parity_with_shared_memory_on_collision(self, tmp_path):
        # Even when many barcodes collide on the first anchor (the case where
        # Step B's permutation table is most loaded), query results must still
        # agree with the SHM indexer byte-for-byte (set of returned barcodes
        # plus their scores).
        shared = "ACTGACTGACTGAC"
        barcodes = [
            shared + "AAAAAAAAAAA",
            shared + "CCCCCCCCCCC",
            shared + "GGGGGGGGGGG",
            shared + "TTTTTTTTTTT",
            "GCGCGCGCGCGCGC" + "ACTGACTGACT",
        ]
        seqs = _bin_seqs(barcodes)
        mapped = MappedSparseAnchorIndexer(
            seqs, kmer_size=14, seq_len=25,
            cache_dir=str(tmp_path), low_mem=True,
        )
        shm = SharedMemorySparseAnchorIndexer(
            seqs, kmer_size=14, seq_len=25, low_mem=True,
        )
        try:
            for q in barcodes:
                mapped_results = mapped.get_occurrences(q, min_kmers=1)
                shm_results = shm.get_occurrences(q, min_kmers=1)
                assert sorted(mapped_results) == sorted(shm_results)
        finally:
            del shm


class TestMappedIndexInfo:

    def test_getstate_setstate(self):
        info = MappedIndexInfo(
            cache_dir="/some/where",
            barcode_count=100,
            kmer_size=14,
            seq_len=25,
            index_size=300,
            layout_version=MappedSparseAnchorIndexer.LAYOUT_VERSION,
        )
        state = info.__getstate__()
        assert isinstance(state, tuple)
        assert len(state) == 6
        new_info = MappedIndexInfo.__new__(MappedIndexInfo)
        new_info.__setstate__(state)
        assert new_info.cache_dir == "/some/where"
        assert new_info.barcode_count == 100
        assert new_info.kmer_size == 14
        assert new_info.seq_len == 25
        assert new_info.index_size == 300
        assert new_info.layout_version == MappedSparseAnchorIndexer.LAYOUT_VERSION


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
