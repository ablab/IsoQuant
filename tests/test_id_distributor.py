############################################################################
# Copyright (c) 2022-2024 University of Helsinki
# Copyright (c) 2019-2022 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import unittest
import threading
from unittest.mock import MagicMock, patch
from src.id_policy import (
    SimpleIDDistributor,
    ExcludingIdDistributor,
    FeatureIdStorage,
    AtomicIDDistributor
)


class TestSimpleIDDistributor(unittest.TestCase):
    def test_get_next_id(self):
        distributor = SimpleIDDistributor()
        self.assertEqual(distributor.increment(), 1)
        self.assertEqual(distributor.increment(), 2)
        self.assertEqual(distributor.increment(), 3)

    def test_initial_value(self):
        distributor = SimpleIDDistributor()
        self.assertEqual(distributor.value, 0)

    def test_sequential_increment(self):
        distributor = SimpleIDDistributor()
        # Generate 100 IDs
        ids = [distributor.increment() for _ in range(100)]
        # Should be sequential 1 to 100
        self.assertEqual(ids, list(range(1, 101)))


class TestExcludingIdDistributor(unittest.TestCase):
    def test_without_genedb(self):
        """Test ExcludingIdDistributor works without genedb."""
        distributor = ExcludingIdDistributor(None, "chr1")
        self.assertEqual(distributor.increment(), 1)
        self.assertEqual(distributor.increment(), 2)
        self.assertEqual(distributor.increment(), 3)

    def test_skips_forbidden_gene_ids(self):
        """Test that forbidden gene IDs from genedb are skipped."""
        # Create mock genedb
        mock_genedb = MagicMock()

        # Mock gene with novel gene prefix (novel_gene_5 per TranscriptNaming)
        mock_gene = MagicMock()
        mock_gene.id = "novel_gene_5"

        # Mock transcript with novel transcript prefix (transcript10.1 per TranscriptNaming)
        mock_transcript = MagicMock()
        mock_transcript.id = "transcript10.1"

        # region() returns genes and transcripts
        # Note: for transcripts, the code passes a tuple ("transcript", "mRNA")
        def mock_region(seqid, start, featuretype):
            if featuretype == "gene":
                return [mock_gene]
            elif featuretype == ("transcript", "mRNA"):
                return [mock_transcript]
            return []

        mock_genedb.region = mock_region

        distributor = ExcludingIdDistributor(mock_genedb, "chr1")

        # Should skip 5 and 10
        self.assertIn(5, distributor.forbidden_ids)
        self.assertIn(10, distributor.forbidden_ids)

        # Generate IDs - should skip 5 and 10
        ids = []
        for _ in range(12):
            ids.append(distributor.increment())

        self.assertNotIn(5, ids)
        self.assertNotIn(10, ids)
        # First few IDs should be 1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14
        self.assertEqual(ids[:4], [1, 2, 3, 4])
        self.assertEqual(ids[4], 6)  # Skipped 5
        self.assertEqual(ids[8], 11)  # Skipped 10

    def test_handles_malformed_gene_ids(self):
        """Test that malformed gene IDs don't crash."""
        mock_genedb = MagicMock()

        mock_gene1 = MagicMock()
        mock_gene1.id = "ISOQUANT_gene_abc"  # Not a number

        mock_gene2 = MagicMock()
        mock_gene2.id = "ISOQUANT_gene_"  # Missing number

        mock_gene3 = MagicMock()
        mock_gene3.id = "other_gene_123"  # Wrong prefix

        def mock_region(seqid, start, featuretype):
            if featuretype == "gene":
                return [mock_gene1, mock_gene2, mock_gene3]
            return []

        mock_genedb.region = mock_region

        # Should not raise exception
        distributor = ExcludingIdDistributor(mock_genedb, "chr1")

        # Forbidden IDs should be empty (all malformed)
        self.assertEqual(len(distributor.forbidden_ids), 0)


class TestFeatureIdStorage(unittest.TestCase):
    def test_without_genedb(self):
        """Test FeatureIdStorage works without genedb."""
        id_distributor = SimpleIDDistributor()
        storage = FeatureIdStorage(id_distributor)

        # First call for new feature returns int (from increment())
        id1 = storage.get_id("chr1", (100, 200), "+")
        self.assertEqual(id1, 1)  # Returns int from increment()

        id2 = storage.get_id("chr1", (300, 400), "+")
        self.assertEqual(id2, 2)  # Returns int from increment()

        # Second call for same feature returns string (from id_dict)
        id3 = storage.get_id("chr1", (100, 200), "+")
        self.assertEqual(id3, "chr1.1")  # Returns stored string

    def test_same_feature_cached_id(self):
        """Test that same feature coordinates return cached ID on second call."""
        id_distributor = SimpleIDDistributor()
        storage = FeatureIdStorage(id_distributor)

        # First call returns int
        id1 = storage.get_id("chr1", (100, 200), "+")
        self.assertEqual(id1, 1)

        # Second call returns cached string
        id2 = storage.get_id("chr1", (100, 200), "+")
        self.assertEqual(id2, "chr1.1")

    def test_different_strand_different_id(self):
        """Test that same coordinates but different strand get different IDs."""
        id_distributor = SimpleIDDistributor()
        storage = FeatureIdStorage(id_distributor)

        id1 = storage.get_id("chr1", (100, 200), "+")
        id2 = storage.get_id("chr1", (100, 200), "-")

        # Both are first calls, so both return int
        self.assertEqual(id1, 1)
        self.assertEqual(id2, 2)
        self.assertNotEqual(id1, id2)

    def test_different_chr_different_id(self):
        """Test that same coordinates on different chromosomes get different IDs."""
        id_distributor = SimpleIDDistributor()
        storage = FeatureIdStorage(id_distributor)

        id1 = storage.get_id("chr1", (100, 200), "+")
        id2 = storage.get_id("chr2", (100, 200), "+")

        # Both are first calls, so both return int
        self.assertEqual(id1, 1)
        self.assertEqual(id2, 2)
        self.assertNotEqual(id1, id2)

    def test_with_genedb_preloaded_exons(self):
        """Test that existing exon IDs from genedb are preserved."""
        mock_genedb = MagicMock()

        # Mock existing exon feature
        mock_exon = MagicMock()
        mock_exon.start = 100
        mock_exon.end = 200
        mock_exon.strand = "+"
        mock_exon.attributes = {"exon_id": ["existing_exon_1"]}

        def mock_region(seqid, start, featuretype):
            if featuretype == "exon":
                return [mock_exon]
            return []

        mock_genedb.region = mock_region

        id_distributor = SimpleIDDistributor()
        storage = FeatureIdStorage(id_distributor, mock_genedb, "chr1", "exon")

        # Same coordinates should return preloaded ID (already in id_dict)
        existing_id = storage.get_id("chr1", (100, 200), "+")
        self.assertEqual(existing_id, "existing_exon_1")

        # New coordinates should generate new ID (first call returns int)
        new_id = storage.get_id("chr1", (300, 400), "+")
        self.assertEqual(new_id, 1)  # Returns int from increment()


class TestAtomicIDDistributor(unittest.TestCase):
    def test_basic_increment(self):
        """Test basic increment functionality."""
        distributor = AtomicIDDistributor()
        self.assertEqual(distributor.increment(), 1)
        self.assertEqual(distributor.increment(), 2)
        self.assertEqual(distributor.increment(), 3)

    def test_thread_safety(self):
        """Test that AtomicIDDistributor is thread-safe."""
        distributor = AtomicIDDistributor()
        results = []
        num_threads = 10
        increments_per_thread = 100

        def worker():
            for _ in range(increments_per_thread):
                results.append(distributor.increment())

        threads = [threading.Thread(target=worker) for _ in range(num_threads)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        # All IDs should be unique
        self.assertEqual(len(results), num_threads * increments_per_thread)
        self.assertEqual(len(set(results)), len(results))

        # Should have generated IDs 1 through 1000
        self.assertEqual(sorted(results), list(range(1, num_threads * increments_per_thread + 1)))


if __name__ == '__main__':
    unittest.main()
