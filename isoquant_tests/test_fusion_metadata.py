import pytest
from unittest.mock import Mock, MagicMock
from isoquant_lib.fusion_metadata import FusionMetadata

class FakeDetector:
    def __init__(self):
        self.fusion_candidates = {}
        self.fusion_breakpoints = {}
        self.fusion_metadata = {}
        self.fusion_assigned_pairs = {}
        self.fusion_read_scores = {}

    def cluster_breakpoints(self, bp_counts, window):
        # Return a deterministic consensus
        # bp_counts = {(c1,p1,c2,p2): weight}
        ((c1, p1, c2, p2), weight) = next(iter(bp_counts.items()))
        return (c1, p1, c2, p2), weight

    def get_gene_biotype(self, gene, chrom=None, pos=None):
        return "protein_coding"

class TestFusionMetadata:
    """Test basic utility methods of FusionMetadata."""

    def test_process_candidate_fails_low_support(self):
        detector = FakeDetector()
        detector.fusion_candidates = {
            "A--B": {"read1"}  # support = 1
        }
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): 1}
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=2)
        assert "A--B" not in detector.fusion_metadata

    def test_initialize_metadata_on_success(self):
        detector = FakeDetector()
        detector.fusion_candidates = {
            "A--B": {"r1", "r2"}
        }
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): 2}
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        meta = detector.fusion_metadata["A--B"]
        assert meta["support"] == 2
        assert meta["consensus_bp"] == ("chr1", 100, "chr2", 200)
        assert meta["left_gene"] == "A"
        assert meta["right_gene"] == "B"

    def test_fusion_key_renaming(self):
        detector = FakeDetector()
        detector.fusion_candidates = {
            "X--Y": {"r1", "r2"}
        }
        detector.fusion_breakpoints = {
            "X--Y": {("chr1", 100, "chr2", 200): 2}
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        assert "A--B" in detector.fusion_metadata or "X--Y" in detector.fusion_metadata

    def test_biotype_assignment(self):
        detector = FakeDetector()
        detector.fusion_candidates = {
            "A--B": {"r1", "r2"}
        }
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): 2}
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        meta = detector.fusion_metadata["A--B"]
        assert meta["left_biotype"] == "protein_coding"
        assert meta["right_biotype"] == "protein_coding"

    def test_average_read_scores(self):
        detector = FakeDetector()
        detector.fusion_candidates = {"A--B": {"r1", "r2"}}
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): 2}
        }
        detector.fusion_read_scores = {
            "A--B": {
                "r1": (0.8, 0.2),
                "r2": (0.6, 0.4),
            }
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        meta = detector.fusion_metadata["A--B"]
        assert meta["left_score"] == pytest.approx(0.7)
        assert meta["right_score"] == pytest.approx(0.3)