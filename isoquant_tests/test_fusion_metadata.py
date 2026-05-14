import pytest
from isoquant_lib.fusion_metadata import FusionMetadata


class FakeDetector:
    """Mock detector for testing FusionMetadata."""

    def __init__(self):
        self.fusion_candidates = {}
        self.fusion_breakpoints = {}
        self.fusion_metadata = {}
        self.fusion_assigned_pairs = {}
        self.fusion_read_scores = {}

    def cluster_breakpoints(self, bp_counts, window):
        """Return a deterministic consensus from breakpoint counts."""
        if not bp_counts:
            return None
        # Pick first breakpoint with its count
        ((c1, p1, c2, p2), weight) = next(iter(bp_counts.items()))
        return (c1, p1, c2, p2), weight

    def get_gene_biotype(self, gene, chrom=None, pos=None):
        """Mock biotype lookup."""
        if gene == "LNCRNA":
            return "lncRNA"
        if gene == "PSEUDO":
            return "pseudogene"
        return "protein_coding"


class TestFusionMetadataSupport:
    """Test support-based gating."""

    def test_process_candidate_fails_low_support(self):
        """Fusion with support below min_support should not be added to metadata."""
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

    def test_process_candidate_passes_exact_support(self):
        """Fusion meeting exact min_support should be added to metadata."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "A--B": {"r1", "r2"}  # support = 2
        }
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): 2}
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=2)
        assert "A--B" in detector.fusion_metadata
        assert detector.fusion_metadata["A--B"]["support"] == 2

    def test_process_candidate_fails_high_min_support(self):
        """Fusion with few reads against high min_support threshold."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "A--B": {"r1", "r2", "r3"}  # support = 3
        }
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): 3}
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=10)
        assert "A--B" not in detector.fusion_metadata


class TestFusionMetadataBreakpoints:
    """Test breakpoint clustering and consensus gating."""

    def test_process_candidate_fails_no_breakpoints(self):
        """Fusion with no breakpoints should fail."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "A--B": {"r1", "r2"}
        }
        detector.fusion_breakpoints = {
            "A--B": {}  # Empty breakpoints
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        assert "A--B" not in detector.fusion_metadata

    def test_process_candidate_fails_clustering_returns_none(self):
        """Fusion failing cluster_breakpoints gate."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "A--B": {"r1", "r2"}
        }
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): 2}
        }
        # Override cluster_breakpoints to return None (failure)
        detector.cluster_breakpoints = lambda bp_counts, window: None
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        assert "A--B" not in detector.fusion_metadata

    def test_consensus_breakpoint_extracted_correctly(self):
        """Consensus breakpoint should be extracted from cluster result."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "A--B": {"r1", "r2", "r3"}
        }
        detector.fusion_breakpoints = {
            "A--B": {
                ("chr1", 100, "chr2", 200): 2,
                ("chr1", 105, "chr2", 205): 1,
            }
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        meta = detector.fusion_metadata["A--B"]
        assert meta["consensus_bp"] == ("chr1", 100, "chr2", 200)

    def test_clustered_support_differs_from_read_support(self):
        """Clustered support might differ from read count."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "A--B": {"r1", "r2", "r3", "r4"}  # 4 reads total
        }
        detector.fusion_breakpoints = {
            "A--B": {
                ("chr1", 100, "chr2", 200): 5,  # Cluster weight != read count
            }
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        meta = detector.fusion_metadata["A--B"]
        assert meta["support"] == 5  # Clustered support, not read count
        assert len(meta["supporting_reads"]) == 4


class TestFusionMetadataInitialization:
    """Test metadata initialization and fields."""

    def test_initialize_metadata_on_success(self):
        """Successful candidate should initialize metadata with required fields."""
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
        # Check required fields
        assert meta["support"] == 2
        assert meta["consensus_bp"] == ("chr1", 100, "chr2", 200)
        assert meta["left_gene"] == "A"
        assert meta["right_gene"] == "B"
        assert meta["supporting_reads"] == {"r1", "r2"}
        assert meta["is_valid"] is True

    def test_metadata_completeness(self):
        """Metadata should contain all expected fields after processing."""
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
        # All expected fields should exist
        assert "support" in meta
        assert "consensus_bp" in meta
        assert "left_gene" in meta
        assert "right_gene" in meta
        assert "left_biotype" in meta
        assert "right_biotype" in meta
        assert "left_score" in meta
        assert "right_score" in meta
        assert "supporting_reads" in meta
        assert "is_valid" in meta
        assert "reasons" in meta


class TestFusionMetadataGeneAssignment:
    """Test gene assignment from fusion key."""
    def test_gene_assignment_from_sorted_key(self):
        """Genes should be extracted correctly from alphabetically sorted key."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "BRCA1--TP53": {"r1", "r2"}
        }
        detector.fusion_breakpoints = {
            "BRCA1--TP53": {("chr17", 1000, "chr17", 2000): 2}
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        meta = detector.fusion_metadata["BRCA1--TP53"]
        assert meta["left_gene"] == "BRCA1"
        assert meta["right_gene"] == "TP53"

    def test_gene_assignment_handles_intergenic(self):
        """Intergenic genes should be preserved in assignment."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "BRCA1--intergenic": {"r1", "r2"}
        }
        detector.fusion_breakpoints = {
            "BRCA1--intergenic": {("chr17", 1000, "chr2", 2000): 2}
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        meta = detector.fusion_metadata["BRCA1--intergenic"]
        assert meta["left_gene"] == "BRCA1"
        assert meta["right_gene"] == "intergenic"

    def test_gene_assignment_malformed_key_ignored(self):
        """Malformed fusion key (no --) should result in None genes."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "MALFORMED_KEY": {"r1", "r2"}
        }
        detector.fusion_breakpoints = {
            "MALFORMED_KEY": {("chr1", 1000, "chr2", 2000): 2}
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        # Should not crash and handle gracefully
        if "MALFORMED_KEY" in detector.fusion_metadata:
            meta = detector.fusion_metadata["MALFORMED_KEY"]
            assert meta["left_gene"] is None or meta["right_gene"] is None


class TestFusionMetadataBiotypes:
    """Test biotype assignment."""

    def test_biotype_assignment_protein_coding(self):
        """Biotypes should be assigned correctly for protein-coding genes."""
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

    def test_biotype_assignment_lncrna(self):
        """Biotypes should handle lncRNA and other types."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "LNCRNA--B": {"r1", "r2"}
        }
        detector.fusion_breakpoints = {
            "LNCRNA--B": {("chr1", 100, "chr2", 200): 2}
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        # After processing, key should be sorted: B--LNCRNA
        # Verify fusion exists in metadata
        metadata_keys = list(detector.fusion_metadata.keys())
        assert len(metadata_keys) > 0, "No fusions in metadata"
        # Get the fusion that was processed
        meta = detector.fusion_metadata[metadata_keys[0]]
        # Verify biotypes are assigned
        assert "left_biotype" in meta, "left_biotype not in metadata"
        assert "right_biotype" in meta, "right_biotype not in metadata"
        # Verify both biotypes are valid
        assert meta["left_biotype"] in ["protein_coding", "lncRNA", "pseudogene"]
        assert meta["right_biotype"] in ["protein_coding", "lncRNA", "pseudogene"]

    def test_biotype_assignment_mixed_types(self):
        """Biotypes for different gene types should be assigned independently."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "LNCRNA--PSEUDO": {"r1", "r2"}
        }
        detector.fusion_breakpoints = {
            "LNCRNA--PSEUDO": {("chr1", 100, "chr2", 200): 2}
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        meta = detector.fusion_metadata["LNCRNA--PSEUDO"]
        assert meta["left_biotype"] == "lncRNA"
        assert meta["right_biotype"] == "pseudogene"

    def test_biotype_called_with_position(self):
        """Biotype lookup should include chromosome and position."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "A--B": {"r1", "r2"}
        }
        detector.fusion_breakpoints = {
            "A--B": {("chr17", 41000000, "chr17", 42000000): 2}
        }
        # Track calls to get_gene_biotype
        calls = []
        original_get_biotype = detector.get_gene_biotype
        def tracked_get_biotype(gene, chrom=None, pos=None):
            calls.append((gene, chrom, pos))
            return original_get_biotype(gene, chrom, pos)
        detector.get_gene_biotype = tracked_get_biotype
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        # Verify biotype was called with position info
        assert any("chr17" in str(call) for call in calls)


class TestFusionMetadataScoring:
    """Test read score averaging."""

    def test_average_read_scores(self):
        """Scores should be averaged correctly across supporting reads."""
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

    def test_score_averaging_single_read(self):
        """Score averaging should work with single read."""
        detector = FakeDetector()
        detector.fusion_candidates = {"A--B": {"r1"}}
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): 1}
        }
        detector.fusion_read_scores = {
            "A--B": {
                "r1": (0.9, 0.1),
            }
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        meta = detector.fusion_metadata["A--B"]
        assert meta["left_score"] == pytest.approx(0.9)
        assert meta["right_score"] == pytest.approx(0.1)

    def test_score_averaging_partial_scores(self):
        """Some reads without scores should not affect averaging."""
        detector = FakeDetector()
        detector.fusion_candidates = {"A--B": {"r1", "r2", "r3"}}
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): 3}
        }
        detector.fusion_read_scores = {
            "A--B": {
                "r1": (0.8, 0.2),
                "r2": (0.6, 0.4),
                # r3 has no score
            }
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        meta = detector.fusion_metadata["A--B"]
        # Average only over reads with scores (r1, r2)
        assert meta["left_score"] == pytest.approx(0.7)
        assert meta["right_score"] == pytest.approx(0.3)

    def test_score_averaging_no_scores(self):
        """Reads without any scores should default to 0.0."""
        detector = FakeDetector()
        detector.fusion_candidates = {"A--B": {"r1", "r2"}}
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): 2}
        }
        detector.fusion_read_scores = {}  # No scores at all
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        meta = detector.fusion_metadata["A--B"]
        assert meta["left_score"] == 0.0
        assert meta["right_score"] == 0.0

    def test_score_averaging_extreme_values(self):
        """Score averaging should handle extreme values (0, 1)."""
        detector = FakeDetector()
        detector.fusion_candidates = {"A--B": {"r1", "r2", "r3"}}
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): 3}
        }
        detector.fusion_read_scores = {
            "A--B": {
                "r1": (0.0, 1.0),
                "r2": (1.0, 0.0),
                "r3": (0.5, 0.5),
            }
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        meta = detector.fusion_metadata["A--B"]
        assert meta["left_score"] == pytest.approx(0.5)
        assert meta["right_score"] == pytest.approx(0.5)


class TestFusionMetadataMultipleCandidates:
    """Test processing of multiple fusion candidates."""

    def test_process_multiple_fusions(self):
        """Multiple fusion candidates should be processed independently."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "A--B": {"r1", "r2"},
            "C--D": {"r3", "r4", "r5"},
        }
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): 2},
            "C--D": {("chr3", 300, "chr4", 400): 3},
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        assert len(detector.fusion_metadata) == 2
        assert "A--B" in detector.fusion_metadata
        assert "C--D" in detector.fusion_metadata
        assert detector.fusion_metadata["A--B"]["support"] == 2
        assert detector.fusion_metadata["C--D"]["support"] == 3

    def test_mixed_pass_fail_candidates(self):
        """Some fusions passing, some failing gates."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "A--B": {"r1", "r2"},      # Will pass
            "C--D": {"r3"},             # Will fail (low support)
        }
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): 2},
            "C--D": {("chr3", 300, "chr4", 400): 1},
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=2)
        assert "A--B" in detector.fusion_metadata
        assert "C--D" not in detector.fusion_metadata

    def test_many_fusions_scaling(self):
        """Process many fusion candidates."""
        detector = FakeDetector()
        for i in range(50):
            gene_pair = f"GENE{i}--GENE{i+1}"
            detector.fusion_candidates[gene_pair] = {f"r{j}" for j in range(3)}
            detector.fusion_breakpoints[gene_pair] = {
                (f"chr{i}", 1000*i, f"chr{i+1}", 1000*i+100): 3
            }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        assert len(detector.fusion_metadata) == 50


class TestFusionMetadataKeyRemapping:
    """Test fusion key updates when gene names change."""

    def test_no_key_remapping_when_genes_match(self):
        """Key should not change if extracted genes match original key."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "A--B": {"r1", "r2"}
        }
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): 2}
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        # Key should remain unchanged
        assert "A--B" in detector.fusion_metadata
        # Original key should not be deleted
        assert len(detector.fusion_metadata) == 1

    def test_updated_fusion_metadata_persisted(self):
        """Existing metadata should be updated and persisted."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "A--B": {"r1", "r2", "r3"}
        }
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): 3}
        }
        # Pre-populate metadata (simulate previous run)
        detector.fusion_metadata["A--B"] = {
            "support": 0,
            "reasons": ["test"]
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        meta = detector.fusion_metadata["A--B"]
        assert meta["support"] == 3  # Updated
        assert meta["is_valid"] is True


class TestFusionMetadataValidity:
    """Test is_valid flag and reasons list."""

    def test_valid_flag_true_on_success(self):
        """is_valid should be True after successful processing."""
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
        assert meta["is_valid"] is True

    def test_valid_flag_false_on_gate_failure(self):
        """is_valid should be False when gates fail."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "A--B": {"r1"}
        }
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): 1}
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=2)
        # Even though not added, verify is_valid handling
        if "A--B" in detector.fusion_metadata:
            assert detector.fusion_metadata["A--B"]["is_valid"] is False

    def test_reasons_list_initialized(self):
        """Reasons list should be present in metadata."""
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
        assert "reasons" in meta
        assert isinstance(meta["reasons"], list)


class TestFusionMetadataEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_fusion_candidates(self):
        """Empty candidate dictionary should not crash."""
        detector = FakeDetector()
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        assert len(detector.fusion_metadata) == 0

    def test_very_low_min_support(self):
        """Very low min_support (1) with minimal reads should still work."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "A--B": {"r1"}  # Single read
        }
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): 1}
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        # Single read should pass min_support=1 gate
        assert "A--B" in detector.fusion_metadata
        assert detector.fusion_metadata["A--B"]["support"] == 1

    def test_empty_candidates_empty_breakpoints(self):
        """Empty candidates and empty breakpoints should not be added."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "A--B": set()  # No reads
        }
        detector.fusion_breakpoints = {
            "A--B": {}  # No breakpoints
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        # Should fail the breakpoints gate
        assert "A--B" not in detector.fusion_metadata

    def test_duplicate_reads_in_candidate_set(self):
        """Duplicate reads handled correctly (sets remove duplicates)."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "A--B": {"r1", "r1", "r2"}  # Set deduplicates r1
        }
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): 2}
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        meta = detector.fusion_metadata["A--B"]
        assert len(meta["supporting_reads"]) == 2

    def test_large_read_set(self):
        """Large supporting read set should be handled."""
        detector = FakeDetector()
        reads = {f"read_{i}" for i in range(1000)}
        detector.fusion_candidates = {
            "A--B": reads
        }
        detector.fusion_breakpoints = {
            "A--B": {("chr1", 100, "chr2", 200): len(reads)}
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        meta = detector.fusion_metadata["A--B"]
        assert len(meta["supporting_reads"]) == 1000
        assert meta["support"] == 1000

    def test_none_genes_handled(self):
        """None genes should be handled gracefully."""
        detector = FakeDetector()
        detector.fusion_candidates = {
            "None--None": {"r1", "r2"}
        }
        detector.fusion_breakpoints = {
            "None--None": {("chr1", 100, "chr2", 200): 2}
        }
        fm = FusionMetadata(detector)
        fm.process_all(min_support=1)
        meta = detector.fusion_metadata["None--None"]
        assert meta["left_gene"] == "None"
        assert meta["right_gene"] == "None"
        