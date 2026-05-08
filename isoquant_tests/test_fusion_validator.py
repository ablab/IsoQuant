import pytest
from unittest.mock import Mock, MagicMock
from isoquant_lib.fusion_validator import FusionValidator


class MockDetector:
    """Mock detector object for testing FusionValidator."""

    def __init__(self):
        self.fusion_metadata = {}
        self.fusion_candidates = {}
        self.fusion_breakpoints = {}
        self.fusion_assigned_pairs = {}

    def get_gene_biotype(self, gene, chrom=None, pos=None):
        """Mock method to get gene biotype."""
        # Simple mock that returns biotype based on gene name
        biotypes = {
            "BRCA1": "protein_coding",
            "BRCA2": "protein_coding",
            "MYC": "protein_coding",
            "PML": "protein_coding",
            "RARA": "protein_coding",
            "RUNX1": "protein_coding",
            "RUNX1T1": "protein_coding",
            "ETV6": "protein_coding",
            "BCR": "protein_coding",
            "ABL1": "protein_coding",
            "RPL5": "ribosomal_RNA",
            "RPS3": "ribosomal_RNA",
            "HIST1H1A": "histone",
            "H2AZ1": "histone",
            "MIR123": "miRNA",
            "INTERGENIC": None,
        }
        return biotypes.get(gene, "protein_coding")


class TestFusionValidatorBasics:
    """Test basic utility methods of FusionValidator."""

    def setup_method(self):
        self.detector = MockDetector()
        self.validator = FusionValidator(self.detector)

    def test_is_allowed_biotype_protein_coding(self):
        """Test that protein_coding is allowed."""
        assert self.validator._is_allowed_biotype("protein_coding") is True

    def test_is_allowed_biotype_ig_genes(self):
        """Test that IG genes are allowed."""
        assert self.validator._is_allowed_biotype("IG_C_gene") is True
        assert self.validator._is_allowed_biotype("IG_V_gene") is True
        assert self.validator._is_allowed_biotype("IG_D_gene") is True
        assert self.validator._is_allowed_biotype("IG_J_gene") is True

    def test_is_allowed_biotype_tr_genes(self):
        """Test that TR genes are allowed."""
        assert self.validator._is_allowed_biotype("TR_V_gene") is True
        assert self.validator._is_allowed_biotype("TR_D_gene") is True
        assert self.validator._is_allowed_biotype("TR_J_gene") is True
        assert self.validator._is_allowed_biotype("TR_C_gene") is True

    def test_is_allowed_biotype_non_coding(self):
        """Test that non-coding biotypes are not allowed."""
        assert self.validator._is_allowed_biotype("miRNA") is False
        assert self.validator._is_allowed_biotype("lncRNA") is False
        assert self.validator._is_allowed_biotype("ribosomal_RNA") is False
        assert self.validator._is_allowed_biotype("histone") is False

    def test_is_allowed_biotype_none(self):
        """Test handling of None biotype."""
        assert self.validator._is_allowed_biotype(None) is False

    def test_is_ribosomal_or_histone_gene_rpl(self):
        """Test detection of ribosomal L genes."""
        assert self.validator._is_ribosomal_or_histone_gene("RPL5") is True
        assert self.validator._is_ribosomal_or_histone_gene("RPL22") is True
        assert self.validator._is_ribosomal_or_histone_gene("rpl5") is True

    def test_is_ribosomal_or_histone_gene_rps(self):
        """Test detection of ribosomal S genes."""
        assert self.validator._is_ribosomal_or_histone_gene("RPS3") is True
        assert self.validator._is_ribosomal_or_histone_gene("RPS6") is True

    def test_is_ribosomal_or_histone_gene_rplp(self):
        """Test detection of ribosomal P genes."""
        assert self.validator._is_ribosomal_or_histone_gene("RPLP1") is True
        assert self.validator._is_ribosomal_or_histone_gene("RPLP2") is True

    def test_is_ribosomal_or_histone_gene_hist(self):
        """Test detection of histone genes."""
        assert self.validator._is_ribosomal_or_histone_gene("HIST1H1A") is True
        assert self.validator._is_ribosomal_or_histone_gene("HIST2H2A") is True
        assert self.validator._is_ribosomal_or_histone_gene("H2AZ1") is True
        assert self.validator._is_ribosomal_or_histone_gene("H3") is True
        assert self.validator._is_ribosomal_or_histone_gene("H4") is True

    def test_is_ribosomal_or_histone_gene_h1(self):
        """Test detection of H1 histone genes."""
        assert self.validator._is_ribosomal_or_histone_gene("H1F0") is True
        assert self.validator._is_ribosomal_or_histone_gene("H1FX") is True

    def test_is_ribosomal_or_histone_gene_non_match(self):
        """Test that non-ribosomal genes don't match."""
        assert self.validator._is_ribosomal_or_histone_gene("BRCA1") is False
        assert self.validator._is_ribosomal_or_histone_gene("MYC") is False
        assert self.validator._is_ribosomal_or_histone_gene("EGFR") is False

    def test_is_ribosomal_or_histone_gene_none(self):
        """Test handling of None input."""
        assert self.validator._is_ribosomal_or_histone_gene(None) is False

    def test_is_ribosomal_or_histone_gene_non_string(self):
        """Test handling of non-string input."""
        assert self.validator._is_ribosomal_or_histone_gene(123) is False
        assert self.validator._is_ribosomal_or_histone_gene([]) is False

    def test_is_driver_gene_true(self):
        """Test that known driver genes are identified."""
        driver_genes = [
            "BCR", "ABL1", "RUNX1", "RUNX1T1", "PML", "RARA", "ETV6",
            "NUP98", "NUP214", "KMT2A", "MLLT3", "EWSR1", "ERG",
            "FGFR1", "RET", "ROS1", "ALK", "TMPRSS2", "FLI1",
            "MYH11", "CBFB"
        ]
        for gene in driver_genes:
            assert self.validator.is_driver_gene(gene) is True

    def test_is_driver_gene_false(self):
        """Test that non-driver genes are not identified as drivers."""
        assert self.validator.is_driver_gene("BRCA1") is False
        assert self.validator.is_driver_gene("MYC") is False
        assert self.validator.is_driver_gene("TP53") is False

    def test_is_driver_gene_none(self):
        """Test handling of None input."""
        assert self.validator.is_driver_gene(None) is False

    def test_is_multicopy_artifact_family_ribosomal(self):
        """Test detection of ribosomal protein families."""
        assert self.validator.is_multicopy_artifact_family("RPL5") is True
        assert self.validator.is_multicopy_artifact_family("RPS3") is True
        assert self.validator.is_multicopy_artifact_family("MRPL1") is True
        assert self.validator.is_multicopy_artifact_family("MRPS7") is True

    def test_is_multicopy_artifact_family_histone(self):
        """Test detection of histone families."""
        assert self.validator.is_multicopy_artifact_family("HIST1H1A") is True
        assert self.validator.is_multicopy_artifact_family("H2AZ1") is True

    def test_is_multicopy_artifact_family_hla(self):
        """Test detection of HLA families."""
        assert self.validator.is_multicopy_artifact_family("HLA-A") is True
        assert self.validator.is_multicopy_artifact_family("HLA-B") is True
        assert self.validator.is_multicopy_artifact_family("MICA") is True
        assert self.validator.is_multicopy_artifact_family("MICB") is True

    def test_is_multicopy_artifact_family_ig_tr(self):
        """Test detection of immunoglobulin and T-cell receptor families."""
        assert self.validator.is_multicopy_artifact_family("IGH") is True
        assert self.validator.is_multicopy_artifact_family("IGK") is True
        assert self.validator.is_multicopy_artifact_family("IGL") is True
        assert self.validator.is_multicopy_artifact_family("TRAV") is True
        assert self.validator.is_multicopy_artifact_family("TRBV") is True

    def test_is_multicopy_artifact_family_other(self):
        """Test detection of other multicopy families."""
        assert self.validator.is_multicopy_artifact_family("KRT1") is True
        assert self.validator.is_multicopy_artifact_family("OR1A1") is True
        assert self.validator.is_multicopy_artifact_family("ZNF100") is True

    def test_is_multicopy_artifact_family_non_match(self):
        """Test that normal genes don't match."""
        assert self.validator.is_multicopy_artifact_family("BRCA1") is False
        assert self.validator.is_multicopy_artifact_family("MYC") is False
        assert self.validator.is_multicopy_artifact_family("EGFR") is False


class TestGetBreakpointCoords:
    """Test the _get_breakpoint_coords method."""

    def setup_method(self):
        self.detector = MockDetector()
        self.validator = FusionValidator(self.detector)

    def test_get_breakpoint_coords_valid(self):
        """Test extraction of valid breakpoint coordinates."""
        fusion_key = "fusion1"
        bp = ("chr1", 1000, "chr2", 2000)
        self.detector.fusion_breakpoints[fusion_key] = {bp: 5}
        left_chr, left_pos, right_chr, right_pos = self.validator._get_breakpoint_coords(fusion_key)
        assert (left_chr, left_pos, right_chr, right_pos) == bp

    def test_get_breakpoint_coords_multiple_bps(self):
        """Test extraction of first breakpoint when multiple exist."""
        fusion_key = "fusion1"
        bp1 = ("chr1", 1000, "chr2", 2000)
        bp2 = ("chr1", 1001, "chr2", 2001)
        self.detector.fusion_breakpoints[fusion_key] = {bp1: 5, bp2: 3}
        left_chr, left_pos, right_chr, right_pos = self.validator._get_breakpoint_coords(fusion_key)
        assert (left_chr, left_pos, right_chr, right_pos) == bp1

    def test_get_breakpoint_coords_missing_key(self):
        """Test handling of missing fusion key."""
        left_chr, left_pos, right_chr, right_pos = self.validator._get_breakpoint_coords("nonexistent")
        assert (left_chr, left_pos, right_chr, right_pos) == (None, None, None, None)

    def test_get_breakpoint_coords_empty_dict(self):
        """Test handling of empty breakpoint dict."""
        fusion_key = "fusion1"
        self.detector.fusion_breakpoints[fusion_key] = {}

        left_chr, left_pos, right_chr, right_pos = self.validator._get_breakpoint_coords(fusion_key)
        assert (left_chr, left_pos, right_chr, right_pos) == (None, None, None, None)


class TestSalvageAndCheckGenePair:
    """Test the _salvage_and_check_gene_pair method."""

    def setup_method(self):
        self.detector = MockDetector()
        self.validator = FusionValidator(self.detector)

    def test_salvage_both_allowed(self):
        """Test checking when both genes are allowed."""
        left_gene, left_biotype, right_gene, right_biotype, has_non_coding, reason = (
            self.validator._salvage_and_check_gene_pair(
                "BRCA1", "BRCA2", "chr1", 1000, "chr2", 2000
            )
        )
        assert left_gene == "BRCA1"
        assert left_biotype == "protein_coding"
        assert right_gene == "BRCA2"
        assert right_biotype == "protein_coding"
        assert has_non_coding is False
        assert reason is None

    def test_salvage_left_non_coding(self):
        """Test detection of non-coding left gene."""
        left_gene, left_biotype, right_gene, right_biotype, has_non_coding, reason = (
            self.validator._salvage_and_check_gene_pair(
                "RPL5", "BRCA2", "chr1", 1000, "chr2", 2000
            )
        )
        assert has_non_coding is True
        assert "RPL5" in reason

    def test_salvage_right_non_coding(self):
        """Test detection of non-coding right gene."""
        left_gene, left_biotype, right_gene, right_biotype, has_non_coding, reason = (
            self.validator._salvage_and_check_gene_pair(
                "BRCA1", "HIST1H1A", "chr1", 1000, "chr2", 2000
            )
        )
        assert has_non_coding is True
        assert "HIST1H1A" in reason

    def test_salvage_both_non_coding(self):
        """Test detection when both genes are non-coding."""
        left_gene, left_biotype, right_gene, right_biotype, has_non_coding, reason = (
            self.validator._salvage_and_check_gene_pair(
                "RPL5", "HIST1H1A", "chr1", 1000, "chr2", 2000
            )
        )
        assert has_non_coding is True
        assert reason is not None


class TestFilterNonCodingGenes:
    """Test the filter_non_coding_genes method."""

    def setup_method(self):
        self.detector = MockDetector()
        self.validator = FusionValidator(self.detector)

    def test_filter_keeps_all_coding(self):
        """Test that fusions with all coding partners are kept."""
        fusion_key = "fusion1"
        self.detector.fusion_metadata[fusion_key] = {
            "left_gene": "BRCA1",
            "right_gene": "BRCA2",
            "consensus_bp": ("chr1", 1000, "chr2", 2000)
        }

        self.validator.filter_non_coding_genes()
        assert fusion_key in self.detector.fusion_metadata

    def test_filter_removes_non_coding_left(self):
        """Test that fusions with non-coding left partner are removed."""
        fusion_key = "fusion1"
        self.detector.fusion_metadata[fusion_key] = {
            "left_gene": "RPL5",
            "right_gene": "BRCA2",
            "consensus_bp": ("chr1", 1000, "chr2", 2000)
        }
        self.validator.filter_non_coding_genes()
        assert fusion_key not in self.detector.fusion_metadata

    def test_filter_removes_non_coding_right(self):
        """Test that fusions with non-coding right partner are removed."""
        fusion_key = "fusion1"
        self.detector.fusion_metadata[fusion_key] = {
            "left_gene": "BRCA1",
            "right_gene": "HIST1H1A",
            "consensus_bp": ("chr1", 1000, "chr2", 2000)
        }
        self.validator.filter_non_coding_genes()
        assert fusion_key not in self.detector.fusion_metadata

    def test_filter_multiple_fusions(self):
        """Test filtering of multiple fusions."""
        self.detector.fusion_metadata = {
            "fusion1": {
                "left_gene": "BRCA1",
                "right_gene": "BRCA2",
                "consensus_bp": ("chr1", 1000, "chr2", 2000)
            },
            "fusion2": {
                "left_gene": "RPL5",
                "right_gene": "BRCA2",
                "consensus_bp": ("chr1", 1000, "chr3", 3000)
            },
            "fusion3": {
                "left_gene": "PML",
                "right_gene": "RARA",
                "consensus_bp": ("chr15", 1000, "chr17", 4000)
            }
        }
        self.validator.filter_non_coding_genes()
        assert "fusion1" in self.detector.fusion_metadata
        assert "fusion2" not in self.detector.fusion_metadata
        assert "fusion3" in self.detector.fusion_metadata


class TestFilterMulticopyArtifactPairs:
    """Test the filter_multicopy_artifact_pairs method."""

    def setup_method(self):
        self.detector = MockDetector()
        self.validator = FusionValidator(self.detector)

    def test_keeps_driver_gene_pairs(self):
        """Test that driver gene pairs are not marked as invalid."""
        fusion_key = "bcr_abl"
        self.detector.fusion_metadata[fusion_key] = {
            "left_gene": "BCR",
            "right_gene": "ABL1",
            "is_valid": True
        }

        self.validator.filter_multicopy_artifact_pairs()
        assert self.detector.fusion_metadata[fusion_key].get("is_valid", True) is True

    def test_marks_invalid_ribosomal_pairs(self):
        """Test that ribosomal pairs are marked as invalid."""
        fusion_key = "rpl_rps"
        self.detector.fusion_metadata[fusion_key] = {
            "left_gene": "RPL5",
            "right_gene": "RPS3",
            "is_valid": True
        }
        self.validator.filter_multicopy_artifact_pairs()
        assert self.detector.fusion_metadata[fusion_key]["is_valid"] is False
        assert "reasons" in self.detector.fusion_metadata[fusion_key]

    def test_marks_invalid_histone_pairs(self):
        """Test that histone pairs are marked as invalid."""
        fusion_key = "hist_hist"
        self.detector.fusion_metadata[fusion_key] = {
            "left_gene": "HIST1H1A",
            "right_gene": "H2AZ1",
            "is_valid": True
        }
        self.validator.filter_multicopy_artifact_pairs()
        assert self.detector.fusion_metadata[fusion_key]["is_valid"] is False

    def test_keeps_mixed_pairs(self):
        """Test that mixed gene pairs (one ribosomal, one normal) are kept."""
        fusion_key = "rpl_brca"
        self.detector.fusion_metadata[fusion_key] = {
            "left_gene": "RPL5",
            "right_gene": "BRCA1",
            "is_valid": True
        }
        self.validator.filter_multicopy_artifact_pairs()
        assert self.detector.fusion_metadata[fusion_key].get("is_valid", True) is True


class TestApplyFrequencyFilters:
    """Test the apply_frequency_filters method."""

    def setup_method(self):
        self.detector = MockDetector()
        self.validator = FusionValidator(self.detector)

    def test_keeps_low_frequency_genes(self):
        """Test that fusions with low-frequency genes are kept."""
        self.detector.fusion_metadata = {
            "fusion1": {
                "left_gene": "BRCA1",
                "right_gene": "BRCA2",
                "is_valid": True
            }
        }
        self.validator.apply_frequency_filters()
        assert self.detector.fusion_metadata["fusion1"].get("is_valid", True) is True

    def test_marks_invalid_high_frequency_left(self):
        """Test that high-frequency left genes trigger invalidity."""
        self.detector.fusion_metadata = {
            "fusion1": {"left_gene": "RPL5", "right_gene": "BRCA2", "is_valid": True},
            "fusion2": {"left_gene": "RPL5", "right_gene": "BRCA1", "is_valid": True},
            "fusion3": {"left_gene": "RPL5", "right_gene": "TP53", "is_valid": True},
            "fusion4": {"left_gene": "RPL5", "right_gene": "MYC", "is_valid": True},
            "fusion5": {"left_gene": "RPL5", "right_gene": "EGFR", "is_valid": True},
        }
        self.validator.apply_frequency_filters()
        # RPL5 appears 5 times, which exceeds threshold of 2 for ribosomal genes
        for fusion_key in self.detector.fusion_metadata:
            assert self.detector.fusion_metadata[fusion_key]["is_valid"] is False

    def test_ribosomal_lower_threshold(self):
        """Test that ribosomal genes have lower threshold (2) than others (4)."""
        # Test ribosomal gene with 3 occurrences (above threshold of 2)
        self.detector.fusion_metadata = {
            "fusion1": {"left_gene": "RPL5", "right_gene": "BRCA2", "is_valid": True},
            "fusion2": {"left_gene": "RPL5", "right_gene": "BRCA1", "is_valid": True},
            "fusion3": {"left_gene": "RPL5", "right_gene": "TP53", "is_valid": True},
        }
        self.validator.apply_frequency_filters()
        for fusion_key in self.detector.fusion_metadata:
            assert self.detector.fusion_metadata[fusion_key]["is_valid"] is False

    def test_normal_genes_higher_threshold(self):
        """Test that normal genes need 5+ occurrences to be marked invalid."""
        # MYC appears 4 times (below threshold of 4, so should be kept)
        self.detector.fusion_metadata = {
            "fusion1": {"left_gene": "MYC", "right_gene": "BRCA2", "is_valid": True},
            "fusion2": {"left_gene": "MYC", "right_gene": "BRCA1", "is_valid": True},
            "fusion3": {"left_gene": "MYC", "right_gene": "TP53", "is_valid": True},
            "fusion4": {"left_gene": "MYC", "right_gene": "EGFR", "is_valid": True},
        }
        self.validator.apply_frequency_filters()
        # All should be kept since count (4) is not > threshold (4)
        for fusion_key in self.detector.fusion_metadata:
            assert self.detector.fusion_metadata[fusion_key].get("is_valid", True) is True

    def test_intergenic_not_counted(self):
        """Test that intergenic genes are not counted."""
        self.detector.fusion_metadata = {
            "fusion1": {"left_gene": "intergenic", "right_gene": "BRCA2", "is_valid": True},
            "fusion2": {"left_gene": "intergenic", "right_gene": "BRCA1", "is_valid": True},
        }
        self.validator.apply_frequency_filters()
        # Should be kept since intergenic is not counted
        for fusion_key in self.detector.fusion_metadata:
            assert self.detector.fusion_metadata[fusion_key].get("is_valid", True) is True


class TestConfidence:
    """Test the confidence method."""

    def setup_method(self):
        self.detector = MockDetector()
        self.validator = FusionValidator(self.detector)

    def test_confidence_minimum(self):
        """Test that confidence is never negative."""
        meta = {
            "support": 0,
            "reconstruction_ok": False,
            "realignment_hits": 0,
            "best_hit_mapq": 0,
        }
        confidence = self.validator.confidence(meta)
        assert confidence >= 0.0

    def test_confidence_maximum(self):
        """Test that confidence is capped at 1.0."""
        meta = {
            "support": 100,
            "reconstruction_ok": True,
            "realignment_hits": 10,
            "best_hit_mapq": 60,
            "left_gene": "BCR",
            "right_gene": "ABL1"
        }
        confidence = self.validator.confidence(meta)
        assert confidence <= 1.0

    def test_confidence_support_component(self):
        """Test support component calculation."""
        meta = {
            "support": 5,
            "reconstruction_ok": False,
            "realignment_hits": 0,
            "best_hit_mapq": 0,
        }
        confidence = self.validator.confidence(meta)
        # support is 5/10 = 0.5
        assert confidence == 0.5

    def test_confidence_reconstruction_bonus(self):
        """Test reconstruction bonus."""
        meta_no_recon = {
            "support": 1,
            "reconstruction_ok": False,
            "realignment_hits": 0,
            "best_hit_mapq": 0,
        }
        meta_with_recon = {
            "support": 1,
            "reconstruction_ok": True,
            "realignment_hits": 0,
            "best_hit_mapq": 0,
        }
        conf_no_recon = self.validator.confidence(meta_no_recon)
        conf_with_recon = self.validator.confidence(meta_with_recon)
        # Reconstruction bonus is 0.2
        assert conf_with_recon - conf_no_recon == pytest.approx(0.2)

    def test_confidence_driver_gene_bonus(self):
        """Test driver gene bonus."""
        meta_normal = {
            "support": 1,
            "reconstruction_ok": False,
            "realignment_hits": 0,
            "best_hit_mapq": 0,
            "left_gene": "BRCA1",
            "right_gene": "BRCA2"
        }
        meta_driver = {
            "support": 1,
            "reconstruction_ok": False,
            "realignment_hits": 0,
            "best_hit_mapq": 0,
            "left_gene": "BCR",
            "right_gene": "ABL1"
        }
        conf_normal = self.validator.confidence(meta_normal)
        conf_driver = self.validator.confidence(meta_driver)
        # Driver gene bonus is 0.2
        assert conf_driver - conf_normal == pytest.approx(0.2)

    def test_confidence_artifact_penalty(self):
        """Test artifact penalty."""
        meta_normal = {
            "support": 5,
            "reconstruction_ok": False,
            "realignment_hits": 0,
            "best_hit_mapq": 0,
            "left_gene": "BRCA1",
            "right_gene": "BRCA2"
        }
        meta_artifact = {
            "support": 5,
            "reconstruction_ok": False,
            "realignment_hits": 0,
            "best_hit_mapq": 0,
            "left_gene": "RPL5",
            "right_gene": "RPS3"
        }
        conf_normal = self.validator.confidence(meta_normal)
        conf_artifact = self.validator.confidence(meta_artifact)
        # Artifact penalty is 0.3
        assert conf_normal - conf_artifact == pytest.approx(0.3)


class TestMergeFusions:
    """Test fusion merging methods."""

    def setup_method(self):
        self.detector = MockDetector()
        self.validator = FusionValidator(self.detector)

    def test_merge_fusion_candidates(self):
        """Test merging of fusion candidates."""
        keep_key = "fusion1"
        discard_key = "fusion2"
        self.detector.fusion_candidates[keep_key] = {"read1", "read2"}
        self.detector.fusion_candidates[discard_key] = {"read3", "read4"}
        self.validator._merge_fusion_candidates(keep_key, discard_key)
        assert self.detector.fusion_candidates[keep_key] == {"read1", "read2", "read3", "read4"}

    def test_merge_fusion_breakpoints(self):
        """Test merging of fusion breakpoints."""
        keep_key = "fusion1"
        discard_key = "fusion2"
        bp1 = ("chr1", 1000, "chr2", 2000)
        bp2 = ("chr1", 1001, "chr2", 2001)
        self.detector.fusion_breakpoints[keep_key] = {bp1: 5}
        self.detector.fusion_breakpoints[discard_key] = {bp2: 3}
        self.validator._merge_fusion_candidates(keep_key, discard_key)
        assert self.detector.fusion_breakpoints[keep_key][bp1] == 5
        assert self.detector.fusion_breakpoints[keep_key][bp2] == 3

    def test_merge_fusion_metadata(self):
        """Test merging of fusion metadata (supporting reads)."""
        keep_key = "fusion1"
        discard_key = "fusion2"
        self.detector.fusion_metadata[keep_key] = {"supporting_reads": {"read1", "read2"}}
        self.detector.fusion_metadata[discard_key] = {"supporting_reads": {"read3", "read4"}}
        self.validator._merge_fusion_candidates(keep_key, discard_key)
        assert self.detector.fusion_metadata[keep_key]["supporting_reads"] == {"read1", "read2", "read3", "read4"}
        assert self.detector.fusion_metadata[keep_key]["support"] == 4

    def test_merge_fully_identical(self):
        """Test merging of fully identical fusions."""
        bp = ("chr1", 1000, "chr2", 2000)
        fusion1 = "fusion1"
        fusion2 = "fusion2"
        self.detector.fusion_metadata = {
            fusion1: {
                "left_gene": "BRCA1",
                "right_gene": "BRCA2",
                "consensus_bp": bp,
                "supporting_reads": {"read1", "read2"}
            },
            fusion2: {
                "left_gene": "BRCA1",
                "right_gene": "BRCA2",
                "consensus_bp": bp,
                "supporting_reads": {"read3"}
            }
        }
        self.detector.fusion_candidates = {
            fusion1: {"read1", "read2"},
            fusion2: {"read3"}
        }
        self.detector.fusion_breakpoints = {
            fusion1: {bp: 2},
            fusion2: {bp: 1}
        }
        fusions_to_discard = self.validator._merge_fully_identical()
        assert fusion2 in fusions_to_discard
        assert len(self.detector.fusion_candidates[fusion1]) == 3

    def test_remove_discarded_fusions(self):
        """Test removal of discarded fusions."""
        fusion_keys = ["fusion1", "fusion2", "fusion3"]

        for key in fusion_keys:
            self.detector.fusion_metadata[key] = {"gene": "test"}
            self.detector.fusion_candidates[key] = {"read1"}
            self.detector.fusion_breakpoints[key] = {("chr1", 1000, "chr2", 2000): 1}
        fusions_to_discard = {"fusion2"}
        self.validator._remove_discarded_fusions(fusions_to_discard)
        assert "fusion1" in self.detector.fusion_metadata
        assert "fusion2" not in self.detector.fusion_metadata
        assert "fusion3" in self.detector.fusion_metadata


class TestFilterEarlyNonCodingGenes:
    """Test the filter_early_non_coding_genes method."""

    def setup_method(self):
        self.detector = MockDetector()
        self.validator = FusionValidator(self.detector)

    def test_filter_keeps_coding(self):
        """Test that coding fusions are kept."""
        fusion_key = "fusion1"
        bp = ("chr1", 1000, "chr2", 2000)
        self.detector.fusion_candidates[fusion_key] = {"read1"}
        self.detector.fusion_assigned_pairs[fusion_key] = {"read1": ("BRCA1", "BRCA2")}
        self.detector.fusion_breakpoints[fusion_key] = {bp: 1}
        self.validator.filter_early_non_coding_genes()
        assert fusion_key in self.detector.fusion_assigned_pairs

    def test_filter_removes_non_coding(self):
        """Test that non-coding fusions are removed early."""
        fusion_key = "fusion1"
        bp = ("chr1", 1000, "chr2", 2000)
        self.detector.fusion_candidates[fusion_key] = {"read1"}
        self.detector.fusion_assigned_pairs[fusion_key] = {"read1": ("RPL5", "BRCA2")}
        self.detector.fusion_breakpoints[fusion_key] = {bp: 1}
        self.validator.filter_early_non_coding_genes()
        # The fusion should be removed from fusion_assigned_pairs
        assert fusion_key not in self.detector.fusion_assigned_pairs


class TestIntegration:
    """Integration tests combining multiple methods."""

    def setup_method(self):
        self.detector = MockDetector()
        self.validator = FusionValidator(self.detector)

    def test_full_pipeline(self):
        """Test the full validation pipeline."""
        # Create a set of test fusions
        self.detector.fusion_metadata = {
            "valid_fusion": {
                "left_gene": "BCR",
                "right_gene": "ABL1",
                "consensus_bp": ("chr9", 133589000, "chr22", 23632000),
                "is_valid": True,
                "support": 10,
                "reconstruction_ok": True,
                "realignment_hits": 2,
                "best_hit_mapq": 30
            },
            "non_coding_fusion": {
                "left_gene": "HIST1H1A",
                "right_gene": "BRCA2",
                "consensus_bp": ("chr1", 1000, "chr2", 2000),
                "is_valid": True,
                "support": 3
            }
        }
        # Apply filters
        self.validator.filter_non_coding_genes()
        # Check results
        assert "valid_fusion" in self.detector.fusion_metadata
        assert self.detector.fusion_metadata["valid_fusion"]["is_valid"] is True
        # The non-coding fusion should be removed
        assert "non_coding_fusion" not in self.detector.fusion_metadata
