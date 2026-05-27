import pytest
import logging
from unittest.mock import MagicMock, patch
from isoquant_lib.fusion_detector import FusionDetector, _CIGAR_CACHE, _ALIGNER_MAP_CACHE

logger = logging.getLogger('IsoQuant')


class TestFusionDetectorInitialization:
    """Test FusionDetector initialization."""

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.mp.Aligner')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_init_with_reference(self, mock_interval_index, mock_aligner, mock_featuredb):
        """Test initialization with reference FASTA file."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta="/path/to/reference.fa"
        )
        assert detector.bam_path == "/path/to/file.bam"
        assert detector.genedb_path == "/path/to/genes.db"
        assert detector.reference_fasta == "/path/to/reference.fa"
        assert detector.aligner is not None
        assert isinstance(detector.fusion_candidates, dict)
        assert isinstance(detector.fusion_breakpoints, dict)

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_init_without_reference(self, mock_interval_index, mock_featuredb):
        """Test initialization without reference FASTA file."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        assert detector.reference_fasta is None
        assert detector.aligner is None

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_init_creates_empty_collections(self, mock_interval_index, mock_featuredb):
        """Test that initialization creates empty collections."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        assert len(detector.fusion_candidates) == 0
        assert len(detector.fusion_breakpoints) == 0
        assert len(detector.fusion_metadata) == 0
        assert len(detector.fusion_assigned_pairs) == 0
        assert len(detector.fusion_read_scores) == 0


class TestCIGARParsing:
    """Test CIGAR string and tuple parsing."""
    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_aligned_len_from_cigartuples_basic(self, mock_interval_index, mock_featuredb):
        """Test aligned length from CIGAR tuples with basic operations."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        # CIGAR: 50M (50 matches)
        cigartuples = [(0, 50)]
        assert detector.aligned_len_from_cigartuples(cigartuples) == 50
        # CIGAR: 30M 10D 20M (50 aligned bases)
        cigartuples = [(0, 30), (2, 10), (0, 20)]
        assert detector.aligned_len_from_cigartuples(cigartuples) == 50
        # CIGAR: 20M 5I 10M (30 aligned bases, insertion doesn't count)
        cigartuples = [(0, 20), (1, 5), (0, 10)]
        assert detector.aligned_len_from_cigartuples(cigartuples) == 30

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_aligned_len_from_cigartuples_empty(self, mock_interval_index, mock_featuredb):
        """Test aligned length from empty CIGAR tuples."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        assert detector.aligned_len_from_cigartuples([]) == 0
        assert detector.aligned_len_from_cigartuples(None) == 0

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_aligned_len_from_cigarstring_basic(self, mock_interval_index, mock_featuredb):
        """Test aligned length from CIGAR strings."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        # Clear cache to ensure clean test
        _CIGAR_CACHE.clear()
        # Simple CIGAR: 50M
        assert detector.aligned_len_from_cigarstring("50M") == 50
        # Complex CIGAR: 30M10D20M
        assert detector.aligned_len_from_cigarstring("30M10D20M") == 50
        # With insertions: 20M5I10M
        assert detector.aligned_len_from_cigarstring("20M5I10M") == 30
        # Sequence match operators: 20=5I10X
        assert detector.aligned_len_from_cigarstring("20=5I10X") == 30

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_aligned_len_from_cigarstring_empty(self, mock_interval_index, mock_featuredb):
        """Test aligned length from empty/None CIGAR strings."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        _CIGAR_CACHE.clear()
        assert detector.aligned_len_from_cigarstring("") == 0
        assert detector.aligned_len_from_cigarstring(None) == 0

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_aligned_len_cigar_caching(self, mock_interval_index, mock_featuredb):
        """Test that CIGAR parsing results are cached."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        _CIGAR_CACHE.clear()
        cigar = "50M10D20M"
        result1 = detector.aligned_len_from_cigarstring(cigar)
        assert cigar in _CIGAR_CACHE
        result2 = detector.aligned_len_from_cigarstring(cigar)
        assert result1 == result2 == 70


class TestSoftClipDetection:
    """Test soft-clip detection from CIGAR."""

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_detect_softclip_from_cigartuples_left(self, mock_interval_index, mock_featuredb):
        """Test detection of left soft-clip from CIGAR tuples."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        # Mock read with left soft-clip: 60S40M
        read = MagicMock()
        read.cigartuples = [(4, 60), (0, 40)]
        side, length = detector.detect_softclip(read, min_len=50)
        assert side == "left"
        assert length == 60

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_detect_softclip_from_cigartuples_right(self, mock_interval_index, mock_featuredb):
        """Test detection of right soft-clip from CIGAR tuples."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        # Mock read with right soft-clip: 40M60S
        read = MagicMock()
        read.cigartuples = [(0, 40), (4, 60)]
        side, length = detector.detect_softclip(read, min_len=50)
        assert side == "right"
        assert length == 60

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_detect_softclip_both_ends_ambiguous(self, mock_interval_index, mock_featuredb):
        """Test that dual soft-clips return ambiguous (None, 0)."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        # Mock read with both ends soft-clipped: 60S40M60S
        read = MagicMock()
        read.cigartuples = [(4, 60), (0, 40), (4, 60)]
        side, length = detector.detect_softclip(read, min_len=50)
        assert side is None
        assert length == 0

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_detect_softclip_below_threshold(self, mock_interval_index, mock_featuredb):
        """Test that soft-clips below min_len threshold are ignored."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        # Mock read with small soft-clip: 20S80M
        read = MagicMock()
        read.cigartuples = [(4, 20), (0, 80)]
        side, length = detector.detect_softclip(read, min_len=50)
        assert side is None
        assert length == 0

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_detect_softclip_fallback_cigarstring(self, mock_interval_index, mock_featuredb):
        """Test soft-clip detection fallback to CIGAR string."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        # Mock read without cigartuples, use cigarstring
        read = MagicMock()
        read.cigartuples = None
        read.cigarstring = "60S40M"
        side, length = detector.detect_softclip(read, min_len=50)
        assert side == "left"
        assert length == 60


class TestSATagParsing:
    """Test SA tag parsing."""

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_parse_sa_entries_single(self, mock_interval_index, mock_featuredb):
        """Test parsing of single SA entry."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        sa_tag = "chr2,1000,+,50M,30,0"
        entries = detector.parse_sa_entries(sa_tag)
        assert len(entries) == 1
        assert entries[0][0] == "chr2"
        assert entries[0][1] == 1000
        assert entries[0][2] == "+"
        assert entries[0][3] == "50M"
        assert entries[0][4] == 30
        assert entries[0][5] == 0

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_parse_sa_entries_multiple(self, mock_interval_index, mock_featuredb):
        """Test parsing of multiple SA entries."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        sa_tag = "chr2,1000,+,50M,30,0;chr3,2000,-,40M,25,1"
        entries = detector.parse_sa_entries(sa_tag)
        assert len(entries) == 2
        assert entries[0][0] == "chr2"
        assert entries[0][1] == 1000
        assert entries[1][0] == "chr3"
        assert entries[1][1] == 2000
        assert entries[1][2] == "-"

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_parse_sa_entries_empty(self, mock_interval_index, mock_featuredb):
        """Test parsing of empty/None SA tag."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        assert detector.parse_sa_entries(None) == []
        assert detector.parse_sa_entries("") == []

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_parse_sa_entries_partial_fields(self, mock_interval_index, mock_featuredb):
        """Test parsing of SA entry with missing fields."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        # SA tag with fewer than 6 fields (should use fallback parsing)
        sa_tag = "chr2,1000,+,50M"
        entries = detector.parse_sa_entries(sa_tag)
        assert len(entries) == 1
        assert entries[0][0] == "chr2"
        assert entries[0][1] == 1000


class TestBreakpointEstimation:
    """Test breakpoint estimation."""

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_estimate_breakpoint_with_right_clip(self, mock_interval_index, mock_featuredb):
        """Test breakpoint estimation with right soft-clip."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        # Mock read on chr1, position 100-150
        read = MagicMock()
        read.reference_name = "chr1"
        read.reference_start = 99  # 0-based
        read.reference_end = 150
        # SA at chr2, position 1000
        left_chr, left_pos, right_chr, right_pos = detector.estimate_breakpoint(
            read, 1000, clip_side="right", sa_cigar="50M"
        )
        assert left_chr == "chr1"
        assert left_pos == 150  # read end
        assert right_chr == "chr1"
        assert right_pos == 1000  # sa_pos start

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_estimate_breakpoint_with_left_clip(self, mock_interval_index, mock_featuredb):
        """Test breakpoint estimation with left soft-clip."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        # Mock read on chr1, position 100-150
        read = MagicMock()
        read.reference_name = "chr1"
        read.reference_start = 99  # 0-based
        read.reference_end = 150
        # SA at chr2, position 1000
        left_chr, left_pos, right_chr, right_pos = detector.estimate_breakpoint(
            read, 1000, clip_side="left", sa_cigar=None
        )
        assert left_chr == "chr1"
        assert left_pos == 100  # read start + 1 (safe_reference_start converts to 1-based)


class TestGeneNameNormalization:
    """Test gene name normalization and resolution."""

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_normalize_gene_label_none_input(self, mock_interval_index, mock_featuredb):
        """Test normalization of None gene label."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        result = detector.normalize_gene_label(None)
        assert result == "intergenic"

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_normalize_gene_label_unresolved_identifier(self, mock_interval_index, mock_featuredb):
        """Test normalization of unresolved identifiers."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        # Mock resolve_gene_name to return an ENSG ID (unresolved)
        detector.resolve_gene_name = MagicMock(return_value="ENSG00000000001")
        detector.canonical_locus_name = MagicMock(return_value="ENSG00000000001")
        result = detector.normalize_gene_label("ENSG00000000001")
        assert result == "intergenic"  # Unresolved ENSGs become intergenic

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_has_antisense_suffix(self, mock_interval_index, mock_featuredb):
        """Test detection of antisense suffix."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        assert detector.has_antisense_suffix("BRCA1-AS1") is True
        assert detector.has_antisense_suffix("BRCA1-DT") is True
        assert detector.has_antisense_suffix("BRCA1-NAT") is True
        assert detector.has_antisense_suffix("BRCA1") is False
        assert detector.has_antisense_suffix(None) is False
        assert detector.has_antisense_suffix("") is False

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_strip_antisense_suffix(self, mock_interval_index, mock_featuredb):
        """Test stripping of antisense suffix."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        assert detector.strip_antisense_suffix("BRCA1-AS1") == "BRCA1"
        assert detector.strip_antisense_suffix("BRCA1-DT") == "BRCA1"
        assert detector.strip_antisense_suffix("BRCA1") == "BRCA1"
        assert detector.strip_antisense_suffix(None) is None
        assert detector.strip_antisense_suffix("") == ""


class TestFusionRecording:
    """Test fusion recording functionality."""

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_record_fusion_basic(self, mock_interval_index, mock_featuredb):
        """Test basic fusion recording."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        # Mock gene label normalization
        detector.normalize_gene_label = MagicMock(side_effect=lambda x: x)
        detector._is_mitochondrial_candidate = MagicMock(return_value=False)
        detector.record_fusion(
            context1="GENE1",
            context2="GENE2",
            read_name="read1",
            chrom1="chr1",
            pos1=1000,
            chrom2="chr2",
            pos2=2000
        )
        # Check fusion candidates recorded
        assert "GENE1--GENE2" in detector.fusion_candidates
        assert "read1" in detector.fusion_candidates["GENE1--GENE2"]
        # Check fusion metadata
        assert "GENE1--GENE2" in detector.fusion_metadata
        assert detector.fusion_metadata["GENE1--GENE2"]["support"] == 1

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_record_fusion_sorted_genes(self, mock_interval_index, mock_featuredb):
        """Test that genes are sorted alphabetically in fusion keys."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        detector.normalize_gene_label = MagicMock(side_effect=lambda x: x)
        detector._is_mitochondrial_candidate = MagicMock(return_value=False)
        # Record fusion with GENE2-GENE1 order
        detector.record_fusion(
            context1="GENE2",
            context2="GENE1",
            read_name="read1",
            chrom1="chr2",
            pos1=2000,
            chrom2="chr1",
            pos2=1000
        )
        # Should be sorted as GENE1--GENE2
        assert "GENE1--GENE2" in detector.fusion_candidates

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_record_fusion_with_scores(self, mock_interval_index, mock_featuredb):
        """Test fusion recording with per-read scores."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        detector.normalize_gene_label = MagicMock(side_effect=lambda x: x)
        detector._is_mitochondrial_candidate = MagicMock(return_value=False)
        detector.record_fusion(
            context1="GENE1",
            context2="GENE2",
            read_name="read1",
            chrom1="chr1",
            pos1=1000,
            chrom2="chr2",
            pos2=2000,
            left_score=0.9,
            right_score=0.8
        )
        assert "GENE1--GENE2" in detector.fusion_read_scores
        assert detector.fusion_read_scores["GENE1--GENE2"]["read1"] == (0.9, 0.8)

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_record_fusion_mitochondrial_filtered(self, mock_interval_index, mock_featuredb):
        """Test that mitochondrial fusions are filtered."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        detector._is_mitochondrial_candidate = MagicMock(return_value=True)
        detector.record_fusion(
            context1="GENE1",
            context2="GENE2",
            read_name="read1",
            chrom1="chrM",
            pos1=1000,
            chrom2="chr2",
            pos2=2000
        )
        # Mitochondrial fusion should not be recorded
        assert len(detector.fusion_candidates) == 0


class TestUtilityFunctions:

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_safe_reference_start(self, mock_interval_index, mock_featuredb):
        """Test safe reference start extraction."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        # Mock read with valid reference_start
        read = MagicMock()
        read.reference_start = 99
        assert detector.safe_reference_start(read) == 100  # Converts to 1-based
        # Mock read with None reference_start
        read.reference_start = None
        assert detector.safe_reference_start(read) is None

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_safe_gene_token(self, mock_interval_index, mock_featuredb):
        """Test safe gene token conversion."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        assert detector._safe_gene_token("BRCA1") == "BRCA1"
        assert detector._safe_gene_token(None) == "intergenic"
        assert detector._safe_gene_token("") == "intergenic"
        assert detector._safe_gene_token("   ") == "intergenic"

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_compute_aligned_length(self, mock_interval_index, mock_featuredb):
        """Test aligned length computation from read."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        # Mock read with cigartuples
        read = MagicMock()
        read.cigartuples = [(0, 50)]
        read.cigarstring = None
        assert detector.compute_aligned_length(read) == 50
        # Mock read with cigarstring
        read.cigartuples = None
        read.cigarstring = "50M"
        _CIGAR_CACHE.clear()
        result = detector.compute_aligned_length(read)
        assert result == 50

class TestStateManagement:
    """Test state management."""

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_clear_state(self, mock_interval_index, mock_featuredb):
        """Test clearing of detector state."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        # Populate with some data
        detector.fusion_candidates["GENE1--GENE2"].add("read1")
        detector.fusion_breakpoints["GENE1--GENE2"][("chr1", 1000, "chr2", 2000)] = 5
        detector.fusion_metadata["GENE1--GENE2"] = {"support": 1}
        # Clear state
        detector.clear_state()
        # Verify all collections are empty
        assert len(detector.fusion_candidates) == 0
        assert len(detector.fusion_breakpoints) == 0
        assert len(detector.fusion_metadata) == 0
        assert len(detector.fusion_assigned_pairs) == 0
        assert len(detector.fusion_read_scores) == 0


class TestClusterBreakpoints:
    """Test breakpoint clustering."""

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_cluster_breakpoints_single(self, mock_interval_index, mock_featuredb):
        """Test clustering of single breakpoint."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        bp_counts = {
            ("chr1", 1000, "chr2", 2000): 10
        }
        result = detector.cluster_breakpoints(bp_counts, window=100)
        assert result is not None
        (c1, p1, c2, p2), total = result
        assert c1 == "chr1"
        assert c2 == "chr2"
        assert total == 10

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_cluster_breakpoints_multiple(self, mock_interval_index, mock_featuredb):
        """Test clustering of multiple nearby breakpoints."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        bp_counts = {
            ("chr1", 1000, "chr2", 2000): 5,
            ("chr1", 1050, "chr2", 2050): 8,  # Close to first
            ("chr1", 2000, "chr2", 3000): 3,  # Far from others
        }
        result = detector.cluster_breakpoints(bp_counts, window=100)
        assert result is not None
        (c1, p1, c2, p2), total = result
        # Should cluster the two closer breakpoints
        assert total == 13

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_cluster_breakpoints_empty(self, mock_interval_index, mock_featuredb):
        """Test clustering of empty breakpoint dictionary."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        result = detector.cluster_breakpoints({}, window=100)
        assert result is None


class TestReadFiltering:
    """Test read filtering."""

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_passes_read_filters_unmapped(self, mock_interval_index, mock_featuredb):
        """Test filtering of unmapped reads."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        read = MagicMock()
        read.is_unmapped = True
        read.is_secondary = False
        read.is_supplementary = False
        passes, clip_side, clip_len = detector._passes_read_filters(read, 50, 20)
        assert passes is False

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_passes_read_filters_short_alignment(self, mock_interval_index, mock_featuredb):
        """Test filtering of short alignments."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        read = MagicMock()
        read.is_unmapped = False
        read.is_secondary = False
        read.is_supplementary = False
        read.cigartuples = [(0, 30)]  # 30bp < 50bp min
        read.cigarstring = None
        read.mapping_quality = 30
        passes, clip_side, clip_len = detector._passes_read_filters(read, 50, 20)
        assert passes is False

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_passes_read_filters_valid(self, mock_interval_index, mock_featuredb):
        """Test valid read passes filters."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        read = MagicMock()
        read.is_unmapped = False
        read.is_secondary = False
        read.is_supplementary = False
        read.cigartuples = [(0, 50)]  # 50bp >= 50bp min
        read.cigarstring = None
        read.mapping_quality = 30
        passes, clip_side, clip_len = detector._passes_read_filters(read, 50, 20)
        assert passes is True


class TestMitochondrialFiltering:
    """Test mitochondrial fusion filtering."""

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_is_mitochondrial_chrM(self, mock_interval_index, mock_featuredb):
        """Test detection of mitochondrial chromosome."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        assert detector._is_mitochondrial_candidate("chrM", "chr1", None, None) is True
        assert detector._is_mitochondrial_candidate("chr1", "MT", None, None) is True

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_is_mitochondrial_gene_name(self, mock_interval_index, mock_featuredb):
        """Test detection of mitochondrial gene names."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        assert detector._is_mitochondrial_candidate("chr1", "chr2", "MT-ND1", None) is True
        assert detector._is_mitochondrial_candidate("chr1", "chr2", None, "MT-TS1") is True

    @patch('isoquant_lib.fusion_detector.gffutils.FeatureDB')
    @patch('isoquant_lib.fusion_detector.GenomicIntervalIndex')
    def test_is_mitochondrial_not_mitochondrial(self, mock_interval_index, mock_featuredb):
        """Test non-mitochondrial fusions."""
        mock_db = MagicMock()
        mock_featuredb.return_value = mock_db
        detector = FusionDetector(
            bam_path="/path/to/file.bam",
            gene_db_path="/path/to/genes.db",
            reference_fasta=None
        )
        assert detector._is_mitochondrial_candidate("chr1", "chr2", "BRCA1", "TP53") is False

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
