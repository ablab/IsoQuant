import logging
import re
from collections import defaultdict
from typing import Any, Optional

import gffutils
import pysam
import mappy as mp
from intervaltree import IntervalTree
from .fusion_validator import FusionValidator
from .genomic_interval_index import GenomicIntervalIndex

logger = logging.getLogger('IsoQuant')
ANTISENSE_SUFFIX_RE = re.compile(r"-(AS\d+|DT|DIVERGENT|NAT)$", re.IGNORECASE)

# Module-level caches
_CONTEXT_CACHE = {}  # (chrom, pos) → (gene_name, region_type)
_GENE_ASSIGNMENT_CACHE = {}  # (chrom, pos) → (gene_name, score)
_ALIGNER_MAP_CACHE = {}  # seq → tuple of hits
_CIGAR_CACHE = {}  # cigar_string → aligned_length

class FusionDetector:
    def __init__(self, bam_path: str, gene_db_path: str, reference_fasta: Optional[str]) -> None:
        if mp is None:
            raise ImportError(
                "mappy is required for fusion detection. "
                "Install with `pip install -r requirements_fusion.txt`."
            )
        self.bam_path = bam_path
        self.genedb_path = gene_db_path
        self.db = gffutils.FeatureDB(gene_db_path, keep_order=True)
        self.fusion_candidates: dict[str, set[str]] = defaultdict(set)
        # store breakpoint counts per fusion key: {(chr1,pos1,chr2,pos2): count}
        self.fusion_breakpoints: dict[str, dict[tuple, int]] = defaultdict(lambda: defaultdict(int))
        # optional reference and aligner for soft-clip realignment
        self.reference_fasta = reference_fasta
        try:
            self.aligner = mp.Aligner(reference_fasta) if reference_fasta else None
        except Exception:
            self.aligner = None
        self.fusion_metadata = {}
        # store per-read assigned raw gene pairs for each fusion key: {fusion_key: {read_name: (left,right)}}
        self.fusion_assigned_pairs = defaultdict(dict)
        # store per-read scores for each fusion key: {fusion_key: {read_name: (left_score, right_score)}}
        self.fusion_read_scores = defaultdict(dict)
        # cache for resolved names: id_or_symbol -> gene_symbol (if found)
        self._resolved_name_cache = {}
        # cache for symbol -> biotype mapping (built on-demand)
        self._symbol_biotype_cache = {}
        # cache for exons by gene_id: {gene_id: [(start, end), ...]}
        self.exon_cache = {}
        # Build interval tree index for fast coordinate-to-gene/exon mapping
        logger.info("Initializing genomic interval index for efficient database queries...")
        if IntervalTree is not None:
            self.interval_index = GenomicIntervalIndex(self.db)
        else:
            logger.warning("intervaltree not available; falling back to gffutils queries")
            self.interval_index = None
        # Build exon cache upfront for efficient assign_fusion_gene lookups
        self._build_exon_cache()
        self.debug = False

    def safe_reference_start(self, read):
        try:
            return read.reference_start + 1 if read.reference_start is not None else None
        except Exception:
            return None

    def compute_aligned_length(self, read):
        try:
            if getattr(read, "cigartuples", None):
                return self.aligned_len_from_cigartuples(read.cigartuples)
            return self.aligned_len_from_cigarstring(getattr(read, "cigarstring", None))
        except Exception:
            return 0

    def clear_state(self) -> None:
        """Clear all fusion detection state for processing a new BAM file."""
        self.fusion_candidates.clear()
        self.fusion_breakpoints.clear()
        self.fusion_metadata.clear()
        self.fusion_assigned_pairs.clear()
        self.fusion_read_scores.clear()

    def _build_exon_cache(self):
        # Cache ordered exon spans for each gene to support exon-boundary
        # and exon-distance calculations. This is used in gene scoring and
        # is distinct from interval-tree–based coordinate lookups

        if self.exon_cache:
            return  # Already built
        logger.debug("Building exon cache from genedb...")
        try:
            # Use interval tree to iterate genes if available (more efficient than db.features_of_type)
            if self.interval_index is not None:
                genes = self.interval_index.get_all_genes()
            else:
                genes = list(self.db.features_of_type('gene'))
            for gene in genes:
                gene_id = getattr(gene, 'id', None)
                if not gene_id:
                    continue
                try:
                    exons = list(self.db.children(gene, featuretype='exon', order_by='start'))
                except Exception:
                    logger.debug(
                        "Failed to fetch exons for gene %s",
                        getattr(gene, "id", gene),
                        exc_info=True
                    )
        except Exception as e:
            logger.warning(f"Failed to build exon cache: {e}")

    def _get_cached_exons(self, gene):
        # Retrieve exons for a gene, preferring cached version.
        gene_id = getattr(gene, 'id', None)
        if not gene_id:
            return []
        # Check cache first
        if gene_id in self.exon_cache:
            return self.exon_cache[gene_id]
        # Fall back to DB if not cached (shouldn't happen after _build_exon_cache)
        try:
            exons = list(self.db.children(gene, featuretype='exon', order_by='start'))
            exon_list = [(int(ex.start), int(ex.end)) for ex in exons]
            # Cache for future use
            self.exon_cache[gene_id] = exon_list
            return exon_list
        except Exception:
            return []

    def _context_query(self, chrom, pos):
        # Query genomic context at position (chrom, pos).
        # Uses module-level cache to avoid holding self reference.
        # Uses interval tree if available for O(log n) performance, otherwise falls back to gffutils.
        cache_key = (chrom, pos)
        if cache_key in _CONTEXT_CACHE:
            return _CONTEXT_CACHE[cache_key]        
        try:
            genes = []
            exons = []
            # Try interval tree first if available
            if self.interval_index is not None:
                genes = self.interval_index.get_genes_at(chrom, pos)
                exons = self.interval_index.get_exons_at(chrom, pos)
            else:
                # Fallback to gffutils
                genes = list(self.db.region(region=(chrom, pos, pos), featuretype='gene'))
                exons = list(self.db.region(region=(chrom, pos, pos), featuretype='exon'))
            if genes:
                gene_name = genes[0].attributes.get('gene_name', [genes[0].id])[0]
                # Are we in an exon of that gene?
                if exons:
                    result = (gene_name, "exonic")
                else:
                    result = (gene_name, "intronic")
            else:
                # Not inside any gene; check if overlapping any exon
                result = (None, ("exonic" if exons else "intergenic"))
        except Exception:
            logger.debug(
                "Context query failed for %s:%s",
                chrom, pos,
                exc_info=True
            )
            result = (None, "unknown")        
        # Cache with size limit
        if len(_CONTEXT_CACHE) > 200000:
            # Remove oldest entries (FIFO)
            for _ in range(20000):
                _CONTEXT_CACHE.pop(next(iter(_CONTEXT_CACHE)))
        _CONTEXT_CACHE[cache_key] = result
        return result

    def _get_parent_gene_symbol(self, feature):
        # Extract gene symbol from a feature, handling gene/transcript/other feature types.
        if feature.featuretype == 'gene':
            # For gene features, prefer gene_name attribute
            return feature.attributes.get('gene_name', [feature.id])[0]
        # For transcript or other features, try to find parent gene
        if 'gene_id' in feature.attributes:
            gid = feature.attributes['gene_id'][0]
            try:
                parent = self.db[gid]
                return parent.attributes.get('gene_name', [parent.id])[0]
            except Exception:
                return gid
        elif 'Parent' in feature.attributes:
            parent_id = feature.attributes['Parent'][0]
            try:
                parent = self.db[parent_id]
                return parent.attributes.get('gene_name', [parent.id])[0]
            except Exception:
                return parent_id
        else:
            # Fallback: try to extract gene_name from feature attributes
            return feature.attributes.get('gene_name', [feature.id])[0] if hasattr(feature, 'attributes') else feature.id

    def _lookup_feature_by_id(self, name_or_id):
        # Try direct DB lookup by ID. Returns (feature, resolved_name) or (None, original_id).
        try:
            feat = self.db[name_or_id]
            resolved = self._get_parent_gene_symbol(feat)
            return feat, resolved
        except Exception:
            return None, name_or_id

    def _fallback_gene_lookup(self, name_or_id: str) -> str:
        # Last-resort resolution of a gene identifier to a gene symbol.
        try:
            genes = list(self.db.features_of_type("gene", id=name_or_id))
            if not genes:
                return name_or_id
            gene = genes[0]
            attrs = getattr(gene, "attributes", {}) or {}
            return attrs.get("gene_name", [gene.id])[0]
        except (AttributeError, KeyError, TypeError):
            logger.debug(
                "Fallback gene lookup failed for identifier %r",
                name_or_id,
                exc_info=True,
            )
            return name_or_id

    def resolve_gene_name(self, name_or_id: Optional[str]) -> Optional[str]:
        """
        Resolve Ensembl or other identifiers (e.g. ENS*, RP11*) to gene symbols
        when possible. Falls back to returning the input unchanged.
        """
        if not name_or_id:
            return name_or_id
        # Cache hit
        cached = self._resolved_name_cache.get(name_or_id)
        if cached is not None:
            return cached
        # Attempt direct lookup
        feat, resolved = self._lookup_feature_by_id(name_or_id)
        if feat is None:
            resolved = self._fallback_gene_lookup(name_or_id)
        # Cache and return
        self._resolved_name_cache[name_or_id] = resolved
        return resolved

    def has_antisense_suffix(self, gene_name: Optional[str]) -> bool:
        """Return True if ``gene_name`` ends with an antisense/regulatory suffix (AS1, DT, etc.)."""
        if not gene_name or not isinstance(gene_name, str):
            return False
        return ANTISENSE_SUFFIX_RE.search(gene_name) is not None

    def strip_antisense_suffix(self, gene_name: Optional[str]) -> Optional[str]:
        """Strip antisense/regulatory suffix from a gene name and return the bare locus name."""
        if not gene_name or not isinstance(gene_name, str):
            return gene_name
        return ANTISENSE_SUFFIX_RE.sub('', gene_name)

    def _lookup_gene_by_name(self, gene_name, chrom, pos):
        # Lookup mode: find gene feature by name at genomic position.
        # Returns (feature, gstart, gend) or (None, None, None) on failure.
        try:
            gene_features = list(self.db.region(
                region=(chrom, max(1, pos - 500), pos + 500),
                featuretype="gene"
            ))
            if not gene_features:
                return None, None, None
            # Find gene by name
            g = None
            for gf in gene_features:
                attrs = getattr(gf, "attributes", {}) or {}
                gname = attrs.get("gene_name", [getattr(gf, "id", None)])[0]
                if gname == gene_name:
                    g = gf
                    break
            if not g:
                # Fallback to closest gene
                g = gene_features[0]
            gstart = int(g.start)
            gend = int(g.end)
            return g, gstart, gend
        except Exception as e:
            logger.debug(f"Error looking up gene {gene_name} at {chrom}:{pos}: {e}")
            return None, None, None

    def _compute_exonic_metrics(self, g, pos):
        # Compute exon-based distance metrics.
        try:
            exons = self._get_cached_exons(g)
            gstart = int(getattr(g, "start", 0))
            gend = int(getattr(g, "end", 0))
            # Determine if position is inside an exon
            exonic_hit = any(es <= pos <= ee for es, ee in exons)
            # Compute distances to exon boundaries
            distances = []
            for es, ee in exons:
                distances.append(abs(pos - es))
                distances.append(abs(pos - ee))
            exon_min_dist = min(distances) if distances else None
            boundary_min_dist = exon_min_dist  # Same metric, used with different thresholds
            # Compute body distance (distance to gene boundaries)
            body_dist = 0
            if pos < gstart:
                body_dist = gstart - pos
            elif pos > gend:
                body_dist = pos - gend
            return exonic_hit, exon_min_dist, boundary_min_dist, body_dist, exons
        except Exception:
            return False, None, None, 0, []

    def _compute_gene_score(self, target, pos, chrom=None, exon_min_dist=None,
                            boundary_min_dist=None, exonic_hit=None, body_dist=None,
                            gstart=None, gend=None, fusion_partner_gene=None,
                            fusion_partner_chrom=None, fusion_partner_pos=None):
        # Simplified gene scoring for long reads: assign gene if breakpoint is inside gene bounds.
        # Returns 1.0 if pos is inside gene, 0.0 if intergenic.
        if isinstance(target, str):
            # Lookup mode: find gene by name at position
            gene_name = target
            if chrom is None:
                raise ValueError("chrom required when target is a gene name")
            g, gstart, gend = self._lookup_gene_by_name(gene_name, chrom, pos)
            if g is None:
                return 0.0
        else:
            # Pre-computed mode: target is already a feature object
            g = target
            if gstart is None or gend is None:
                raise ValueError("gstart and gend required for pre-computed mode")
        # Simple rule: assign score based on whether position is inside gene bounds
        if gstart <= pos <= gend:
            return 1.0  # Position is inside gene
        else:
            return 0.0  # Position is intergenic

    @staticmethod
    def reverse_complement(seq: str) -> str:
        """Return the reverse complement of a DNA sequence (ACGTN alphabet, case preserved on N)."""
        comp = str.maketrans("ACGTN", "TGCAN")
        return seq.translate(comp)[::-1]

    def _get_genes_at(self, chrom, pos, window):
        # Query genes at location using interval tree or gffutils.
        try:
            if self.interval_index is None:
                logger.warning("IntervalTree not available; assign_fusion_gene requires it for efficiency")
                return None
            return self.interval_index.get_genes_at(chrom, pos, window=window)
        except Exception:
            return None

    def _collect_exonic_genes(self, genes, pos):
        # Return a list of genes for which ``pos`` falls inside at least one exon.
        exonic_genes = []
        for gene in genes:
            attrs = getattr(gene, "attributes", {}) or {}
            gene_name = attrs.get("gene_name", [getattr(gene, "id", None)])[0]
            gene_type = (
                attrs.get("gene_type", [None])[0]
                or attrs.get("gene_biotype", [None])[0]
                or attrs.get("transcript_biotype", [None])[0]
            )
            try:
                gene_start = int(gene.start)
                gene_end = int(gene.end)
                exons = self._get_cached_exons(gene)
            except (AttributeError, TypeError, ValueError):
                logger.debug(
                    "Failed to extract gene or exon information for %r",
                    gene,
                    exc_info=True,
                )
                continue
            if any(start <= pos <= end for start, end in exons):
                exonic_genes.append((gene, gene_name, gene_type, gene_start, gene_end, exons))
        return exonic_genes


    def _compute_gene_distances(self, g, pos):
        # Compute exon and body distances for a gene relative to pos.
        exon_min_dist = None
        boundary_min_dist = None
        exonic_hit = False
        try:
            exons = self._get_cached_exons(g)
            for es, ee in exons:
                # Distance to exon span
                if es <= pos <= ee:
                    d_span = 0
                    exonic_hit = True
                else:
                    d_span = min(abs(pos - es), abs(pos - ee))
                exon_min_dist = d_span if exon_min_dist is None else min(exon_min_dist, d_span)
                # Distance to exon boundary
                d_bound = min(abs(pos - es), abs(pos - ee))
                boundary_min_dist = d_bound if boundary_min_dist is None else min(boundary_min_dist, d_bound)
        except Exception:
            exons = []
        # Body distance (distance to gene boundaries)
        gstart = int(getattr(g, "start", 0))
        gend = int(getattr(g, "end", 0))
        if gstart <= pos <= gend:
            body_dist = 0
        else:
            body_dist = min(abs(pos - gstart), abs(pos - gend))
        return exon_min_dist, boundary_min_dist, exonic_hit, body_dist, exons

    def _build_gene_key(self, gtype, score, exonic_hit, boundary_min_dist, exon_min_dist, body_dist, glen, gname):
        # Build a comparison tuple for gene scoring. Higher tuple = better gene."""
        INF = 10**9
        bnd = boundary_min_dist if boundary_min_dist is not None else INF
        exd = exon_min_dist if exon_min_dist is not None else INF
        bod = body_dist if body_dist else 0
        return (
            1 if gtype == "protein_coding" else 0,  # 1) PRIORITY: protein-coding
            score,                                   # 2) main score
            #bool(exonic_hit),                        # 3) prefer exonic
            -bnd,                                    # 4) nearer exon boundary
            -exd,                                    # 5) nearer exon span
            -bod,                                    # 6) nearer gene body
            -min(glen, 500000),                      # 7) gentle preference for longer locus
            gname or ""                              # 8) deterministic
        )

    def _score_exonic_genes(self, exonic_genes, pos, chrom=None, fusion_partner_gene=None,
                            fusion_partner_chrom=None, fusion_partner_pos=None):
        # Score multiple exonic genes and return the best one, or single exonic gene.
        # Returns (gene_name, score)
        if len(exonic_genes) == 1:
            g = exonic_genes[0][0]
            gname = exonic_genes[0][1]
            # Compute score for fast path
            boundary_distances = []
            exons = exonic_genes[0][5]
            for start, end in exons:
                boundary_distances.append(min(abs(pos - start), abs(pos - end)))
            boundary_min_dist = min(boundary_distances) if boundary_distances else 0
            score = self._compute_gene_score(
                g, pos, chrom=chrom, exon_min_dist=0, boundary_min_dist=boundary_min_dist,
                exonic_hit=True, body_dist=0, gstart=exonic_genes[0][3], gend=exonic_genes[0][4],
                fusion_partner_gene=fusion_partner_gene,
                fusion_partner_chrom=fusion_partner_chrom,
                fusion_partner_pos=fusion_partner_pos
            )
            return gname, score
        best_gene = None
        best_score = None
        best_key = None
        for g, gname, gtype, gstart, gend, exons in exonic_genes:
            boundary_distances = [min(abs(pos - start), abs(pos - end)) for start, end in exons]
            boundary_min_dist = min(boundary_distances) if boundary_distances else 0
            score = self._compute_gene_score(
                g, pos, chrom=chrom, exon_min_dist=0, boundary_min_dist=boundary_min_dist,
                exonic_hit=True, body_dist=0, gstart=gstart, gend=gend,
                fusion_partner_gene=fusion_partner_gene,
                fusion_partner_chrom=fusion_partner_chrom,
                fusion_partner_pos=fusion_partner_pos
            )
            glen = max(1, gend - gstart + 1)
            key = self._build_gene_key(gtype, score, True, boundary_min_dist, 0, 0, glen, gname)
            if best_key is None or key > best_key:
                best_key = key
                best_gene = gname or getattr(g, "id", None)
                best_score = score
        return best_gene, best_score if best_score is not None else 0.0

    def _score_intronic_genes(self, genes, pos, chrom=None, fusion_partner_gene=None,
                              fusion_partner_chrom=None, fusion_partner_pos=None):
        # Score intronic/intergenic genes and return the best one.
        # Returns (gene_name, score)
        best_gene = None
        best_score = None
        best_key = None
        for g in genes:
            attrs = getattr(g, "attributes", {}) or {}
            gname = attrs.get("gene_name", [getattr(g, "id", None)])[0]
            gtype = (attrs.get("gene_type", [None])[0]
                    or attrs.get("gene_biotype", [None])[0]
                    or attrs.get("transcript_biotype", [None])[0])
            gstart = int(getattr(g, "start", 0))
            gend = int(getattr(g, "end", 0))
            exon_min_dist, boundary_min_dist, exonic_hit, body_dist, exons = self._compute_gene_distances(g, pos)
            score = self._compute_gene_score(
                g, pos, chrom=chrom, exon_min_dist=exon_min_dist, boundary_min_dist=boundary_min_dist,
                exonic_hit=exonic_hit, body_dist=body_dist, gstart=gstart, gend=gend,
                fusion_partner_gene=fusion_partner_gene,
                fusion_partner_chrom=fusion_partner_chrom,
                fusion_partner_pos=fusion_partner_pos
            )
            glen = max(1, gend - gstart + 1)
            key = self._build_gene_key(gtype, score, exonic_hit, boundary_min_dist, exon_min_dist, body_dist, glen, gname)
            if best_key is None or key > best_key:
                best_key = key
                best_gene = gname or getattr(g, "id", None)
                best_score = score
        return best_gene, best_score if best_score is not None else 0.0

    def assign_fusion_gene(self, chrom: str, pos: int, window: int = 1000,
                            fusion_partner_gene: Optional[str] = None,
                            fusion_partner_chrom: Optional[str] = None,
                            fusion_partner_pos: Optional[int] = None) -> tuple[Optional[str], float]:
        """Assign the best gene at a genomic location using exonic priority and scoring.

        Returns ``(gene_name, score)`` or ``(None, 0.0)`` if no genes are found.
        """
        genes = self._get_genes_at(chrom, pos, window)
        if not genes:
            return None, 0.0
        # FAST PATH: Try exonic genes first (position inside an exon)
        exonic_genes = self._collect_exonic_genes(genes, pos)
        if exonic_genes:
            return self._score_exonic_genes(exonic_genes, pos, chrom=chrom,
                                           fusion_partner_gene=fusion_partner_gene,
                                           fusion_partner_chrom=fusion_partner_chrom,
                                           fusion_partner_pos=fusion_partner_pos)
        # SLOW PATH: No exonic hits; use intronic/intergenic scoring
        return self._score_intronic_genes(genes, pos, chrom=chrom,
                                         fusion_partner_gene=fusion_partner_gene,
                                         fusion_partner_chrom=fusion_partner_chrom,
                                         fusion_partner_pos=fusion_partner_pos)

    def assign_fusion_gene_cached(self, chrom, pos):
        # Cache gene assignment results at (chrom, pos) using module-level dict
        cache_key = (chrom, pos)
        if cache_key in _GENE_ASSIGNMENT_CACHE:
            return _GENE_ASSIGNMENT_CACHE[cache_key]
        result = self.assign_fusion_gene(chrom, pos)
        # Cache with size limit to prevent unbounded growth
        if len(_GENE_ASSIGNMENT_CACHE) > 500000:
            # Remove oldest entries (FIFO)
            for _ in range(50000):
                _GENE_ASSIGNMENT_CACHE.pop(next(iter(_GENE_ASSIGNMENT_CACHE)))
        _GENE_ASSIGNMENT_CACHE[cache_key] = result
        return result

    def _passes_read_filters(self, read, min_al_len_primary, min_sa_mapq):
        # Check if read passes basic quality filters.
        # Returns (passes, clip_side, clip_len).
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            return False, None, 0
        prim_al_len = self.compute_aligned_length(read)
        if prim_al_len < min_al_len_primary:
            return False, None, 0
        prim_mapq = getattr(read, "mapping_quality", 0)
        if prim_mapq < min_sa_mapq:
            return False, None, 0
        clip_side, clip_len = self.detect_softclip(read)
        return True, clip_side, clip_len

    def _filter_sa_entries(self, sa_entries, min_al_len_sa, min_sa_mapq):
        # Filter & prioritize SA entries by mapq and aligned length.
        # Returns list of (chr, pos, strand, cigar, mapq, nm, al_len) tuples.
        filtered = []
        for (sa_chr, sa_pos, sa_strand, sa_cigar, sa_mapq, sa_nm) in sa_entries:
            if sa_pos is None or sa_cigar is None:
                continue
            sa_al_len = self.aligned_len_from_cigarstring(sa_cigar)
            if sa_al_len < min_al_len_sa or sa_mapq < min_sa_mapq:
                continue
            filtered.append((sa_chr, sa_pos, sa_strand, sa_cigar, sa_mapq, sa_nm, sa_al_len))
        # Keep top-N by sa_mapq or aligned length to limit work
        if len(filtered) > 3:
            filtered.sort(key=lambda x: (x[4], x[6]), reverse=True)
            filtered = filtered[:3]
        return filtered

    def _process_sa_entries(self, read, sa_entries, clip_side, jitter_window, min_sa_mapq):
        # Process all SA entries from a read and record fusion candidates.
        # Returns True if at least one fusion was recorded, False otherwise.
        seen_pairs = set()
        sa_used = False
        for (sa_chr, sa_pos, sa_strand, sa_cigar, sa_mapq, sa_nm, sa_al_len) in sa_entries:
            # Estimate breakpoints with clip guidance
            left_chr, left_pos, right_chr, right_pos = self.estimate_breakpoint(
                read, sa_pos, clip_side=clip_side, sa_cigar=sa_cigar
            )
            if left_pos is None or right_pos is None:
                continue
            key = (left_chr, (left_pos // jitter_window),
                sa_chr, (right_pos // jitter_window))
            if key in seen_pairs:
                continue
            seen_pairs.add(key)
            left_gene, left_score = self.assign_fusion_gene_cached(left_chr, left_pos)
            right_gene, right_score = self.assign_fusion_gene_cached(sa_chr, right_pos)
            if left_gene is None or right_gene is None:
                continue
            if left_gene != right_gene:
                # Check soft-clip orientation consistency
                # if clip_side and not self._is_softclip_orientation_valid(read, clip_side, sa_chr, sa_pos):
                    # continue
                self.record_fusion(
                    left_gene, right_gene, read.query_name,
                    left_chr, left_pos, sa_chr, right_pos,
                    left_score=left_score, right_score=right_score
                )
                sa_used = True
        return sa_used

    def detect_softclip(self, read, min_len: int = 50) -> tuple[Optional[str], int]:
        """Detect a single large terminal soft-clip; returns ``(side, length)`` or ``(None, 0)``."""
        cigartuples = getattr(read, "cigartuples", None)
        if cigartuples:
            left_clip = (
                cigartuples[0][0] == 4 and cigartuples[0][1] >= min_len
            )
            right_clip = (
                cigartuples[-1][0] == 4 and cigartuples[-1][1] >= min_len
            )
            # Ambiguous: soft-clips on both ends
            if left_clip and right_clip:
                return None, 0
            if left_clip:
                return "left", cigartuples[0][1]
            if right_clip:
                return "right", cigartuples[-1][1]
            return None, 0
        # Fallback to cigarstring (rare in pysam, but safe)
        cigar = getattr(read, "cigarstring", "") or ""
        m_left = re.match(r'^(\d+)S', cigar)
        m_right = re.search(r'(\d+)S$', cigar)
        left_len = int(m_left.group(1)) if m_left else 0
        right_len = int(m_right.group(1)) if m_right else 0
        # Ambiguous dual clip
        if left_len >= min_len and right_len >= min_len:
            return None, 0
        if left_len >= min_len:
            return "left", left_len
        if right_len >= min_len:
            return "right", right_len
        return None, 0

    def aligned_len_from_cigartuples(self, cigartuples) -> int:
        """Sum the read-aligned operations (M, =, X) from a list of CIGAR tuples."""
        if not cigartuples:
            return 0
        return sum(l for op, l in cigartuples if op in (0, 7, 8))
    
    def aligned_len_from_cigarstring(self, cigar: Optional[str]) -> int:
        """Sum the read-aligned operations (M, =, X) from a CIGAR string."""
        return self._cached_aligned_len_from_cigarstring(cigar)

    @staticmethod
    def _cached_aligned_len_from_cigarstring(cigar):
        if not cigar:
            return 0
        cached = _CIGAR_CACHE.get(cigar)
        if cached is not None:
            return cached
        total = 0
        for m in re.finditer(r'(\d+)([MIDNSHP=XB])', cigar):
            length = int(m.group(1))
            op = m.group(2)
            if op in ("M", "=", "X"):
                total += length
        _CIGAR_CACHE[cigar] = total
        return total

    def _cached_aligner_map(self, seq):
        # Cache aligner.map() results for identical sequences using module-level dict
        if seq in _ALIGNER_MAP_CACHE:
            return _ALIGNER_MAP_CACHE[seq]
        try:
            hits = tuple(self.aligner.map(seq))
            result = hits
        except Exception:
            result = ()
        # Cache with size limit
        if len(_ALIGNER_MAP_CACHE) > 10000:
            # Remove oldest entries (FIFO)
            for _ in range(1000):
                _ALIGNER_MAP_CACHE.pop(next(iter(_ALIGNER_MAP_CACHE)))
        _ALIGNER_MAP_CACHE[seq] = result
        return result

    def _is_softclip_orientation_valid(self, read, clip_side, sa_chr, sa_pos, tol=10):
        if clip_side is None:
            return True
        if read.reference_name != sa_chr:
            return True
        prim_start = self.safe_reference_start(read)
        prim_end = read.reference_end
        if prim_start is None or prim_end is None:
            return True
        if not read.is_reverse:
            if clip_side == "right":
                return sa_pos > prim_end - tol
            elif clip_side == "left":
                return sa_pos < prim_start + tol
        else:
            if clip_side == "left":
                return sa_pos > prim_end - tol
            elif clip_side == "right":
                return sa_pos < prim_start + tol
        return True

    def realign_clipped_seq(self, seq):
        # Realign a clipped sequence using the aligner.
        if not seq or not self.aligner:
            return None
        try:
            hits = self._cached_aligner_map(seq)
            for hit in hits:
                hit_len = (getattr(hit, "r_en", None) - getattr(hit, "r_st", None)) if getattr(hit, "r_en", None) else getattr(hit, "mlen", None) or getattr(hit, "alen", None) or 0
                if hit_len and hit_len >= 50:
                    return hit  # first good hit
            return None
        except Exception:
            return None

    def _safe_gene_token(self, g) -> str:
        """Convert ``None`` or empty/whitespace strings to ``'intergenic'``."""
        if not g or (isinstance(g, str) and not g.strip()):
            return "intergenic"
        return str(g)

    def record_fusion(self, context1: Optional[str], context2: Optional[str],
                       read_name: str, chrom1: str, pos1: int, chrom2: str, pos2: int,
                       left_score: Optional[float] = None,
                       right_score: Optional[float] = None) -> None:
        """Record a fusion observation: register supporting read, breakpoint, and per-read scores."""
        if self._is_mitochondrial_candidate(chrom1, chrom2, context1, context2):
            return
        left_gene = self.normalize_gene_label(context1)
        right_gene = self.normalize_gene_label(context2)
        # Sort genes alphabetically to avoid duplicate entries in both directions
        sorted_genes = sorted([left_gene, right_gene])
        fusion_key = f"{sorted_genes[0]}--{sorted_genes[1]}"
        self.fusion_assigned_pairs[fusion_key][read_name] = (left_gene, right_gene)
        self.fusion_candidates[fusion_key].add(read_name)
        # Store breakpoint without strand information: (chrom1, pos1, chrom2, pos2)
        bp = (chrom1, int(pos1), chrom2, int(pos2))
        self.fusion_breakpoints[fusion_key][bp] += 1
        # Store per-read scores for later averaging
        if left_score is not None and right_score is not None:
            self.fusion_read_scores[fusion_key][read_name] = (left_score, right_score)
        meta = self.fusion_metadata.setdefault(
            fusion_key,
            {"supporting_reads": set(), "consensus_bp": None,
            "left_gene": None, "right_gene": None, "support": 0}
        )
        meta["supporting_reads"].add(read_name)
        meta["support"] = len(meta["supporting_reads"])

    def normalize_gene_label(self, gene: Optional[str]) -> str:
        """Normalize any gene identifier to a canonical gene symbol string (or ``intergenic``)."""
        if not gene:
            return "intergenic"
        # Resolve ENSG / IDs → symbol if possible
        gene = self.resolve_gene_name(gene)
        # Collapse antisense / divergent
        gene = self.canonical_locus_name(gene)
        # Safety net: exclude unresolved identifiers -> mark as intergenic
        # This catches ENSG*, RP11*, and other unresolved gene identifiers
        if gene is None or gene.startswith(("ENSG", "RP11")):
            return "intergenic"
        return gene

    def build_metadata(self, min_support: int = 1) -> None:
        """Apply early non-coding filtering, then run the FusionMetadata pipeline."""
        validator = FusionValidator(self)
        validator.filter_early_non_coding_genes()

        try:
            from .fusion_metadata import FusionMetadata
            FusionMetadata(self).process_all(min_support=min_support)
        except Exception as e:
            logger.error(f"Fusion metadata processing failed: {str(e)}")
            logger.debug("Traceback:", exc_info=True)
            # Fallback: if processing fails, keep original behavior minimal (no-op)
            return

    def cluster_breakpoints(self, bp_counts: dict, window: int):
        """Cluster breakpoint counts within ``window`` and return ``((c1,p1,c2,p2), total)`` or None."""
        if not bp_counts:
            return None
        # Group by chrom pair and collect strand information
        by_pair = defaultdict(list)
        # Handle format (c1, p1, c2, p2)
        for bp_tuple, w in bp_counts.items():
            c1, p1, c2, p2 = bp_tuple
            by_pair[(c1, c2)].append((p1, p2, w))
        def best_window(items, w):
            # Sort by p1; we’ll apply a sliding window on p1, and inside it on p2
            items.sort(key=lambda x: x[0])  # sort by p1
            best_total = 0
            best_p1 = 0
            best_p2 = 0
            n = len(items)
            j = 0
            k = 0
            sum_w = 0
            sum_p1w = 0
            sum_p2w = 0
            # Sort by (p1, p2) once: O(n log n). Then use two-pointer for both windows: O(n).
            items.sort(key=lambda x: (x[0], x[1]))
            # Two-pointer sliding window on both p1 and p2 dimensions
            for i in range(n):
                p1_i, p2_i, wi = items[i]
                # Shrink p1-window from left if needed
                while j <= i and p1_i - items[j][0] > w:
                    sum_w -= items[j][2]
                    sum_p1w -= items[j][0] * items[j][2]
                    sum_p2w -= items[j][1] * items[j][2]
                    j += 1
                # Expand current item into window
                sum_w += wi
                sum_p1w += p1_i * wi
                sum_p2w += p2_i * wi
                # k must stay within p1-window
                k = max(k, j)
                # Shrink p2-window from left if needed (keeping only items within p2-window relative to i)
                while k <= i and p2_i - items[k][1] > w:
                    sum_w -= items[k][2]
                    sum_p1w -= items[k][0] * items[k][2]
                    sum_p2w -= items[k][1] * items[k][2]
                    k += 1
                # Update best window if current is better
                if sum_w > best_total:
                    best_total = sum_w
                    best_p1 = sum_p1w // sum_w if sum_w > 0 else 0
                    best_p2 = sum_p2w // sum_w if sum_w > 0 else 0
            return best_p1, best_p2, best_total
        best_pair = None
        best = (None, None, 0)
        for pair, items in by_pair.items():
            c_p1, c_p2, total = best_window(items, window)
            if total > best[2]:
                best_pair = pair
                best = (c_p1, c_p2, total)

        if best_pair is None:
            return None
        (c1, c2) = best_pair
        (p1, p2, total) = best
        return (c1, p1, c2, p2), total

    def parse_sa_entries(self, sa_tag: Optional[str]) -> list[tuple]:
        """Parse a BAM ``SA`` tag string into a list of ``(chrom, pos, strand, cigar, mapq, nm)`` tuples."""
        entries: list[tuple] = []
        if not sa_tag:
            return entries
        for sa in sa_tag.split(";"):
            if not sa:
                continue
            fields = sa.split(",")
            if len(fields) < 6:
                # fall back to partial parse
                chrom = fields[0]
                pos = int(fields[1]) if len(fields) > 1 and fields[1].isdigit() else None
                strand = fields[2] if len(fields) > 2 else "+"
                cigar = fields[3] if len(fields) > 3 else None
                mapq = int(fields[4]) if len(fields) > 4 and fields[4].isdigit() else 0
                nm = int(fields[5]) if len(fields) > 5 and fields[5].isdigit() else 0
            else:
                chrom = fields[0]
                pos = int(fields[1])
                strand = fields[2]
                cigar = fields[3]
                mapq = int(fields[4])
                nm = int(fields[5])
            entries.append((chrom, pos, strand, cigar, mapq, nm))
        return entries

    def estimate_breakpoint(self, read, sa_pos: Optional[int],
                             clip_side: Optional[str] = None,
                             sa_cigar: Optional[str] = None) -> tuple:
        """Estimate left/right fusion breakpoints from a read and its SA mate."""
        left_chr = read.reference_name
        prim_start = self.safe_reference_start(read)
        prim_end = (read.reference_end if read.reference_end is not None else None)
        if prim_end is not None:
            prim_end = int(prim_end)
        if clip_side == "right" and prim_end is not None:
            left_pos = prim_end
        elif clip_side == "left" and prim_start is not None:
            left_pos = prim_start
        else:
            # No strong clip cue -> choose the side with larger aligned block; default to end
            left_pos = prim_end if prim_end is not None else prim_start
        right_chr = read.reference_name
        right_pos = sa_pos
        if sa_cigar:
            sa_al_len = self.aligned_len_from_cigarstring(sa_cigar)
            # If primary clip was 'left', SA likely anchors the right side near its end; otherwise near its start
            if clip_side == "left" and sa_pos is not None and sa_al_len:
                right_pos = sa_pos + sa_al_len  # approximate SA end (1-based)
            else:
                right_pos = sa_pos  # SA start

        return (left_chr, int(left_pos) if left_pos is not None else None,
                right_chr, int(right_pos) if right_pos is not None else None)

    def realign_softclip(self, read, clip_side, clip_len,
                                seen_pairs, jitter_window=25,
                                min_sa_mapq=20):
        seq = read.query_sequence
        if not seq:
            return
        clipped_seq = seq[:clip_len] if clip_side == "left" else seq[-clip_len:]
        best_hit = self.realign_clipped_seq(clipped_seq)
        if not best_hit:
            return
        hit_mapq = getattr(best_hit, "mapq", 0)
        if hit_mapq < min_sa_mapq:
            return
        sa_chr = getattr(best_hit, "ctg", None)
        sa_pos = getattr(best_hit, "r_st", None)
        if sa_chr is None or sa_pos is None:
            return
        sa_pos = int(sa_pos) + 1  # convert to 1-based
        left_chr, left_pos, right_chr, right_pos = self.estimate_breakpoint(
            read, sa_pos, clip_side=clip_side, sa_cigar=None
        )
        if left_pos is None or right_pos is None:
            return
        key = (left_chr, (left_pos // jitter_window),
            sa_chr, (right_pos // jitter_window))
        if key in seen_pairs:
            return
        left_gene, left_score = self.assign_fusion_gene_cached(left_chr, left_pos)
        right_gene, right_score = self.assign_fusion_gene_cached(sa_chr, right_pos)
        if left_gene is not None and right_gene is not None and left_gene != right_gene:
            self.record_fusion(left_gene, right_gene, read.query_name, left_chr, left_pos, sa_chr, right_pos,
                             left_score=left_score, right_score=right_score)

    def detect_fusions(self,
                       min_al_len_primary: int = 50,
                       min_al_len_sa: int = 40,
                       min_sa_mapq: int = 10,
                       min_softclip_len: int = 30,
                       jitter_window: int = 50) -> None:
        """Main fusion detection entry point. Processes all reads in BAM file."""
        logger.info(f"Starting fusion detection on {self.bam_path}")
        with pysam.AlignmentFile(self.bam_path, "rb") as bam:
            for read in bam:
                passes, clip_side, clip_len = self._passes_read_filters(
                    read, min_al_len_primary, min_sa_mapq
                )
                if not passes:
                    continue
                sa_entries = []
                if read.has_tag("SA"):
                    sa_entries = self.parse_sa_entries(read.get_tag("SA"))
                if sa_entries:
                    filtered = self._filter_sa_entries(
                        sa_entries, min_al_len_sa, min_sa_mapq
                    )
                    self._process_sa_entries(
                        read, filtered, clip_side, jitter_window, min_sa_mapq
                    )
                if clip_side and clip_len >= min_softclip_len:
                    seen_pairs: set = set()
                    self.realign_softclip(
                        read, clip_side, clip_len,
                        seen_pairs,
                        jitter_window=jitter_window,
                        min_sa_mapq=min_sa_mapq
                    )

    def reconstruct_fusion_transcript(self, c1: str, p1: int, c2: str, p2: int,
        exon_padding: int = 100,
    ):
        """
        Attempt to reconstruct a fusion transcript sequence by concatenating
        exon sequences in the vicinity of two breakpoints.
        Returns:
            (sequence, left_exons, right_exons) on success,
            (None, None, None) on failure.
        """
        try:
            left_exons = list(
                self.db.region(
                    region=(c1, max(1, p1 - exon_padding), p1 + exon_padding),
                    featuretype="exon",
                )
            )
            if not left_exons:
                logger.debug("No left exons found near %s:%s", c1, p1)
                return None, None, None
            right_exons = list(
                self.db.region(
                    region=(c2, max(1, p2 - exon_padding), p2 + exon_padding),
                    featuretype="exon",
                )
            )
            if not right_exons:
                logger.debug("No right exons found near %s:%s", c2, p2)
                return None, None, None

            left_coords = sorted((int(e.start), int(e.end)) for e in left_exons)
            right_coords = sorted((int(e.start), int(e.end)) for e in right_exons)

            # Open FASTA lazily
            if not hasattr(self, "_ref_fasta"):
                self._ref_fasta = pysam.FastaFile(self.reference_fasta)
            ref = self._ref_fasta
            left_seq = "".join(ref.fetch(c1, s - 1, e) for s, e in left_coords)
            right_seq = "".join(ref.fetch(c2, s - 1, e) for s, e in right_coords)
            if not left_seq or not right_seq:
                logger.debug("Empty sequence during fusion reconstruction at %s:%s--%s:%s", c1, p1, c2, p2)
                return None, None, None
            return left_seq + right_seq, left_coords, right_coords

        except (ValueError, KeyError, AttributeError, OSError):
            logger.debug(
                "Failed to reconstruct fusion transcript at %s:%s--%s:%s",
                c1, p1, c2, p2,
                exc_info=True,
            )
            return None, None, None

    def realign_fusion_transcript(self, fusion_seq, min_match_len=30, min_mapq=10):
        #Returns list of hits or empty list if alignment fails or is weak.
        if not self.aligner or not fusion_seq or len(fusion_seq) < min_match_len:
            logger.debug("Realign skipped: missing aligner or sequence too short")
            return []
        try:
            hits = list(self.aligner.map(fusion_seq))
            if not hits:
                logger.debug("No raw alignment hits for fusion transcript")
                return []
            good_hits = []
            for h in hits:
                # Support various mappy attribute names (q_st/q_en or qstart/qend or query_start/query_end)
                q_st = getattr(h, 'q_st', None)
                q_en = getattr(h, 'q_en', None)
                if q_st is None or q_en is None:
                    q_st = getattr(h, 'qstart', None) or getattr(h, 'query_start', None)
                    q_en = getattr(h, 'qend', None) or getattr(h, 'query_end', None)
                # Fallback to reported match/aligned lengths
                match_len = None
                if q_st is not None and q_en is not None:
                    try:
                        match_len = int(q_en) - int(q_st)
                    except Exception:
                        match_len = None
                if not match_len:
                    match_len = getattr(h, 'mlen', None) or getattr(h, 'alen', None) or 0
                mapq = getattr(h, 'mapq', 0)
                logger.debug(f"Raw hit: q_st={q_st} q_en={q_en} match_len={match_len} mapq={mapq} ctg={getattr(h,'ctg',None)} r_st={getattr(h,'r_st',None)} r_en={getattr(h,'r_en',None)}")
                if match_len >= min_match_len and mapq >= min_mapq:
                    good_hits.append(h)
            logger.debug(f"Filtered {len(good_hits)} good hits out of {len(hits)} raw hits for fusion transcript")
            return good_hits
        except Exception as e:
            logger.debug(f"Failed to realign fusion transcript: {str(e)}")
            return []

    def _is_mitochondrial_candidate(self, c1: str, c2: str,
                                     left_gene: Optional[str],
                                     right_gene: Optional[str]) -> bool:
        """Return True if either side of a fusion candidate appears mitochondrial."""
        mito_chrs = {"chrM", "MT", "M", "chrMT", "mitochondrion"}
        if (c1 in mito_chrs) or (c2 in mito_chrs):
            return True
        # gene names like MT-TS1, MT-ND1 etc. treat as mitochondrial
        for g in (left_gene, right_gene):
            if not g:
                continue
            if isinstance(g, str) and g.upper().startswith("MT-"):
                return True
            if g.upper() in ("MT", "MTDNA"):
                return True
        return False

    def _apply_classification_and_filters(self, meta, flags, c1, p1, c2, p2,
                                          g1_name, r1, g2_name, r2,
                                          require_gene_names,
                                          max_intra_chr_distance):
        # classify
        if g1_name and g2_name:
            # compare resolved gene symbols for intragenic check
            g1 = meta["left_gene"]
            g2 = meta["right_gene"]
            if g1 == g2:
                flags["class"] = "intragenic"
                flags["is_valid"] = False
                flags["reasons"].append(f"Intralocus (antisense/divergent) event at {g1}")
            else:
                if c1 == c2 and max_intra_chr_distance is not None:
                    dist = abs(p2 - p1)
                    if dist <= max_intra_chr_distance:
                        flags["class"] = "cis-SAGe"
        else:
            if (r1 == "intergenic") or (r2 == "intergenic"):
                flags["class"] = "intergenic"
        # gene-name requirement
        if require_gene_names and not (g1_name and g2_name):
            flags["is_valid"] = False
            flags["reasons"].append("Missing gene name on one side")

    def _attempt_reconstruction_and_realignment(self, meta, flags, c1, p1, c2, p2):
        # Try to reconstruct fusion transcript and realign; update meta/flags accordingly.
        if not (self.reference_fasta and self.aligner):
            meta.setdefault("notes", []).append("Realignment skipped: reference or aligner not available")
            return
        fusion_seq, left_exons, right_exons = self.reconstruct_fusion_transcript(c1, p1, c2, p2, exon_padding=100)
        if not fusion_seq:
            flags["is_valid"] = False
            flags["reasons"].append("Failed to reconstruct fusion transcript")
            meta["reconstruction_ok"] = False
            return
        # mark that reconstruction produced a transcript
        meta["reconstruction_ok"] = True
        hits = self.realign_fusion_transcript(fusion_seq)
        if not hits:
            flags["is_valid"] = False
            flags["reasons"].append("Fusion transcript does not realign cleanly to genome")
        else:
            # Store realignment info for debugging
            meta["realignment_hits"] = len(hits)
            try:
                meta["best_hit_mapq"] = max([h.mapq for h in hits]) if hits else 0
            except Exception:
                meta["best_hit_mapq"] = 0

    def canonical_locus_name(self, gene_name: Optional[str]) -> Optional[str]:
        """Collapse antisense/divergent transcript names to the canonical locus name."""
        if not gene_name:
            return gene_name
        # Strip common antisense / divergent suffixes
        collapsed = ANTISENSE_SUFFIX_RE.sub("", gene_name)
        return collapsed

    def _build_symbol_biotype_index(self):
        # One-time index construction: map HGNC symbols to biotypes using all genes in db
        # This is called lazily when first needed
        if self._symbol_biotype_cache:
            return  # already built
        logger.debug("Building symbol → biotype index from genedb...")
        try:
            # Use interval tree to iterate genes if available (more efficient)
            if self.interval_index is not None:
                genes = self.interval_index.get_all_genes()
            else:
                genes = list(self.db.features_of_type('gene'))
            for gene in genes:
                attrs = getattr(gene, 'attributes', {}) or {}
                gene_name = attrs.get('gene_name', [None])[0]
                if not gene_name:
                    continue
                # Extract biotype
                biotype = (
                    attrs.get('gene_type', [None])[0]
                    or attrs.get('gene_biotype', [None])[0]
                    or attrs.get('transcript_biotype', [None])[0]
                )
                if biotype:
                    self._symbol_biotype_cache[gene_name] = biotype
        except Exception as e:
            logger.warning(f"Failed to build symbol→biotype index: {e}")

    def _get_gene_biotype_by_symbol_or_coords(self, gene_symbol, chrom=None, pos=None):
        # Resolve HGNC symbol to biotype using two strategies:
        # 1) Direct lookup in symbol cache (built on-demand)
        # 2) If cache miss, query genes at (chrom, pos) and match by gene_name attribute
        if not gene_symbol:
            return None
        # Strategy 1: Check symbol cache first
        if not self._symbol_biotype_cache:
            self._build_symbol_biotype_index()
        if gene_symbol in self._symbol_biotype_cache:
            biotype = self._symbol_biotype_cache[gene_symbol]
            logger.debug(f"Found biotype for {gene_symbol} in cache: {biotype}")
            return biotype
        # Strategy 2: Use consensus coordinates if available
        if chrom is not None and pos is not None:
            try:
                if self.interval_index is not None:
                    genes = self.interval_index.get_genes_at(chrom, pos, window=2000)
                else:
                    genes = list(self.db.region(
                        region=(chrom, max(1, pos - 2000), pos + 2000),
                        featuretype='gene'
                    ))
                for g in genes:
                    attrs = getattr(g, 'attributes', {}) or {}
                    g_name = attrs.get('gene_name', [None])[0]
                    # Match by symbol
                    if g_name and g_name.upper() == gene_symbol.upper():
                        biotype = (
                            attrs.get('gene_type', [None])[0]
                            or attrs.get('gene_biotype', [None])[0]
                            or attrs.get('transcript_biotype', [None])[0]
                        )
                        if biotype:
                            logger.debug(f"Found biotype for {gene_symbol} by coordinates {chrom}:{pos}: {biotype}")
                            # Cache for future lookups
                            self._symbol_biotype_cache[gene_symbol] = biotype
                            return biotype
            except Exception as e:
                logger.debug(f"Failed to query genes at {chrom}:{pos} for symbol {gene_symbol}: {e}")
        logger.debug(f"Could not find biotype for {gene_symbol} (chrom={chrom}, pos={pos}, cache_size={len(self._symbol_biotype_cache)})")
        return None

    def get_gene_biotype(self, gene_name: Optional[str],
                          chrom: Optional[str] = None,
                          pos: Optional[int] = None) -> Optional[str]:
        """Retrieve the biotype for a gene name, optionally disambiguating by ``(chrom, pos)``."""
        if not gene_name:
            return None
        gene_name = self.normalize_gene_label(gene_name)
        return self._get_gene_biotype_by_symbol_or_coords(gene_name, chrom=chrom, pos=pos)

    def validate_candidates(self, min_support: int = 2, window: int = 25,
                             require_gene_names: bool = True,
                             require_mapq: int = 10,
                             allow_cis_sage: bool = True,
                             require_exon_boundary: bool = True,
                             max_intra_chr_distance: Optional[int] = None) -> None:
        """Run the validation pipeline: gating, classification, reconstruction, and scoring."""
        self.build_metadata(min_support=min_support)
        self.fusion_assigned_pairs.clear()
        self.fusion_read_scores.clear()
        # Delegate validation and filtering to FusionValidator
        validator = FusionValidator(self)
        # validator.filter_unresolved_genes()
        validator.filter_non_coding_genes()
        validator.filter_multicopy_artifact_pairs()
        validator.apply_frequency_filters()
        for fusion_key, meta in list(self.fusion_metadata.items()):
            # Check if already marked invalid by frequency filter
            if meta.get("is_valid", True) is False:
                # Keep existing reasons, skip to next
                continue
            support = meta.get("support", 0)
            flags = {"is_valid": True, "reasons": [], "class": "canonical"}
            if support < min_support:
                flags["is_valid"] = False
                flags["reasons"].append(f"Low support ({support} < {min_support})")
                meta.update(flags)
                continue
            # Consensus breakpoint already computed in build_metadata - no need to re-cluster
            consensus_bp = meta.get("consensus_bp")
            if not consensus_bp:
                flags["is_valid"] = False
                flags["reasons"].append("No stable breakpoint cluster")
                meta.update(flags)
                continue

            # perform classification and basic filtering
            c1, p1, c2, p2 = consensus_bp
            g1_name, r1 = self._context_query(c1, p1)
            g2_name, r2 = self._context_query(c2, p2)
            self._apply_classification_and_filters(meta, flags, c1, p1, c2, p2,
                                                   g1_name, r1, g2_name, r2,
                                                   require_gene_names, max_intra_chr_distance)
            # skip reconstruction for mitochondrial candidates and invalid ones
            if flags["is_valid"] and not self._is_mitochondrial_candidate(c1, c2, meta.get("left_gene"), meta.get("right_gene")):
                self._attempt_reconstruction_and_realignment(meta, flags, c1, p1, c2, p2)
            else:
                if flags["is_valid"] and self._is_mitochondrial_candidate(c1, c2, meta.get("left_gene"), meta.get("right_gene")):
                    flags["is_valid"] = False
                    flags["reasons"].append("Mitochondrial fusion candidate filtered out")

            # cis-SAGe policy enforcement
            if flags["class"] == "cis-SAGe" and not allow_cis_sage:
                flags["is_valid"] = False
                flags["reasons"].append("cis-SAGe disallowed by policy")

            # Compute confidence
            conf = validator.confidence(meta, flags)
            meta["confidence"] = conf
            # convert very low-confidence to invalid
            if conf < 0.30 and meta.get("support", 0) < 2 and flags["is_valid"]:
                flags["is_valid"] = False
                flags["reasons"].append("Low confidence")
            meta.update(flags)

    def report(self, output_path: str = "fusion_candidates.tsv",
                min_support: int = 2,
                include_classes: tuple = ("canonical", "cis-SAGe"),
                min_confidence: float = 0.3,
                only_valid: bool = False) -> None:
            """Validate candidates and write the fusion report TSV."""
            self.validate_candidates(min_support=min_support)
            # Merge only fully identical fusions (exact duplicates)
            validator = FusionValidator(self)
            validator._merge_fully_identical()
            with open(output_path, "w") as f:
                f.write("LeftGene\tLeftBiotype\tLeftScore\tLeftChromosome\tLeftBreakpoint\t"
                        "RightGene\tRightBiotype\tRightScore\tRightChromosome\tRightBreakpoint\t"
                        "SupportingReads\tFusionName\tClass\tValid\tConfidence\tReasons\n")
                for fusion_key, meta in sorted(self.fusion_metadata.items(), key=lambda x: -x[1].get("support", 0)):
                    if meta.get("confidence", 0) < min_confidence:
                        continue
                    if only_valid and not meta.get("is_valid", True):
                        continue
                    if meta.get("class") not in include_classes:
                        continue
                    if not meta.get("consensus_bp"):
                        continue
                    # Handle both old format (c1, p1, c2, p2) and new format (c1, p1, s1, c2, p2, s2)
                    consensus_bp = meta["consensus_bp"]
                    left_chr, left_pos, right_chr, right_pos = consensus_bp
                    left_gene  = meta["left_gene"]
                    right_gene = meta["right_gene"]
                    # Skip if collapse makes them equal
                    if left_gene == right_gene:
                        continue
                    # Skip if either gene has antisense/regulatory suffix
                    if self.has_antisense_suffix(left_gene) or self.has_antisense_suffix(right_gene):
                        continue
                    fusion_name = f"{left_gene}-{right_gene}"
                    reasons = ";".join(meta.get("reasons", []))
                    left_score = meta.get("left_score", 0.0)
                    right_score = meta.get("right_score", 0.0)
                    left_biotype = meta.get("left_biotype", "unknown")
                    right_biotype = meta.get("right_biotype", "unknown")
                    f.write(f"{left_gene}\t{left_biotype}\t{left_score:.2f}\t{left_chr}\t{left_pos}\t"
                            f"{right_gene}\t{right_biotype}\t{right_score:.2f}\t{right_chr}\t{right_pos}\t"
                            f"{meta.get('support', 0)}\t{fusion_name}\t{meta.get('class')}\t"
                            f"{meta.get('is_valid')}\t{meta.get('confidence')}\t{reasons}\n")

