import re
import pysam
import gffutils
from collections import defaultdict
import mappy as mp
import logging
from functools import lru_cache
from intervaltree import IntervalTree
from .genomic_interval_index import GenomicIntervalIndex
from .fusion_validator import FusionValidator

logger = logging.getLogger('IsoQuant')
Antisense_suffix = re.compile(r"-(AS\d+|DT|DIVERGENT|NAT)$", re.IGNORECASE)

class FusionDetector:
    def __init__(self, bam_path, gene_db_path, reference_fasta):
        self.bam_path = bam_path
        self.genedb_path = gene_db_path
        self.db = gffutils.FeatureDB(gene_db_path, keep_order=True)
        self.fusion_candidates = defaultdict(set)
        # store breakpoint counts per fusion key: {(chr1,pos1,chr2,pos2): count}
        self.fusion_breakpoints = defaultdict(lambda: defaultdict(int))
        # optional reference and aligner for soft-clip realignment
        self.reference_fasta = reference_fasta
        try:
            self.aligner = mp.Aligner(reference_fasta) if reference_fasta else None
        except Exception:
            self.aligner = None
        self.fusion_metadata = {}
        # store per-read assigned raw gene pairs for each fusion key: {fusion_key: {read_name: (left,right)}}
        self.fusion_assigned_pairs = defaultdict(dict)
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

    def _build_exon_cache(self):
        # Pre-build exon cache for all genes to avoid repeated DB queries during fusion gene assignment.
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
                    if exons:
                        self.exon_cache[gene_id] = [(int(ex.start), int(ex.end)) for ex in exons]
                except Exception:
                    pass
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

    @lru_cache(maxsize=200000)
    def _context_query(self, chrom, pos):
        # Query genomic context at position (chrom, pos).
        # Uses interval tree if available for O(log n) performance, otherwise falls back to gffutils.
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
                    return gene_name, "exonic"
                else:
                    return gene_name, "intronic"
            else:
                # Not inside any gene; check if overlapping any exon 
                return None, "exonic" if exons else "intergenic"
        except Exception:
            return None, "unknown"

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

    def _fallback_gene_lookup(self, name_or_id):
        # Last resort: search for a gene feature with matching ID.
        try:
            genes = list(self.db.features_of_type('gene', id=name_or_id))
            if genes:
                return genes[0].attributes.get('gene_name', [genes[0].id])[0]
        except Exception:
            pass
        return name_or_id
    
    def sanitize_raw_gene(self, gene, chrom, pos):
        # Keep original if protein coding
        biotype = self._get_gene_biotype_by_symbol_or_coords(gene, chrom, pos)
        if biotype == "protein_coding":
            return gene

        # Try collapsing antisense/pseudogene naming like HMGA1P3 → HMGA1
        collapsed = self.canonical_locus_name(gene)

        # If collapsed is different (i.e. suffix removed), try using that
        if collapsed != gene:
            return collapsed

        # Rescue: find nearby coding gene
        rescue = self._find_nearby_protein_coding(chrom, pos, window=2000)
        if rescue:
            return rescue

        # If still nothing, leave gene unchanged (do not turn into intergenic!)
        return gene

    def resolve_gene_name(self, name_or_id):
        # Resolve Ensembl or other IDs (e.g. ENS*, RP11*) to gene symbols when possible.
        if not name_or_id:
            return name_or_id
        # Quick cache hit
        if name_or_id in self._resolved_name_cache:
            return self._resolved_name_cache[name_or_id]
        # Try direct lookup first
        resolved = name_or_id
        try:
            feat, resolved = self._lookup_feature_by_id(name_or_id)
            if feat is None:
                # If direct lookup failed, try fallback
                resolved = self._fallback_gene_lookup(name_or_id)
        except Exception:
            resolved = name_or_id
        # Cache and return
        self._resolved_name_cache[name_or_id] = resolved
        return resolved

    def has_antisense_suffix(self, gene_name):
        # Check if gene name ends with antisense/regulatory suffixes (AS1, DT, etc.)
        if not gene_name or not isinstance(gene_name, str):
            return False
        return Antisense_suffix.search(gene_name) is not None

    def strip_antisense_suffix(self, gene_name):
        # Remove antisense/regulatory suffix and return stripped name
        if not gene_name or not isinstance(gene_name, str):
            return gene_name
        return Antisense_suffix.sub('', gene_name)

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

    def _score_hard_evidence(self, pos, gstart, gend, exonic_hit, boundary_min_dist, exon_min_dist):
        # Score based on proximity and exonic hits (hard evidence).
        # Returns score contribution from location-based metrics.
        score = 0.0
        if gstart - 3000 <= pos <= gend + 3000:
            score += 50  # near edge bonus
        if exonic_hit:
            score += 200.0
        if boundary_min_dist is not None:
            score += max(0.0, 80.0 - min(boundary_min_dist, 200))
        if exon_min_dist is not None:
            score += max(0.0, 40.0 - min(exon_min_dist, 500))
        # Body proximity (weak signal)
        if pos < gstart or pos > gend:
            body_dist = min(abs(pos - gstart), abs(pos - gend))
            score += max(0.0, 20.0 - min(body_dist, 2000) / 100.0)
        return score

    def _score_biotype(self, gtype):
        # Score based on gene biotype.
        if gtype == "protein_coding":
            return 60.0
        # Penalize pseudogenes heavily
        # if gtype and "pseudogene" in gtype.lower():
            #return -90.0
        # Mild penalty for non-coding/ambiguous biotypes
        if gtype in ("antisense", "lncRNA", "lincRNA", "processed_transcript",
                     "sense_intronic", "sense_overlapping", "transcribed_unprocessed_pseudogene"):
            return -15.0
        return 0.0

    def _score_exon_span_bonus(self, g):
        # Score bonus for longer exon span (more likely real coding gene).
        try:
            exons = self._get_cached_exons(g)
            exon_len_sum = sum(ee - es + 1 for es, ee in exons)
            return min(exon_len_sum / 1e4, 10.0)  # cap at +10
        except Exception:
            return 0.0

    def _compute_gene_score(self, target, pos, chrom=None, exon_min_dist=None, 
                            boundary_min_dist=None, exonic_hit=None, body_dist=None, 
                            gstart=None, gend=None):
        # Compute composite gene score incorporating multiple factors.
        # Determine mode and resolve feature + metrics
        if isinstance(target, str):
            # Lookup mode: find gene by name at position
            gene_name = target
            if chrom is None:
                raise ValueError("chrom required when target is a gene name")
            g, gstart, gend = self._lookup_gene_by_name(gene_name, chrom, pos)
            if g is None:
                return 0.0
            # Compute metrics for this gene
            exonic_hit, exon_min_dist, boundary_min_dist, body_dist, _ = self._compute_exonic_metrics(g, pos)
        else:
            # Pre-computed mode: target is already a feature object
            g = target
            if gstart is None or gend is None:
                raise ValueError("gstart and gend required for pre-computed mode")
        # Extract biotype
        attrs = getattr(g, "attributes", {}) or {}
        gtype = (attrs.get("gene_type", [None])[0]
                or attrs.get("gene_biotype", [None])[0]
                or attrs.get("transcript_biotype", [None])[0])
        # Aggregate scores from all components
        score = 0.0
        score += self._score_hard_evidence(pos, gstart, gend, exonic_hit, boundary_min_dist, exon_min_dist)
        score += self._score_biotype(gtype)
        score += self._score_exon_span_bonus(g)
        return score

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
        # Return list of genes where pos falls inside an exon.
        exonic_genes = []
        for g in genes:
            attrs = getattr(g, "attributes", {}) or {}
            gname = attrs.get("gene_name", [getattr(g, "id", None)])[0]
            gtype = (attrs.get("gene_type", [None])[0]
                    or attrs.get("gene_biotype", [None])[0]
                    or attrs.get("transcript_biotype", [None])[0])
            gstart = int(getattr(g, "start", 0))
            gend = int(getattr(g, "end", 0))
            try:
                exons = self._get_cached_exons(g)
                exonic_hit = any(start <= pos <= end for start, end in exons)
                if exonic_hit:
                    exonic_genes.append((g, gname, gtype, gstart, gend, exons))
            except Exception:
                pass
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
            bool(exonic_hit),                        # 3) prefer exonic
            -bnd,                                    # 4) nearer exon boundary
            -exd,                                    # 5) nearer exon span
            -bod,                                    # 6) nearer gene body
            -min(glen, 500000),                      # 7) gentle preference for longer locus
            gname or ""                              # 8) deterministic
        )

    def _score_exonic_genes(self, exonic_genes, pos):
        # Score multiple exonic genes and return the best one, or single exonic gene.
        if len(exonic_genes) == 1:
            return exonic_genes[0][1]  # Fast return for single exonic gene
        best_gene = None
        best_key = None
        for g, gname, gtype, gstart, gend, exons in exonic_genes:
            boundary_distances = [min(abs(pos - start), abs(pos - end)) for start, end in exons]
            boundary_min_dist = min(boundary_distances) if boundary_distances else 0
            score = self._compute_gene_score(
                g, pos, exon_min_dist=0, boundary_min_dist=boundary_min_dist,
                exonic_hit=True, body_dist=0, gstart=gstart, gend=gend
            )
            glen = max(1, gend - gstart + 1)
            key = self._build_gene_key(gtype, score, True, boundary_min_dist, 0, 0, glen, gname)
            if best_key is None or key > best_key:
                best_key = key
                best_gene = gname or getattr(g, "id", None)
        return best_gene

    def _score_intronic_genes(self, genes, pos):
        # Score intronic/intergenic genes and return the best one.
        best_gene = None
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
                g, pos, exon_min_dist=exon_min_dist, boundary_min_dist=boundary_min_dist,
                exonic_hit=exonic_hit, body_dist=body_dist, gstart=gstart, gend=gend
            )
            glen = max(1, gend - gstart + 1)
            key = self._build_gene_key(gtype, score, exonic_hit, boundary_min_dist, exon_min_dist, body_dist, glen, gname)
            if best_key is None or key > best_key:
                best_key = key
                best_gene = gname or getattr(g, "id", None)
        return best_gene

    def assign_fusion_gene(self, chrom, pos, window=1000):
        # Assign the best gene at a genomic location using exonic priority and scoring.
        genes = self._get_genes_at(chrom, pos, window)
        if not genes:
            return None
        # FAST PATH: Try exonic genes first (position inside an exon)
        exonic_genes = self._collect_exonic_genes(genes, pos)
        if exonic_genes:
            return self._score_exonic_genes(exonic_genes, pos)
        # SLOW PATH: No exonic hits; use intronic/intergenic scoring
        return self._score_intronic_genes(genes, pos)

    @lru_cache(maxsize=500000)
    def assign_fusion_gene_cached(self, chrom, pos):
        return self.assign_fusion_gene(chrom, pos)

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
            # Store raw assigned genes (no normalization yet - defer to build_metadata)
            left_gene = self.assign_fusion_gene_cached(left_chr, left_pos)
            right_gene = self.assign_fusion_gene_cached(sa_chr, right_pos)
            if left_gene is None or right_gene is None:
                continue
            if left_gene != right_gene:
                # Check soft-clip orientation consistency
                # if self._is_softclip_orientation_valid(read, clip_side, sa_chr, sa_pos):
                self.record_fusion(
                    left_gene, right_gene, read.query_name,
                    left_chr, left_pos, sa_chr, right_pos
                )
                sa_used = True
        return sa_used

    def detect_softclip(self, read):
        # Detect large soft-clips (>=50 bp) on either end of the read.
        # Returns (clip_side, clip_len) or (None, 0) if no significant clip.
        cigartuples = getattr(read, "cigartuples", None)
        clip_side, clip_len = None, 0
        if cigartuples:
            # Leading soft-clip
            if len(cigartuples) > 0 and cigartuples[0][0] == 4 and cigartuples[0][1] >= 50:
                clip_side, clip_len = "left", cigartuples[0][1]
            # Trailing soft-clip
            elif len(cigartuples) > 0 and cigartuples[-1][0] == 4 and cigartuples[-1][1] >= 50:
                clip_side, clip_len = "right", cigartuples[-1][1]
        else:
            # Fallback: parse cigarstring
            cigarstr = getattr(read, "cigarstring", "") or ""
            m1 = re.match(r'^(\d+)S', cigarstr)
            m2 = re.search(r'(\d+)S$', cigarstr)
            if m1 and int(m1.group(1)) >= 50:
                clip_side, clip_len = "left", int(m1.group(1))
            elif m2 and int(m2.group(1)) >= 50:
                clip_side, clip_len = "right", int(m2.group(1))
        return clip_side, clip_len



    def aligned_len_from_cigartuples(self, cigartuples):
        # count read-aligned operations: M (=0), = (=7), X (=8)
        if not cigartuples:
            return 0
        return sum(l for op, l in cigartuples if op in (0, 7, 8))

    def aligned_len_from_cigarstring(self, cigar):
        # Use cached parser to avoid re-running regex for identical CIGAR strings
        return self._cached_aligned_len_from_cigarstring(cigar)

    @staticmethod
    @lru_cache(maxsize=50000)
    def _cached_aligned_len_from_cigarstring(cigar):
        if not cigar:
            return 0
        total = 0
        for m in re.finditer(r'(\d+)([MIDNSHP=XB])', cigar):
            length = int(m.group(1)); op = m.group(2)
            if op in ('M', '=', 'X'):
                total += length
        return total

    @lru_cache(maxsize=10000)
    def _cached_aligner_map(self, seq):
        # Cache aligner.map() results for identical sequences to avoid expensive re-alignment
        try:
            hits = tuple(self.aligner.map(seq))
            return hits
        except Exception:
            return ()

    def _is_softclip_orientation_valid(self, read, clip_side, sa_chr, sa_pos):
        # Returns True if soft-clip direction is consistent with an outward-facing
        # Same chromosome required for directional check
        if read.reference_name != sa_chr:
            return True
        prim_start = self.safe_reference_start(read)
        prim_end   = read.reference_end
        if prim_start is None or prim_end is None:
            return True  # can't evaluate; don't block the event
        read_is_reverse = read.is_reverse
        # Forward-strand logic
        if not read_is_reverse:
            if clip_side == "right":
                # clipping at 3' end → expect SA alignment downstream
                return sa_pos > prim_end
            elif clip_side == "left":
                # clipping at 5' end → expect SA alignment upstream
                return sa_pos < prim_start
            else:
                return True
        # Reverse-strand logic
        else:
            if clip_side == "left":
                # left clip on reverse strand is 3' end → downstream
                return sa_pos > prim_end
            elif clip_side == "right":
                # right clip on reverse strand is 5' end → upstream
                return sa_pos < prim_start
            else:
                return True

    def realign_clipped_seq(self, seq):
        # Realign a clipped sequence using the aligner.
        if not seq or not self.aligner:
            return None
        try:
            hits = self._cached_aligner_map(seq)
            for hit in hits:
                hit_len = self._get_hit_length(hit)
                if hit_len and hit_len >= 50:
                    return hit  # first good hit
            return None
        except Exception:
            return None

    def _get_hit_length(self, hit):
        # Calculate the aligned length from a hit object.
        r_en = getattr(hit, "r_en", None)
        if r_en is not None:
            r_st = getattr(hit, "r_st", None)
            return r_en - r_st
        else:
            mlen = getattr(hit, "mlen", None)
            if mlen is not None:
                return mlen
            else:
                alen = getattr(hit, "alen", None)
                if alen is not None:
                    return alen
                else:
                    return 0

    def _safe_gene_token(self, g):
        # Convert None/empty to 'intergenic'
        if not g or (isinstance(g, str) and not g.strip()):
            return "intergenic"
        return str(g)

    def record_fusion(self, context1, context2, read_name, chrom1, pos1, chrom2, pos2):
        if self._is_mitochondrial_candidate(chrom1, chrom2, context1, context2):
            return
        left_clean  = self.sanitize_raw_gene(context1, chrom1, pos1)
        right_clean = self.sanitize_raw_gene(context2, chrom2, pos2)

        # 2. Now normalize AFTER sanitization
        left_symbol  = self.normalize_gene_label(left_clean)
        right_symbol = self.normalize_gene_label(right_clean)

        # 3. Use the same values for raw assignment
        left_raw  = left_clean
        right_raw = right_clean

        # 4. Build fusion key from the sanitized & normalized names
        fusion_key = f"{left_symbol}--{right_symbol}"

        # 5. Now store the raw genes
        self.fusion_assigned_pairs[fusion_key][read_name] = (left_raw, right_raw)

        self.fusion_candidates[fusion_key].add(read_name)
        # Store breakpoint without strand information: (chrom1, pos1, chrom2, pos2)
        bp = (chrom1, int(pos1), chrom2, int(pos2))
        self.fusion_breakpoints[fusion_key][bp] += 1

        meta = self.fusion_metadata.setdefault(
            fusion_key,
            {"supporting_reads": set(), "consensus_bp": None,
            "left_gene": None, "right_gene": None, "support": 0}
        )
        meta["supporting_reads"].add(read_name)
        meta["support"] = len(meta["supporting_reads"])

    def normalize_gene_label(self, gene):
        if not gene or gene == "intergenic":
            return "intergenic"

            # Collapse antisense suffixes
        gene = self.canonical_locus_name(gene)

        # If it's an ENSG ID, resolve it to symbol.
        if gene.startswith("ENSG"):
            resolved = self.resolve_gene_name(gene)
            return resolved or gene

        # Otherwise, assume it's already an HGNC gene symbol.
        # DO NOT resolve again. DO NOT turn into intergenic.
        return gene

    def build_metadata(self, min_support=1):
        # Apply early filtering to drop non-protein-coding fusions BEFORE breakpoint clustering
        validator = FusionValidator(self)
        validator.filter_raw_non_coding_genes()

        try:
            from .fusion_metadata import FusionMetadata
            FusionMetadata(self).process_all(min_support=min_support)
        except Exception as e:
            logger.error(f"Fusion metadata processing failed: {str(e)}")
            logger.debug("Traceback:", exc_info=True)
            # Fallback: if processing fails, keep original behavior minimal (no-op)
            return
        # After metadata is finalized, consolidate duplicate fusions (swapped partners)
        self.consolidate_duplicate_fusions()

    def cluster_breakpoints(self, bp_counts, window):
        if not bp_counts:
            return None
        # Group by chrom pair and collect strand information
        from collections import defaultdict
        by_pair = defaultdict(list)
        # Handle format (c1, p1, c2, p2)
        for bp_tuple, w in bp_counts.items():
            c1, p1, c2, p2 = bp_tuple
            by_pair[(c1, c2)].append((p1, p2, w))
        best_pair = None
        best = (None, None, 0)
        for pair, items in by_pair.items():
            c_p1, c_p2, total = self._best_window_for_breakpoints(items, window)
            if total > best[2]:
                best_pair = pair
                best = (c_p1, c_p2, total)

        if best_pair is None:
            return None
        (c1, c2) = best_pair
        (p1, p2, total) = best
        return (c1, p1, c2, p2), total

    def _best_window_for_breakpoints(self, items, w):
        # Find the best weighted window of breakpoints using two-pointer sliding window on both dimensions.
        # Sort by p1; we'll apply a sliding window on p1, and inside it on p2
        items.sort(key=lambda x: x[0])  # sort by p1
        best_total = 0
        best_p1 = 0
        best_p2 = 0
        n = len(items)
        j = 0
        sum_w = 0
        # Maintain a multiset (sorted list) on p2 within current p1-window using two-pointer
        # For speed, we keep an array slice [j..i] and then slide a second pointer k on p2.
        for i in range(n):
            p1_i, p2_i, wi = items[i]
            # Expand p1-window
            sum_w += wi
            # Shrink from left until p1-range <= w
            while p1_i - items[j][0] > w:
                sum_w -= items[j][2]
                j += 1
            # Now consider only slice items[j:i+1]; build p2-sorted view (small slice)
            slice_view = items[j:i+1]
            slice_view.sort(key=lambda x: x[1])  # sort by p2
            # Secondary sliding window on p2
            k = 0
            cur_sum = 0
            for t in range(len(slice_view)):
                cur_sum += slice_view[t][2]
                while slice_view[t][1] - slice_view[k][1] > w:
                    cur_sum -= slice_view[k][2]
                    k += 1
                if cur_sum > best_total:
                    # Weighted mean p1/p2 for current p2-window [k..t]
                    s_w = 0
                    s_p1 = 0
                    s_p2 = 0
                    for x in slice_view[k:t+1]:
                        s_w += x[2]
                        s_p1 += x[0] * x[2]
                        s_p2 += x[1] * x[2]
                    best_total = s_w
                    best_p1 = s_p1 // s_w
                    best_p2 = s_p2 // s_w
        return best_p1, best_p2, best_total

    def _compute_canonical_fusion_key(self, left_gene, right_gene):
        # Compute a canonical fusion key from final gene partners.
        if not left_gene or not right_gene:
            return None
        genes = sorted([left_gene, right_gene])
        return f"{genes[0]}--{genes[1]}"

    def _find_duplicate_fusion_key(self, canonical_key):
        # Find if a canonical fusion key already exists in fusion_metadata.
        for existing_key in self.fusion_metadata.keys():
            existing_canonical = self._compute_canonical_fusion_key(
                self.fusion_metadata[existing_key].get("left_gene"),
                self.fusion_metadata[existing_key].get("right_gene")
            )
            if existing_canonical == canonical_key:
                return existing_key
        return None

    def _rename_fusion_key(self, old_key, new_key):
        # Rename a fusion key across all internal structures (single key, no merging).
        if old_key == new_key:
            return  # No-op if keys are identical
        if old_key not in self.fusion_candidates and old_key not in self.fusion_metadata:
            return  # Key doesn't exist, nothing to do
        logger.debug(f"Renaming fusion key: '{old_key}' → '{new_key}'")
        # Rename in fusion_candidates
        if old_key in self.fusion_candidates:
            self.fusion_candidates[new_key] = self.fusion_candidates.pop(old_key)
        # Rename in fusion_breakpoints
        if old_key in self.fusion_breakpoints:
            self.fusion_breakpoints[new_key] = self.fusion_breakpoints.pop(old_key)
        # Rename in fusion_metadata
        if old_key in self.fusion_metadata:
            self.fusion_metadata[new_key] = self.fusion_metadata.pop(old_key)
        # Rename in fusion_assigned_pairs
        if old_key in self.fusion_assigned_pairs:
            self.fusion_assigned_pairs[new_key] = self.fusion_assigned_pairs.pop(old_key)

    def _merge_fusion_structures(self, keep_key, discard_key):
        # Merge all internal structures from discard_key into keep_key atomically.
        logger.info(f"Consolidating duplicate fusion: '{keep_key}' ← '{discard_key}'")
        # Merge fusion_candidates: union of read sets
        if discard_key in self.fusion_candidates and keep_key in self.fusion_candidates:
            self.fusion_candidates[keep_key].update(self.fusion_candidates[discard_key])
        elif discard_key in self.fusion_candidates:
            self.fusion_candidates[keep_key] = self.fusion_candidates[discard_key].copy()
        # Merge fusion_breakpoints: sum of breakpoint counts
        if discard_key in self.fusion_breakpoints and keep_key in self.fusion_breakpoints:
            for bp, count in self.fusion_breakpoints[discard_key].items():
                self.fusion_breakpoints[keep_key][bp] = self.fusion_breakpoints[keep_key].get(bp, 0) + count
        elif discard_key in self.fusion_breakpoints:
            self.fusion_breakpoints[keep_key] = self.fusion_breakpoints[discard_key].copy()
        # Merge fusion_metadata: combine supporting reads and update support count
        if discard_key in self.fusion_metadata and keep_key in self.fusion_metadata:
            discard_meta = self.fusion_metadata[discard_key]
            keep_meta = self.fusion_metadata[keep_key]
            if "supporting_reads" in discard_meta:
                keep_meta["supporting_reads"].update(discard_meta["supporting_reads"])
            keep_meta["support"] = len(keep_meta.get("supporting_reads", set()))
        elif discard_key in self.fusion_metadata:
            self.fusion_metadata[keep_key] = self.fusion_metadata[discard_key]
        # Merge fusion_assigned_pairs: combine read-to-gene-pair mappings
        if discard_key in self.fusion_assigned_pairs and keep_key in self.fusion_assigned_pairs:
            self.fusion_assigned_pairs[keep_key].update(self.fusion_assigned_pairs[discard_key])
        elif discard_key in self.fusion_assigned_pairs:
            self.fusion_assigned_pairs[keep_key] = self.fusion_assigned_pairs[discard_key].copy()
        # Remove discard_key from all structures
        self.fusion_candidates.pop(discard_key, None)
        self.fusion_breakpoints.pop(discard_key, None)
        self.fusion_metadata.pop(discard_key, None)
        self.fusion_assigned_pairs.pop(discard_key, None)

    def consolidate_duplicate_fusions(self):
        # handles cases where A--B and B--A are kept as separate fusions.
        # Map canonical key to original fusion key found in metadata
        canonical_to_original = {}
        fusions_to_merge = {}  # canonical_key → list of original keys with this canonical
        # First pass: identify all canonical keys and find duplicates
        for fusion_key in list(self.fusion_metadata.keys()):
            meta = self.fusion_metadata[fusion_key]
            left_gene = meta.get("left_gene")
            right_gene = meta.get("right_gene")
            if not left_gene or not right_gene:
                # Skip fusions without both genes
                continue
            canonical_key = self._compute_canonical_fusion_key(left_gene, right_gene)
            if canonical_key is None:
                continue
            if canonical_key not in canonical_to_original:
                canonical_to_original[canonical_key] = fusion_key
                fusions_to_merge[canonical_key] = [fusion_key]
            else:
                fusions_to_merge[canonical_key].append(fusion_key)
        # Second pass: merge all duplicates with same canonical key
        merged_count = 0
        for canonical_key, fusion_keys in fusions_to_merge.items():
            if len(fusion_keys) <= 1:
                continue  # No duplicates
            # Keep first, merge all others into it
            keep_key = fusion_keys[0]
            for discard_key in fusion_keys[1:]:
                self._merge_fusion_structures(keep_key, discard_key)
                merged_count += 1
        if merged_count > 0:
            logger.info(f"Consolidated {merged_count} duplicate fusion(s) (swapped partners)")

    def parse_sa_entries(self, sa_tag):
        entries = []
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

    def estimate_breakpoint(self, read, sa_pos, clip_side=None, sa_cigar=None):
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
        left_gene = self.assign_fusion_gene_cached(left_chr, left_pos)
        right_gene = self.assign_fusion_gene_cached(sa_chr, right_pos)
        # Store raw genes without normalization - defer to build_metadata
        if left_gene is not None and right_gene is not None and left_gene != right_gene:
            self.record_fusion(left_gene, right_gene, read.query_name, left_chr, left_pos, sa_chr, right_pos)

    def detect_fusions(self,
                    min_al_len_primary=50,
                    min_al_len_sa=40,
                    min_sa_mapq=10,
                    min_softclip_len=30,
                    jitter_window=50):
        # Main fusion detection entry point. Processes all reads in BAM file.
        logger.info(f"Starting fusion detection on {self.bam_path}")
        bam = pysam.AlignmentFile(self.bam_path, "rb")
        processed_reads = 0
        for read in bam:
            processed_reads += 1
            # Filter reads by quality and alignment length
            passes, clip_side, clip_len = self._passes_read_filters(
                read, min_al_len_primary, min_sa_mapq
            )
            if not passes:
                continue
            # Extract and filter SA entries
            sa_entries = []
            if read.has_tag("SA"):
                sa_entries = self.parse_sa_entries(read.get_tag("SA"))
            # Process SA entries to find fusions
            sa_used = False
            if sa_entries:
                filtered = self._filter_sa_entries(
                    sa_entries, min_al_len_sa, min_sa_mapq
                )
                sa_used = self._process_sa_entries(
                    read, filtered, clip_side, jitter_window, min_sa_mapq
                )
            if clip_side and clip_len >= min_softclip_len:
                seen_pairs = set()  # Local scope for this read's softclip realignment
                self.realign_softclip(
                    read, clip_side, clip_len,
                    seen_pairs,
                    jitter_window=jitter_window,
                    min_sa_mapq=min_sa_mapq
                )
        bam.close()
        '''
        logger.info(f"Fusion detection complete: processed {processed_reads} reads, found {len(self.fusion_candidates)} fusion keys")
        logger.info(f"Fusion keys detected: {list(self.fusion_candidates.keys())}")
        '''

    def reconstruct_fusion_transcript(self, c1, p1, c2, p2, exon_padding=100):
        # Returns (sequence, left_exons, right_exons) or (None, None, None) on failure.
        try:
            # Get exons around left breakpoint (c1:p1)
            left_exons = list(self.db.region(region=(c1, max(1, p1 - exon_padding), p1 + exon_padding), featuretype='exon'))
            if not left_exons:
                return None, None, None
            left_exons = sorted([(int(e.start), int(e.end)) for e in left_exons])
            # Get exons around right breakpoint (c2:p2)
            right_exons = list(self.db.region(region=(c2, max(1, p2 - exon_padding), p2 + exon_padding), featuretype='exon'))
            if not right_exons:
                return None, None, None
            right_exons = sorted([(int(e.start), int(e.end)) for e in right_exons])
            # Extract sequences
            try:
                ref = pysam.FastaFile(self.reference_fasta)
            except Exception:
                return None, None, None
            left_seq = ""
            for s, e in left_exons:
                seq = ref.fetch(c1, s - 1, e)
                if seq:
                    left_seq += seq
            right_seq = ""
            for s, e in right_exons:
                seq = ref.fetch(c2, s - 1, e)
                if seq:
                    right_seq += seq
            ref.close()
            if not left_seq or not right_seq:
                return None, None, None
            # Construct fusion transcript: left portion + right portion
            fusion_transcript = left_seq + right_seq
            return fusion_transcript, left_exons, right_exons
        except Exception as e:
            logger.debug(f"Failed to reconstruct fusion transcript at {c1}:{p1}--{c2}:{p2}: {str(e)}")
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
                mapq = getattr(h, 'mapq', 0) or getattr(h, 'mapq', 0)
                logger.debug(f"Raw hit: q_st={q_st} q_en={q_en} match_len={match_len} mapq={mapq} ctg={getattr(h,'ctg',None)} r_st={getattr(h,'r_st',None)} r_en={getattr(h,'r_en',None)}")
                if match_len >= min_match_len and mapq >= min_mapq:
                    good_hits.append(h)
            logger.debug(f"Filtered {len(good_hits)} good hits out of {len(hits)} raw hits for fusion transcript")
            return good_hits
        except Exception as e:
            logger.debug(f"Failed to realign fusion transcript: {str(e)}")
            return []

    def _is_mitochondrial_candidate(self, c1, c2, left_gene, right_gene):
        # Return True if either side appears mitochondrial by chromosome or gene name.
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
        try:
            if isinstance(meta.get("left_gene"), str) and (meta["left_gene"].startswith("RP11") or meta["left_gene"].upper().startswith("ENS")):
                nearby = self._find_nearby_protein_coding(c1, p1, window=1000)
                if nearby:
                    meta["left_gene"] = nearby
                    meta.setdefault("notes", []).append(f"Left gene resolved to nearby coding {nearby}")
            if isinstance(meta.get("right_gene"), str) and (meta["right_gene"].startswith("RP11") or meta["right_gene"].upper().startswith("ENS")):
                nearby = self._find_nearby_protein_coding(c2, p2, window=1000)
                if nearby:
                    meta["right_gene"] = nearby
                    meta.setdefault("notes", []).append(f"Right gene resolved to nearby coding {nearby}")
        except Exception:
            pass
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

    @lru_cache(maxsize=50000)
    def _find_nearby_protein_coding(self, chrom, pos, window):
        # Search for a nearby protein-coding gene within +/- window and return its gene_name if found.
        try:
            # Use interval tree if available for faster lookups
            if self.interval_index is not None:
                genes = self.interval_index.get_genes_at(chrom, pos, window=window)
            else:
                start = max(1, pos - window)
                end = pos + window
                genes = list(self.db.region(region=(chrom, start, end), featuretype='gene'))
            for g in genes:
                # check common biotype attribute keys
                attrs = g.attributes if hasattr(g, 'attributes') else {}
                biotype = None
                for key in ('gene_type', 'gene_biotype', 'transcript_type', 'transcript_biotype'):
                    if key in attrs:
                        biotype = attrs.get(key, [None])[0]
                        break
                if biotype and biotype == 'protein_coding':
                    return attrs.get('gene_name', [g.id])[0]
            return None
        except Exception:
            return None

    def _compute_early_confidence(self, support_count):
        # Early confidence based only on support (before reconstruction/validation)
        # Used in build_metadata to decide gene reassignment strategy
        return min(support_count, 10) / 10.0

    def canonical_locus_name(self, gene_name):
        # Collapse antisense / divergent transcript names and map pseudogenes to parent genes.
        if not gene_name:
            return gene_name
        # Strip common antisense / divergent suffixes
        collapsed = Antisense_suffix.sub("", gene_name)
        # Try to map pseudogene to parent (pattern: ends with P followed by digits)
        # E.g. RPL23P6 → RPL23, SAE1P1 → SAE1
        match = re.match(r"^(.+?)P\d+$", collapsed)
        if match:
            parent_candidate = match.group(1)
            # Verify parent gene exists and is protein-coding
            try:
                parent_biotype = self.get_gene_biotype(parent_candidate)
                if parent_biotype == "protein_coding":
                    logger.debug(f"canonical_locus_name: Mapped pseudogene {gene_name} → parent {parent_candidate}")
                    return parent_candidate
            except Exception:
                pass
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

    def get_gene_biotype(self, gene_name, chrom=None, pos=None):
        # Retrieve the biotype (gene_type or gene_biotype) for a given gene name.
        # If gene_name is an HGNC symbol (not found as direct key), uses coordinates to match.
        # Returns biotype string or None if not found.
        if not gene_name:
            return None
        # First normalize via canonical_locus_name to handle pseudogenes and antisense
        canonical_name = self.canonical_locus_name(gene_name)
        # Then further normalize via normalize_gene_label
        normalized = self.normalize_gene_label(canonical_name)
        if normalized == "intergenic":
            # If normalization fails completely, try original gene_name too
            return self._get_gene_biotype_by_symbol_or_coords(gene_name, chrom=chrom, pos=pos)
        return self._get_gene_biotype_by_symbol_or_coords(normalized, chrom=chrom, pos=pos)

    def get_gene_strand(self, gene_name, chrom=None, pos=None):
        # Retrieve the strand (+/-) for a given gene name.
        if not gene_name:
            return None
        try:
            # Try direct lookup by gene ID first
            try:
                gene = self.db[gene_name]
                if gene.featuretype == "gene":
                    return gene.strand
            except Exception:
                pass
            # If not found and we have coordinates, look up genes at position
            if chrom and pos:
                if self.interval_index is not None:
                    genes = self.interval_index.get_genes_at(chrom, pos, window=500)
                else:
                    genes = list(self.db.region(region=(chrom, max(1, pos - 500), pos + 500), featuretype='gene'))
                for g in genes:
                    attrs = getattr(g, 'attributes', {}) or {}
                    gname = attrs.get('gene_name', [None])[0]
                    if gname == gene_name:
                        return g.strand
            # Last resort: search through all genes
            if self.interval_index is not None:
                genes = self.interval_index.find_genes_by_name(gene_name)
                if genes:
                    return genes[0].strand
        except Exception:
            pass
        return None

    def validate_candidates(self, min_support=2, window=25, require_gene_names=True,
                        require_mapq=10, allow_cis_sage=True, require_exon_boundary=True,
                        max_intra_chr_distance=None):
        self.build_metadata(min_support=min_support)
        self.fusion_assigned_pairs.clear()
        # Delegate validation and filtering to FusionValidator
        validator = FusionValidator(self)
        validator.filter_non_coding_genes()
        validator.apply_frequency_filters()
        for fusion_key, meta in self.fusion_metadata.items():
            # Check if already marked invalid by frequency filter
            if meta.get("is_valid") == False:
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

            # Compute full confidence (includes reconstruction/realignment bonuses)
            conf = validator.confidence(meta, flags)
            meta["confidence"] = conf
            if conf < 0.20 and meta.get("support", 0) < 2 and flags["is_valid"]:
                flags["is_valid"] = False
                flags["reasons"].append("Low confidence")
            meta.update(flags)

    def report(self, output_path="fusion_candidates.tsv", min_support=2,
                include_classes=("canonical","cis-SAGe"),
                min_confidence=0.3, only_valid=False):
            self.validate_candidates(min_support=min_support)
            validator = FusionValidator(self)
            validator.merge_nearly_identical()
            with open(output_path, "w") as f:
                f.write("LeftGene\tLeftBiotype\tRawLeftGene\tLeftChromosome\tLeftBreakpoint\t"
                        "RightGene\tRightBiotype\tRawRightGene\tRightChromosome\tRightBreakpoint\t"
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
                    raw_left = meta.get("raw_left_gene")
                    raw_right = meta.get("raw_right_gene")
                    left_biotype = meta.get("left_biotype", "unknown")
                    right_biotype = meta.get("right_biotype", "unknown")
                    f.write(f"{left_gene}\t{left_biotype}\t{raw_left}\t{left_chr}\t{left_pos}\t{right_gene}\t{right_biotype}\t"
                            f"{raw_right}\t{right_chr}\t{right_pos}\t{meta.get('support', 0)}\t{fusion_name}\t{meta.get('class')}\t"
                            f"{meta.get('is_valid')}\t{meta.get('confidence')}\t{reasons}\n")