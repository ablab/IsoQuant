import re
import pysam
import gffutils
from collections import defaultdict
import mappy as mp
import logging
from functools import lru_cache
from intervaltree import IntervalTree
from .genomic_interval_index import GenomicIntervalIndex

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
        # Build interval tree index for fast coordinate-to-gene/exon mapping
        logger.info("Initializing genomic interval index for efficient database queries...")
        if IntervalTree is not None:
            self.interval_index = GenomicIntervalIndex(self.db)
        else:
            logger.warning("intervaltree not available; falling back to gffutils queries")
            self.interval_index = None
        self.debug = False
        self.debug_sink = None

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

    def resolve_gene_name(self, name_or_id):
        # Resolve Ensembl or other IDs (e.g. ENS*, RP11*) to gene_name when possible.
        # If `name_or_id` is already a gene symbol or cannot be resolved, return it unchanged.
        if not name_or_id:
            return name_or_id
        # quick cache hit
        if name_or_id in self._resolved_name_cache:
            return self._resolved_name_cache[name_or_id]
        # if looks like an Ensembl or RP11 identifier, try to resolve via genedb
        resolved = name_or_id
        try:
            # Try direct lookup by id
            try:
                feat = self.db[name_or_id]
            except Exception:
                feat = None
            if feat is not None:
                # If it's a gene, prefer gene_name attribute
                if feat.featuretype == 'gene':
                    resolved = feat.attributes.get('gene_name', [feat.id])[0]
                else:
                    # For transcript or other features, try to get parent gene id
                    if 'gene_id' in feat.attributes:
                        gid = feat.attributes['gene_id'][0]
                        try:
                            g = self.db[gid]
                            resolved = g.attributes.get('gene_name', [g.id])[0]
                        except Exception:
                            resolved = gid
                    elif 'Parent' in feat.attributes:
                        gid = feat.attributes['Parent'][0]
                        try:
                            g = self.db[gid]
                            resolved = g.attributes.get('gene_name', [g.id])[0]
                        except Exception:
                            resolved = gid
                    else:
                        # fallback to feature id
                        resolved = feat.attributes.get('gene_name', [feat.id])[0] if hasattr(feat, 'attributes') else feat.id
            else:
                # As a last resort, try to find a gene feature whose id equals the provided id
                try:
                    genes = list(self.db.features_of_type('gene', id=name_or_id))
                    if genes:
                        resolved = genes[0].attributes.get('gene_name', [genes[0].id])[0]
                except Exception:
                    pass
        except Exception:
            resolved = name_or_id
        # cache and return
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

    def normalize_and_resolve_gene_name(self, chrom, pos, original_name, window=500):
        # Normalize gene name: strip antisense suffixes and try to find parent gene
        if not original_name:
            return original_name
        # If has antisense suffix, try to find the parent gene
        if self.has_antisense_suffix(original_name):
            stripped = self.strip_antisense_suffix(original_name)
            # Try to find parent gene without the suffix nearby
            try:
                # Use interval tree if available
                if self.interval_index is not None:
                    genes = self.interval_index.get_genes_at(chrom, pos, window=window)
                else:
                    start = max(1, pos - window)
                    end = pos + window
                    genes = list(self.db.region(region=(chrom, start, end), featuretype='gene'))
                for g in genes:
                    gene_name = g.attributes.get('gene_name', [g.id])[0]
                    if gene_name and gene_name.upper() == stripped.upper():
                        return gene_name
            except Exception:
                pass
            # If we can't find parent, return the stripped version
            return stripped
        return original_name

    def _compute_fusion_gene_score(self, g, pos, exon_min_dist, boundary_min_dist, exonic_hit, body_dist, gstart, gend):
        attrs  = getattr(g, "attributes", {}) or {}
        gtype  = (attrs.get("gene_type", [None])[0]
                or attrs.get("gene_biotype", [None])[0]
                or attrs.get("transcript_biotype", [None])[0])
        score = 0.0
        # 1) hard evidence first
        if gstart - 3000 <= pos <= gend + 3000:
            score += 50  # near edge bonus
        if exonic_hit:
            score += 200.0
        if boundary_min_dist is not None:
            score += max(0.0, 80.0 - min(boundary_min_dist, 200))  # closer boundary → higher
        if exon_min_dist is not None:
            score += max(0.0, 40.0 - min(exon_min_dist, 500))
        # 2) gene body proximity (weak)
        if body_dist:
            score += max(0.0, 20.0 - min(body_dist, 2000) / 100.0)
        # 3) biotype weighting
        if gtype == "protein_coding":
            score += 60.0
        else:
            # Penalize any pseudogene annotations 
            if gtype and "pseudogene" in gtype.lower():
                score -= 90.0
            # Mild penalty for various non-coding/ambiguous biotypes
            elif gtype in ("antisense", "lncRNA", "lincRNA", "processed_transcript",
                        "sense_intronic", "sense_overlapping", "transcribed_unprocessed_pseudogene"):
                score -= 15.0
        # 4) gentle preference for longer exon span (more likely real coding gene)
        try:
            # sum exon lengths (per gene) cheaply
            exon_len_sum = 0
            exons = self.interval_index.get_exons_of_gene(g) if self.interval_index else list(self.db.children(g, featuretype="exon", order_by="start"))
            for ex in exons:
                exon_len_sum += int(ex.end) - int(ex.start) + 1
            score += min(exon_len_sum / 1e4, 10.0)  # cap at +10
        except Exception:
            pass
        return score
    
    def assign_fusion_gene(self, chrom, pos, window=2000):
        try:
            if self.interval_index is not None:
                genes = self.interval_index.get_genes_at(chrom, pos, window=window)
            else:
                genes = list(self.db.region(
                    region=(chrom, max(1, pos - window), pos + window),
                    featuretype="gene"
                ))
        except Exception:
            return None
        if not genes:
            return None
        INF = 10**9
        best_gene = None
        best_key  = None  # tuple key for lexicographic comparison
        for g in genes:
            attrs  = getattr(g, "attributes", {}) or {}
            gname  = attrs.get("gene_name", [getattr(g, "id", None)])[0]
            gtype  = (attrs.get("gene_type", [None])[0]
                or attrs.get("gene_biotype", [None])[0]
                or attrs.get("transcript_biotype", [None])[0])
            gstart = int(getattr(g, "start", 0))
            gend   = int(getattr(g, "end",   0))
            # gstr = getattr(g, "strand", ".")  # currently unused here

            # --- exon distances ---
            exon_min_dist = None
            boundary_min_dist = None
            exonic_hit = False
            try:
                exons = (self.interval_index.get_exons_of_gene(g)
                        if self.interval_index
                        else list(self.db.children(g, featuretype="exon", order_by="start")))
                for ex in exons:
                    es = int(ex.start); ee = int(ex.end)
                    # span distance
                    if es <= pos <= ee:
                        d_span = 0
                        exonic_hit = True
                    else:
                        d_span = min(abs(pos-es), abs(pos-ee))
                    exon_min_dist = d_span if exon_min_dist is None else min(exon_min_dist, d_span)
                    # boundary distance
                    d_bound = min(abs(pos-es), abs(pos-ee))
                    boundary_min_dist = d_bound if boundary_min_dist is None else min(boundary_min_dist, d_bound)
            except Exception:
                exon_min_dist = None
                boundary_min_dist = None
                exonic_hit = False
            # body distance if outside gene
            if gstart <= pos <= gend:
                body_dist = 0
            else:
                body_dist = min(abs(pos - gstart), abs(pos - gend))
            # compute score with your refactored function
            score = self._compute_fusion_gene_score(
                g, pos, exon_min_dist, boundary_min_dist, exonic_hit, body_dist, gstart, gend
            )
            # normalize None distances
            bnd = boundary_min_dist if boundary_min_dist is not None else INF
            exd = exon_min_dist     if exon_min_dist     is not None else INF
            bod = body_dist if body_dist else 0
            # optional: gene length as a gentle preference (cap to avoid domination)
            glen = max(1, gend - gstart + 1)
            # build a stable key; higher is better for earlier elements
            # CRITICAL: protein_coding must come EARLY to prevent pseudogenes from winning based on proximity
            key = (
                1 if gtype == "protein_coding" else 0,  # 1) PRIORITY: protein-coding must come first
                score,                           # 2) main score (including biotype penalties)
                bool(exonic_hit),                # 3) prefer exonic
                -bnd,                            # 4) nearer exon boundary
                -exd,                            # 5) nearer exon span
                -bod,                            # 6) nearer gene body
                -min(glen, 500000),              # 7) gentle preference for longer locus
                gname or getattr(g, "id", "")    # 8) deterministic
            )
            if best_key is None or key > best_key:
                best_key  = key
                best_gene = gname or getattr(g, "id", None)
        return best_gene

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
            left_gene = self.assign_fusion_gene(left_chr, left_pos)
            right_gene = self.assign_fusion_gene(sa_chr, right_pos)
            if left_gene is None or right_gene is None:
                continue
            if left_gene != right_gene:
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

    def _safe_gene_token(self, g):
        # Convert None/empty to 'intergenic'
        if not g or (isinstance(g, str) and not g.strip()):
            return "intergenic"
        return str(g)

    def record_fusion(self, context1, context2, read_name, chrom1, pos1, chrom2, pos2):
        # Skip fusions involving antisense/regulatory genes (AS1, DT, etc.)
        if self.has_antisense_suffix(context1) or self.has_antisense_suffix(context2):
            return
        # Skip fusions with intergenic partners
        if context1 == "intergenic" or context2 == "intergenic":
            return
        left = self._safe_gene_token(context1)
        right = self._safe_gene_token(context2)
        fusion_key = "--".join(sorted([left, right]))

        self.fusion_candidates[fusion_key].add(read_name)
        bp = (chrom1, int(pos1), chrom2, int(pos2))
        self.fusion_breakpoints[fusion_key][bp] += 1

        meta = self.fusion_metadata.setdefault(
            fusion_key,
            {"supporting_reads": set(), "consensus_bp": None,
            "left_gene": None, "right_gene": None, "support": 0}
        )
        meta["supporting_reads"].add(read_name)
        meta["support"] = len(meta["supporting_reads"])
        # store the raw assigned pair for this read so we don't need to re-run assignment later
        try:
            self.fusion_assigned_pairs[fusion_key][read_name] = (context1, context2)
        except Exception:
            pass

    def normalize_gene_label(self, gene):
        # Normalize any gene identifier to a canonical gene symbol string.
        if not gene:
            return "intergenic"
        # Resolve ENSG / IDs → symbol if possible
        gene = self.resolve_gene_name(gene)
        # Collapse antisense / divergent
        gene = self.canonical_locus_name(gene)
        # Safety net
        if gene is None or gene.startswith("ENSG"):
            return "intergenic"
        return gene

    def build_metadata(self, min_support=1):
        # Delegate the heavy lifting to the new FusionMetadata helper class
        try:
            from .fusion_metadata import FusionMetadata
            FusionMetadata(self).process_all(min_support=min_support)
        except Exception:
            # Fallback: if import fails, keep original behavior minimal (no-op)
            return

    def cluster_breakpoints(self, bp_counts, window):
        # bp_counts: dict[(chr1,pos1,chr2,pos2)] -> count
        # Group by chr pairs, then cluster positions by window
        from collections import defaultdict
        clusters = defaultdict(list)
        for (c1,p1,c2,p2), cnt in bp_counts.items():
            key = (c1, c2)
            clusters[key].append((p1, p2, cnt))
        def cluster_positions(items):
            items = sorted(items)
            current = [items[0]]
            clustered = []
            for itm in items[1:]:
                if abs(itm[0]-current[-1][0]) <= window and abs(itm[1]-current[-1][1]) <= window:
                    current.append(itm)
                else:
                    clustered.append(current)
                    current = [itm]
            clustered.append(current)
            # pick the cluster with largest total count
            best = max(clustered, key=lambda cl: sum(x[2] for x in cl))
            # consensus: weighted median
            import numpy as np
            p1s = np.array([x[0] for x in best]); w = np.array([x[2] for x in best])
            p2s = np.array([x[1] for x in best])
            cons_p1 = int(np.average(p1s, weights=w))
            cons_p2 = int(np.average(p2s, weights=w))
            total = int(np.sum(w))
            return cons_p1, cons_p2, total
        # choose the chr pair with max clustered support
        best_pair, best_cons = None, (None,None,0)
        for key, items in clusters.items():
            c_p1, c_p2, total = cluster_positions(items)
            if total > best_cons[2]:
                best_pair, best_cons = key, (c_p1, c_p2, total)
        if best_pair is None:
            return None
        (c1, c2), (p1, p2, total) = best_pair, best_cons
        return (c1, p1, c2, p2), total

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
        left_gene = self.assign_fusion_gene(left_chr, left_pos)
        right_gene = self.assign_fusion_gene(sa_chr, right_pos)
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
        bam = pysam.AlignmentFile(self.bam_path, "rb")
        for read in bam:
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

    def check_gene_exon_boundary(self, gene_name, chrom, pos, delta=50):
        # find the gene feature by name, get its exons, and check |pos - exon boundary| <= delta
        try:
            # Resolve to a gene feature
            g = None
            try:
                g = self.db[gene_name]  # if name==id
            except Exception:
                # fallback: search genes in small window and match by gene_name attribute
                # Use interval tree if available
                genes = []
                if self.interval_index is not None:
                    genes = self.interval_index.get_genes_at(chrom, pos, window=2000)
                else:
                    for cand in self.db.region(region=(chrom, max(1, pos-2000), pos+2000), featuretype="gene"):
                        genes.append(cand)
                for cand in genes:
                    if cand.attributes.get("gene_name", [""])[0] == gene_name:
                        g = cand; break
            if g is None:  # give benefit of doubt
                return True, ""
            exons = self.interval_index.get_exons_of_gene(g) if self.interval_index else list(self.db.children(g, featuretype="exon", order_by="start"))
            for ex in exons:
                if abs(pos - int(ex.start)) <= delta or abs(pos - int(ex.end)) <= delta:
                    return True, ""
            return False, f"{gene_name} breakpoint {chrom}:{pos} >±{delta}bp from its exon boundary"
        except Exception as e:
            return False, f"Exon boundary check failed: {e}"

    # validate_candidates moved to FusionMetadata

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

    def confidence(self, meta, flags=None):
        # Full confidence computation including reconstruction/realignment bonuses
        # Start from clustered support normalized
        support = min(meta.get("support", 0), 10) / 10.0
        # Reconstruction bonus
        recon = 0.2 if meta.get("reconstruction_ok") else 0.0
        # Realignment bonus: scale by number of hits and MAPQ
        hits = meta.get("realignment_hits", 0)
        mapq = meta.get("best_hit_mapq", 0)
        realign = min(hits, 3) * 0.1 + min(mapq, 30) / 300.0  # up to ~0.2
        # Penalties for missing gene names and intergenic
        penalties = 0.0
        if "Missing gene name on one side" in meta.get("reasons", []):
            penalties += 0.1
        if meta.get("class") == "intergenic":
            penalties += 0.1
        conf = max(0.0, min(1.0, support + recon + realign - penalties))
        return conf

    def canonical_locus_name(self, gene_name):
        # Collapse antisense / divergent transcript names to the canonical locus name.
        if not gene_name:
            return gene_name
        # Strip common antisense / divergent suffixes
        collapsed = Antisense_suffix.sub("", gene_name)
        return collapsed

    def _meta_has_pseudogene_partner(self, meta):
        # Return True if either consensus-assigned partner is a pseudogene/ncRNA
        consensus = meta.get("consensus_bp")
        if not consensus:
            return False
        c1, p1, c2, p2 = consensus
        left = meta.get("left_gene")
        right = meta.get("right_gene")
        try:
            if self.is_bad_gene(left, c1, p1):
                return True
            if self.is_bad_gene(right, c2, p2):
                return True
        except Exception:
            return False
        return False

    def _exon_boundary_issues(self, meta, c1, p1, c2, p2):
        # Return list of human-readable exon-boundary issue strings (empty if none)
        issues = []
        left_gene = meta.get("left_gene")
        right_gene = meta.get("right_gene")
        try:
            okL, whyL = self.check_gene_exon_boundary(left_gene, c1, p1, delta=50)
            okR, whyR = self.check_gene_exon_boundary(right_gene, c2, p2, delta=50)
            if not okL and whyL:
                issues.append(whyL)
            if not okR and whyR:
                issues.append(whyR)
        except Exception:
            issues.append("Exon boundary check failed")
        return issues

    def _build_symbol_biotype_index(self):
        # One-time index construction: map HGNC symbols to biotypes using all genes in db
        # This is called lazily when first needed
        if self._symbol_biotype_cache:
            return  # already built
        logger.debug("Building symbol → biotype index from genedb...")
        try:
            for gene in self.db.features_of_type('gene'):
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
            return self._symbol_biotype_cache[gene_symbol]
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
                            # Cache for future lookups
                            self._symbol_biotype_cache[gene_symbol] = biotype
                            return biotype
            except Exception as e:
                logger.debug(f"Failed to query genes at {chrom}:{pos} for symbol {gene_symbol}: {e}")
        return None

    def get_gene_biotype(self, gene_name, chrom=None, pos=None):
        # Retrieve the biotype (gene_type or gene_biotype) for a given gene name.
        # If gene_name is an HGNC symbol (not found as direct key), uses coordinates to match.
        # Returns biotype string or None if not found.
        if not gene_name:
            return None
        return self._get_gene_biotype_by_symbol_or_coords(gene_name, chrom=chrom, pos=pos)

    def is_bad_gene(self, gene_name, chrom=None, pos=None):
        # Check if a gene is "bad" (pseudogene, ncRNA, etc.) based on its biotype.
        # Uses HGNC symbol with coordinate-based fallback if needed.
        if not gene_name:
            return False
        biotype = self._get_gene_biotype_by_symbol_or_coords(gene_name, chrom=chrom, pos=pos)
        if biotype is None:
            return False
        biotype_lower = biotype.lower()
        bad_types = (
            'pseudogene', 'processed_pseudogene', 'unprocessed_pseudogene',
            'transcribed_unitary_pseudogene', 'polymorphic_pseudogene', 'unitary_pseudogene',
            'lncrna', 'lincrna', 'antisense', 'sense_intronic', 'sense_overlapping',
            'snrna', 'mirna'
        )
        return any(bt in biotype_lower for bt in bad_types)

    def validate_candidates(self, min_support=2, window=25, require_gene_names=True,
                        require_mapq=10, allow_cis_sage=True, require_exon_boundary=True,
                        max_intra_chr_distance=None):
        self.build_metadata(min_support=min_support)
        for fusion_key, meta in self.fusion_metadata.items():
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
            if flags["is_valid"] and require_exon_boundary:
                left_gene  = meta.get("left_gene")
                right_gene = meta.get("right_gene")
                okL, whyL = self.check_gene_exon_boundary(left_gene,  c1, p1, delta=50)
                okR, whyR = self.check_gene_exon_boundary(right_gene, c2, p2, delta=50)
                if not okL or not okR:
                    flags["is_valid"] = False
                    # Append specific reasons (only non-empty strings)
                    reasons = []
                    if whyL: reasons.append(whyL)
                    if whyR: reasons.append(whyR)
                    if reasons:
                        flags["reasons"].append("; ".join(reasons))

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
            conf = self.confidence(meta, flags)
            meta["confidence"] = conf
            # convert very low-confidence to invalid
            if conf < 0.20 and meta.get("support", 0) < 2 and flags["is_valid"]:
                flags["is_valid"] = False
                flags["reasons"].append("Low confidence")
            meta.update(flags)

    def report(self, output_path="fusion_candidates.tsv", min_support=2,
            include_classes=("canonical","cis-SAGe"),
            min_confidence=0.3, only_valid=False):
        self.validate_candidates(min_support=min_support)
        with open(output_path, "w") as f:
            f.write("LeftGene\tLeftBiotype\tRawLeftGene\tLeftChromosome\tLeftBreakpoint\t"
                    "RightGene\tRightBiotype\tRawRightGene\tRightChromosome\tRightBreakpoint\t"
                    "SupportingReads\tFusionName\tClass\tValid\tReasons\n")
            for fusion_key, meta in sorted(self.fusion_metadata.items(), key=lambda x: -x[1].get("support", 0)):
                if meta.get("confidence", 0) < min_confidence:
                    continue
                if only_valid and not meta.get("is_valid", True):
                    continue
                if meta.get("class") not in include_classes:
                    continue
                if not meta.get("consensus_bp"):
                    continue
                left_chr, left_pos, right_chr, right_pos = meta["consensus_bp"]
                left_gene  = meta["left_gene"]
                right_gene = meta["right_gene"]
                # Exclude intergenic fusions (no true positives)
                if left_gene == "intergenic" or right_gene == "intergenic":
                    continue
                # Exclude mitochondrial fusion candidates from report
                if self._is_mitochondrial_candidate(left_chr, right_chr, left_gene, right_gene):
                    continue
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
                if self.is_bad_gene(left_gene,  left_chr,  left_pos) or self.is_bad_gene(right_gene, right_chr, right_pos):
                    continue
                f.write(f"{left_gene}\t{left_biotype}\t{raw_left}\t{left_chr}\t{left_pos}\t{right_gene}\t{right_biotype}\t"
                        "{raw_right}\t{right_chr}\t{right_pos}\t{meta.get('support', 0)}\t{fusion_name}\t{meta.get('class')}\t"
                        "{meta.get('is_valid')}\t{reasons}\n")