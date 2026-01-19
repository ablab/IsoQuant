import re
import pysam
import gffutils
from collections import defaultdict
import mappy as mp
import logging
from functools import lru_cache
from sklearn.cluster import DBSCAN
import numpy as np


logger = logging.getLogger('IsoQuant')
ANTISENSE_SUFFIX_RE = re.compile(r"-(AS\d+|DT|DIVERGENT|NAT)$", re.IGNORECASE)

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
        # cache for resolved names: id_or_symbol -> gene_symbol (if found)
        self._resolved_name_cache = {}

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
        # Minimal gffutils calls; returns (gene_name, region_type)
        try:
            genes = list(self.db.region(region=(chrom, pos, pos), featuretype='gene'))
            if genes:
                gene_name = genes[0].attributes.get('gene_name', [genes[0].id])[0]
                # Are we in an exon of that gene?
                exons = list(self.db.region(region=(chrom, pos, pos), featuretype='exon'))
                if exons:
                    return gene_name, "exonic"
                else:
                    return gene_name, "intronic"
            else:
                # Not inside any gene; check if overlapping any exon 
                exons = list(self.db.region(region=(chrom, pos, pos), featuretype='exon'))
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

    def get_context(self, chrom, pos):
        gene_name, region_type = self._context_query(chrom, pos)
        return gene_name if gene_name else region_type

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
        if not cigar:
            return 0
        total = 0
        for m in re.finditer(r'(\d+)([MIDNSHP=XB])', cigar):
            length = int(m.group(1)); op = m.group(2)
            if op in ('M', '=', 'X'):
                total += length
        return total

    def realign_clipped_seq(self, seq):
        try:
            for hit in self.aligner.map(seq):
                hit_len = (getattr(hit, "r_en", None) - getattr(hit, "r_st", None)) if getattr(hit, "r_en", None) else getattr(hit, "mlen", None) or getattr(hit, "alen", None) or 0
                if hit_len and hit_len >= 50:
                    return hit  # first good hit
            return None
        except Exception:
            return None

    def record_fusion(self, context1, context2, read_name, chrom1, pos1, chrom2, pos2):        
        c1 = self.canonical_locus_name(context1)
        c2 = self.canonical_locus_name(context2)
        fusion_key = "--".join(sorted([c1, c2]))
        self.fusion_candidates[fusion_key].add(read_name)
        bp = (chrom1, int(pos1), chrom2, int(pos2))
        self.fusion_breakpoints[fusion_key][bp] += 1
        # also maintain supporting_reads set for quick access
        meta = self.fusion_metadata.setdefault(fusion_key, {"supporting_reads": set(), "consensus_bp": None,
                                                            "left_gene": None, "right_gene": None, "support": 0})
        meta["supporting_reads"].add(read_name)
        meta["support"] = len(meta["supporting_reads"])
        # consensus_bp will be calculated later in build_metadata()

    def build_metadata(self, min_support=1):
        # Compute consensus breakpoint and gene names for all fusion keys,
        for fusion_key, reads in self.fusion_candidates.items():
            support = len(reads)
            meta = self.fusion_metadata.setdefault(fusion_key, {"supporting_reads": set(), "consensus_bp": None,
                                                                "left_gene": None, "right_gene": None, "support": 0})
            meta["supporting_reads"].update(reads)
            meta["support"] = len(meta["supporting_reads"])
            if meta["support"] < min_support:
                continue
            bp_counts = self.fusion_breakpoints.get(fusion_key, {})
            if not bp_counts:
                continue
            consensus_bp, _ = max(bp_counts.items(), key=lambda item: item[1])
            meta["consensus_bp"] = consensus_bp
            left_chr, left_pos, right_chr, right_pos = consensus_bp
            # Normalize possible Ensembl / RP11 ids to gene symbols when possible
            left_ctx = self.get_context(left_chr, left_pos)
            right_ctx = self.get_context(right_chr, right_pos)
            meta["left_gene"] = self.resolve_gene_name(left_ctx)
            meta["right_gene"] = self.resolve_gene_name(right_ctx)

    def cluster_breakpoints(self, bp_counts, window=100):
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
        left_context = self.get_context(left_chr, left_pos)
        right_context = self.get_context(sa_chr, right_pos)
        left_context  = self.canonical_locus_name(left_context)
        right_context = self.canonical_locus_name(right_context)
        if left_context != right_context:
            self.record_fusion(left_context, right_context, read.query_name,
                            left_chr, left_pos, sa_chr, right_pos)
            seen_pairs.add(key)

    def detect_fusions(self,
                    min_al_len_primary=50,
                    min_al_len_sa=40,
                    min_sa_mapq=10,
                    min_softclip_len=40,
                    jitter_window=25,
                    realign_when_sa_present=True):
        bam = pysam.AlignmentFile(self.bam_path, "rb")
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            prim_al_len = self.compute_aligned_length(read)
            if prim_al_len < min_al_len_primary:
                continue
            prim_mapq = getattr(read, "mapping_quality", 0)
            if prim_mapq < min_sa_mapq:
                continue
            clip_side, clip_len = self.detect_softclip(read)
            sa_entries = []
            if read.has_tag("SA"):
                sa_entries = self.parse_sa_entries(read.get_tag("SA"))
            seen_pairs = set()
            sa_used = False
            if sa_entries:
                # Filter & prioritize SA entries by mapq and aligned length
                filtered = []
                for (sa_chr, sa_pos, sa_strand, sa_cigar, sa_mapq, sa_nm) in sa_entries:
                    if sa_pos is None or sa_cigar is None:
                        continue
                    sa_al_len = self.aligned_len_from_cigarstring(sa_cigar)
                    if sa_al_len < min_al_len_sa or sa_mapq < min_sa_mapq:
                        continue
                    filtered.append((sa_chr, sa_pos, sa_strand, sa_cigar, sa_mapq, sa_nm, sa_al_len))
                # keep top-N by sa_mapq or aligned length to limit work
                if len(filtered) > 3:
                    filtered.sort(key=lambda x: (x[4], x[6]), reverse=True)
                    filtered = filtered[:3]

                for (sa_chr, sa_pos, sa_strand, sa_cigar, sa_mapq, sa_nm, sa_al_len) in filtered:
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
                    left_context = self.get_context(left_chr, left_pos)
                    right_context = self.get_context(sa_chr, right_pos)
                    if left_context != right_context:
                        self.record_fusion(left_context, right_context, read.query_name,
                                        left_chr, left_pos, sa_chr, right_pos)
                        sa_used = True
            do_realign = False
            if clip_side and clip_len >= min_softclip_len:
                if not sa_used:
                    do_realign = True
                elif realign_when_sa_present:
                    do_realign = True
            if do_realign:
                self.realign_softclip(read, clip_side, clip_len,
                                            seen_pairs,
                                            jitter_window=jitter_window,
                                            min_sa_mapq=min_sa_mapq)
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

    @lru_cache(maxsize=200000)
    def _get_nearby_exons(self, chrom, start, end):
        # Cached query for exons in a region
        try:
            return tuple(list(self.db.region(region=(chrom, start, end), featuretype='exon')))
        except Exception:
            return ()

    def check_exon_boundary_proximity(self, c1, p1, c2, p2, delta=100):
        # Check if both breakpoints are within delta bp of exon boundaries.
        try:
            # Check left breakpoint proximity to exons
            left_exons = self._get_nearby_exons(c1, max(1, p1 - delta), p1 + delta)
            left_valid = False
            for exon in left_exons:
                # Check if breakpoint is near exon start or end
                if abs(p1 - int(exon.start)) <= delta or abs(p1 - int(exon.end)) <= delta:
                    left_valid = True
                    break
            # Check right breakpoint proximity to exons
            right_exons = self._get_nearby_exons(c2, max(1, p2 - delta), p2 + delta)
            right_valid = False
            for exon in right_exons:
                # Check if breakpoint is near exon start or end
                if abs(p2 - int(exon.start)) <= delta or abs(p2 - int(exon.end)) <= delta:
                    right_valid = True
                    break
            if not left_valid and not right_valid:
                return False, f"Both breakpoints >±{delta}bp from exon boundaries"
            elif not left_valid:
                return False, f"Left breakpoint {c1}:{p1} >±{delta}bp from exon boundary"
            elif not right_valid:
                return False, f"Right breakpoint {c2}:{p2} >±{delta}bp from exon boundary"
            return True, ""
        except Exception as e:
            logger.debug(f"Exon boundary check failed at {c1}:{p1}--{c2}:{p2}: {str(e)}")
            return False, f"Exon boundary check failed: {str(e)}"

    def validate_candidates(self, min_support=2, window=25, require_gene_names=True,
                        require_mapq=10, allow_cis_sage=True, require_exon_boundary=True,
                        max_intra_chr_distance=None):
        self.build_metadata(min_support=min_support)
        for fusion_key, meta in self.fusion_metadata.items():
            support = meta.get("support", 0)
            bp_counts = self.fusion_breakpoints.get(fusion_key, {})
            flags = {"is_valid": True, "reasons": [], "class": "canonical"}
            if support < min_support or not bp_counts:
                flags["is_valid"] = False
                flags["reasons"].append(f"Low support ({support} < {min_support})")
                meta.update(flags)
                continue

            # consensus via clustering
            cons = self.cluster_breakpoints(bp_counts, window=window)
            if not cons:
                flags["is_valid"] = False
                flags["reasons"].append("No stable breakpoint cluster")
                meta.update(flags)
                continue
            consensus_bp, clustered_support = cons
            meta["consensus_bp"] = consensus_bp
            meta["support"] = clustered_support

            # perform classification and basic filtering
            c1, p1, c2, p2 = consensus_bp
            g1_name, r1 = self._context_query(c1, p1)
            g2_name, r2 = self._context_query(c2, p2)
            self._apply_classification_and_filters(meta, flags, c1, p1, c2, p2,
                                                   g1_name, r1, g2_name, r2,
                                                   require_gene_names, require_exon_boundary,
                                                   allow_cis_sage, max_intra_chr_distance)

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

            conf = self.confidence(meta, flags)
            meta["confidence"] = conf
            # convert very low-confidence to invalid
            if conf < 0.20 and meta.get("support", 0) < 2 and flags["is_valid"]:
                flags["is_valid"] = False
                flags["reasons"].append("Low confidence")
            meta.update(flags)

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
                                          require_gene_names, require_exon_boundary,
                                          allow_cis_sage, max_intra_chr_distance):
        # Apply classification rules and basic filters, update meta and flags.
        # Resolve Ensembl/RP11 ids to gene symbols when possible, but keep region type fallback
        resolved_left = self.resolve_gene_name(g1_name) if g1_name else r1
        resolved_right = self.resolve_gene_name(g2_name) if g2_name else r2
        meta["left_gene"] = resolved_left
        meta["right_gene"] = resolved_right

        # Prefer nearby protein-coding genes if the chosen gene looks like RP11/ENSG/non-coding
        try:
            if isinstance(meta.get("left_gene"), str) and (meta["left_gene"].startswith("RP11") or meta["left_gene"].upper().startswith("ENS")):
                nearby = self._find_nearby_protein_coding(c1, p1, window=500)
                if nearby:
                    meta["left_gene"] = nearby
                    meta.setdefault("notes", []).append(f"Left gene resolved to nearby coding {nearby}")
            if isinstance(meta.get("right_gene"), str) and (meta["right_gene"].startswith("RP11") or meta["right_gene"].upper().startswith("ENS")):
                nearby = self._find_nearby_protein_coding(c2, p2, window=500)
                if nearby:
                    meta["right_gene"] = nearby
                    meta.setdefault("notes", []).append(f"Right gene resolved to nearby coding {nearby}")
        except Exception:
            pass

        # classify
        if g1_name and g2_name:
            # compare resolved gene symbols for intragenic check
            g1 = self.canonical_locus_name(self.resolve_gene_name(g1_name))
            g2 = self.canonical_locus_name(self.resolve_gene_name(g2_name))
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
        # exon boundary proximity filter: ±30 bp of exon junctions
        if require_exon_boundary:
            exon_valid, exon_reason = self.check_exon_boundary_proximity(c1, p1, c2, p2, delta=30)
            if not exon_valid:
                flags["is_valid"] = False
                flags["reasons"].append(exon_reason)

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
    def _find_nearby_protein_coding(self, chrom, pos, window=500):
        # Search for a nearby protein-coding gene within +/- window and return its gene_name if found.
        try:
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

    def confidence(self, meta, flags=None):
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
        collapsed = ANTISENSE_SUFFIX_RE.sub("", gene_name)
        return collapsed


    def report(self, output_path="fusion_candidates.tsv", min_support=2,
            include_classes=("canonical","cis-SAGe"),
            min_confidence=0.3, only_valid=False):
        self.validate_candidates(min_support=min_support)
        with open(output_path, "w") as f:
            f.write("LeftGene\tLeftChromosome\tLeftBreakpoint\tRightGene\tRightChromosome\tRightBreakpoint\tSupportingReads\tFusionName\tClass\tValid\tReasons\n")
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
                left_gene = self.resolve_gene_name(meta.get("left_gene") or self.get_context(left_chr, left_pos))
                right_gene = self.resolve_gene_name(meta.get("right_gene") or self.get_context(right_chr, right_pos))
                # Skip intragenic fusions (same gene on both sides)
                if left_gene == right_gene:
                    continue
                # Create fusion name using resolved gene names
                fusion_name = f"{left_gene}-{right_gene}"
                reasons = ";".join(meta.get("reasons", []))
                f.write(f"{left_gene}\t{left_chr}\t{left_pos}\t{right_gene}\t{right_chr}\t{right_pos}\t{meta.get('support', 0)}\t{fusion_name}\t{meta.get('class')}\t{meta.get('is_valid')}\t{reasons}\n")