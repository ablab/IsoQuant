from collections import defaultdict
import logging

logger = logging.getLogger('IsoQuant')

class FusionMetadata:
    def __init__(self, detector):
        self.detector = detector

    def process_all(self, min_support=1):
        detector = self.detector
        logger.info(
            f"FusionMetadata: Processing {len(detector.fusion_candidates)} fusion keys "
            f"with min_support={min_support}"
        )
        # Iterate through fusion candidates created earlier in record_fusion()
        for fusion_key, reads in list(detector.fusion_candidates.items()):
            self._process_fusion_candidate(fusion_key, reads, min_support)

    def _process_fusion_candidate(self, fusion_key, reads, min_support):
        # Process a single fusion candidate through the metadata pipeline."""
        detector = self.detector
        # 1. Apply early support + breakpoint-clustering gates
        gate_result = self._passes_all_gates(fusion_key, reads, min_support)
        if not gate_result:
            meta = detector.fusion_metadata.get(fusion_key)
            if meta is not None:
                meta["is_valid"] = False
                meta.setdefault("reasons", []).append("Failed gating (support or clustering)")
            return
        support, consensus_bp, clustered_support = gate_result
        left_chr, left_pos, right_chr, right_pos = consensus_bp
        # 2. Get or initialize metadata
        meta = self._initialize_or_get_metadata(fusion_key, reads, clustered_support, consensus_bp)
        # 3. Set raw genes and compute raw biotypes
        self._compute_raw_genes_and_biotypes(fusion_key, meta, left_chr, left_pos, right_chr, right_pos)
        # 4. Determine final genes based on raw gene scores (threshold: 230)
        final_left_gene, final_right_gene = self._compute_final_genes(
            meta, left_chr, left_pos, right_chr, right_pos
        )
        # 5. Update fusion_key mappings if genes changed
        original_key = fusion_key
        new_key = self._update_fusion_key_mappings(original_key, final_left_gene, final_right_gene, meta)
        # 6. Update metadata with final genes and biotypes
        meta["left_gene"] = final_left_gene
        meta["right_gene"] = final_right_gene
        self._compute_final_biotypes(meta, final_left_gene, final_right_gene, left_chr, left_pos, right_chr, right_pos)
        # 7. Mark fusion valid
        meta["is_valid"] = True
        meta.setdefault("reasons", [])

    def _initialize_or_get_metadata(self, fusion_key, reads, clustered_support, consensus_bp):
        # Get existing metadata or create minimal metadata if missing.
        detector = self.detector
        meta = detector.fusion_metadata.get(fusion_key)
        if meta is None:
            meta = {
                "supporting_reads": set(),
                "left_gene": None,
                "right_gene": None,
            }
            detector.fusion_metadata[fusion_key] = meta
        # Update support information
        meta["supporting_reads"] = set(reads)
        meta["support"] = clustered_support
        meta["consensus_bp"] = consensus_bp
        return meta

    def _compute_raw_genes_and_biotypes(self, fusion_key, meta, left_chr, left_pos, right_chr, right_pos):
        # Assign raw genes from fusion_assigned_pairs and compute their biotypes.
        detector = self.detector
        # Populate raw genes from read-level assignments
        self._assign_raw_genes(fusion_key, meta, left_chr, left_pos, right_chr, right_pos)
        # Compute biotypes for raw genes
        if meta.get("raw_left_gene"):
            meta["raw_left_biotype"] = detector.get_gene_biotype(
                meta["raw_left_gene"], chrom=left_chr, pos=left_pos
            )
        if meta.get("raw_right_gene"):
            meta["raw_right_biotype"] = detector.get_gene_biotype(
                meta["raw_right_gene"], chrom=right_chr, pos=right_pos
            )
        logger.debug(
            f"Raw assignments for {fusion_key}: "
            f"left={meta.get('raw_left_gene')}({meta.get('raw_left_biotype')}), "
            f"right={meta.get('raw_right_gene')}({meta.get('raw_right_biotype')})"
        )

    def _compute_final_genes(self, meta, left_chr, left_pos, right_chr, right_pos):
        # Determine final genes based on raw gene scores (threshold: 230).
        # If raw score >= 230, use raw gene; otherwise re-assign at consensus breakpoint.
        detector = self.detector
        SCORE_THRESHOLD = 230
        # Left gene decision
        raw_left_score = meta.get("raw_left_score", 0.0)
        if raw_left_score >= SCORE_THRESHOLD:
            # Score is good enough; use raw gene
            final_left_gene = meta.get("raw_left_gene")
            meta["final_left_score"] = raw_left_score
            logger.debug(f"Left gene score {raw_left_score:.1f} >= {SCORE_THRESHOLD}: Using raw gene {final_left_gene}")
        else:
            # Score is low; re-assign at consensus breakpoint
            final_left_gene, final_left_score = detector.assign_fusion_gene(left_chr, left_pos)
            meta["final_left_score"] = final_left_score
            logger.debug(
                f"Left gene score {raw_left_score:.1f} < {SCORE_THRESHOLD}: "
                f"Re-assigned to {final_left_gene} with score {final_left_score:.1f}"
            )
        # Right gene decision
        raw_right_score = meta.get("raw_right_score", 0.0)
        if raw_right_score >= SCORE_THRESHOLD:
            # Score is good enough; use raw gene
            final_right_gene = meta.get("raw_right_gene")
            meta["final_right_score"] = raw_right_score
            logger.debug(f"Right gene score {raw_right_score:.1f} >= {SCORE_THRESHOLD}: Using raw gene {final_right_gene}")
        else:
            # Score is low; re-assign at consensus breakpoint
            final_right_gene, final_right_score = detector.assign_fusion_gene(right_chr, right_pos)
            meta["final_right_score"] = final_right_score
            logger.debug(
                f"Right gene score {raw_right_score:.1f} < {SCORE_THRESHOLD}: "
                f"Re-assigned to {final_right_gene} with score {final_right_score:.1f}"
            )
        return final_left_gene, final_right_gene

    def _update_fusion_key_mappings(self, original_key, final_left_gene, final_right_gene, meta):
        # Update fusion_key mappings if gene names changed during re-assignment.
        detector = self.detector
        new_left = final_left_gene or "intergenic"
        new_right = final_right_gene or "intergenic"
        new_key = f"{new_left}--{new_right}"
        # If gene names changed, update all mappings
        if new_key != original_key:
            logger.info(f"Gene assignment changed: {original_key} → {new_key}")
            detector.fusion_metadata[new_key] = meta
            if original_key in detector.fusion_metadata:
                del detector.fusion_metadata[original_key]
            # Update fusion_candidates mapping
            detector.fusion_candidates[new_key] = detector.fusion_candidates.pop(original_key, set())
            # Update fusion_breakpoints mapping
            if original_key in detector.fusion_breakpoints:
                detector.fusion_breakpoints[new_key] = detector.fusion_breakpoints.pop(original_key)
            # Update fusion_assigned_pairs mapping
            if original_key in detector.fusion_assigned_pairs:
                detector.fusion_assigned_pairs[new_key] = detector.fusion_assigned_pairs.pop(original_key)
        return new_key

    def _compute_final_biotypes(self, meta, final_left_gene, final_right_gene, left_chr, left_pos, right_chr, right_pos):
        detector = self.detector
        if final_left_gene:
            meta["left_biotype"] = detector.get_gene_biotype(
                final_left_gene, chrom=left_chr, pos=left_pos
            )
        if final_right_gene:
            meta["right_biotype"] = detector.get_gene_biotype(
                final_right_gene, chrom=right_chr, pos=right_pos
            )

    def _passes_all_gates(self, fusion_key, reads, min_support):
        #  Validate fusion against support, breakpoint, and consensus gates.
        # Returns (support, consensus_bp, clustered_support) if all gates pass, None otherwise.
        # Gate 1: Support threshold
        support = len(reads)
        if support < min_support:
            return None
        # Gate 2: Breakpoint counts exist
        bp_counts = self.detector.fusion_breakpoints.get(fusion_key, {})
        if not bp_counts:
            return None
        # Gate 3: Consensus clustering succeeds
        consensus_result = self.detector.cluster_breakpoints(bp_counts, window=2000)
        if not consensus_result:
            return None
        consensus_bp, clustered_support = consensus_result
        return support, consensus_bp, clustered_support

    def _assign_raw_genes(self, fusion_key, meta, left_chr, left_pos, right_chr, right_pos):
        # Collect raw gene assignments from read-level assignments in fusion_assigned_pairs.
        # Directly populate meta["raw_left_gene"], meta["raw_right_gene"] and their scores
        assigned = self.detector.fusion_assigned_pairs.get(fusion_key, {})
        raw_left, raw_right = self._collect_raw_assignments(
            meta["supporting_reads"], assigned, left_chr, left_pos, right_chr, right_pos
        )
        # Fall back to coordinate-based assignment if needed
        if raw_left is None:
            raw_left, left_score = self.detector.assign_fusion_gene(left_chr, left_pos)
            meta["raw_left_score"] = left_score
        else:
            # Compute score for read-consensus gene too
            _, left_score = self.detector.assign_fusion_gene(left_chr, left_pos)
            meta["raw_left_score"] = left_score
        if raw_right is None:
            raw_right, right_score = self.detector.assign_fusion_gene(right_chr, right_pos)
            meta["raw_right_score"] = right_score
        else:
            # Compute score for read-consensus gene too
            _, right_score = self.detector.assign_fusion_gene(right_chr, right_pos)
            meta["raw_right_score"] = right_score
        # Populate metadata with raw genes
        meta["raw_left_gene"] = raw_left
        meta["raw_right_gene"] = raw_right

    def _collect_raw_assignments(self, supporting_reads, assigned, left_chr, left_pos, right_chr, right_pos):
        left_counts = defaultdict(int)
        right_counts = defaultdict(int)
        for r in supporting_reads:
            if r in assigned:
                lpair, rpair = assigned[r]
                if lpair:
                    left_counts[lpair] += 1
                if rpair:
                    right_counts[rpair] += 1

        raw_left = self._pick_most_common(left_counts)
        raw_right = self._pick_most_common(right_counts)
        return raw_left, raw_right

    @staticmethod
    def _pick_most_common(counts_dict):
        if not counts_dict:
            return None
        return sorted(counts_dict.items(), key=lambda x: (-x[1], x[0]))[0][0]

    def _is_bad_gene(self, gene_name, chrom=None, pos=None):
        # Check if a gene is "bad" (pseudogene, ncRNA, etc.) using detector's method
        # which handles HGNC symbol resolution with coordinate fallback
        if not gene_name:
            return False
        try:
            return self.detector.is_bad_gene(gene_name, chrom=chrom, pos=pos)
        except Exception:
            return False

    def _replace_bad_gene_if_needed(self, gene, chrom, pos, meta):
        if self._is_bad_gene(gene, chrom, pos):
            replacement = self.detector._find_nearby_protein_coding(chrom, pos, window=500)
            if replacement:
                meta.setdefault('notes', []).append(f"Replaced pseudogene/ncRNA {gene} → {replacement} at consensus")
