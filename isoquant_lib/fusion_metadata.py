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
        # 3. Assign genes and compute biotypes
        self._assign_genes(fusion_key, meta, left_chr, left_pos, right_chr, right_pos)
        # 4. Extract assigned genes
        left_gene = meta.get("left_gene")
        right_gene = meta.get("right_gene")
        # 5. Update fusion_key mappings if genes changed
        original_key = fusion_key
        new_key = self._update_fusion_key_mappings(original_key, left_gene, right_gene, meta)
        # 6. Update metadata with biotypes
        self._compute_biotypes(meta, left_gene, right_gene, left_chr, left_pos, right_chr, right_pos)
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

    def _compute_biotypes(self, meta, left_gene, right_gene, left_chr, left_pos, right_chr, right_pos):
        # Compute biotypes for the assigned genes.
        detector = self.detector
        if left_gene:
            meta["left_biotype"] = detector.get_gene_biotype(
                left_gene, chrom=left_chr, pos=left_pos
            )
        if right_gene:
            meta["right_biotype"] = detector.get_gene_biotype(
                right_gene, chrom=right_chr, pos=right_pos
            )
        logger.debug(
            f"Gene assignments for {meta.get('supporting_reads')}: "
            f"left={left_gene}({meta.get('left_biotype')}), "
            f"right={right_gene}({meta.get('right_biotype')})"
        )

    def _update_fusion_key_mappings(self, original_key, final_left_gene, final_right_gene, meta):
        # Update fusion_key mappings if gene names changed during re-assignment.
        detector = self.detector
        new_left = final_left_gene or "intergenic"
        new_right = final_right_gene or "intergenic"
        # Sort genes alphabetically to match the ordering used in record_fusion
        sorted_genes = sorted([new_left, new_right])
        new_key = f"{sorted_genes[0]}--{sorted_genes[1]}"
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

    def _assign_genes(self, fusion_key, meta, left_chr, left_pos, right_chr, right_pos):
        # Extract genes directly from fusion_key (format: GENE1--GENE2, alphabetically sorted).
        # This avoids majority voting corruption where reads can assign genes to different sides.
        parts = fusion_key.split("--")
        if len(parts) == 2:
            left_gene, right_gene = parts[0], parts[1]
        else:
            # Fallback to empty if key format is unexpected
            left_gene, right_gene = None, None
        # Populate metadata with genes
        meta["left_gene"] = left_gene
        meta["right_gene"] = right_gene
        # Compute average scores from per-read scores
        read_scores = self.detector.fusion_read_scores.get(fusion_key, {})
        left_scores = []
        right_scores = []
        for read_name in meta["supporting_reads"]:
            if read_name in read_scores:
                left_score, right_score = read_scores[read_name]
                left_scores.append(left_score)
                right_scores.append(right_score)
        # Average the scores, or fallback to 0.0 if no scores available
        meta["left_score"] = sum(left_scores) / len(left_scores) if left_scores else 0.0
        meta["right_score"] = sum(right_scores) / len(right_scores) if right_scores else 0.0

