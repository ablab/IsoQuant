import logging
from collections import defaultdict

logger = logging.getLogger('IsoQuant')
import re

class FusionValidator:

    def __init__(self, detector):
        self.detector = detector

    def _get_breakpoint_coords(self, fusion_key):
        # Extract chromosome and position coordinates from the first breakpoint of a fusion.
        left_chr, left_pos, right_chr, right_pos = None, None, None, None
        bp_counts = self.detector.fusion_breakpoints.get(fusion_key, {})
        if bp_counts:
            first_bp = next(iter(bp_counts.keys()), None)
            if first_bp and len(first_bp) == 4:
                left_chr, left_pos, right_chr, right_pos = first_bp
                logger.debug(f"Got breakpoint for {fusion_key}: {left_chr}:{left_pos} - {right_chr}:{right_pos}")
        return left_chr, left_pos, right_chr, right_pos

    def _is_allowed_biotype(self, biotype):
        # Check if a biotype is allowed for fusion detection.
        # Includes protein-coding genes and immune receptor genes (IG and TR).
        allowed_biotypes = {
            "protein_coding",
            "IG_C_gene",
            "IG_V_gene",
            "IG_D_gene",
            "IG_J_gene",
            "TR_V_gene",
            "TR_D_gene",
            "TR_J_gene",
            "TR_C_gene",
            "GBA3"
        }
        return biotype in allowed_biotypes

    def _salvage_and_check_gene_pair(self, left_gene, right_gene, left_chr, left_pos, right_chr, right_pos):
        # Check biotypes for both genes against the allowed whitelist.
        left_biotype = self.detector.get_gene_biotype(left_gene, chrom=left_chr, pos=left_pos)
        right_biotype = self.detector.get_gene_biotype(right_gene, chrom=right_chr, pos=right_pos)
        logger.debug(f"Gene pair biotypes: {left_gene}={left_biotype}, {right_gene}={right_biotype}")
        # Check if either biotype is not in the allowed whitelist
        has_non_coding = False
        reason = None
        if left_gene and left_biotype and not self._is_allowed_biotype(left_biotype):
            has_non_coding = True
            reason = f"{left_gene}={left_biotype}"
        elif right_gene and right_biotype and not self._is_allowed_biotype(right_biotype):
            has_non_coding = True
            reason = f"{right_gene}={right_biotype}"
        return left_gene, left_biotype, right_gene, right_biotype, has_non_coding, reason

    def filter_raw_non_coding_genes(self):
        # Drop fusions with non-protein-coding genes at the raw assignment stage.
        fusions_to_discard = set()
        logger.info(f"filter_raw_non_coding_genes: Processing {len(self.detector.fusion_assigned_pairs)} fusion keys")
        for fusion_key, read_assignments in self.detector.fusion_assigned_pairs.items():
            if fusion_key not in self.detector.fusion_candidates:
                continue
            # Extract breakpoint coordinates
            left_chr, left_pos, right_chr, right_pos = self._get_breakpoint_coords(fusion_key)
            # Check all raw assigned gene pairs for this fusion
            has_non_coding = False
            non_coding_reason = None
            for read_name, (left_gene, right_gene) in read_assignments.items():
                # Check biotypes for both genes
                left_gene, left_biotype, right_gene, right_biotype, is_non_coding, reason = (
                    self._salvage_and_check_gene_pair(left_gene, right_gene, left_chr, left_pos, right_chr, right_pos)
                )
                # Update assignment with salvaged genes
                if self.detector.fusion_assigned_pairs.get(fusion_key, {}).get(read_name):
                    self.detector.fusion_assigned_pairs[fusion_key][read_name] = (left_gene, right_gene)
                if is_non_coding:
                    has_non_coding = True
                    reasons = []
                    if not self._is_allowed_biotype(left_biotype):
                        reasons.append(f"Left {left_gene}={left_biotype}")
                    if not self._is_allowed_biotype(right_biotype):
                        reasons.append(f"Right {right_gene}={right_biotype}")
                    non_coding_reason = "; ".join(reasons)
                    break
            if has_non_coding:
                logger.info(f"Dropping early non-coding fusion: {fusion_key} - {non_coding_reason}")
                fusions_to_discard.add(fusion_key)
        # Remove all discarded fusions from data structures
        self._remove_discarded_fusions_internal(fusions_to_discard)

    def _remove_discarded_fusions_internal(self, fusions_to_discard):
        # Remove fusions from all internal data structures.
        for fusion_key in fusions_to_discard:
            if fusion_key in self.detector.fusion_metadata:
                del self.detector.fusion_metadata[fusion_key]
            if fusion_key in self.detector.fusion_candidates:
                del self.detector.fusion_candidates[fusion_key]
            if fusion_key in self.detector.fusion_breakpoints:
                del self.detector.fusion_breakpoints[fusion_key]
            if fusion_key in self.detector.fusion_assigned_pairs:
                del self.detector.fusion_assigned_pairs[fusion_key]

    def filter_non_coding_genes(self):
        # Removes fusions with partners not in the allowed biotype whitelist from all data structures.
        fusions_to_discard = []
        for fusion_key, meta in list(self.detector.fusion_metadata.items()):
            left_gene = meta.get("left_gene")
            right_gene = meta.get("right_gene")
            left_chr, left_pos = None, None
            right_chr, right_pos = None, None
            consensus_bp = meta.get("consensus_bp")
            if consensus_bp:
                left_chr, left_pos, right_chr, right_pos = consensus_bp
            # Get biotypes for both genes
            left_biotype = self.detector.get_gene_biotype(left_gene, chrom=left_chr, pos=left_pos)
            right_biotype = self.detector.get_gene_biotype(right_gene, chrom=right_chr, pos=right_pos)
            # Discard if either partner has a biotype not in the allowed whitelist
            left_is_allowed = self._is_allowed_biotype(left_biotype)
            right_is_allowed = self._is_allowed_biotype(right_biotype)
            if not left_is_allowed or not right_is_allowed:
                reason = []
                if not left_is_allowed:
                    reason.append(f"Left {left_gene}={left_biotype}")
                if not right_is_allowed:
                    reason.append(f"Right {right_gene}={right_biotype}")
                # logger.info(f"Discarding fusion with non-allowed biotype: {fusion_key} - {'; '.join(reason)}")
                fusions_to_discard.append(fusion_key)
        # Remove all discarded fusions from data structures
        self._remove_discarded_fusions_internal(fusions_to_discard)

    def filter_multicopy_artifact_pairs(self):
        # Filter out fusions where both partners belong to multicopy artifact families.
        # These are biologically non-relevant fusions like TRAJ17-TRAV1-2, ZNF124-ZNF670, H2AC13-H2BC13.
        if not self.detector.fusion_metadata:
            return
        for fusion_key, meta in self.detector.fusion_metadata.items():
            if meta.get("is_valid") == False:
                continue
            left_gene = meta.get("left_gene")
            right_gene = meta.get("right_gene")
            # Check if both genes belong to multicopy artifact families
            if (left_gene and right_gene and 
                self.is_multicopy_artifact_family(left_gene) and 
                self.is_multicopy_artifact_family(right_gene)) or (left_gene and right_gene and 
                self._is_ribosomal_or_histone_gene(left_gene) and 
                self._is_ribosomal_or_histone_gene(right_gene)):
                meta["is_valid"] = False
                if "reasons" not in meta:
                    meta["reasons"] = []
                meta["reasons"].append(
                    f"Both genes belong to multicopy artifact families: {left_gene} - {right_gene}"
                )
                logger.info(f"Filtering multicopy artifact pair fusion: {fusion_key} ({left_gene} - {right_gene})")

    def apply_frequency_filters(self):
        # Filter out multicopy artifacts based on gene frequency within the sample.
        if not self.detector.fusion_metadata:
            return
        # Count gene frequencies across all fusions
        left_gene_count = defaultdict(int)
        right_gene_count = defaultdict(int)
        for fusion_key, meta in self.detector.fusion_metadata.items():
            left_gene = meta.get("left_gene")
            right_gene = meta.get("right_gene")
            if left_gene and left_gene != "intergenic":
                left_gene_count[left_gene] += 1
            if right_gene and right_gene != "intergenic":
                right_gene_count[right_gene] += 1
        # Mark fusions with high-frequency genes as artifacts
        for fusion_key, meta in self.detector.fusion_metadata.items():
            if meta.get("is_valid") == False:
                continue
            left_gene = meta.get("left_gene")
            right_gene = meta.get("right_gene")
            artifact_reasons = []
            # Check left gene frequency
            if left_gene and left_gene != "intergenic":
                is_rib_hist = self._is_ribosomal_or_histone_gene(left_gene)
                # lower threshold for ribosomal protein artifacts
                threshold = 2 if is_rib_hist else 4
                count = left_gene_count.get(left_gene, 0)
                if count > threshold:
                    artifact_reasons.append(
                        f"Left gene '{left_gene}' appears {count} times "
                        f"(threshold={threshold} for {'ribosomal/histone' if is_rib_hist else 'other'} genes)"
                    )
            # Check right gene frequency
            if right_gene and right_gene != "intergenic":
                is_rib_hist = self._is_ribosomal_or_histone_gene(right_gene)
                threshold = 2 if is_rib_hist else 4
                count = right_gene_count.get(right_gene, 0)
                if count > threshold:
                    artifact_reasons.append(
                        f"Right gene '{right_gene}' appears {count} times "
                        f"(threshold={threshold} for {'ribosomal/histone' if is_rib_hist else 'other'} genes)"
                    )
            # Mark as artifact if any gene exceeds threshold
            if artifact_reasons:
                meta["is_valid"] = False
                if "reasons" not in meta:
                    meta["reasons"] = []
                meta["reasons"].extend(artifact_reasons)

    def _is_ribosomal_or_histone_gene(self, gene_name):
        # Check if a gene is ribosomal or histone-related based on name patterns.
        if not gene_name or not isinstance(gene_name, str):
            return False
        gene_upper = gene_name.upper()
        rib_hist_patterns = ('RPL', 'RPS', 'RPLP', 'HIST', 'H1', 'H2A', 'H2B', 'H3', 'H4')
        return any(gene_upper.startswith(p) for p in rib_hist_patterns)

    def is_driver_gene(self, gene):
        driver_genes = {
            "BCR", "ABL1", "RUNX1", "RUNX1T1", "PML", "RARA", "ETV6",
            "NUP98", "NUP214", "KMT2A", "MLLT3", "EWSR1", "ERG",
            "FGFR1", "RET", "ROS1", "ALK", "TMPRSS2", "FLI1",
            "MYH11", "CBFB"
        }
        return gene in driver_genes

    def is_multicopy_artifact_family(self, gene):
        # Check if a gene belongs to a multicopy artifact family (ribosomal, histone, etc.)
        prefixes = ("RPL", "RPS", "MRPL", "MRPS", "H1-", "H2A", "H2B", "H3-", "H4-",
                    "HIST1", "HIST2", "HIST3", "HLA-", "MICA", "MICB", "KRT", "KRT",
                    "OR", "ZNF", "DEFA", "DEFB", "IGH", "IGK", "IGL", "TRAV", "TRBV", "TRGV", "TRDV",
                    "MUC", "AMY", "CYP")
        g = gene.upper()
        return g.startswith(prefixes)

    def confidence(self, meta, flags=None):
        # Compute full confidence score including reconstruction/realignment bonuses.
        # Start from clustered support normalized
        support = min(meta.get("support", 0), 10) / 10.0
        # Reconstruction bonus
        recon = 0.2 if meta.get("reconstruction_ok") else 0.0
        # Realignment bonus: scale by number of hits and MAPQ
        hits = meta.get("realignment_hits", 0)
        mapq = meta.get("best_hit_mapq", 0)
        realign = min(hits, 3) * 0.1 + min(mapq, 30) / 300.0
        # Driver gene bonus vs artifact penalty
        priors = 0.0
        left = meta.get("left_gene")
        right = meta.get("right_gene")
        if left and right:
            if self.is_multicopy_artifact_family(left) or self.is_multicopy_artifact_family(right):
                priors -= 0.30
            if self.is_driver_gene(left) or self.is_driver_gene(right):
                priors += 0.20
        conf = max(0.0, min(1.0, support + recon + realign + priors))
        return conf

    def _merge_fully_identical(self):
        # Detect and merge fully identical fusions (same genes AND same consensus breakpoints).
        # This handles exact duplicates that result from processing.
        bp_to_fusion_keys = defaultdict(list)
        for fusion_key, meta in self.detector.fusion_metadata.items():
            consensus_bp = meta.get("consensus_bp")
            left_gene = meta.get("left_gene")
            right_gene = meta.get("right_gene")
            if not consensus_bp or not left_gene or not right_gene:
                continue
            # Create a signature: (left_gene, right_gene, consensus_bp)
            signature = (left_gene, right_gene, consensus_bp)
            bp_to_fusion_keys[signature].append(fusion_key)
        fusions_to_discard = set()
        for signature, fusion_keys in bp_to_fusion_keys.items():
            if len(fusion_keys) < 2:
                continue
            # Multiple identical fusions found; keep first, merge rest into it
            keep_key = fusion_keys[0]
            for discard_key in fusion_keys[1:]:
                logger.info(f"Merging identical duplicates: {keep_key} ← {discard_key}")
                self._merge_fusion_candidates(keep_key, discard_key)
                fusions_to_discard.add(discard_key)
        return fusions_to_discard

    def _remove_discarded_fusions(self, fusions_to_discard):
        # Remove discarded fusions from all internal data structures.
        for fusion_key in fusions_to_discard:
            logger.info(f"Discarding fusion candidate: {fusion_key}")
            if fusion_key in self.detector.fusion_metadata:
                del self.detector.fusion_metadata[fusion_key]
            if fusion_key in self.detector.fusion_candidates:
                del self.detector.fusion_candidates[fusion_key]
            if fusion_key in self.detector.fusion_breakpoints:
                del self.detector.fusion_breakpoints[fusion_key]

    def _merge_fusion_candidates(self, keep_fusion_key, discard_fusion_key):
        # Merge two fusion candidates: add supporting reads from discard_fusion_key to keep_fusion_key
        if discard_fusion_key in self.detector.fusion_candidates:
            keep_reads = self.detector.fusion_candidates.get(keep_fusion_key, set())
            discard_reads = self.detector.fusion_candidates[discard_fusion_key]
            keep_reads.update(discard_reads)
            self.detector.fusion_candidates[keep_fusion_key] = keep_reads
        # Merge breakpoint counts
        if discard_fusion_key in self.detector.fusion_breakpoints:
            keep_bps = self.detector.fusion_breakpoints.get(keep_fusion_key, {})
            discard_bps = self.detector.fusion_breakpoints[discard_fusion_key]
            for bp, count in discard_bps.items():
                keep_bps[bp] = keep_bps.get(bp, 0) + count
            self.detector.fusion_breakpoints[keep_fusion_key] = keep_bps
        # Update metadata with merged read count
        if keep_fusion_key in self.detector.fusion_metadata:
            meta = self.detector.fusion_metadata[keep_fusion_key]
            if discard_fusion_key in self.detector.fusion_metadata:
                discard_meta = self.detector.fusion_metadata[discard_fusion_key]
                if "supporting_reads" in meta and "supporting_reads" in discard_meta:
                    meta["supporting_reads"].update(discard_meta["supporting_reads"])
                    meta["support"] = len(meta["supporting_reads"])
