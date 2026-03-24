import logging
import re
from collections import defaultdict

logger = logging.getLogger('IsoQuant')

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

    def _map_pseudogene_to_parent(self, gene, chrom=None, pos=None):
        # Attempt to map a pseudogene to its parent gene.
        # Pattern: Remove suffix like P1, P2, etc. or other numeric pseudogene suffixes.
        # Examples: RPL23P6 → RPL23, SAE1P1 → SAE1, ZFYVE9P1 → ZFYVE9
        if not gene or not isinstance(gene, str):
            return None
        # Match pattern: ends with P followed by one or more digits (e.g. P1, P6, P123)
        # Also handles patterns like P followed by other suffixes
        match = re.match(r"^(.+?)P\d+$", gene)
        if not match:
            return None
        parent_candidate = match.group(1)
        # Verify parent gene exists and is protein-coding
        parent_biotype = self.detector.get_gene_biotype(parent_candidate, chrom=chrom, pos=pos)
        if parent_biotype == "protein_coding":
            logger.debug(f"Mapped pseudogene {gene} → parent gene {parent_candidate}")
            return parent_candidate
        return None

    def _salvage_and_check_gene_pair(self, left_gene, right_gene, left_chr, left_pos, right_chr, right_pos):
        # Salvage non-coding genes by finding nearby protein-coding alternatives.
        def _salvage_to_nearby_coding(gene, chrom, pos):
            # Only attempt rescue if we have coordinates and the gene is non-coding/unknown
            if not gene or chrom is None or pos is None:
                return gene, None
            biotype = self.detector.get_gene_biotype(gene, chrom=chrom, pos=pos)
            if biotype and biotype == "protein_coding":
                return gene, biotype
            # Try to map pseudogene to parent gene first
            parent_gene = self._map_pseudogene_to_parent(gene, chrom=chrom, pos=pos)
            if parent_gene:
                return parent_gene, "protein_coding"
            # Try a small window to avoid spurious swaps
            nearby = self.detector._find_nearby_protein_coding(chrom, pos, window=500)
            if nearby:
                return nearby, "protein_coding"
            return gene, biotype
        left_gene, left_biotype = _salvage_to_nearby_coding(left_gene, left_chr, left_pos)
        right_gene, right_biotype = _salvage_to_nearby_coding(right_gene, right_chr, right_pos)
        logger.debug(f"Gene pair after salvage: {left_gene}={left_biotype}, {right_gene}={right_biotype}")
        # Check if either is still non-protein-coding
        has_non_coding = False
        reason = None
        if left_gene and left_biotype and left_biotype != "protein_coding":
            has_non_coding = True
            reason = f"{left_gene}={left_biotype}"
        elif right_gene and right_biotype and right_biotype != "protein_coding":
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
                # Salvage genes and check biotypes
                left_gene, left_biotype, right_gene, right_biotype, is_non_coding, reason = (
                    self._salvage_and_check_gene_pair(left_gene, right_gene, left_chr, left_pos, right_chr, right_pos)
                )
                # Update assignment with salvaged genes
                if self.detector.fusion_assigned_pairs.get(fusion_key, {}).get(read_name):
                    self.detector.fusion_assigned_pairs[fusion_key][read_name] = (left_gene, right_gene)
                if is_non_coding:
                    has_non_coding = True
                    non_coding_reason = reason
                    break
            if has_non_coding:
                logger.info(f"Dropping early non-coding fusion: {fusion_key} - {non_coding_reason}")
                fusions_to_discard.add(fusion_key)
        # Remove all discarded fusions from data structures
        self._remove_discarded_fusions_internal(fusions_to_discard)
        
        # Consolidate fusions with pseudogene/parent gene variants
        self._consolidate_pseudogene_variants()

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

    def _consolidate_pseudogene_variants(self):
        # Merge fusions that differ only by pseudogene/parent gene variants.
        # E.g., merge "RPL23P6--SOMETHING" into "RPL23--SOMETHING"
        def canonicalize_gene_in_key(gene):
            # Try to map pseudogene to parent
            parent = self._map_pseudogene_to_parent(gene)
            return parent if parent else gene
        # Build a canonical key mapping
        canonical_keys = {}
        for fusion_key in list(self.detector.fusion_candidates.keys()):
            parts = fusion_key.split("--")
            if len(parts) != 2:
                continue
            left_gene, right_gene = parts
            canonical_left = canonicalize_gene_in_key(left_gene)
            canonical_right = canonicalize_gene_in_key(right_gene)
            canonical_key = f"{canonical_left}--{canonical_right}"
            canonical_keys[fusion_key] = canonical_key
        # Merge fusions with their canonical keys
        fusions_to_merge = {}
        for fusion_key, canonical_key in canonical_keys.items():
            if canonical_key != fusion_key:
                if canonical_key not in fusions_to_merge:
                    fusions_to_merge[canonical_key] = []
                fusions_to_merge[canonical_key].append(fusion_key)
        # Execute merges: for each canonical key, merge all variant keys into it
        for canonical_key, variant_keys in fusions_to_merge.items():
            # Ensure the canonical key exists
            if canonical_key not in self.detector.fusion_candidates:
                self.detector.fusion_candidates[canonical_key] = set()
                self.detector.fusion_breakpoints[canonical_key] = defaultdict(int)
                self.detector.fusion_assigned_pairs[canonical_key] = {}
            # Merge all variants into the canonical key
            for variant_key in variant_keys:
                # Merge fusion_candidates
                if variant_key in self.detector.fusion_candidates:
                    self.detector.fusion_candidates[canonical_key].update(
                        self.detector.fusion_candidates[variant_key]
                    )
                    del self.detector.fusion_candidates[variant_key]
                # Merge fusion_breakpoints
                if variant_key in self.detector.fusion_breakpoints:
                    for bp, count in self.detector.fusion_breakpoints[variant_key].items():
                        self.detector.fusion_breakpoints[canonical_key][bp] += count
                    del self.detector.fusion_breakpoints[variant_key]
                # Merge fusion_assigned_pairs
                if variant_key in self.detector.fusion_assigned_pairs:
                    self.detector.fusion_assigned_pairs[canonical_key].update(
                        self.detector.fusion_assigned_pairs[variant_key]
                    )
                    del self.detector.fusion_assigned_pairs[variant_key]
                # Merge fusion_metadata if it exists
                if variant_key in self.detector.fusion_metadata:
                    if canonical_key not in self.detector.fusion_metadata:
                        self.detector.fusion_metadata[canonical_key] = self.detector.fusion_metadata[variant_key]
                    else:
                        # Merge supporting reads
                        if "supporting_reads" in self.detector.fusion_metadata[variant_key]:
                            if "supporting_reads" not in self.detector.fusion_metadata[canonical_key]:
                                self.detector.fusion_metadata[canonical_key]["supporting_reads"] = set()
                            self.detector.fusion_metadata[canonical_key]["supporting_reads"].update(
                                self.detector.fusion_metadata[variant_key]["supporting_reads"]
                            )
                    del self.detector.fusion_metadata[variant_key]
                logger.info(f"Consolidated pseudogene variant: '{variant_key}' → '{canonical_key}'")


    def filter_non_coding_genes(self):
        # Removes fusions with non-protein-coding partners completely from all data structures.
        # NOTE: Do NOT discard fusions where biotype is None (unknown) - only discard if explicitly non-coding
        fusions_to_discard = []
        discard_reasons = {}  # Track reason for each discarded fusion
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
            # Only discard if BOTH conditions are met:
            # 1. Biotype is not None (i.e., we actually found it)
            # 2. Biotype is not "protein_coding"
            # This way we KEEP fusions where biotype is unknown (None) or protein_coding
            left_is_non_coding = left_biotype is not None and left_biotype != "protein_coding"
            right_is_non_coding = right_biotype is not None and right_biotype != "protein_coding"
            
            if left_is_non_coding or right_is_non_coding:
                reason = []
                if left_is_non_coding:
                    reason.append(f"Left {left_gene}={left_biotype}")
                if right_is_non_coding:
                    reason.append(f"Right {right_gene}={right_biotype}")
                fusions_to_discard.append(fusion_key)
                discard_reasons[fusion_key] = "; ".join(reason)
        # Log discarded fusions with details
        if fusions_to_discard:
            logger.info(f"filter_non_coding_genes: Discarding {len(fusions_to_discard)} fusion(s)")
            for fusion_key in fusions_to_discard:
                logger.info(f"  Discard: {fusion_key} - {discard_reasons[fusion_key]}")
        else:
            logger.info("filter_non_coding_genes: No fusions discarded (all have protein-coding partners or unknown biotypes)")
        # Remove all discarded fusions from data structures
        self._remove_discarded_fusions_internal(fusions_to_discard)

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

    def calculate_exon_boundary_bonus(self, meta, candidate):
        # Calculate exon boundary proximity bonus for confidence scoring.
        bL = meta.get("left_min_exon_boundary_delta")
        bR = meta.get("right_min_exon_boundary_delta")
        boundary_bonus = 0.0
        if isinstance(bL, int):
            boundary_bonus += max(0.0, min(0.20, (200 - min(bL, 200)) / 200.0 * 0.20))
        if isinstance(bR, int):
            boundary_bonus += max(0.0, min(0.20, (200 - min(bR, 200)) / 200.0 * 0.20))
        boundary_bonus = min(boundary_bonus, 0.35)
        return boundary_bonus

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
        # Exon boundary bonus
        boundary_bonus = self.calculate_exon_boundary_bonus(meta, None)
        # Driver gene bonus vs artifact penalty
        priors = 0.0
        left = meta.get("left_gene")
        right = meta.get("right_gene")
        if left and right:
            if self.is_multicopy_artifact_family(left) or self.is_multicopy_artifact_family(right):
                priors -= 0.30
            if self.is_driver_gene(left) or self.is_driver_gene(right):
                priors += 0.20
        conf = max(0.0, min(1.0, support + recon + realign + boundary_bonus + priors))
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

    def _merge_similar_candidates_by_gene(self, shared_gene_index, index_type="left",
                                           distance_threshold=10000):
        # Merge similar fusion candidates sharing a gene on one side.
        fusions_to_discard = set()
        for shared_gene, entries in shared_gene_index.items():
            if len(entries) < 2:
                continue
            for i in range(len(entries)):
                for j in range(i + 1, len(entries)):
                    if fusions_to_discard & {entries[i][0], entries[j][0]}:
                        continue  # Skip if either fusion already marked for discard
                    fusion_key_i, other_gene_i, other_chr_i, other_pos_i = entries[i]
                    fusion_key_j, other_gene_j, other_chr_j, other_pos_j = entries[j]
                    # Same chromosome and within distance threshold?
                    if other_chr_i != other_chr_j:
                        continue
                    distance = abs(other_pos_i - other_pos_j)
                    if distance > distance_threshold:
                        continue
                    # If genes are identical and very close (< 1kb), definitely merge as duplicates
                    if other_gene_i == other_gene_j and distance < 1000:
                        logger.info(f"Merging near-identical duplicates: {fusion_key_i} ← {fusion_key_j} "
                                   f"(shared {index_type}={shared_gene}, distance={distance}bp)")
                        self._merge_fusion_candidates(fusion_key_i, fusion_key_j)
                        fusions_to_discard.add(fusion_key_j)
                        continue
                    # Skip if genes are identical and far apart (different events)
                    if other_gene_i == other_gene_j:
                        continue
                    # Score both genes at their respective positions
                    side = "right" if index_type == "left" else "left"
                    logger.debug(f"Merging candidates: {fusion_key_i} vs {fusion_key_j} "
                               f"(shared {index_type}={shared_gene}, {side} distance={distance}bp)")
                    score_i = self.detector._compute_gene_score(other_gene_i, other_pos_i, chrom=other_chr_i)
                    score_j = self.detector._compute_gene_score(other_gene_j, other_pos_j, chrom=other_chr_j)
                    # Keep the one with higher score, discard the other
                    if score_i >= score_j:
                        self._merge_fusion_candidates(fusion_key_i, fusion_key_j)
                        fusions_to_discard.add(fusion_key_j)
                    else:
                        self._merge_fusion_candidates(fusion_key_j, fusion_key_i)
                        fusions_to_discard.add(fusion_key_i)
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

    def merge_nearly_identical(self, distance_threshold=10000):
        # Merge fully identical and nearly identical fusion candidates.
        # Step 1: Merge fully identical fusions
        fully_identical_discards = self._merge_fully_identical()
        # Step 2: Build indexes for nearly identical fusion detection
        left_gene_index = defaultdict(list)
        right_gene_index = defaultdict(list)
        for fusion_key, meta in self.detector.fusion_metadata.items():
            if fusion_key in fully_identical_discards:
                continue  # Skip already-discarded fusions
            consensus_bp = meta.get("consensus_bp")
            if not consensus_bp:
                continue
            # Handle format (c1, p1, c2, p2)
            left_chr, left_pos, right_chr, right_pos = consensus_bp
            left_gene = meta.get("left_gene")
            right_gene = meta.get("right_gene")
            if not left_gene or not right_gene:
                continue
            # Index by left_gene: store (fusion_key, right_gene, right_chr, right_pos)
            left_gene_index[left_gene].append((fusion_key, right_gene, right_chr, right_pos))
            # Index by right_gene: store (fusion_key, left_gene, left_chr, left_pos)
            right_gene_index[right_gene].append((fusion_key, left_gene, left_chr, left_pos))
        # Step 3: Merge nearly identical by shared left gene
        left_gene_discards = self._merge_similar_candidates_by_gene(
            left_gene_index, index_type="left", distance_threshold=distance_threshold
        )
        # Step 4: Merge nearly identical by shared right gene
        right_gene_discards = self._merge_similar_candidates_by_gene(
            right_gene_index, index_type="right", distance_threshold=distance_threshold
        )
        # Step 5: Remove all discarded fusions
        all_discards = fully_identical_discards | left_gene_discards | right_gene_discards
        self._remove_discarded_fusions(all_discards)

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
