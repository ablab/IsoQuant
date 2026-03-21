from collections import defaultdict
import logging

logger = logging.getLogger('IsoQuant')

class FusionMetadata:
    def __init__(self, detector):
        self.detector = detector

    def process_all(self, min_support=1):
        # Process all fusion candidates and generate metadata for valid fusions."""
        logger.info(f"FusionMetadata: Processing {len(self.detector.fusion_candidates)} fusion keys with min_support={min_support}")
        for fusion_key, reads in self.detector.fusion_candidates.items():
            # Validate all filtering gates
            gate_result = self._passes_all_gates(fusion_key, reads, min_support)
            if not gate_result:
                continue
            support, consensus_bp, clustered_support = gate_result
            left_chr, left_pos, right_chr, right_pos = consensus_bp
            # Initialize metadata with consensus information
            meta = self._initialize_metadata(fusion_key, reads, consensus_bp, clustered_support)
            # Assign raw genes from read assignments or consensus coordinates
            raw_left, raw_right = self._assign_raw_genes(fusion_key, meta, left_chr, left_pos, right_chr, right_pos)
            meta["raw_left_gene"] = raw_left
            meta["raw_right_gene"] = raw_right
            # Select final genes based on confidence and bad gene status
            left_gene, right_gene = self._select_final_genes(
                meta, raw_left, raw_right, left_chr, left_pos, right_chr, right_pos
            )
            # Clean up and normalize gene labels
            left_gene = self.detector.normalize_gene_label(left_gene) if left_gene else "intergenic"
            right_gene = self.detector.normalize_gene_label(right_gene) if right_gene else "intergenic"
            # Finalize metadata with biotypes and validation
            self._finalize_metadata(meta, left_gene, right_gene, raw_left, raw_right, left_chr, left_pos, right_chr, right_pos)

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

    def _initialize_metadata(self, fusion_key, reads, consensus_bp, clustered_support):
        # Create and initialize metadata record with consensus information.
        meta = self._ensure_meta(fusion_key)
        meta["supporting_reads"].update(reads)
        meta["consensus_bp"] = consensus_bp
        meta["support"] = clustered_support
        meta["confidence"] = self.detector._compute_early_confidence(clustered_support)
        return meta

    def _assign_raw_genes(self, fusion_key, meta, left_chr, left_pos, right_chr, right_pos):
        # Collect raw gene assignments from read-level assignments.
        assigned = self.detector.fusion_assigned_pairs.get(fusion_key, {})
        raw_left, raw_right = self._collect_raw_assignments(
            meta["supporting_reads"], assigned, left_chr, left_pos, right_chr, right_pos
        )
        # Fall back to coordinate-based assignment if needed
        if raw_left is None:
            raw_left = self.detector.assign_fusion_gene(left_chr, left_pos)
        if raw_right is None:
            raw_right = self.detector.assign_fusion_gene(right_chr, right_pos)
        return raw_left, raw_right

    def _select_final_genes(self, meta, raw_left, raw_right, left_chr, left_pos, right_chr, right_pos):
        #Select final genes based on confidence level with bad-gene override.
        confidence = meta["confidence"]
        raw_left_is_bad = self._is_bad_gene(raw_left, left_chr, left_pos)
        raw_right_is_bad = self._is_bad_gene(raw_right, right_chr, right_pos)
        if confidence > 0.6:
            # High confidence: prefer raw genes, but replace bad ones
            left_gene = self.detector.assign_fusion_gene(left_chr, left_pos) if raw_left_is_bad else raw_left
            right_gene = self.detector.assign_fusion_gene(right_chr, right_pos) if raw_right_is_bad else raw_right
        else:
            # Low confidence: always use consensus assignment
            left_gene = self.detector.assign_fusion_gene(left_chr, left_pos)
            right_gene = self.detector.assign_fusion_gene(right_chr, right_pos)
        # Handle any remaining bad genes from consensus assignment
        left_gene = self._replace_bad_gene_if_needed(left_gene, left_chr, left_pos, meta)
        right_gene = self._replace_bad_gene_if_needed(right_gene, right_chr, right_pos, meta)
        return left_gene, right_gene

    def _finalize_metadata(self, meta, left_gene, right_gene, raw_left, raw_right, left_chr, left_pos, right_chr, right_pos):
        # Store final gene labels, biotypes, and perform final validation.
        meta["left_gene"] = left_gene
        meta["right_gene"] = right_gene
        # Retrieve and store biotypes
        left_biotype = self.detector.get_gene_biotype(raw_left, chrom=left_chr, pos=left_pos)
        right_biotype = self.detector.get_gene_biotype(raw_right, chrom=right_chr, pos=right_pos)
        meta["left_biotype"] = left_biotype
        meta["right_biotype"] = right_biotype
        # Final validation: reject fusions with bad genes
        if self._is_bad_gene(left_gene, left_chr, left_pos) or self._is_bad_gene(right_gene, right_chr, right_pos):
            meta["is_valid"] = False
            meta["reason_invalid"] = f"Left or right gene is pseudogene/ncRNA: {left_gene} / {right_gene}"

    def _ensure_meta(self, fusion_key):
        return self.detector.fusion_metadata.setdefault(
            fusion_key,
            {
                "supporting_reads": set(),
                "consensus_bp": None,
                "left_gene": None,
                "right_gene": None,
                "support": 0,
                "raw_left_gene": None,
                "raw_right_gene": None,
                "confidence": 0.0,
            }
        )

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
                return replacement
        return gene

    def build_exon_boundaries_for_genes(self, gene_names):
        # Build a dict mapping gene names to their exon boundary positions.
        boundaries_dict = {}
        for gene_name in gene_names:
            if not gene_name or gene_name == "intergenic":
                continue
            # Try to find the gene using interval tree (if available)
            g = None
            try:
                # First, try direct DB lookup by ID
                g = self.detector.db[gene_name]
            except Exception:
                # If not found by ID, search by gene_name attribute via interval tree
                if self.detector.interval_index is not None:
                    matching = self.detector.interval_index.find_genes_by_name(gene_name)
                    if matching:
                        g = matching[0]  # Use first match
                else:
                    # Fallback to DB search if no interval tree
                    try:
                        genes = list(self.detector.db.features_of_type('gene'))
                        for gene in genes:
                            if gene.attributes.get('gene_name', [None])[0] == gene_name:
                                g = gene
                                break
                    except Exception:
                        pass
            if g is None:
                continue
            # Get exons and extract boundary positions
            try:
                exons = self.detector._get_cached_exons(g)
                boundaries = []
                for ex_start, ex_end in exons:
                    boundaries.append(ex_start)
                    boundaries.append(ex_end)
                if boundaries:
                    boundaries_dict[gene_name] = boundaries
            except Exception:
                pass
        return boundaries_dict

