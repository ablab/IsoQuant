from collections import defaultdict

class FusionMetadata:
    def __init__(self, detector):
        self.detector = detector

    def process_all(self, min_support=1):
        for fusion_key, reads in self.detector.fusion_candidates.items():
            meta = self._ensure_meta(fusion_key)
            meta["supporting_reads"].update(reads)
            meta["support"] = len(meta["supporting_reads"])
            if meta["support"] < min_support:
                continue
            bp_counts = self.detector.fusion_breakpoints.get(fusion_key, {})
            if not bp_counts:
                continue
            consensus_result = self.detector.cluster_breakpoints(bp_counts, window=2000)
            if not consensus_result:
                continue
            consensus_bp, clustered_support = consensus_result
            meta["consensus_bp"] = consensus_bp
            meta["support"] = clustered_support
            
            # Handle new consensus_bp format with strands: (c1, p1, s1, c2, p2, s2)
            # or old format without strands: (c1, p1, c2, p2) for backward compatibility
            if len(consensus_bp) == 6:
                left_chr, left_pos, left_strand, right_chr, right_pos, right_strand = consensus_bp
            else:
                left_chr, left_pos, right_chr, right_pos = consensus_bp
                left_strand, right_strand = "+", "+"
            
            meta["read_strand1"] = left_strand
            meta["read_strand2"] = right_strand
            early_confidence = self.detector._compute_early_confidence(meta["support"])
            meta["confidence"] = early_confidence

            assigned = self.detector.fusion_assigned_pairs.get(fusion_key, {})
            raw_left, raw_right = self._collect_raw_assignments(meta["supporting_reads"], assigned, left_chr, left_pos, right_chr, right_pos)
            if raw_left is None:
                raw_left = self.detector.assign_fusion_gene(left_chr, left_pos)
            if raw_right is None:
                raw_right = self.detector.assign_fusion_gene(right_chr, right_pos)
            meta["raw_left_gene"] = raw_left
            meta["raw_right_gene"] = raw_right

            # Check if raw genes are "bad" (pseudogene, ncRNA, etc.) BEFORE confidence decision
            # If a raw gene is bad, skip it in favor of consensus assignment, even if confidence > 0.6
            raw_left_is_bad = self._is_bad_gene(raw_left, left_chr, left_pos)
            raw_right_is_bad = self._is_bad_gene(raw_right, right_chr, right_pos)

            # Confidence-based gene selection WITH bad-gene override
            if early_confidence > 0.6:
                # High confidence: prefer raw genes, BUT replace if they are bad
                left_gene = self.detector.assign_fusion_gene(left_chr, left_pos) if raw_left_is_bad else raw_left
                right_gene = self.detector.assign_fusion_gene(right_chr, right_pos) if raw_right_is_bad else raw_right
            else:
                # Low confidence: always use consensus assignment
                left_gene = self.detector.assign_fusion_gene(left_chr, left_pos)
                right_gene = self.detector.assign_fusion_gene(right_chr, right_pos)

            # Final normalization: handle any remaining bad genes from consensus assignment
            left_gene = self._replace_bad_gene_if_needed(left_gene, left_chr, left_pos, meta)
            right_gene = self._replace_bad_gene_if_needed(right_gene, right_chr, right_pos, meta)

            left_gene = self.detector.normalize_gene_label(left_gene) if left_gene else "intergenic"
            right_gene = self.detector.normalize_gene_label(right_gene) if right_gene else "intergenic"

            # Store biotypes and strands for reporting, passing coordinates for symbol resolution
            left_biotype = self.detector.get_gene_biotype(left_gene, chrom=left_chr, pos=left_pos)
            right_biotype = self.detector.get_gene_biotype(right_gene, chrom=right_chr, pos=right_pos)
            left_gene_strand = self.detector.get_gene_strand(left_gene, chrom=left_chr, pos=left_pos)
            right_gene_strand = self.detector.get_gene_strand(right_gene, chrom=right_chr, pos=right_pos)
            
            meta["left_biotype"] = left_biotype
            meta["right_biotype"] = right_biotype
            meta["left_gene_strand"] = left_gene_strand
            meta["right_gene_strand"] = right_gene_strand
            
            # Determine functional orientation and reorder if necessary
            # to ensure genes are in 5' to 3' direction relative to the fusion
            final_left_gene, final_right_gene, orientation_flipped = self._determine_functional_order(
                left_gene, right_gene, left_strand, right_strand, left_gene_strand, right_gene_strand
            )
            
            if orientation_flipped:
                # Swap all coordinate and strand pairs
                final_left_chr, final_left_pos, final_right_chr, final_right_pos = right_chr, right_pos, left_chr, left_pos
                final_left_strand, final_right_strand = right_strand, left_strand
                final_left_biotype, final_right_biotype = right_biotype, left_biotype
                final_left_gene_strand, final_right_gene_strand = right_gene_strand, left_gene_strand
                meta["consensus_bp"] = (final_left_chr, final_left_pos, final_left_strand, final_right_chr, final_right_pos, final_right_strand)
            else:
                final_left_chr, final_left_pos, final_right_chr, final_right_pos = left_chr, left_pos, right_chr, right_pos
                final_left_strand, final_right_strand = left_strand, right_strand
                final_left_biotype, final_right_biotype = left_biotype, right_biotype
                final_left_gene_strand, final_right_gene_strand = left_gene_strand, right_gene_strand
            
            meta["left_gene"] = final_left_gene
            meta["right_gene"] = final_right_gene
            meta["left_biotype"] = final_left_biotype
            meta["right_biotype"] = final_right_biotype
            meta["left_gene_strand"] = final_left_gene_strand
            meta["right_gene_strand"] = final_right_gene_strand
            meta["orientation_flipped"] = orientation_flipped

            # FINAL VALIDATION: If either final gene is still bad (pseudogene, ncRNA, etc.),
            # mark the entire fusion as invalid to prevent false-alt duplicates in output
            if self._is_bad_gene(final_left_gene, final_left_chr, final_left_pos) or self._is_bad_gene(final_right_gene, final_right_chr, final_right_pos):
                meta["is_valid"] = False
                meta["reason_invalid"] = f"Left or right gene is pseudogene/ncRNA: {final_left_gene} / {final_right_gene}"

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

    def _determine_functional_order(self, left_gene, right_gene, left_read_strand, right_read_strand, 
                                    left_gene_strand, right_gene_strand):
        # Determine if genes should be reordered to maintain 5' to 3' direction.

        # Default: no reordering
        if left_gene == "intergenic" or right_gene == "intergenic":
            return left_gene, right_gene, False
        # Primary decision based on gene strand annotations
        # For a canonical fusion pattern: left_gene on + strand, right_gene on - strand
        # This represents 5' → 3' direction conceptually
        if left_gene_strand == "+" and right_gene_strand == "-":
            # Classic canonical pattern: keep as is
            return left_gene, right_gene, False
        elif left_gene_strand == "-" and right_gene_strand == "+":
            # Reverse of canonical pattern: swap to get canonical orientation
            return right_gene, left_gene, True
        elif left_gene_strand == right_gene_strand and left_gene_strand is not None:
            # Both genes on same strand - use read strand patterns as tie-breaker
            # If reads also show opposite strands, that's consistent with the gene strands
            # suggesting that one read maps to + and other to - (which is expected in alignment)
            # In this case, keep the order since it's consistent
            return left_gene, right_gene, False
        
        # Default: keep current order if we don't have enough info
        return left_gene, right_gene, False

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

