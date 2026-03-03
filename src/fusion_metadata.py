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
            left_chr, left_pos, right_chr, right_pos = consensus_bp
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
            meta["left_gene"] = left_gene
            meta["right_gene"] = right_gene

            # Store biotypes for reporting, passing coordinates for symbol resolution
            meta["left_biotype"] = self.detector.get_gene_biotype(left_gene, chrom=left_chr, pos=left_pos)
            meta["right_biotype"] = self.detector.get_gene_biotype(right_gene, chrom=right_chr, pos=right_pos)

            # FINAL VALIDATION: If either final gene is still bad (pseudogene, ncRNA, etc.),
            # mark the entire fusion as invalid to prevent false-alt duplicates in output
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
            # Try to find the gene
            try:
                g = self.detector.db[gene_name]
            except Exception:
                # Try to find by gene_name attribute
                try:
                    genes = list(self.detector.db.features_of_type('gene'))
                    g = None
                    for gene in genes:
                        if gene.attributes.get('gene_name', [None])[0] == gene_name:
                            g = gene
                            break
                    if g is None:
                        continue
                except Exception:
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

