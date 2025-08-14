"""
Enhanced gene ranking utility for long‑read sequencing experiments without replicates.

This module refines the original SimpleGeneRanker by incorporating best practices
from recent long‑read isoform switching studies.  It defines "interesting" genes
as those that are highly expressed, exhibit bona‑fide isoform switching (at
least two isoforms changing in opposite directions), and potentially show
functional consequences (e.g., gain or loss of coding potential).  Genes with
extreme overall expression changes or very complex isoform architectures are
penalised to reduce false positives.

Key features:

1. **Isoform count filter** – Genes with fewer than two transcripts are
   excluded, as isoform usage cannot change.  Genes with excessive numbers of
   isoforms can be down‑weighted via a complexity penalty.
2. **Bidirectional isoform switching** – For a gene to be considered a
   candidate isoform switcher, at least one transcript must increase in usage
   while another decreases between reference and target conditions.  This helps
   distinguish true isoform switches from uniform scaling of all isoforms.
3. **Functional impact assessment** – When transcript annotations include
   attributes such as coding status, ORF length or predicted NMD sensitivity,
   the ranker rewards genes where isoform switches change these properties.
   Lacking such annotations, this component defaults to zero influence.
4. **Adaptive thresholds** – Gating thresholds for expression level, fold
   change and usage delta are derived from quantiles of the observed
   distributions, making the algorithm robust across datasets with different
   scales.
5. **Extreme change penalty** – Genes with very large gene‑level fold changes
   are down‑weighted to prioritise isoform regulation over conventional
   differential expression.
6. **Categorised output** – The ranker labels the top genes according to
   whether they are isoform switchers, high expressers or conventional DEGs.

Example usage:

    from enhanced_gene_ranker import EnhancedGeneRanker
    ranker = EnhancedGeneRanker(output_dir="./out",
                                ref_conditions=["ref"],
                                target_conditions=["tgt"],
                                updated_gene_dict=gene_dict)
    top_genes = ranker.rank(top_n=50)
    # top_genes is a list of gene names

Note: This implementation assumes that `updated_gene_dict` follows the same
structure as in the original SimpleGeneRanker.  Transcript annotations may
contain keys such as ``coding`` (bool), ``orf_length`` (int) or
``functional_consequence`` (str).  Missing annotations are handled gracefully.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

# Configure module‑level logger
logger = logging.getLogger("IsoQuant.visualization.enhanced_ranker")
logger.setLevel(logging.INFO)


class SimpleGeneRanker:
    """Rank genes by integrating expression level, fold change, isoform switching and functional impact.

    Parameters
    ----------
    output_dir : str or Path
        Directory where intermediate results could be written (not used here but kept
        for compatibility).
    ref_conditions : List[str]
        List of keys in ``updated_gene_dict`` corresponding to reference
        conditions.
    target_conditions : List[str]
        List of keys in ``updated_gene_dict`` corresponding to target
        conditions.
    ref_only : bool, optional
        If True, only reference conditions will be considered (fold change
        computation disabled).  Defaults to False.
    updated_gene_dict : Dict, optional
        Nested dictionary with expression and transcript information.  See
        SimpleGeneRanker for expected format.
    """

    def __init__(
        self,
        output_dir: str | Path,
        ref_conditions: List[str],
        target_conditions: List[str],
        ref_only: bool = False,
        updated_gene_dict: Dict | None = None,
    ) -> None:
        self.output_dir = Path(output_dir)
        self.ref_conditions = list(ref_conditions)
        self.target_conditions = list(target_conditions)
        self.ref_only = ref_only
        self.updated_gene_dict: Dict[str, Dict] = updated_gene_dict or {}

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def rank(self, top_n: int = 100) -> List[str]:
        """Return a ranked list of gene names based on the enhanced scoring algorithm.

        Parameters
        ----------
        top_n : int
            Maximum number of genes to return.  If fewer genes meet the
            thresholds, the returned list may be shorter.

        Returns
        -------
        List[str]
            List of gene names (uppercase) ranked by decreasing composite score.
        """
        logger.info("Running EnhancedGeneRanker")

        if not self.updated_gene_dict:
            logger.warning("No updated_gene_dict provided; returning empty list")
            return []

        # Scan transcript ID patterns in updated_gene_dict
        self._scan_transcript_id_patterns()

        # 1. Extract gene‑level TPMs
        gene_expr_ref, gene_expr_tgt = self._extract_gene_tpms_from_dict()
        if gene_expr_ref.empty or gene_expr_tgt.empty:
            logger.warning("No gene expression data found; returning empty list")
            return []

        # 2. Compute gene list intersection
        common_genes = gene_expr_ref.index.intersection(gene_expr_tgt.index)
        gene_expr_ref = gene_expr_ref.loc[common_genes]
        gene_expr_tgt = gene_expr_tgt.loc[common_genes]
        logger.info(f"{len(common_genes)} genes present in both conditions")

        # 3. Compute isoform usage deltas and switching flags
        usage_delta, switch_flags, func_flags, isoform_counts = self._compute_isoform_usage_metrics(common_genes)

        # 4. Normalize features
        abs_log2fc = np.abs(np.log2(gene_expr_tgt + 1) - np.log2(gene_expr_ref + 1))
        geom_expr = np.sqrt(gene_expr_ref * gene_expr_tgt)

        norm_expr = self._normalize_feature(geom_expr, name="Expression")
        norm_change = self._normalize_feature(abs_log2fc, name="FoldChange")
        norm_usage = self._normalize_feature(usage_delta, name="UsageDelta")
        # Functional impact does not need normalization (0/1), but we convert to series
        func_series = pd.Series(func_flags, index=common_genes, dtype=float)

        # 5. Derive adaptive thresholds
        expr_gate = norm_expr.quantile(0.5)  # median
        change_gate = norm_change.quantile(0.75)  # upper quartile
        usage_gate = norm_usage.quantile(0.75)
        logger.info(
            f"Adaptive gates – expression: {expr_gate:.3f}, fold change: {change_gate:.3f}, usage delta: {usage_gate:.3f}"
        )

        # 6. Compute composite scores
        scores = pd.Series(0.0, index=common_genes)
        categories = {}
        for gene in common_genes:
            expr_val = norm_expr.at[gene]
            change_val = norm_change.at[gene]
            usage_val = norm_usage.at[gene]
            is_switch = switch_flags.get(gene, False)
            func_val = func_series.at[gene]
            iso_count = isoform_counts.get(gene, 0)

            # Handle single-transcript genes differently
            if iso_count == 1:
                # Single-transcript genes: focus on expression change, require higher thresholds
                single_transcript_expr_gate = norm_expr.quantile(0.7)  # Higher than multi-isoform
                single_transcript_change_gate = norm_change.quantile(0.85)  # Much higher fold change required
                
                passes_expr = expr_val > single_transcript_expr_gate
                passes_change = change_val > single_transcript_change_gate
                
                # Single-transcript penalty: they need to work harder to compete
                single_transcript_penalty = 0.6
                
                # Score based only on expression and fold change (no isoform switching possible)
                base = 0.0
                if passes_expr and passes_change:
                    base = expr_val * 0.4 + change_val * 0.6  # Weight fold change more heavily
                    if func_val > 0:
                        base += 0.3  # Smaller functional bonus than multi-isoform
                
                # Apply single-transcript penalty
                scores.at[gene] = base * single_transcript_penalty
                
                # Assign category
                if base == 0:
                    categories[gene] = "LOW_EXPR"
                else:
                    categories[gene] = "SINGLE_TRANSCRIPT_DE"  # New category
                    
            elif iso_count < 2:
                # Skip genes with 0 transcripts (shouldn't happen but safety check)
                continue
            else:
                # Multi-transcript genes: original logic with isoform switching
                passes_expr = expr_val > expr_gate
                passes_change = change_val > change_gate
                passes_usage = usage_val > usage_gate and is_switch

                # Penalise extreme expression changes (>90th percentile)
                penalty = 1.0
                if change_val > norm_change.quantile(0.9):
                    penalty *= 0.5

                # Complexity penalty: down‑weight genes with many isoforms (top 10%)
                if iso_count > np.quantile(list(isoform_counts.values()), 0.9):
                    penalty *= 0.7

                # Compute base score; weight usage more heavily when switching
                base = 0.0
                if passes_expr and (passes_change or passes_usage):
                    base = expr_val * 0.3 + change_val * 0.3 + usage_val * 1.2
                    if func_val > 0:
                        base += 0.5  # Functional impact bonus

                # Final score after penalty
                scores.at[gene] = base * penalty

                # Assign category for top genes later
                if base == 0:
                    categories[gene] = "LOW_EXPR"
                elif passes_usage and is_switch:
                    categories[gene] = "ISOFORM_SWITCHER"
                elif passes_expr and passes_change:
                    categories[gene] = "DIFFERENTIAL_EXPRESSION"
                else:
                    categories[gene] = "HIGH_EXPRESSION"

        # 7. Sort and select top genes
        ranked = scores[scores > 0].sort_values(ascending=False)
        ranked_gene_ids = ranked.head(top_n).index.tolist()

        # 8. Map gene IDs to names using updated_gene_dict
        ranked_gene_names: List[str] = []
        for gene_id in ranked_gene_ids:
            gene_name = None
            for cond in self.updated_gene_dict.values():
                if gene_id in cond:
                    gene_info = cond[gene_id]
                    gene_name = gene_info.get("name")
                    break
            ranked_gene_names.append(gene_name.upper() if gene_name else gene_id)

        # Log top entries with categories
        for gene_id in ranked_gene_ids[:10]:
            cat = categories.get(gene_id, "UNKNOWN")
            score = scores.at[gene_id]
            isoform_count = isoform_counts.get(gene_id, 0)
            logger.info(f"Top gene {gene_id}: score={score:.3f}, category={cat}, isoforms={isoform_count}")

        # Log single-transcript gene statistics
        self._log_single_transcript_statistics(categories, scores, isoform_counts)

        # Add lncRNA-specific statistics and biotype distribution
        self._log_biotype_distribution()
        self._log_lncrna_statistics(ranked_gene_ids, categories, scores)

        return ranked_gene_names

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _scan_transcript_id_patterns(self) -> None:
        """Scan updated_gene_dict for transcript ID patterns and report findings."""
        transcript_generic = 0  # Count of IDs starting with "transcript"
        transcript_ensembl = 0  # Count of proper Ensembl IDs (ENSMUST, ENST, etc.)
        transcript_other = 0    # Count of other patterns
        
        generic_examples = []
        ensembl_examples = []
        other_examples = []
        
        # Scan all conditions and genes
        for condition, genes in self.updated_gene_dict.items():
            for gene_id, gene_info in genes.items():
                transcripts_dict = gene_info.get("transcripts", {})
                for tx_id in transcripts_dict.keys():
                    if tx_id.lower().startswith("transcript"):
                        transcript_generic += 1
                        if len(generic_examples) < 5:
                            generic_examples.append(f"{gene_id}:{tx_id}")
                    elif tx_id.startswith("ENSMUST") or tx_id.startswith("ENST") or tx_id.startswith("ENS"):
                        transcript_ensembl += 1
                        if len(ensembl_examples) < 5:
                            ensembl_examples.append(f"{gene_id}:{tx_id}")
                    else:
                        transcript_other += 1
                        if len(other_examples) < 5:
                            other_examples.append(f"{gene_id}:{tx_id}")
        
        total_transcripts = transcript_generic + transcript_ensembl + transcript_other
        
        logger.info(f"TRANSCRIPT_SCAN: Found {total_transcripts} total transcript entries")
        
        if total_transcripts > 0:
            logger.info(f"TRANSCRIPT_SCAN: {transcript_generic} transcripts start with 'transcript' ({transcript_generic/total_transcripts*100:.1f}%)")
            logger.info(f"TRANSCRIPT_SCAN: {transcript_ensembl} transcripts are Ensembl IDs ({transcript_ensembl/total_transcripts*100:.1f}%)")
            logger.info(f"TRANSCRIPT_SCAN: {transcript_other} transcripts have other patterns ({transcript_other/total_transcripts*100:.1f}%)")
        else:
            logger.warning("TRANSCRIPT_SCAN: No transcripts found in updated_gene_dict")
        
        if generic_examples:
            logger.info(f"TRANSCRIPT_SCAN: Generic transcript examples: {generic_examples}")
        if ensembl_examples:
            logger.info(f"TRANSCRIPT_SCAN: Ensembl transcript examples: {ensembl_examples}")
        if other_examples:
            logger.info(f"TRANSCRIPT_SCAN: Other transcript examples: {other_examples}")
        
        # Warn if too many generic transcripts
        if transcript_generic > transcript_ensembl:
            logger.warning(f"TRANSCRIPT_SCAN: More generic 'transcript' IDs ({transcript_generic}) than Ensembl IDs ({transcript_ensembl}) - this may indicate annotation issues")

    def _normalize_feature(self, feature: pd.Series, name: str) -> pd.Series:
        """Normalize a numeric Series to the 0–1 range.  If all values are equal,
        return zeros.

        Parameters
        ----------
        feature : pd.Series
            The numeric data to normalize.
        name : str
            Name of the feature for logging.

        Returns
        -------
        pd.Series
            Normalized series with the same index.
        """
        if feature.empty:
            return pd.Series(dtype=float)
        fmin, fmax = feature.min(), feature.max()
        if fmax == fmin:
            logger.warning(f"{name} values are identical; returning zeros")
            return pd.Series(0.0, index=feature.index)
        norm = (feature - fmin) / (fmax - fmin)
        return norm.fillna(0.0)

    def _extract_gene_tpms_from_dict(self) -> Tuple[pd.Series, pd.Series]:
        """Extract average gene TPMs for reference and target conditions.

        This method sums transcript TPMs within each gene and averages across
        multiple conditions.  It mirrors the corresponding method in the
        SimpleGeneRanker but returns pandas Series indexed by gene ID.
        """
        ref_gene_tpms: Dict[str, float] = {}
        tgt_gene_tpms: Dict[str, float] = {}
        ref_count = 0
        tgt_count = 0

        # Process reference conditions
        for cond in self.ref_conditions:
            if cond not in self.updated_gene_dict:
                continue
            ref_count += 1
            genes = self.updated_gene_dict[cond]
            for gene_id, gene_info in genes.items():
                tpm = 0.0
                for tx_info in gene_info.get("transcripts", {}).values():
                    if isinstance(tx_info, dict) and "value" in tx_info:
                        tpm += tx_info["value"]
                ref_gene_tpms[gene_id] = ref_gene_tpms.get(gene_id, 0.0) + tpm

        # Average across reference conditions
        for gene_id in ref_gene_tpms:
            if ref_count > 0:
                ref_gene_tpms[gene_id] /= ref_count

        # Process target conditions
        for cond in self.target_conditions:
            if cond not in self.updated_gene_dict:
                continue
            tgt_count += 1
            genes = self.updated_gene_dict[cond]
            for gene_id, gene_info in genes.items():
                tpm = 0.0
                for tx_info in gene_info.get("transcripts", {}).values():
                    if isinstance(tx_info, dict) and "value" in tx_info:
                        tpm += tx_info["value"]
                tgt_gene_tpms[gene_id] = tgt_gene_tpms.get(gene_id, 0.0) + tpm

        # Average across target conditions
        for gene_id in tgt_gene_tpms:
            if tgt_count > 0:
                tgt_gene_tpms[gene_id] /= tgt_count

        ref_series = pd.Series(ref_gene_tpms, name="ref_tpm")
        tgt_series = pd.Series(tgt_gene_tpms, name="tgt_tpm")
        return ref_series, tgt_series

    def _compute_isoform_usage_metrics(self, genes: List[str]) -> Tuple[pd.Series, Dict[str, bool], Dict[str, bool], Dict[str, int]]:
        """Compute isoform usage delta, switching flag, functional impact flag and isoform count.

        Parameters
        ----------
        genes : list of str
            Gene identifiers to process.

        Returns
        -------
        usage_delta : pd.Series
            Maximum absolute change in isoform usage per gene.
        switch_flags : dict
            Dictionary mapping gene IDs to True if at least one transcript
            increases and another decreases in usage between reference and target.
        func_flags : dict
            Dictionary mapping gene IDs to True if the isoform switching implies
            a change in coding status or functional consequence.
        isoform_counts : dict
            Dictionary mapping gene IDs to the number of isoforms detected.
        """
        usage_delta = pd.Series(0.0, index=genes)
        switch_flags: Dict[str, bool] = {}
        func_flags: Dict[str, bool] = {}
        isoform_counts: Dict[str, int] = {}

        # Build per‑condition transcript TPM dictionaries
        ref_tx = self._aggregate_transcript_tpms(self.ref_conditions)
        tgt_tx = self._aggregate_transcript_tpms(self.target_conditions)

        # For each gene, compute usage change
        for gene in genes:
            # Collect transcripts and counts
            tx_ids = set()
            for cond in self.ref_conditions:
                cond_dict = self.updated_gene_dict.get(cond, {})
                if gene in cond_dict:
                    tx_ids.update(cond_dict[gene].get("transcripts", {}).keys())
            for cond in self.target_conditions:
                cond_dict = self.updated_gene_dict.get(cond, {})
                if gene in cond_dict:
                    tx_ids.update(cond_dict[gene].get("transcripts", {}).keys())
            isoform_counts[gene] = len(tx_ids)
            if len(tx_ids) < 2:
                switch_flags[gene] = False
                func_flags[gene] = False
                usage_delta.at[gene] = 0.0
                continue

            # Compute usage per condition
            ref_total = 0.0
            tgt_total = 0.0
            ref_usages: Dict[str, float] = {}
            tgt_usages: Dict[str, float] = {}
            for tx_id in tx_ids:
                r_tpm = ref_tx.get(tx_id, 0.0)
                t_tpm = tgt_tx.get(tx_id, 0.0)
                ref_total += r_tpm
                tgt_total += t_tpm
                ref_usages[tx_id] = r_tpm
                tgt_usages[tx_id] = t_tpm
            ref_total += 1e-6  # avoid zero division
            tgt_total += 1e-6

            # Compute usage fractions
            deltas = []
            directions = []
            func_change = False
            for tx_id in tx_ids:
                ref_u = ref_usages[tx_id] / ref_total
                tgt_u = tgt_usages[tx_id] / tgt_total
                delta = tgt_u - ref_u
                deltas.append(abs(delta))
                directions.append(np.sign(delta))

                # Assess functional impact if annotation exists
                # We check across any condition; assume annotation consistent
                for cond in self.updated_gene_dict:
                    cond_dict = self.updated_gene_dict[cond]
                    if gene in cond_dict:
                        tx_info = cond_dict[gene].get("transcripts", {}).get(tx_id, {})
                        # Compare coding status and ORF length relative to other transcripts
                        coding = tx_info.get("coding")
                        orf_len = tx_info.get("orf_length")
                        func = tx_info.get("functional_consequence")
                        break
                # Simple heuristic: if any transcript has non‑zero functional_consequence
                if func is not None:
                    func_change = True
            # Determine switching: at least one positive and one negative change
            switch_flags[gene] = (1 in directions) and (-1 in directions)
            func_flags[gene] = func_change
            usage_delta.at[gene] = max(deltas) if deltas else 0.0

        return usage_delta, switch_flags, func_flags, isoform_counts

    def _aggregate_transcript_tpms(self, conditions: List[str]) -> Dict[str, float]:
        """Aggregate transcript TPMs across a list of conditions, averaging across
        conditions.

        Parameters
        ----------
        conditions : list of str
            Conditions to aggregate.

        Returns
        -------
        Dict[str, float]
            Mapping from transcript ID to averaged TPM value.
        """
        tx_totals: Dict[str, float] = {}
        count = 0
        for cond in conditions:
            if cond not in self.updated_gene_dict:
                continue
            count += 1
            cond_dict = self.updated_gene_dict[cond]
            for gene_info in cond_dict.values():
                for tx_id, tx_info in gene_info.get("transcripts", {}).items():
                    if isinstance(tx_info, dict) and "value" in tx_info:
                        tx_totals[tx_id] = tx_totals.get(tx_id, 0.0) + tx_info["value"]
        # Average
        if count > 0:
            for tx in tx_totals:
                tx_totals[tx] /= count
        return tx_totals

    def _log_single_transcript_statistics(self, categories: Dict[str, str], scores: pd.Series, isoform_counts: Dict[str, int]) -> None:
        """Log statistics about single-transcript genes and how they performed."""
        single_transcript_genes = [gene_id for gene_id, count in isoform_counts.items() if count == 1]
        multi_transcript_genes = [gene_id for gene_id, count in isoform_counts.items() if count > 1]
        
        single_transcript_scored = [gene_id for gene_id in single_transcript_genes if scores.get(gene_id, 0) > 0]
        single_transcript_in_categories = [gene_id for gene_id in single_transcript_genes if categories.get(gene_id) == "SINGLE_TRANSCRIPT_DE"]
        
        # Get top single-transcript genes by score
        single_transcript_scores = {gene_id: scores.get(gene_id, 0) for gene_id in single_transcript_genes}
        top_single_transcript = sorted(single_transcript_scores.items(), key=lambda x: x[1], reverse=True)[:5]
        
        # Get expression info for top single-transcript genes
        top_single_examples = []
        for gene_id, score in top_single_transcript[:3]:
            if score > 0:
                # Find gene name and expression values
                gene_name = gene_id
                ref_tpm = 0
                tgt_tpm = 0
                
                for condition, genes in self.updated_gene_dict.items():
                    if gene_id in genes:
                        gene_info = genes[gene_id]
                        gene_name = gene_info.get("name", gene_id)
                        
                        # Get the single transcript's TPM
                        transcripts = gene_info.get("transcripts", {})
                        if transcripts:
                            tx_id, tx_info = next(iter(transcripts.items()))
                            tpm = tx_info.get("value", 0)
                            
                            if condition in self.ref_conditions:
                                ref_tpm += tpm / len(self.ref_conditions)
                            elif condition in self.target_conditions:
                                tgt_tpm += tpm / len(self.target_conditions)
                
                fold_change = (tgt_tpm + 0.1) / (ref_tpm + 0.1)  # Add pseudocount
                top_single_examples.append({
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "score": score,
                    "ref_tpm": ref_tpm,
                    "tgt_tpm": tgt_tpm,
                    "fold_change": fold_change
                })
        
        logger.info("=== SINGLE-TRANSCRIPT GENE ANALYSIS ===")
        logger.info(f"Total single-transcript genes: {len(single_transcript_genes)}")
        logger.info(f"Total multi-transcript genes: {len(multi_transcript_genes)}")
        logger.info(f"Single-transcript genes with scores > 0: {len(single_transcript_scored)}")
        logger.info(f"Single-transcript genes passing high thresholds: {len(single_transcript_in_categories)}")
        
        if len(single_transcript_genes) > 0:
            pass_rate = (len(single_transcript_in_categories) / len(single_transcript_genes)) * 100
            logger.info(f"Single-transcript gene pass rate: {pass_rate:.1f}% (requires 85th percentile fold change)")
        
        if top_single_examples:
            logger.info("Top single-transcript genes by score:")
            for i, example in enumerate(top_single_examples, 1):
                logger.info(f"  {i}. {example['gene_name']} ({example['gene_id']}): "
                           f"score={example['score']:.3f}, "
                           f"ref_TPM={example['ref_tpm']:.1f}, "
                           f"tgt_TPM={example['tgt_tpm']:.1f}, "
                           f"FC={example['fold_change']:.2f}x")
        
        # Compare to multi-transcript genes
        multi_transcript_scored = [gene_id for gene_id in multi_transcript_genes if scores.get(gene_id, 0) > 0]
        if len(multi_transcript_genes) > 0:
            multi_pass_rate = (len(multi_transcript_scored) / len(multi_transcript_genes)) * 100
            logger.info(f"Multi-transcript gene pass rate: {multi_pass_rate:.1f}% (for comparison)")

    def _log_biotype_distribution(self) -> None:
        """Log the distribution of gene biotypes in the dataset."""
        biotype_counts = {}
        total_genes = 0
        
        # Count biotypes across all genes (avoiding duplicates by using first condition)
        for condition, genes in self.updated_gene_dict.items():
            for gene_id, gene_info in genes.items():
                biotype = gene_info.get("biotype", "unknown")
                biotype_counts[biotype] = biotype_counts.get(biotype, 0) + 1
                total_genes += 1
            break  # Only count from one condition to avoid duplicates
        
        if total_genes == 0:
            logger.warning("No genes found for biotype distribution analysis")
            return
        
        # Sort biotypes by count
        sorted_biotypes = sorted(biotype_counts.items(), key=lambda x: x[1], reverse=True)
        
        logger.info("=== GENE BIOTYPE DISTRIBUTION ===")
        logger.info(f"Total genes analyzed: {total_genes}")
        
        for biotype, count in sorted_biotypes:
            percentage = (count / total_genes) * 100
            logger.info(f"{biotype}: {count} genes ({percentage:.1f}%)")
        
        # Highlight key biotypes
        lncrna_count = biotype_counts.get("lncRNA", 0) + biotype_counts.get("long_noncoding_rna", 0)
        protein_coding_count = biotype_counts.get("protein_coding", 0)
        pseudogene_counts = sum(count for biotype, count in biotype_counts.items() 
                               if "pseudogene" in biotype.lower())
        
        logger.info("=== KEY BIOTYPE SUMMARY ===")
        logger.info(f"Protein coding genes: {protein_coding_count} ({(protein_coding_count/total_genes)*100:.1f}%)")
        logger.info(f"lncRNA genes: {lncrna_count} ({(lncrna_count/total_genes)*100:.1f}%)")
        logger.info(f"Pseudogenes (all types): {pseudogene_counts} ({(pseudogene_counts/total_genes)*100:.1f}%)")

    def _log_lncrna_statistics(self, ranked_gene_ids: List[str], categories: Dict[str, str], scores: pd.Series) -> None:
        """Analyze and log statistics about lncRNAs in the ranked gene list."""
        lncrna_stats = {
            "total_lncrnas": 0,
            "top_50_lncrnas": 0,
            "top_10_lncrnas": 0,
            "isoform_switcher_lncrnas": 0,
            "high_expr_lncrnas": 0,
            "lncrna_examples": []
        }
        
        all_lncrnas = []
        
        # Scan all genes to find lncRNAs
        for condition, genes in self.updated_gene_dict.items():
            for gene_id, gene_info in genes.items():
                gene_biotype = gene_info.get("biotype", "").lower()
                if gene_biotype in ["lncrna", "long_noncoding_rna", "lincrna"]:
                    lncrna_stats["total_lncrnas"] += 1
                    gene_name = gene_info.get("name", gene_id)
                    score = scores.get(gene_id, 0.0)
                    category = categories.get(gene_id, "LOW_EXPR")
                    
                    all_lncrnas.append({
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                        "score": score,
                        "category": category,
                        "biotype": gene_biotype
                    })
                    
                    # Count categories
                    if category == "ISOFORM_SWITCHER":
                        lncrna_stats["isoform_switcher_lncrnas"] += 1
                    elif category == "HIGH_EXPRESSION":
                        lncrna_stats["high_expr_lncrnas"] += 1
                    
                    # Check if in top rankings
                    if gene_id in ranked_gene_ids[:50]:
                        lncrna_stats["top_50_lncrnas"] += 1
                    if gene_id in ranked_gene_ids[:10]:
                        lncrna_stats["top_10_lncrnas"] += 1
            break  # Only check one condition to avoid duplicates
        
        # Sort lncRNAs by score and get top examples
        all_lncrnas.sort(key=lambda x: x["score"], reverse=True)
        lncrna_stats["lncrna_examples"] = all_lncrnas[:5]
        
        # Log comprehensive lncRNA statistics
        logger.info("=== lncRNA ANALYSIS ===")
        logger.info(f"Total lncRNAs detected: {lncrna_stats['total_lncrnas']}")
        logger.info(f"lncRNAs in top 50 genes: {lncrna_stats['top_50_lncrnas']}")
        logger.info(f"lncRNAs in top 10 genes: {lncrna_stats['top_10_lncrnas']}")
        logger.info(f"lncRNA isoform switchers: {lncrna_stats['isoform_switcher_lncrnas']}")
        logger.info(f"High-expression lncRNAs: {lncrna_stats['high_expr_lncrnas']}")
        
        if lncrna_stats["total_lncrnas"] > 0:
            top_50_pct = (lncrna_stats["top_50_lncrnas"] / lncrna_stats["total_lncrnas"]) * 100
            logger.info(f"Percentage of lncRNAs in top 50: {top_50_pct:.1f}%")
        
        # Log top lncRNA examples
        if lncrna_stats["lncrna_examples"]:
            logger.info("Top scoring lncRNAs:")
            for i, lncrna in enumerate(lncrna_stats["lncrna_examples"], 1):
                logger.info(f"  {i}. {lncrna['gene_name']} ({lncrna['gene_id']}): "
                           f"score={lncrna['score']:.3f}, category={lncrna['category']}")
        
        # Analyze transcript complexity for top lncRNAs
        self._analyze_lncrna_transcript_complexity(lncrna_stats["lncrna_examples"][:3])
    
    def _analyze_lncrna_transcript_complexity(self, top_lncrnas: List[Dict]) -> None:
        """Analyze transcript complexity and biotype diversity for top lncRNAs."""
        if not top_lncrnas:
            return
            
        logger.info("=== lncRNA TRANSCRIPT COMPLEXITY ===")
        
        for lncrna in top_lncrnas:
            gene_id = lncrna["gene_id"]
            gene_name = lncrna["gene_name"]
            
            # Find transcript information across conditions
            transcript_info = {}
            for condition, genes in self.updated_gene_dict.items():
                if gene_id in genes:
                    transcripts = genes[gene_id].get("transcripts", {})
                    for tx_id, tx_info in transcripts.items():
                        tx_biotype = tx_info.get("biotype", "unknown")
                        tx_value = tx_info.get("value", 0.0)
                        tx_name = tx_info.get("name", tx_id)
                        
                        if tx_id not in transcript_info:
                            transcript_info[tx_id] = {
                                "name": tx_name,
                                "biotype": tx_biotype,
                                "values": []
                            }
                        transcript_info[tx_id]["values"].append(tx_value)
            
            # Calculate transcript statistics
            transcript_count = len(transcript_info)
            biotype_counts = {}
            active_transcripts = 0
            
            for tx_id, tx_data in transcript_info.items():
                biotype = tx_data["biotype"]
                biotype_counts[biotype] = biotype_counts.get(biotype, 0) + 1
                
                avg_value = sum(tx_data["values"]) / len(tx_data["values"]) if tx_data["values"] else 0
                if avg_value > 1.0:  # Consider active if TPM > 1
                    active_transcripts += 1
            
            logger.info(f"lncRNA {gene_name} ({gene_id}):")
            logger.info(f"  Total transcripts: {transcript_count}")
            logger.info(f"  Active transcripts (TPM>1): {active_transcripts}")
            logger.info(f"  Transcript biotypes: {dict(biotype_counts)}")
            
            # Show top 3 most expressed transcripts
            if transcript_info:
                sorted_transcripts = sorted(
                    transcript_info.items(),
                    key=lambda x: sum(x[1]["values"]) / len(x[1]["values"]) if x[1]["values"] else 0,
                    reverse=True
                )
                logger.info("  Top transcripts by expression:")
                for i, (tx_id, tx_data) in enumerate(sorted_transcripts[:3], 1):
                    avg_expr = sum(tx_data["values"]) / len(tx_data["values"]) if tx_data["values"] else 0
                    logger.info(f"    {i}. {tx_data['name']}: {avg_expr:.2f} TPM ({tx_data['biotype']})")
