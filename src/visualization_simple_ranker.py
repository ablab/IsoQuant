#!/usr/bin/env python3
"""A lightweight gene ranking utility for experiments lacking biological replicates.

This module implements a fallback algorithm used when each experimental
condition has fewer than two biological replicates and, therefore, formal
statistics with DESeq2 are inappropriate.

The scoring heuristic combines two intuitive effect-size metrics:
    1. Absolute log2 fold-change of gene-level TPM between target and
       reference groups.
    2. Maximum change in isoform usage for any transcript belonging to a
       gene. Usage is defined as transcript TPM / total gene TPM.

The final score is:
    score = |log2FC_gene| + max_delta_isoform_usage

Genes are ranked by the score in descending order.

The implementation is designed to mirror the interface expected by
visualize.py – namely, a ``rank(top_n)`` method that returns a list of gene
names (mapped from gene IDs using the updated gene dictionary).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Dict

import numpy as np
import pandas as pd

logger = logging.getLogger("IsoQuant.visualization.simple_ranker")
logger.setLevel(logging.INFO)


class SimpleGeneRanker:
    """Rank genes by combined gene-expression and isoform-usage change."""

    def __init__(
        self,
        output_dir: str | Path,
        ref_conditions: List[str],
        target_conditions: List[str],
        ref_only: bool = False,
        updated_gene_dict: Dict = None,
    ) -> None:
        self.output_dir = Path(output_dir)
        self.ref_conditions = list(ref_conditions)
        self.target_conditions = list(target_conditions)
        self.ref_only = ref_only
        self.updated_gene_dict = updated_gene_dict or {}

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def rank(self, top_n: int = 100) -> List[str]:
        logger.info("Running SimpleGeneRanker (no replicates detected)")
        
        # Debug: Log updated_gene_dict structure
        if self.updated_gene_dict:
            logger.info(f"updated_gene_dict has {len(self.updated_gene_dict)} conditions: {list(self.updated_gene_dict.keys())}")
            # Show sample genes from each condition
            for condition, genes in self.updated_gene_dict.items():
                sample_gene_ids = list(genes.keys())[:3]
                logger.info(f"Condition '{condition}' has {len(genes)} genes. Sample gene IDs: {sample_gene_ids}")
                # Show gene structure for first gene
                if sample_gene_ids:
                    first_gene_id = sample_gene_ids[0]
                    first_gene_info = genes[first_gene_id]
                    logger.info(f"  Sample gene '{first_gene_id}' structure: name='{first_gene_info.get('name', 'MISSING')}', keys={list(first_gene_info.keys())}")
        else:
            logger.warning("No updated_gene_dict provided!")

        # 1. Load gene-level TPM aggregated by sample.
        logger.info(f"SIMPLE_RANKER: Loading TPM data for reference conditions: {self.ref_conditions}")
        gene_expr_ref = self._aggregate_gene_tpm(self.ref_conditions)
        logger.info(f"SIMPLE_RANKER: Reference TPM loaded - shape: {gene_expr_ref.shape}")
        
        logger.info(f"SIMPLE_RANKER: Loading TPM data for target conditions: {self.target_conditions}")
        gene_expr_tgt = self._aggregate_gene_tpm(self.target_conditions)
        logger.info(f"SIMPLE_RANKER: Target TPM loaded - shape: {gene_expr_tgt.shape}")

        # Genes common to both groups.
        common_genes = gene_expr_ref.index.intersection(gene_expr_tgt.index)
        logger.info(f"SIMPLE_RANKER: Found {len(common_genes)} genes common to both ref and target groups")
        logger.info(f"SIMPLE_RANKER: Sample common genes: {list(common_genes)[:10]}")
        
        gene_expr_ref = gene_expr_ref.loc[common_genes]
        gene_expr_tgt = gene_expr_tgt.loc[common_genes]
        logger.info(f"SIMPLE_RANKER: After filtering to common genes - ref shape: {gene_expr_ref.shape}, target shape: {gene_expr_tgt.shape}")

        # Filter to only genes present in updated_gene_dict (like differential expression analysis does)
        if self.updated_gene_dict:
            available_genes = set()
            for condition_dict in self.updated_gene_dict.values():
                available_genes.update(condition_dict.keys())
            
            # Log sample of available genes from updated_gene_dict
            sample_genes = list(available_genes)[:10]
            logger.info(f"Sample genes available in updated_gene_dict: {sample_genes}")
            
            # Keep only genes that are in the updated_gene_dict
            filtered_common_genes = [g for g in common_genes if g in available_genes]
            logger.info(f"Filtered genes to {len(filtered_common_genes)} from {len(common_genes)} based on updated_gene_dict availability")
            
            if not filtered_common_genes:
                logger.warning("No genes remain after filtering by updated_gene_dict. Returning empty list.")
                return []
            
            gene_expr_ref = gene_expr_ref.loc[filtered_common_genes]
            gene_expr_tgt = gene_expr_tgt.loc[filtered_common_genes]
            common_genes = filtered_common_genes

        # 2. Compute log2 fold-change (add pseudocount of 1).
        logger.info("SIMPLE_RANKER: Computing log2 fold-change...")
        log2fc = np.log2(gene_expr_tgt + 1) - np.log2(gene_expr_ref + 1)
        abs_log2fc = log2fc.abs()
        logger.info(f"SIMPLE_RANKER: Log2FC stats - min: {log2fc.min():.3f}, max: {log2fc.max():.3f}, mean: {log2fc.mean():.3f}")
        logger.info(f"SIMPLE_RANKER: Abs Log2FC stats - min: {abs_log2fc.min():.3f}, max: {abs_log2fc.max():.3f}, mean: {abs_log2fc.mean():.3f}")
        
        # Show top log2FC examples
        top_log2fc_genes = abs_log2fc.nlargest(5)
        logger.info(f"SIMPLE_RANKER: Top 5 genes by abs log2FC:")
        for gene_id, score in top_log2fc_genes.items():
            ref_val = gene_expr_ref[gene_id]
            tgt_val = gene_expr_tgt[gene_id]
            logger.info(f"  {gene_id}: ref_tpm={ref_val:.3f}, tgt_tpm={tgt_val:.3f}, abs_log2fc={score:.3f}")

        # 3. Compute isoform-usage change per gene.
        logger.info("SIMPLE_RANKER: Computing isoform usage delta...")
        delta_usage = self._compute_isoform_usage_delta(common_genes)
        logger.info(f"SIMPLE_RANKER: Delta usage stats - min: {delta_usage.min():.3f}, max: {delta_usage.max():.3f}, mean: {delta_usage.mean():.3f}")

        # 4. Combined score.
        logger.info("SIMPLE_RANKER: Computing combined score = |log2FC| + max_delta_isoform_usage...")
        combined_score = abs_log2fc + delta_usage
        combined_score.name = "score"
        logger.info(f"SIMPLE_RANKER: Combined score stats - min: {combined_score.min():.3f}, max: {combined_score.max():.3f}, mean: {combined_score.mean():.3f}")

        # Show detailed scoring for top genes
        top_combined_genes = combined_score.nlargest(10)
        logger.info(f"SIMPLE_RANKER: Top 10 genes by combined score:")
        for gene_id, score in top_combined_genes.items():
            log2fc_contrib = abs_log2fc[gene_id]
            usage_contrib = delta_usage[gene_id]
            logger.info(f"  {gene_id}: total_score={score:.3f} (log2fc={log2fc_contrib:.3f} + usage={usage_contrib:.3f})")

        # 5. Rank and get top N gene IDs.
        ranked_gene_ids = combined_score.sort_values(ascending=False).head(top_n).index.tolist()
        logger.info(f"SimpleGeneRanker selected {len(ranked_gene_ids)} genes (top {top_n}) by score.")
        logger.info(f"Top 10 ranked gene IDs: {ranked_gene_ids[:10]}")
        
        # Show the final scores for the top ranked genes
        final_scores = combined_score.loc[ranked_gene_ids[:10]]
        logger.info(f"Final scores for top 10 genes: {final_scores.to_dict()}")

        # 6. Map gene IDs to gene names directly from updated_gene_dict to ensure exact compatibility with plotter
        if self.updated_gene_dict:
            ranked_gene_names = []
            mapped_count = 0
            mapping_details = []  # For detailed logging
            
            for gene_id in ranked_gene_ids:
                gene_name_found = None
                # Look for this gene_id in updated_gene_dict to get the exact gene name
                for condition_dict in self.updated_gene_dict.values():
                    if gene_id in condition_dict:
                        gene_info = condition_dict[gene_id]
                        if "name" in gene_info and gene_info["name"]:
                            gene_name_found = gene_info["name"]
                            mapped_count += 1
                            mapping_details.append(f"{gene_id} -> {gene_name_found}")
                        else:
                            mapping_details.append(f"{gene_id} -> NO_NAME (name field: {gene_info.get('name', 'MISSING')})")
                        break
                else:
                    mapping_details.append(f"{gene_id} -> NOT_FOUND_IN_DICT")
                
                # Use the found gene name (uppercase to match plotter expectations) or fallback to gene_id
                ranked_gene_names.append(gene_name_found.upper() if gene_name_found else gene_id)
            
            # Log mapping details for first 10 genes
            logger.info(f"Gene mapping details (first 10):")
            for detail in mapping_details[:10]:
                logger.info(f"  {detail}")
            
            logger.info(f"SimpleGeneRanker mapped {mapped_count}/{len(ranked_gene_ids)} gene IDs to gene names from updated_gene_dict.")
            logger.info(f"Final gene names (first 10): {ranked_gene_names[:10]}")
            return ranked_gene_names
        else:
            logger.warning("No updated_gene_dict provided, returning raw gene IDs.")
            return ranked_gene_ids

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    # Cached root-level matrices to avoid re-reading
    _root_gene_tpm: pd.DataFrame | None = None
    _root_transcript_tpm: pd.DataFrame | None = None

    def _aggregate_gene_tpm(self, conditions: List[str]) -> pd.Series:
        """Return mean TPM vector across all samples of the given conditions.

        Works with two possible layouts:
        1. YAML / multi-condition run – counts live under <output_dir>/<condition>/.
        2. Single-sample run – a single *gene_grouped_tpm.tsv in <output_dir>/.
        """
        logger.info(f"SIMPLE_RANKER: Aggregating gene TPM for conditions: {conditions}")
        tpm_values: list[pd.Series] = []
        
        for cond in conditions:
            cond_dir = self.output_dir / cond
            files = list(cond_dir.glob("*gene_grouped_tpm.tsv"))
            logger.info(f"SIMPLE_RANKER: Condition '{cond}' - found {len(files)} gene TPM files in {cond_dir}")
            
            if files:
                for fp in files:
                    logger.info(f"SIMPLE_RANKER: Reading gene TPM from: {fp}")
                    df = self._read_tpm(fp)
                    logger.info(f"SIMPLE_RANKER: Gene TPM file shape: {df.shape}, columns: {list(df.columns)[:5]}...")
                    logger.info(f"SIMPLE_RANKER: Sample gene IDs: {list(df.index)[:5]}...")
                    tpm_series = df.sum(axis=1)
                    logger.info(f"SIMPLE_RANKER: Summed TPM series shape: {tpm_series.shape}, sample values: {tpm_series.head()}")
                    tpm_values.append(tpm_series)
                continue  # next condition

            # Fallback: root-level file present
            logger.info(f"SIMPLE_RANKER: No condition-specific files found for '{cond}', trying root-level gene TPM...")
            root_df = self._get_root_gene_tpm()
            logger.info(f"SIMPLE_RANKER: Root gene TPM shape: {root_df.shape}, columns: {list(root_df.columns)}")
            cond_cols = [c for c in root_df.columns if c == cond or c.startswith(f"{cond}__")]
            logger.info(f"SIMPLE_RANKER: Found {len(cond_cols)} columns for condition '{cond}': {cond_cols}")
            if not cond_cols:
                logger.warning("Condition '%s' columns not found in root gene_grouped_tpm.tsv; treating as missing.", cond)
                continue
            cond_tpm = root_df[cond_cols].mean(axis=1)
            logger.info(f"SIMPLE_RANKER: Condition TPM series shape: {cond_tpm.shape}, sample values: {cond_tpm.head()}")
            tpm_values.append(cond_tpm)

        if not tpm_values:
            logger.error("SIMPLE_RANKER: No TPM values found for provided conditions.")
            raise FileNotFoundError("No TPM values found for provided conditions.")
        
        logger.info(f"SIMPLE_RANKER: Collected {len(tpm_values)} TPM series for conditions {conditions}")
        stacked = pd.concat(tpm_values, axis=1)
        logger.info(f"SIMPLE_RANKER: Stacked TPM shape: {stacked.shape}")
        aggregated = stacked.mean(axis=1)
        logger.info(f"SIMPLE_RANKER: Final aggregated TPM shape: {aggregated.shape}, sample values: {aggregated.head()}")
        logger.info(f"SIMPLE_RANKER: TPM stats - min: {aggregated.min():.3f}, max: {aggregated.max():.3f}, mean: {aggregated.mean():.3f}")
        return aggregated

    def _get_root_gene_tpm(self) -> pd.DataFrame:
        if self._root_gene_tpm is not None:
            return self._root_gene_tpm
        files = list(self.output_dir.glob("*gene_grouped_tpm.tsv"))
        if not files:
            raise FileNotFoundError("Root-level gene_grouped_tpm.tsv not found in output directory.")
        df = self._read_tpm(files[0])
        self._root_gene_tpm = df
        return df

    def _read_tpm(self, fp: Path) -> pd.DataFrame:
        df = pd.read_csv(fp, sep="\t")
        first_col = df.columns[0]
        if first_col.startswith("#"):
            df.rename(columns={first_col: "feature_id"}, inplace=True)
        return df.set_index("feature_id")

    def _compute_isoform_usage_delta(self, gene_list: List[str]) -> pd.Series:
        """Return a Series of maximal isoform-usage change for each gene."""
        # Load transcript TPMs for each group and compute usage.
        usage_diff = pd.Series(0.0, index=gene_list)

        # Attempt to locate transcript TPM files.
        trans_ref = self._aggregate_transcript_tpm(self.ref_conditions)
        trans_tgt = self._aggregate_transcript_tpm(self.target_conditions)

        if trans_ref.empty or trans_tgt.empty:
            logger.warning("Transcript TPM files missing – setting isoform usage component to 0.")
            return usage_diff

        # Common transcripts.
        common_tx = trans_ref.index.intersection(trans_tgt.index)
        trans_ref = trans_ref.loc[common_tx]
        trans_tgt = trans_tgt.loc[common_tx]

        # Map transcripts to genes by simple split (before '.')
        gene_ids = trans_ref.index.to_series().str.split(".").str[0]
        trans_ref_grouped = trans_ref.groupby(gene_ids).sum()
        trans_tgt_grouped = trans_tgt.groupby(gene_ids).sum()

        # Compute per-gene usage change.
        for gene in gene_list:
            if gene not in trans_ref_grouped.index or gene not in trans_tgt_grouped.index:
                continue
            # Filter transcripts belonging to this gene.
            mask = gene_ids == gene
            gene_tx_ref = trans_ref[mask]
            gene_tx_tgt = trans_tgt[mask]

            # Gene totals (add 1e-6 to avoid divide-by-zero)
            ref_total = gene_tx_ref.sum() + 1e-6
            tgt_total = gene_tx_tgt.sum() + 1e-6

            ref_usage = gene_tx_ref / ref_total
            tgt_usage = gene_tx_tgt / tgt_total
            max_delta = (ref_usage - tgt_usage).abs().max()
            usage_diff.at[gene] = max_delta
        return usage_diff

    def _aggregate_transcript_tpm(self, conditions: List[str]) -> pd.Series:
        """Return mean TPM per transcript across all samples in conditions (handles both layouts)."""
        tpm_values: list[pd.Series] = []
        pattern_default = "*transcript_grouped_tpm.tsv"
        pattern_model = "*transcript_model_grouped_tpm.tsv"

        for cond in conditions:
            cond_dir = self.output_dir / cond
            pattern = pattern_model if self.ref_only else pattern_default
            files = list(cond_dir.glob(pattern))
            if files:
                for fp in files:
                    df = self._read_tpm(fp)
                    tpm_values.append(df.sum(axis=1))
                continue

            # Fallback: root-level transcript TPM
            root_df = self._get_root_transcript_tpm()
            cond_cols = [c for c in root_df.columns if c == cond or c.startswith(f"{cond}__")]
            if not cond_cols:
                continue
            tpm_values.append(root_df[cond_cols].mean(axis=1))

        if not tpm_values:
            return pd.Series(dtype=float)
        stacked = pd.concat(tpm_values, axis=1)
        return stacked.mean(axis=1)

    def _get_root_transcript_tpm(self) -> pd.DataFrame:
        if self._root_transcript_tpm is not None:
            return self._root_transcript_tpm
        pattern = "*transcript_model_grouped_tpm.tsv" if self.ref_only else "*transcript_grouped_tpm.tsv"
        files = list(self.output_dir.glob(pattern))
        if not files:
            return pd.DataFrame()
        df = self._read_tpm(files[0])
        self._root_transcript_tpm = df
        return df
