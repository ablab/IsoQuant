"""
Visualization and summary module for differential expression analysis results.
"""

from pathlib import Path
import logging
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List


class ExpressionVisualizer:
    def __init__(self, output_path: Path):
        """
        Initialize visualizer with output directory.

        Args:
            output_path: Path to output directory
        """
        self.output_path = Path(output_path)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def create_volcano_plot(
        self,
        df: pd.DataFrame,
        target_label: str,
        reference_label: str,
        padj_threshold: float = 0.05,
        lfc_threshold: float = 1,
        top_n: int = 10,
    ) -> None:
        """Create volcano plot from differential expression results."""
        plt.figure(figsize=(10, 8))

        # Prepare data
        df["padj"] = df["padj"].replace(0, 1e-300)
        df = df[df["padj"] > 0]
        df = df.copy()  # Create a copy to avoid the warning
        df.loc[:, "-log10(padj)"] = -np.log10(df["padj"])

        # Define significant genes
        significant = (df["padj"] < padj_threshold) & (
            abs(df["log2FoldChange"]) > lfc_threshold
        )
        up_regulated = significant & (df["log2FoldChange"] > lfc_threshold)
        down_regulated = significant & (df["log2FoldChange"] < -lfc_threshold)

        # Plot points
        plt.scatter(
            df.loc[~significant, "log2FoldChange"],
            df.loc[~significant, "-log10(padj)"],
            color="grey",
            alpha=0.5,
            label="Not Significant",
        )
        plt.scatter(
            df.loc[up_regulated, "log2FoldChange"],
            df.loc[up_regulated, "-log10(padj)"],
            color="red",
            alpha=0.7,
            label=f"Up-regulated in ({target_label})",
        )
        plt.scatter(
            df.loc[down_regulated, "log2FoldChange"],
            df.loc[down_regulated, "-log10(padj)"],
            color="blue",
            alpha=0.7,
            label=f"Up-regulated in ({reference_label})",
        )

        # Add threshold lines and labels
        plt.axhline(-np.log10(padj_threshold), color="grey", linestyle="--")
        plt.axvline(lfc_threshold, color="grey", linestyle="--")
        plt.axvline(-lfc_threshold, color="grey", linestyle="--")

        plt.xlabel("log2 Fold Change")
        plt.ylabel("-log10(adjusted p-value)")
        plt.title(f"Volcano Plot: {target_label} vs {reference_label}")
        plt.legend()

        # Add labels for top significant features
        sig_df = df.loc[significant].nsmallest(top_n, "padj")
        for _, row in sig_df.iterrows():
            symbol = row["symbol"] if pd.notnull(row["symbol"]) else row["feature_id"]
            plt.text(
                row["log2FoldChange"],
                row["-log10(padj)"],
                symbol,
                fontsize=8,
                ha="center",
                va="bottom",
            )

        plt.tight_layout()
        plot_path = self.output_path / "volcano_plot.png"
        plt.savefig(str(plot_path))
        plt.close()
        logging.info(f"Volcano plot saved to {plot_path}")

    def create_ma_plot(
        self, df: pd.DataFrame, target_label: str, reference_label: str
    ) -> None:
        """Create MA plot from differential expression results."""
        plt.figure(figsize=(10, 8))

        # Prepare data
        df = df[df["baseMean"] > 0]
        df["log10(baseMean)"] = np.log10(df["baseMean"])

        # Create plot
        plt.scatter(
            df["log10(baseMean)"], df["log2FoldChange"], alpha=0.5, color="grey"
        )
        plt.axhline(y=0, color="red", linestyle="--")

        plt.xlabel("log10(Base Mean)")
        plt.ylabel("log2 Fold Change")
        plt.title(f"MA Plot: {target_label} vs {reference_label}")

        plt.tight_layout()
        plot_path = self.output_path / "ma_plot.png"
        plt.savefig(str(plot_path))
        plt.close()
        logging.info(f"MA plot saved to {plot_path}")

    def create_summary(
        self,
        res_df: pd.DataFrame,
        target_label: str,
        reference_label: str,
        min_count: int,
        feature_type: str,
    ) -> None:
        """
        Create and save analysis summary.

        Args:
            res_df: Results DataFrame
            target_label: Target condition label
            reference_label: Reference condition label
            min_count: Minimum count threshold used in filtering
            feature_type: Type of features analyzed ("genes" or "transcripts")
        """
        total_features = len(res_df)
        sig_features = (
            (res_df["padj"] < 0.05) & (res_df["log2FoldChange"].abs() > 1)
        ).sum()
        up_regulated = ((res_df["padj"] < 0.05) & (res_df["log2FoldChange"] > 1)).sum()
        down_regulated = (
            (res_df["padj"] < 0.05) & (res_df["log2FoldChange"] < -1)
        ).sum()

        summary_path = self.output_path / "analysis_summary.txt"
        with summary_path.open("w") as f:
            f.write(f"Analysis Summary: {target_label} vs {reference_label}\n")
            f.write("================================\n")
            f.write(
                f"{feature_type.capitalize()} after filtering "
                f"(mean count >= {min_count} in both groups): {total_features}\n"
            )
            f.write(f"Significantly differential {feature_type}: {sig_features}\n")
            f.write(f"Up-regulated {feature_type}: {up_regulated}\n")
            f.write(f"Down-regulated {feature_type}: {down_regulated}\n")
        logging.info(f"Analysis summary saved to {summary_path}")

    def visualize_results(
        self,
        results: pd.DataFrame,
        target_label: str,
        reference_label: str,
        min_count: int,
        feature_type: str,
    ) -> None:
        """
        Create all visualizations and summary for the analysis results.

        Args:
            results: DataFrame containing differential expression results
            target_label: Target condition label
            reference_label: Reference condition label
            min_count: Minimum count threshold used in filtering
            feature_type: Type of features analyzed ("genes" or "transcripts")
        """
        try:
            self.create_volcano_plot(results, target_label, reference_label)
            self.create_ma_plot(results, target_label, reference_label)
            self.create_summary(
                results, target_label, reference_label, min_count, feature_type
            )
        except Exception as e:
            logging.exception("Failed to create visualizations")
            raise
