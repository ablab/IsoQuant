import logging
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rpy2 import robjects
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
from rpy2.rinterface_lib import callbacks


class GSEAAnalysis:
    def __init__(self, output_path: Path):
        """
        Initialize GSEA analysis.

        Args:
            output_path: Path to save GSEA results
        """
        self.output_path = Path(output_path) / "gsea_results"
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # Configure R to be quiet

        
        def quiet_cb(x):
            pass
        
        callbacks.logger.setLevel(logging.WARNING)
        callbacks.consolewrite_print = quiet_cb
        callbacks.consolewrite_warnerror = quiet_cb

    def run_gsea_analysis(self, results: pd.DataFrame, target_label: str) -> None:
        """
        Run GSEA analysis using DESeq2 stat value as ranking metric.
        Creates visualizations for top enriched pathways in each GO category.
        
        Args:
            results: DataFrame containing DESeq2 results
            target_label: Label indicating the comparison being made
        """
        if results is None or results.empty:
            logging.error("No DESeq2 results provided for GSEA analysis")
            return
            
        logging.info("Starting GSEA analysis...")
        logging.debug(f"Full DE results shape: {results.shape}")
        logging.debug(f"DE results columns: {results.columns.tolist()}")

        # Filter for significant DE genes
        sig_genes = results.dropna(subset=["padj"])
        logging.debug(f"After dropping genes with NaN padj: {sig_genes.shape}")
        sig_genes = sig_genes[sig_genes["padj"] < 0.05]
        logging.debug(f"Significantly DE genes (padj<0.05): {sig_genes.shape}")
        if sig_genes.empty:
            logging.info("No significantly DE genes found for GSEA.")
            return

        sig_genes = sig_genes[sig_genes["pvalue"] > 0]
        logging.debug(f"Significantly DE genes with pvalue>0: {sig_genes.shape}")
        if sig_genes.empty:
            logging.info("No genes with valid p-values for GSEA.")
            return

        # Use gene_name instead of symbol
        gene_symbols_final = sig_genes["gene_name"].values

        ranked_genes = pd.Series(
            sig_genes["stat"].values, index=gene_symbols_final
        ).dropna()
        ranked_genes = ranked_genes[~ranked_genes.index.duplicated(keep="first")]
        logging.debug(f"Final ranked genes count: {len(ranked_genes)}")

        if ranked_genes.empty:
            logging.info("No valid ranked genes after processing.")
            return

        # Save the ranked genes
        ranked_outfile = self.output_path / "ranked_genes.csv"
        ranked_genes_df = pd.DataFrame(
            {"gene": ranked_genes.index, "rank": ranked_genes.values}
        )
        ranked_genes_df.to_csv(ranked_outfile, index=False)
        logging.info(f"Ranked genes saved to {ranked_outfile}")

        # Import required R packages
        clusterProfiler = importr("clusterProfiler")
        r("library(org.Hs.eg.db)")

        with localconverter(robjects.default_converter + pandas2ri.converter):
            r_ranked_genes = pandas2ri.py2rpy(ranked_genes.sort_values(ascending=False))

        def plot_pathways(df: pd.DataFrame, direction: str, ont: str):
            if df.empty:
                logging.info(f"No {direction} pathways to plot.")
                return

            df["label"] = df["ID"] + ": " + df["Description"]
            df["-log10(p.adjust)"] = -np.log10(df["p.adjust"])
            values = df["-log10(p.adjust)"]
            norm = plt.Normalize(vmin=values.min(), vmax=values.max())
            cmap = plt.cm.get_cmap("viridis")
            colors_for_bars = [cmap(norm(v)) for v in values]

            plt.figure(figsize=(12, 8))
            plt.barh(
                df["label"].iloc[::-1],
                df["-log10(p.adjust)"].iloc[::-1],
                color=colors_for_bars[::-1],
            )
            plt.xlabel("-log10(adjusted p-value)")

            # Split target label into reference and target parts
            target_parts = target_label.split("_vs_")
            target_condition = target_parts[0]
            ref_condition = target_parts[1]

            # Create title based on direction
            if direction == "up":
                condition_str = f"Pathways enriched in {target_condition}\nvs {ref_condition} - {ont}"
            else:
                condition_str = f"Pathways enriched in {ref_condition}\nvs {target_condition} - {ont}"

            plt.title(condition_str, fontsize=10)
            plt.tight_layout()
            plot_path = self.output_path / f"GSEA_top_pathways_{direction}_{ont}.png"
            plt.savefig(plot_path)
            plt.close()
            logging.info(f"GSEA {direction} pathways plot saved to {plot_path}")

        # Run GO analysis for each ontology
        ontologies = ["BP", "MF", "CC"]
        for ont in ontologies:
            logging.debug(f"Running gseGO for {ont}...")

            gsea_res = clusterProfiler.gseGO(
                geneList=r_ranked_genes,
                OrgDb="org.Hs.eg.db",
                keyType="SYMBOL",
                ont=ont,
                minGSSize=5,
                maxGSSize=1000,
                pvalueCutoff=1,
                verbose=True,
                nPermSimple=10000,
            )

            gsea_table = r("data.frame")(gsea_res)
            with localconverter(robjects.default_converter + pandas2ri.converter):
                gsea_df = pandas2ri.rpy2py(gsea_table)

            # Log detailed results
            logging.debug(f"GSEA results for {ont}:")
            logging.debug(f"  Total pathways tested: {len(gsea_df)}")
            if not gsea_df.empty:
                logging.debug(
                    f"  P-value range: {gsea_df['pvalue'].min():.2e} - {gsea_df['pvalue'].max():.2e}"
                )
                logging.debug(
                    f"  Adjusted p-value range: {gsea_df['p.adjust'].min():.2e} - {gsea_df['p.adjust'].max():.2e}"
                )
                logging.debug(
                    f"  NES range: {gsea_df['NES'].min():.2f} - {gsea_df['NES'].max():.2f}"
                )
                logging.debug(
                    f"  Pathways with adj.P<0.1: {len(gsea_df[gsea_df['p.adjust'] < 0.1])}"
                )
                logging.debug(
                    f"  Pathways with adj.P<0.05: {len(gsea_df[gsea_df['p.adjust'] < 0.05])}"
                )

                # Save all results
                gsea_outfile = self.output_path / f"GSEA_results_{ont}.csv"
                gsea_df.to_csv(gsea_outfile, index=False)
                logging.info(f"Complete GSEA results for {ont} saved to {gsea_outfile}")

                # Plot significant pathways
                sig_gsea_df = gsea_df[
                    gsea_df["p.adjust"] < 0.05
                ].copy()  # Using 0.05 threshold
                if not sig_gsea_df.empty:
                    up_pathways = sig_gsea_df[sig_gsea_df["NES"] > 0].nsmallest(
                        10, "p.adjust"
                    )
                    down_pathways = sig_gsea_df[sig_gsea_df["NES"] < 0].nsmallest(
                        10, "p.adjust"
                    )

                    if not up_pathways.empty:
                        plot_pathways(up_pathways, "up", ont)
                    if not down_pathways.empty:
                        plot_pathways(down_pathways, "down", ont)
                else:
                    logging.info(f"No pathways with adj.P<0.05 found for {ont}")

        logging.info("GSEA analysis completed.")
