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
from matplotlib.patches import Patch


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

        # Don't filter for significant genes - use ALL genes with valid statistics
        # Just remove NaN values
        valid_genes = results.dropna(subset=["stat", "gene_name"])
        logging.debug(f"Genes with valid statistics: {valid_genes.shape}")
        
        if valid_genes.empty:
            logging.info("No genes with valid statistics found for GSEA.")
            return

        # Use gene_name instead of symbol
        gene_symbols = valid_genes["gene_name"].values

        # Create ranked list using ALL genes (not just significant ones)
        ranked_genes = pd.Series(
            valid_genes["stat"].values, index=gene_symbols
        ).dropna()
        ranked_genes = ranked_genes[~ranked_genes.index.duplicated(keep="first")]
        logging.debug(f"Final ranked genes count: {len(ranked_genes)}")

        if ranked_genes.empty or len(ranked_genes) < 50:  # Ensure we have enough genes
            logging.info(f"Not enough valid ranked genes for GSEA: {len(ranked_genes)}")
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
            
            # Sort by NES value - for up-regulated, highest NES first; for down-regulated, lowest NES first
            if direction == "up":
                df = df.sort_values(by="NES", ascending=False)
                plot_values = df["NES"]  # Use NES directly for up-regulated
            else:  # down
                df = df.sort_values(by="NES", ascending=True)
                plot_values = df["NES"].abs()  # Use absolute NES for down-regulated
            
            # Use NES for bar length but -log10(p.adjust) for color
            values = df["-log10(p.adjust)"]
            
            # Use the data's own range for each direction
            vmin = values.min()
            vmax = values.max()
            
            norm = plt.Normalize(vmin=vmin, vmax=vmax)
            cmap = plt.cm.get_cmap("viridis")
            colors_for_bars = [cmap(norm(v)) for v in values]

            plt.figure(figsize=(14, 8))  # Wider figure to accommodate legend
            
            # Use NES for bar length (absolute value for down-regulated)
            plt.barh(
                df["label"].iloc[::-1],
                plot_values.iloc[::-1],  # Use appropriate values based on direction
                color=colors_for_bars[::-1],  # Still color by significance
            )
            
            # Add a colorbar to show the significance scale
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm)
            cbar.set_label("-log10(adjusted p-value)")
            
            # Add legend explaining the visualization
            legend_elements = [
                Patch(facecolor='gray', alpha=0.5, 
                      label='Bar length: Normalized Enrichment Score (NES)'),
                Patch(facecolor=cmap(0.25), alpha=0.8, 
                      label='Bar color: Statistical significance'),
            ]
            # Move legend much further to the right
            plt.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.25, 1))
            
            # Set x-axis label based on direction
            if direction == "up":
                plt.xlabel("Normalized Enrichment Score (NES)")
            else:
                plt.xlabel("Absolute Normalized Enrichment Score (|NES|)")

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
            
            # Adjust layout to make room for legend
            plt.tight_layout()
            # Save with extra space for the legend
            plot_path = self.output_path / f"GSEA_top_pathways_{direction}_{ont}.pdf"
            plt.savefig(plot_path, format="pdf", bbox_inches="tight", dpi=300)
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
                        15, "p.adjust"
                    )
                    down_pathways = sig_gsea_df[sig_gsea_df["NES"] < 0].nsmallest(
                        15, "p.adjust"
                    )

                    # Use separate color scales for each direction
                    if not up_pathways.empty:
                        plot_pathways(up_pathways, "up", ont)
                    if not down_pathways.empty:
                        plot_pathways(down_pathways, "down", ont)
                else:
                    logging.info(f"No pathways with adj.P<0.05 found for {ont}")

        logging.info("GSEA analysis completed.")
