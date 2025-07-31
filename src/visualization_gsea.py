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

        def plot_pathways(up_df: pd.DataFrame, down_df: pd.DataFrame, ont: str):
            if up_df.empty and down_df.empty:
                logging.info(f"No pathways to plot for {ont}.")
                return
                
            # Process DataFrame if not empty
            if not up_df.empty:
                # Remove GO IDs from labels, keeping only the description
                up_df["label"] = up_df["Description"]
                up_df["-log10(p.adjust)"] = -np.log10(up_df["p.adjust"])
                up_df = up_df.sort_values(by="NES", ascending=False)
            
            if not down_df.empty:
                # Remove GO IDs from labels, keeping only the description
                down_df["label"] = down_df["Description"]
                down_df["-log10(p.adjust)"] = -np.log10(down_df["p.adjust"])
                down_df = down_df.sort_values(by="NES", ascending=True)
            
            # Find the global min and max for -log10(p.adjust) for consistent coloring
            all_pvals = []
            if not up_df.empty:
                all_pvals.extend(up_df["-log10(p.adjust)"].tolist())
            if not down_df.empty:
                all_pvals.extend(down_df["-log10(p.adjust)"].tolist())
                
            if not all_pvals:
                return  # Skip if no values
                
            # Get global min and max p-values
            global_vmin = min(all_pvals)
            global_vmax = max(all_pvals)
            
            # Adjust the maximum value to prevent saturation of highly significant pathways
            # Use either actual max or a higher percentile value, whichever is higher
            # This prevents all highly significant pathways from appearing with the same color
            if len(all_pvals) > 1:
                # Calculate 90th percentile of p-values
                percentile_90 = np.percentile(all_pvals, 90)
                
                # If max is much larger than 90th percentile, use an intermediate value
                if global_vmax > 2 * percentile_90:
                    adjusted_vmax = percentile_90 + (global_vmax - percentile_90) / 3
                    # But ensure we don't lower the max too much
                    global_vmax = max(adjusted_vmax, global_vmax * 0.7)
                    
                # Log the adjustment for debugging
                logging.debug(f"P-value color scale: original max={max(all_pvals):.2f}, adjusted max={global_vmax:.2f}")
            
            # Create a consistent color normalization across both plots
            norm = plt.Normalize(vmin=global_vmin, vmax=global_vmax)
            cmap = plt.cm.get_cmap("viridis")
            
            # Split target label into reference and target parts
            target_parts = target_label.split("_vs_")
            target_condition = target_parts[0]
            ref_condition = target_parts[1]
            
            # Plot UP-regulated pathways
            if not up_df.empty:
                up_values = up_df["-log10(p.adjust)"]
                up_colors = [cmap(norm(v)) for v in up_values]
                
                # Adjust figure size - no need for extra space for legend
                plt.figure(figsize=(12, 10))
                
                # Create horizontal bar plot
                bars = plt.barh(
                    up_df["label"].iloc[::-1],
                    up_df["NES"].iloc[::-1],
                    color=up_colors[::-1],
                )
                
                # Remove the p-value text labels
                
                # Add colorbar with the global scale
                sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
                sm.set_array([])
                cbar = plt.colorbar(sm)
                cbar.set_label("-log10(adjusted p-value)", fontsize=12)
                
                plt.xlabel("Normalized Enrichment Score (NES)", fontsize=12)
                condition_str = f"Pathways enriched in {target_condition}\nvs {ref_condition} - {ont}"
                plt.title(condition_str, fontsize=14)
                
                # Ensure y-axis labels are fully visible
                plt.tight_layout()
                plt.subplots_adjust(left=0.3)  # Add more space on the left for labels
                
                plot_path = self.output_path / f"GSEA_top_pathways_up_{ont}.pdf"
                plt.savefig(plot_path, format="pdf", bbox_inches="tight", dpi=600)
                plt.close()
                logging.info(f"GSEA up-regulated pathways plot saved to {plot_path} with high resolution")
            
            # Plot DOWN-regulated pathways
            if not down_df.empty:
                down_values = down_df["-log10(p.adjust)"]
                down_colors = [cmap(norm(v)) for v in down_values]
                
                # Adjust figure size - no need for extra space for legend
                plt.figure(figsize=(12, 10))
                
                # Create horizontal bar plot
                bars = plt.barh(
                    down_df["label"].iloc[::-1],
                    down_df["NES"].abs().iloc[::-1],  # Use absolute NES for down-regulated
                    color=down_colors[::-1],
                )
                
                # Add colorbar with the global scale
                sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
                sm.set_array([])
                cbar = plt.colorbar(sm)
                cbar.set_label("-log10(adjusted p-value)", fontsize=12)
                
                plt.xlabel("Absolute Normalized Enrichment Score (|NES|)", fontsize=12)
                condition_str = f"Pathways enriched in {ref_condition}\nvs {target_condition} - {ont}"
                plt.title(condition_str, fontsize=14)
                
                # Ensure y-axis labels are fully visible
                plt.tight_layout()
                plt.subplots_adjust(left=0.3)  # Add more space on the left for labels
                
                plot_path = self.output_path / f"GSEA_top_pathways_down_{ont}.pdf"
                plt.savefig(plot_path, format="pdf", bbox_inches="tight", dpi=600)
                plt.close()
                logging.info(f"GSEA down-regulated pathways plot saved to {plot_path} with high resolution")

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

                # Process significant pathways
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

                    # Use consistent color scales across both plots
                    plot_pathways(up_pathways, down_pathways, ont)
                else:
                    logging.info(f"No pathways with adj.P<0.05 found for {ont}")

        logging.info("GSEA analysis completed.")
