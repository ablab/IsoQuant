import os
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import logging
import pandas as pd
import matplotlib.patches as patches
import seaborn as sns


class PlotOutput:
    def __init__(
        self,
        updated_gene_dict,
        gene_names,
        gene_visualizations_dir,
        read_assignments_dir,
        reads_and_class=None,
        filter_transcripts=None,
        conditions=False,
        use_counts=False,
    ):
        self.updated_gene_dict = updated_gene_dict
        self.gene_names = gene_names
        self.gene_visualizations_dir = gene_visualizations_dir
        self.read_assignments_dir = read_assignments_dir
        self.reads_and_class = reads_and_class
        self.filter_transcripts = filter_transcripts
        self.conditions = conditions
        self.use_counts = use_counts

        # Ensure output directories exist
        if self.gene_visualizations_dir:
            os.makedirs(self.gene_visualizations_dir, exist_ok=True)
        os.makedirs(self.read_assignments_dir, exist_ok=True)

    def plot_transcript_map(self):
        """Plot transcript structure with different colors for reference and novel exons."""
        if not self.gene_visualizations_dir:
            logging.warning("No gene_visualizations_dir provided. Skipping transcript map plotting.")
            return

        for gene_name in self.gene_names:
            gene_data = {}
            for condition, genes in self.updated_gene_dict.items():
                if gene_name in genes:
                    gene_data = genes[gene_name]
                    break

            if not gene_data:
                logging.warning(f"Gene {gene_name} not found in the data.")
                continue

            # Get chromosome info and calculate buffer
            chromosome = gene_data.get("chromosome", "Unknown")
            start = gene_data.get("start", 0)
            end = gene_data.get("end", 0)
            
            # Calculate buffer (5% of total width)
            width = end - start
            buffer = width * 0.05
            plot_start = start - buffer
            plot_end = end + buffer

            plot_height = max(8, len(gene_data["transcripts"]) * 0.4)
            logging.debug(f"Creating transcript map for gene '{gene_name}' with {len(gene_data['transcripts'])} transcripts")

            # Collect all reference exon coordinates from reference transcripts
            reference_exons = set()
            for transcript_id, transcript_info in gene_data["transcripts"].items():
                if transcript_id.startswith("ENST"):
                    for exon in transcript_info["exons"]:
                        # Store exon coordinates as tuple for easy comparison
                        reference_exons.add((exon["start"], exon["end"]))
            
            logging.debug(f"Found {len(reference_exons)} reference exons for gene '{gene_name}'")

            fig, ax = plt.subplots(figsize=(12, plot_height))
            
            # Add legend handles
            legend_elements = [
                patches.Patch(facecolor='skyblue', label='Reference Exon'),
                patches.Patch(facecolor='red', alpha=0.6, label='Novel Exon')
            ]

            # Plot each transcript
            y_ticks = []
            y_labels = []
            for i, (transcript_id, transcript_info) in enumerate(gene_data["transcripts"].items()):
                # Plot direction marker
                direction_marker = ">" if gene_data["strand"] == "+" else "<"
                marker_pos = (
                    transcript_info["end"] + 100
                    if gene_data["strand"] == "+"
                    else transcript_info["start"] - 100
                )
                ax.plot(
                    marker_pos, i, marker=direction_marker, markersize=5, color="blue"
                )

                # Draw the line for the whole transcript
                ax.plot(
                    [transcript_info["start"], transcript_info["end"]],
                    [i, i],
                    color="grey",
                    linewidth=2,
                )

                # Exon blocks with color based on reference status
                for exon in transcript_info["exons"]:
                    exon_length = exon["end"] - exon["start"]
                    # Check if this exon's coordinates match any reference exon
                    is_reference_exon = (exon["start"], exon["end"]) in reference_exons
                    exon_color = "skyblue" if is_reference_exon else "red"
                    exon_alpha = 1.0 if is_reference_exon else 0.6
                    
                    ax.add_patch(
                        plt.Rectangle(
                            (exon["start"], i - 0.4),
                            exon_length,
                            0.8,
                            color=exon_color,
                            alpha=exon_alpha
                        )
                    )

                if not any((exon["start"], exon["end"]) in reference_exons for exon in transcript_info["exons"]):
                    logging.debug(f"Transcript {transcript_id} in gene {gene_name} contains all novel exons")

                # Store y-axis label information
                y_ticks.append(i)
                # Get transcript name with fallback options
                transcript_name = (transcript_info.get("name") or 
                                 transcript_info.get("transcript_id") or 
                                 transcript_id)
                y_labels.append(transcript_name)

            # Set up the plot formatting with just chromosome
            if self.filter_transcripts:
                title = f"Transcript Structure - {gene_name} (Chromosome {chromosome}) (Count > {self.filter_transcripts})"
            else:
                title = f"Transcript Structure - {gene_name} (Chromosome {chromosome})"
            
            ax.set_title(title, pad=20)  # Increase padding to move title up
            ax.set_xlabel("Chromosomal position")
            ax.set_ylabel("Transcripts")
            
            # Set y-axis ticks and labels
            ax.set_yticks(y_ticks)
            ax.set_yticklabels(y_labels)

            # Add legend in upper right corner
            ax.legend(handles=legend_elements, loc='upper right')

            # Set plot limits with buffer
            ax.set_xlim(plot_start, plot_end)
            ax.invert_yaxis()  # First transcript at the top

            # Add grid lines
            ax.grid(True, axis='y', linestyle='--', alpha=0.3)

            plt.tight_layout()
            plot_path = os.path.join(
                self.gene_visualizations_dir, f"{gene_name}_splicing.png"
            )
            plt.savefig(plot_path, bbox_inches='tight', dpi=300)
            plt.close(fig)
            logging.debug(f"Saved transcript map for gene '{gene_name}' at: {plot_path}")

   
    def plot_transcript_usage(self):
        """Visualize transcript usage for each gene in gene_names across different conditions."""
        if not self.gene_visualizations_dir:
            logging.warning("No gene_visualizations_dir provided. Skipping transcript usage plotting.")
            return

        for gene_name in self.gene_names:
            gene_data = {}
            for condition, genes in self.updated_gene_dict.items():
                if gene_name in genes:
                    gene_data[condition] = genes[gene_name]["transcripts"]

            if not gene_data:
                logging.warning(f"Gene {gene_name} not found in the data.")
                continue

            conditions = list(gene_data.keys())
            n_bars = len(conditions)

            fig, ax = plt.subplots(figsize=(12, 8))
            index = np.arange(n_bars)
            bar_width = 0.35
            opacity = 0.8
            max_transcripts = max(len(gene_data[condition]) for condition in conditions)
            colors = plt.cm.plasma(np.linspace(0, 1, num=max_transcripts))

            bottom_val = np.zeros(n_bars)
            for i, condition in enumerate(conditions):
                transcripts = gene_data[condition]
                for j, (transcript_id, transcript_info) in enumerate(transcripts.items()):
                    color = colors[j % len(colors)]
                    value = transcript_info["value"]
                    # Get transcript name with fallback options
                    transcript_name = (transcript_info.get("name") or 
                                     transcript_info.get("transcript_id") or 
                                     transcript_id)
                    ax.bar(
                        i,
                        float(value),
                        bar_width,
                        bottom=bottom_val[i],
                        alpha=opacity,
                        color=color,
                        label=transcript_name if i == 0 else "",
                    )
                    bottom_val[i] += float(value)

            ax.set_xlabel("Sample Type")
            ax.set_ylabel("Transcript Usage (TPM)")
            ax.set_title(f"Transcript Usage for {gene_name} by Sample Type")
            ax.set_xticks(index)
            ax.set_xticklabels(conditions)
            ax.legend(
                title="Transcript IDs",
                bbox_to_anchor=(1.05, 1),
                loc="upper left",
                fontsize=8,
            )

            plt.tight_layout()
            plot_path = os.path.join(
                self.gene_visualizations_dir,
                f"{gene_name}_transcript_usage_by_sample_type.png",
            )
            plt.savefig(plot_path)
            plt.close(fig)

    def make_pie_charts(self):
        """
        Create pie charts for transcript alignment classifications and read assignment consistency.
        Handles both combined and separate sample data structures.
        """

        titles = ["Transcript Alignment Classifications", "Read Assignment Consistency"]

        for title, data in zip(titles, self.reads_and_class):
            if isinstance(data, dict):
                if any(isinstance(v, dict) for v in data.values()):
                    # Separate 'Mutants' and 'WildType' case
                    for sample_name, sample_data in data.items():
                        self._create_pie_chart(f"{title} - {sample_name}", sample_data)
                else:
                    # Combined data case
                    self._create_pie_chart(title, data)
            else:
                print(f"Skipping unexpected data type for {title}: {type(data)}")

    def _create_pie_chart(self, title, data):
        """
        Helper method to create a single pie chart.
        """
        labels = list(data.keys())
        sizes = list(data.values())
        total = sum(sizes)

        # Generate a file-friendly title
        file_title = title.lower().replace(" ", "_").replace("-", "_")

        plt.figure(figsize=(12, 8))
        wedges, texts, autotexts = plt.pie(
            sizes,
            labels=labels,
            autopct=lambda pct: f"{pct:.1f}%\n({int(pct/100.*total):d})",
            startangle=140,
            textprops=dict(color="w"),
        )
        plt.setp(autotexts, size=8, weight="bold")
        plt.setp(texts, size=7)

        plt.axis("equal")  # Equal aspect ratio ensures that pie is drawn as a circle.
        plt.title(f"{title}\nTotal: {total}")

        plt.legend(
            wedges,
            labels,
            title="Categories",
            loc="center left",
            bbox_to_anchor=(1, 0, 0.5, 1),
            fontsize=8,
        )
        # Save pie charts in the read_assignments directory
        plot_path = os.path.join(
            self.read_assignments_dir, f"{file_title}_pie_chart.png"
        )
        plt.savefig(plot_path, bbox_inches="tight", dpi=300)
        plt.close()


class ExpressionVisualizer:
    def __init__(self, output_path):
        """Initialize with output path for plots."""
        self.output_path = Path(output_path)
        self.output_path.mkdir(parents=True, exist_ok=True)
        self.logger = logging.getLogger(__name__)  # Logger for this class
        # Suppress matplotlib font debug messages
        logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)

    def create_volcano_plot(
        self,
        df: pd.DataFrame,
        target_label: str,
        reference_label: str,
        padj_threshold: float = 0.05,
        lfc_threshold: float = 1,
        top_n: int = 10,
        feature_type: str = "genes",
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
            if feature_type == "genes":
                symbol = row["gene_name"] if pd.notnull(row["gene_name"]) else row["feature_id"]
            elif feature_type == "transcripts":
                symbol = row["transcript_symbol"] if pd.notnull(row["transcript_symbol"]) else row["feature_id"]
            else: # Fallback to feature_id if feature_type is not recognized
                symbol = row["feature_id"]
            plt.text(
                row["log2FoldChange"],
                row["-log10(padj)"],
                symbol,
                fontsize=8,
                ha="center",
                va="bottom",
            )

        plt.tight_layout()
        plot_path = (
            self.output_path / f"volcano_plot_{feature_type}.png"
        )  # Modified line
        plt.savefig(str(plot_path))
        plt.close()
        logging.info(f"Volcano plot saved to {plot_path}")

    def create_ma_plot(
        self,
        df: pd.DataFrame,
        target_label: str,
        reference_label: str,
        feature_type: str = "genes",
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
        plot_path = self.output_path / f"ma_plot_{feature_type}.png"  # Modified line
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
        Create and save analysis summary with correct filtering criteria reporting.

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

        # Incorporate feature_type into the summary filename
        summary_filename = f"analysis_summary_{feature_type}.txt"
        summary_path = self.output_path / summary_filename

        with summary_path.open("w") as f:
            f.write(f"Analysis Summary: {target_label} vs {reference_label}\n")
            f.write("================================\n")
            
            # Different filtering description based on feature type
            if feature_type == "genes":
                f.write(
                    f"{feature_type.capitalize()} after filtering "
                    f"(mean count >= {min_count} in either condition group): {total_features}\n"
                )
            else:  # transcripts
                f.write(
                    f"{feature_type.capitalize()} after filtering "
                    f"(count >= {min_count} in at least half of all samples): {total_features}\n"
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
            self.create_volcano_plot(
                results, target_label, reference_label, feature_type=feature_type
            )
            self.create_ma_plot(
                results, target_label, reference_label, feature_type=feature_type
            )
            self.create_summary(
                results,
                target_label,
                reference_label,
                min_count,
                feature_type=feature_type,
            )
        except Exception as e:
            logging.exception("Failed to create visualizations")
            raise

    
    def plot_pca(self, pca_df: pd.DataFrame, title: str, output_prefix: str) -> Path:
        """Plot PCA scatter plot."""
        plt.figure(figsize=(8, 6))
        
        # Extract variance info from title for axis labels only
        pc1_var = title.split("PC1 (")[1].split("%)")[0]
        pc2_var = title.split("PC2 (")[1].split("%)")[0]
        
        # Get clean title without PCs and variance - using string literal instead of \n
        base_title = title.split(' Level PCA: ')[0]
        comparison = title.split(': ')[1].split('PC1')[0].strip()
        clean_title = f"{base_title} Level PCA: {comparison}"
        
        # Update group labels in the DataFrame
        condition_mapping = {'Target': title.split(": ")[1].split(" vs ")[0],
                            'Reference': title.split(" vs ")[1].split("PC1")[0].strip()}
        pca_df['group'] = pca_df['group'].map(condition_mapping)
        
        # Create plot with updated labels
        sns.scatterplot(x='PC1', y='PC2', hue='group', data=pca_df, s=100)
        plt.xlabel(f'PC1 ({pc1_var}%)')
        plt.ylabel(f'PC2 ({pc2_var}%)')
        plt.title(clean_title)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.tight_layout()

        output_path = self.output_path / f"{output_prefix}_pca.png"
        plt.savefig(output_path)
        plt.close()
        return output_path

