import os
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import logging
import pandas as pd
import matplotlib.patches as patches
import seaborn as sns
from typing import List
from matplotlib.colors import Normalize
import matplotlib.cm as cm

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
        ref_only=False,
        ref_conditions=None,
        target_conditions=None,
    ):
        self.updated_gene_dict = updated_gene_dict
        self.gene_names = gene_names
        self.gene_visualizations_dir = gene_visualizations_dir
        self.read_assignments_dir = read_assignments_dir
        self.reads_and_class = reads_and_class
        self.conditions = conditions
        self.ref_only = ref_only
        self.min_tpm_threshold = filter_transcripts
        
        # Explicitly set reference and target conditions
        self.ref_conditions = ref_conditions if ref_conditions else []
        self.target_conditions = target_conditions if target_conditions else []
        
        # Log conditions for debugging
        if self.ref_conditions and self.target_conditions:
            logging.info(f"Filtering plots to include only ref conditions: {self.ref_conditions} and target conditions: {self.target_conditions}")
        else:
            logging.warning("No ref_conditions or target_conditions set, filtering may not work correctly")
        
        # Log TPM threshold if set
        if self.min_tpm_threshold:
            logging.info(f"Filtering transcripts with TPM value < {self.min_tpm_threshold}")

        # Ensure output directories exist
        if self.gene_visualizations_dir:
            os.makedirs(self.gene_visualizations_dir, exist_ok=True)
        os.makedirs(self.read_assignments_dir, exist_ok=True)

    def plot_transcript_map(self):
        """Plot transcript structure with different colors for reference and novel exons."""
        if not self.gene_visualizations_dir:
            logging.warning("No gene_visualizations_dir provided. Skipping transcript map plotting.")
            return

        # Check if reference and target conditions are defined
        has_specific_conditions = (hasattr(self, 'ref_conditions') and hasattr(self, 'target_conditions') and 
                                 self.ref_conditions and self.target_conditions)
        
        if has_specific_conditions:
            logging.info(f"Filtering transcript map to include only ref conditions: {self.ref_conditions} and target conditions: {self.target_conditions}")
            # Define all allowed conditions
            allowed_conditions = set(self.ref_conditions + self.target_conditions)

        for gene_name_or_id in self.gene_names:  # gene_names list contains gene names (symbols)
            gene_data = None  # Initialize gene_data to None
            found_condition = None  # Track which condition we found the gene in
            
            # First pass: Try to find the gene in allowed conditions only
            if has_specific_conditions:
                # Search only in allowed conditions
                for condition in allowed_conditions:
                    if condition not in self.updated_gene_dict:
                        continue
                        
                    genes = self.updated_gene_dict[condition]
                    for gene_id, gene_info in genes.items():
                        if "name" in gene_info and gene_info["name"] == gene_name_or_id.upper():  # Compare gene names (case-insensitive)
                            gene_data = gene_info
                            found_condition = condition
                            break
                    if gene_data:
                        break  # Found gene, stop searching
            
            # Second pass: If not found and we're allowing fallback, try all conditions
            if not gene_data:
                for condition, genes in self.updated_gene_dict.items():
                    # Skip conditions we already checked if using specific conditions
                    if has_specific_conditions and condition in allowed_conditions:
                        continue
                        
                    for gene_id, gene_info in genes.items():
                        if "name" in gene_info and gene_info["name"] == gene_name_or_id.upper():
                            gene_data = gene_info
                            found_condition = condition
                            break
                    if gene_data:
                        break  # Found gene, stop searching

            if gene_data:
                if has_specific_conditions and found_condition in allowed_conditions:
                    logging.debug(f"Gene {gene_name_or_id} found in prioritized condition: {found_condition}")
                else:
                    logging.debug(f"Gene {gene_name_or_id} found in fallback condition: {found_condition}")
            else:
                logging.warning(f"Gene {gene_name_or_id} not found in any condition.")
                continue  # Skip to the next gene if not found

            # Get chromosome info and calculate buffer
            chromosome = gene_data.get("chromosome", "Unknown")
            start = gene_data.get("start", 0)
            end = gene_data.get("end", 0)
            
            # Calculate buffer (5% of total width)
            width = end - start
            buffer = width * 0.05
            plot_start = start - buffer
            plot_end = end + buffer

            # NEW APPROACH: If we have ref/target conditions AND TPM filtering,
            # we need to consider transcript expression across ALL relevant conditions
            if has_specific_conditions and self.min_tpm_threshold is not None:
                # First, collect the max TPM for each transcript across all ref/target conditions
                transcript_max_tpm = {}
                
                # Collect all transcripts from the current condition first
                for transcript_id, transcript_info in gene_data["transcripts"].items():
                    value = float(transcript_info.get("value", 0))
                    transcript_max_tpm[transcript_id] = value
                
                # Check other ref/target conditions for the same gene to find max TPM values
                for condition in allowed_conditions:
                    if condition == found_condition or condition not in self.updated_gene_dict:
                        continue  # Skip the condition we already processed
                        
                    genes = self.updated_gene_dict[condition]
                    for gene_id, gene_info in genes.items():
                        # Check if this is the same gene in another condition
                        if "name" in gene_info and gene_info["name"] == gene_name_or_id.upper():
                            # Found the same gene in another condition, check transcript TPM values
                            for transcript_id, transcript_info in gene_info["transcripts"].items():
                                value = float(transcript_info.get("value", 0))
                                # Update max TPM if higher in this condition
                                if transcript_id in transcript_max_tpm:
                                    transcript_max_tpm[transcript_id] = max(transcript_max_tpm[transcript_id], value)
                                else:
                                    transcript_max_tpm[transcript_id] = value
                            break  # Found the gene in this condition, no need to check other genes
                
                # Now filter transcripts based on their max TPM across ref/target conditions
                filtered_transcripts = {}
                for transcript_id, transcript_info in gene_data["transcripts"].items():
                    max_tpm = transcript_max_tpm.get(transcript_id, 0)
                    if max_tpm >= self.min_tpm_threshold:
                        filtered_transcripts[transcript_id] = transcript_info
                
                # Log filtering results
                total_transcripts = len(gene_data["transcripts"])
                filtered_count = len(filtered_transcripts)
                logging.debug(f"Cross-condition TPM filtering: {filtered_count} of {total_transcripts} transcripts have TPM >= {self.min_tpm_threshold} in any ref/target condition for gene {gene_name_or_id}")
            else:
                # Original filtering approach for single condition
                filtered_transcripts = {}
                total_transcripts = len(gene_data["transcripts"])
                filtered_count = 0
                
                for transcript_id, transcript_info in gene_data["transcripts"].items():
                    # Apply TPM filtering if threshold is set
                    if self.min_tpm_threshold is not None:
                        value = float(transcript_info.get("value", 0))
                        if value >= self.min_tpm_threshold:
                            filtered_transcripts[transcript_id] = transcript_info
                            filtered_count += 1
                    else:
                        # No filtering, include all transcripts
                        filtered_transcripts[transcript_id] = transcript_info
                
                if self.min_tpm_threshold is not None:
                    logging.debug(f"Single-condition TPM filtering: {filtered_count} of {total_transcripts} transcripts have TPM >= {self.min_tpm_threshold} for gene {gene_name_or_id}")
            
            # Skip plotting if no transcripts pass the filter
            if not filtered_transcripts:
                logging.warning(f"No transcripts for gene {gene_name_or_id} pass the TPM threshold of {self.min_tpm_threshold}. Skipping plot.")
                continue
                
            # Calculate plot height based on number of filtered transcripts
            num_transcripts = len(filtered_transcripts)
            plot_height = max(8, num_transcripts * 0.4)
            #logging.debug(f"Creating transcript map for gene '{gene_name_or_id}' with {num_transcripts} transcripts from {found_condition}")

            fig, ax = plt.subplots(figsize=(12, plot_height))
            
            # Add legend handles
            legend_elements = [
                patches.Patch(facecolor='skyblue', label='Exon'),
            ]
            if not self.ref_only:
                legend_elements.append(patches.Patch(facecolor='red', alpha=0.6, label='Novel Exon'))

            # Plot each transcript
            y_ticks = []
            y_labels = []
            for i, (transcript_id, transcript_info) in enumerate(filtered_transcripts.items()):
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
                    if self.ref_only: # Check ref_only flag
                        exon_color = "skyblue" # If ref_only, always treat as reference
                        exon_alpha = 1.0
                    else:
                        is_reference_exon = exon["exon_id"].startswith("ENSE") # Original logic
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

                if not any(exon["exon_id"].startswith("ENSE") for exon in transcript_info["exons"]):
                    logging.debug(f"Transcript {transcript_id} in gene {gene_name_or_id} contains NO reference exons (based on ENSEMBL IDs)")
                    #log the exon_ids
                    #logging.debug(f"Exon IDs: {[exon['exon_id'] for exon in transcript_info['exons']]}")
                else:
                    #logging.debug(f"Transcript {transcript_id} in gene {gene_name_or_id} contains at least one reference exon (based on ENSEMBL IDs)")
                    pass  # Added explicit pass statement for the empty block
                
                # Store y-axis label information
                y_ticks.append(i)
                # Get transcript name with fallback options
                transcript_name = (transcript_info.get("name") or 
                                  transcript_info.get("transcript_id") or 
                                  transcript_id)
                value = float(transcript_info.get("value", 0))
                
                # If cross-condition filtering was used, show the max TPM value in the label
                if has_specific_conditions and self.min_tpm_threshold is not None:
                    max_tpm = transcript_max_tpm.get(transcript_id, 0)
                    y_labels.append(f"{transcript_name}")
                else:
                    y_labels.append(f"{transcript_name}")

            # Set up the plot formatting with just chromosome
            gene_display_name = gene_data.get("name", gene_name_or_id) # Fallback to ID if no name
            
            # Update title to include TPM threshold if applied
            if self.min_tpm_threshold is not None:
                if has_specific_conditions:
                    title = f"Transcript Structure - {gene_display_name} (Chromosome {chromosome}) (TPM >= {self.min_tpm_threshold} in any ref/target condition)"
                else:
                    title = f"Transcript Structure - {gene_display_name} (Chromosome {chromosome}) (TPM >= {self.min_tpm_threshold})"
            else:
                title = f"Transcript Structure - {gene_display_name} (Chromosome {chromosome})"
            
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
                self.gene_visualizations_dir, f"{gene_name_or_id}_splicing.pdf" # Changed from .png to .pdf
            )
            plt.savefig(plot_path, bbox_inches='tight', dpi=300)
            plt.close(fig)
            logging.debug(f"Saved transcript map for gene '{gene_name_or_id}' at: {plot_path}")

   
    def plot_transcript_usage(self):
        """Visualize transcript usage for each gene in gene_names across different conditions."""
        if not self.gene_visualizations_dir:
            logging.warning("No gene_visualizations_dir provided. Skipping transcript usage plotting.")
            return
        
        # Add this section near the beginning of the method
        logging.info("=== SPECIAL DEBUG FOR YBX1 TRANSCRIPTS ===")
        for condition in self.updated_gene_dict:
            for gene_id, gene_info in self.updated_gene_dict[condition].items():
                if gene_info.get("name") == "YBX1" or gene_id == "YBX1":
                    logging.info(f"Found YBX1 in condition {condition}")
                    logging.info(f"Total transcripts before filtering: {len(gene_info.get('transcripts', {}))}")
                    for transcript_id, transcript_info in gene_info.get('transcripts', {}).items():
                        value = transcript_info.get('value', 0)
                        logging.info(f"  Transcript {transcript_id}: TPM = {value:.2f}")
        
        # Check if reference and target conditions are defined
        has_specific_conditions = (hasattr(self, 'ref_conditions') and hasattr(self, 'target_conditions') and 
                                 self.ref_conditions and self.target_conditions)
        
        if has_specific_conditions:
            logging.debug(f"Filtering transcript usage plot to include only ref conditions: {self.ref_conditions} and target conditions: {self.target_conditions}")
            # Define all allowed conditions
            allowed_conditions = set(self.ref_conditions + self.target_conditions)
        
        for gene_name_or_id in self.gene_names:  # gene_names list contains gene names (symbols)
            gene_data_per_condition = {}  # Store gene data per condition
            found_gene_any_condition = False  # Flag if gene found in any condition
            
            # Only process allowed conditions if specific conditions are defined
            for condition, genes in self.updated_gene_dict.items():
                # Skip conditions that aren't in ref or target if we have those defined
                if has_specific_conditions and condition not in allowed_conditions:
                    #logging.debug(f"Skipping condition {condition} for gene {gene_name_or_id} (not in allowed conditions)")
                    continue
                    
                condition_gene_data = None
                for gene_id, gene_info in genes.items():
                    if "name" in gene_info and gene_info["name"] == gene_name_or_id.upper():  # Compare gene names (case-insensitive)
                        condition_gene_data = gene_info["transcripts"]  # Only need transcripts for usage plot
                        found_gene_any_condition = True
                        #logging.debug(f"Found gene {gene_name_or_id} in condition {condition}")
                        break  # Found gene in this condition, break inner loop
                if condition_gene_data:
                    gene_data_per_condition[condition] = condition_gene_data  # Store transcripts for this condition
            
            if not found_gene_any_condition:
                logging.warning(f"Gene {gene_name_or_id} not found in any allowed condition.")
                continue  # Skip to the next gene if not found
            
            # NEW APPROACH: If TPM threshold is set, collect max TPM for each transcript across all conditions
            if self.min_tpm_threshold is not None:
                # Create list of all allowed conditions
                allowed_conditions = set(self.ref_conditions + self.target_conditions) if has_specific_conditions else set(gene_data_per_condition.keys())
                
                # First, collect the max TPM for each transcript ONLY in allowed conditions
                transcript_max_tpm = {}
                
                for condition, transcripts in gene_data_per_condition.items():
                    # Skip conditions that aren't in ref or target if we have those defined
                    if has_specific_conditions and condition not in allowed_conditions:
                        continue
                        
                    for transcript_id, transcript_info in transcripts.items():
                        value = float(transcript_info.get("value", 0))
                        if transcript_id not in transcript_max_tpm:
                            transcript_max_tpm[transcript_id] = value
                        else:
                            transcript_max_tpm[transcript_id] = max(transcript_max_tpm[transcript_id], value)
                
                # Filter out transcripts that don\'t meet threshold in ANY allowed condition
                valid_transcripts = {t_id: t_val for t_id, t_val in transcript_max_tpm.items() 
                                    if t_val >= self.min_tpm_threshold}
                
                # Now keep ALL instances of valid transcripts in allowed conditions, even if below threshold
                filtered_gene_data_per_condition = {}
                for condition, transcripts in gene_data_per_condition.items():
                    # Skip conditions that aren\'t in ref or target if we have those defined
                    if has_specific_conditions and condition not in allowed_conditions:
                        continue
                        
                    filtered_transcripts = {}
                    for transcript_id, transcript_info in transcripts.items():
                        # Include this transcript if it\'s in our valid list
                        if transcript_id in valid_transcripts:
                            filtered_transcripts[transcript_id] = transcript_info
                    
                    # Only include conditions with transcripts
                    if filtered_transcripts:
                        filtered_gene_data_per_condition[condition] = filtered_transcripts
                
                # Replace original data with filtered data
                gene_data_per_condition = filtered_gene_data_per_condition
                
                # Log filtering results
                total_unique_transcripts = len(transcript_max_tpm)
                kept_transcripts = len(valid_transcripts)
                logging.debug(f"TPM filtering for {gene_name_or_id}: {kept_transcripts} of {total_unique_transcripts} unique transcripts have TPM >= {self.min_tpm_threshold} in at least one ref/target condition")
            
            if not gene_data_per_condition:
                logging.warning(f"No data available for gene {gene_name_or_id} after filtering. Skipping plot.")
                continue
            
            conditions = list(gene_data_per_condition.keys())
            n_bars = len(conditions)
            
            if n_bars == 0:
                logging.warning(f"No conditions found for gene {gene_name_or_id} after filtering.")
                continue
            
            
            fig, ax = plt.subplots(figsize=(12, 8))
            index = np.arange(n_bars)
            bar_width = 0.35
            opacity = 0.8
            max_transcripts = max(len(gene_data_per_condition[condition]) for condition in conditions)
            colors = plt.cm.plasma(np.linspace(0, 1, num=max_transcripts))

            bottom_val = np.zeros(n_bars)
            for i, condition in enumerate(conditions):
                transcripts = gene_data_per_condition[condition]
                if not transcripts:  # Skip if no transcript data for this condition
                    continue
                    
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

            ax.set_xlabel("Condition")
            ax.set_ylabel("Transcript Usage (TPM)")
            gene_display_name = list(gene_data_per_condition.values())[0].get("name", gene_name_or_id)  # Fallback to ID if no name
            
            # Update title to include TPM threshold if applied
            if self.min_tpm_threshold is not None:
                ax.set_title(f"Transcript Usage for {gene_display_name} by Condition (TPM >= {self.min_tpm_threshold} in any ref/target condition)")
            else:
                ax.set_title(f"Transcript Usage for {gene_display_name} by Condition")
                
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
                f"{gene_name_or_id}_transcript_usage_by_sample_type.pdf",  # Changed from .png to .pdf
            )
            plt.savefig(plot_path)
            plt.close(fig)

    def make_pie_charts(self):
        """
        Create pie charts for transcript alignment classifications and read assignment consistency.
        Handles both combined and separate sample data structures.
        """
        # Skip if reads_and_class is not provided
        if not self.reads_and_class:
            logging.warning("No reads_and_class data provided. Skipping pie chart creation.")
            return

        titles = ["Transcript Alignment Classifications", "Read Assignment Consistency"]
        
        # Check if reference and target conditions are defined
        has_specific_conditions = hasattr(self, 'ref_conditions') and hasattr(self, 'target_conditions')
        
        for title, data in zip(titles, self.reads_and_class):
            if isinstance(data, dict):
                if any(isinstance(v, dict) for v in data.values()):
                    # Separate sample data case (e.g. 'Mutants' and 'WildType')
                    for sample_name, sample_data in data.items():
                        # Skip conditions that aren't in ref or target if we have those defined
                        if has_specific_conditions and sample_name not in self.ref_conditions and sample_name not in self.target_conditions:
                            logging.debug(f"Skipping pie chart for {sample_name} (not in ref/target conditions)")
                            continue
                        self._create_pie_chart(f"{title} - {sample_name}", sample_data)
                else:
                    # Combined data case - always create this as it's an overall summary
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
            self.read_assignments_dir, f"{file_title}_pie_chart.pdf"  # Changed from .png to .pdf
        )
        plt.savefig(plot_path, bbox_inches="tight", dpi=300)
        plt.close()

    def plot_novel_transcript_contribution(self):
        """
        Creates a plot showing the percentage of expression from novel transcripts.
        - Y-axis: Percentage of expression from novel transcripts (combined across conditions)
        - X-axis: Expression log2 fold change between conditions
        - Point size: Overall expression level
        - Color: Red (target) to Blue (reference) indicating which condition contributes more to novel transcript expression
        """
        logging.info("Creating novel transcript contribution plot")
        
        # Skip if we don't have reference vs target conditions defined
        if not hasattr(self, 'ref_conditions') or not hasattr(self, 'target_conditions'):
            logging.warning("Cannot create novel transcript plot: missing reference or target conditions")
            return
        
        # Get actual condition labels
        ref_label = "+".join(self.ref_conditions)
        target_label = "+".join(self.target_conditions)
        
        # Set TPM threshold for transcript inclusion
        min_tpm_threshold = 10
        
        # Track all unique genes across all conditions
        all_genes = {}  # Dictionary to track gene_id -> gene_info mapping across conditions
        
        # First, collect all genes from all conditions
        for condition, genes in self.updated_gene_dict.items():
            for gene_id, gene_info in genes.items():
                gene_name = gene_info.get('name', gene_id)
                if gene_id not in all_genes:
                    all_genes[gene_id] = {'name': gene_name, 'conditions': {}}
                
                # Store condition-specific data
                all_genes[gene_id]['conditions'][condition] = gene_info
        
        logging.info(f"Total unique genes found across all conditions: {len(all_genes)}")
        
        # First, let's investigate the discrepancy between our filtering and the GTF filtering
        logging.info("Investigating transcript filtering discrepancy")
        
        # Track unique transcript IDs that pass our threshold (to match GTF filtering method)
        unique_transcripts_above_threshold = set()
        transcript_values = {}  # To store max values for debugging
        
        # First pass: identify all unique transcripts with TPM >= threshold
        for condition, genes in self.updated_gene_dict.items():
            # Check if this condition is in our ref or target groups
            is_relevant_condition = condition in self.ref_conditions or condition in self.target_conditions
            if not is_relevant_condition:
                continue
                
            for gene_id, gene_info in genes.items():
                transcripts = gene_info.get('transcripts', {})
                for transcript_id, transcript_info in transcripts.items():
                    value = float(transcript_info.get("value", 0))
                    
                    # Track max value for this transcript across all conditions
                    if transcript_id not in transcript_values:
                        transcript_values[transcript_id] = value
                    else:
                        transcript_values[transcript_id] = max(transcript_values[transcript_id], value)
                    
                    # Check if transcript meets threshold
                    if value >= min_tpm_threshold:
                        unique_transcripts_above_threshold.add(transcript_id)
        
        logging.info(f"Filtering comparison: Found {len(unique_transcripts_above_threshold)} unique transcripts with TPM >= {min_tpm_threshold}")
        logging.info(f"Filtering comparison: This compares to 8,230 transcripts reported by GTF filter")
        
        # Analyze distribution of TPM values to understand filtering
        tpm_value_counts = {
            "0-1": 0,
            "1-5": 0,
            "5-10": 0,
            "10-20": 0,
            "20-50": 0,
            "50-100": 0,
            "100+": 0
        }
        
        for transcript_id, max_value in transcript_values.items():
            if max_value < 1:
                tpm_value_counts["0-1"] += 1
            elif max_value < 5:
                tpm_value_counts["1-5"] += 1
            elif max_value < 10:
                tpm_value_counts["5-10"] += 1
            elif max_value < 20:
                tpm_value_counts["10-20"] += 1
            elif max_value < 50:
                tpm_value_counts["20-50"] += 1
            elif max_value < 100:
                tpm_value_counts["50-100"] += 1
            else:
                tpm_value_counts["100+"] += 1
        
        logging.info(f"TPM value distribution across transcripts: {tpm_value_counts}")
        
        # Check if there are any TTN transcripts in the unique set
        ttn_transcripts = [t for t in unique_transcripts_above_threshold if "TTN" in t.upper()]
        if ttn_transcripts:
            logging.info(f"Found {len(ttn_transcripts)} TTN transcripts in high TPM set: {ttn_transcripts}")
        
        # Prepare data storage for the main plot
        plot_data = []  # Re-initialize plot_data here
        
        # Track transcripts that pass TPM threshold
        total_transcripts = 0
        transcripts_above_threshold = 0
        total_genes = 0
        genes_with_high_expr_transcripts = 0
        
        # Process each gene from all_genes
        for gene_id, gene_data in all_genes.items():
            total_genes += 1
            gene_name = gene_data['name']
            conditions_data = gene_data['conditions']
            
            # Calculate expression for each condition group
            ref_total_exp = {cond: 0 for cond in self.ref_conditions}
            target_total_exp = {cond: 0 for cond in self.target_conditions}
            ref_novel_exp = {cond: 0 for cond in self.ref_conditions}
            target_novel_exp = {cond: 0 for cond in self.target_conditions}
            
            # Track if this gene has any high-expression transcripts
            gene_has_high_expr_transcript = False
            
            # Process each condition
            for condition, gene_info in conditions_data.items():
                transcripts = gene_info.get('transcripts', {})
                if not transcripts:
                    continue
                    
                # Check if this condition is in our condition groups
                is_ref = condition in self.ref_conditions
                is_target = condition in self.target_conditions
                
                if not (is_ref or is_target):
                    continue  # Skip conditions that aren't in ref or target groups
                    
                for transcript_id, transcript_info in transcripts.items():
                    total_transcripts += 1
                    
                    # Improved novel transcript identification - transcript is novel if not from Ensembl
                    transcript_is_reference = transcript_id.startswith("ENST")
                    is_novel = not transcript_is_reference
                    
                    value = float(transcript_info.get("value", 0))
                    
                    # Filter by TPM threshold - only count transcripts with TPM >= threshold
                    # We only check TPM threshold in ref and target conditions (not other conditions)
                    if value >= min_tpm_threshold:
                        transcripts_above_threshold += 1
                        gene_has_high_expr_transcript = True
                        
                        if is_ref:
                            ref_total_exp[condition] += value
                            if is_novel:
                                ref_novel_exp[condition] += value
                        
                        if is_target:
                            target_total_exp[condition] += value
                            if is_novel:
                                target_novel_exp[condition] += value
                    elif gene_name == "TTN" and condition in self.ref_conditions:
                        pass # Add pass to avoid empty block
            
            # Count genes with high-expression transcripts
            if gene_has_high_expr_transcript:
                genes_with_high_expr_transcripts += 1
                
                # Calculate average expression for each condition group
                ref_novel_pct = 0
                target_novel_pct = 0
                ref_expr_total = 0
                target_expr_total = 0
                ref_novel_expr_total = 0
                target_novel_expr_total = 0
                
                # Sum up expression values across conditions
                for cond in self.ref_conditions:
                    ref_expr_total += ref_total_exp.get(cond, 0)
                    ref_novel_expr_total += ref_novel_exp.get(cond, 0)
                    
                    # Also calculate percentages per condition for color coding
                    if ref_total_exp.get(cond, 0) > 0:
                        ref_novel_pct += (ref_novel_exp.get(cond, 0) / ref_total_exp[cond]) * 100
                
                for cond in self.target_conditions:
                    target_expr_total += target_total_exp.get(cond, 0)
                    target_novel_expr_total += target_novel_exp.get(cond, 0)
                    
                    # Also calculate percentages per condition for color coding
                    if target_total_exp.get(cond, 0) > 0:
                        target_novel_pct += (target_novel_exp.get(cond, 0) / target_total_exp[cond]) * 100
                
                # Average the condition-specific percentages (for color coding only)
                ref_novel_pct /= len([c for c in self.ref_conditions if c in ref_total_exp and ref_total_exp[c] > 0]) or 1
                target_novel_pct /= len([c for c in self.target_conditions if c in target_total_exp and target_total_exp[c] > 0]) or 1
                
                # Calculate overall novel percentage (for y-axis)
                combined_expr_total = ref_expr_total + target_expr_total
                combined_novel_expr_total = ref_novel_expr_total + target_novel_expr_total
                
                # Calculate log2 fold change using the total expression values
                if ref_expr_total > 0 and target_expr_total > 0:
                    log2fc = np.log2(target_expr_total / ref_expr_total)
                    
                    # Calculate novel transcript contribution difference (for color)
                    novel_pct_diff = target_novel_pct - ref_novel_pct
                    
                    # Calculate overall novel percentage (for y-axis)
                    if combined_expr_total > 0:
                        overall_novel_pct = (combined_novel_expr_total / combined_expr_total) * 100
                        
                        # Add data point
                        plot_data.append({
                            'gene_id': gene_id,
                            'gene_name': gene_name,
                            'ref_novel_pct': ref_novel_pct,
                            'target_novel_pct': target_novel_pct,
                            'novel_pct_diff': novel_pct_diff,
                            'overall_novel_pct': overall_novel_pct,
                            'log2fc': log2fc,
                            'total_expr': combined_expr_total
                        })
        
        # Report filtering results
        logging.info(f"TPM filtering: {transcripts_above_threshold} of {total_transcripts} transcripts have TPM >= {min_tpm_threshold} in ref or target conditions")
        logging.info(f"TPM filtering: {genes_with_high_expr_transcripts} of {total_genes} genes have at least one transcript with TPM >= {min_tpm_threshold} in ref or target conditions")
            
        # Create dataframe
        df = pd.DataFrame(plot_data)
        # Get the parent directory of gene_visualizations_dir
        parent_dir = os.path.dirname(self.gene_visualizations_dir)
        
        # Save the CSV to parent directory instead of gene_visualizations_dir
        df.to_csv(os.path.join(parent_dir, "novel_transcript_expression_data.csv"), index=False)
        
        # Log the number of genes used in the plot
        logging.info(f"Number of genes used in novel transcript plot after transcript-level TPM filtering: {len(df)}")
        
        if df.empty:
            logging.warning("No data available for novel transcript plot after transcript-level TPM filtering")
            return
            
        # Create the plot with more space on right for legend
        plt.figure(figsize=(16, 10))  # Increased width from 14 to 16
        
        # Define red-blue colormap
        norm = Normalize(vmin=-50, vmax=50)  # Normalize based on difference range
        cmap = cm.get_cmap('coolwarm')  # Red-Blue colormap
        
        # More dramatic scaling for point sizes
        min_size = 30
        max_size = 800  # Much larger maximum size
        
        # Use np.power for more dramatic scaling differences
        expression_values = df['total_expr'].values
        max_expr = expression_values.max()
        min_expr = expression_values.min()
        
        # Log the actual min and max expression values for reference
        logging.debug(f"Expression range in data: min={min_expr}, max={max_expr}")
        
        # Define the scaling function that will be used for both data points and legend
        def scale_point_size(expr_value, min_expr, max_expr, min_size, max_size, power=0.3):
            """Scale expression values to point sizes using the same formula for data and legend"""
            # Normalize the expression value to [0,1] range
            if max_expr == min_expr:  # Avoid division by zero
                normalized = 0
            else:
                normalized = (expr_value - min_expr) / (max_expr - min_expr)
            # Apply power scaling and convert to point size
            return min_size + (max_size - min_size) * (normalized ** power)
        
        # Apply scaling to actual data points
        scaled_sizes = [scale_point_size(val, min_expr, max_expr, min_size, max_size) for val in expression_values]
        
        # Plot points with scaled sizes
        sc = plt.scatter(df['log2fc'], df['overall_novel_pct'],
                        s=scaled_sizes, 
                        c=df['novel_pct_diff'], 
                        cmap=cmap, 
                        norm=norm,
                        alpha=0.8,
                        edgecolors='black')
        
        # Add color legend on the right
        cbar = plt.colorbar(sc, orientation='vertical', pad=0.02)
        cbar.set_label('Novel transcript usage difference (%)', size=12)
        cbar.ax.tick_params(labelsize=10)
        
        # Use red and blue blocks to explain the colormap
        plt.figtext(0.92, 0.72, f'Blue = higher (%) in {self.ref_conditions}', fontsize=12, ha='center')
        plt.figtext(0.92, 0.75, f'Red = higher (%) in {self.target_conditions}', fontsize=12, ha='center')
        
        # Add size legend directly to the plot
        # Create legend elements for different sizes with new values: 50, 500, 5000
        size_legend_values = [50, 500, 5000]
        size_legend_elements = []
        
        # Calculate sizes for legend using EXACTLY the same scaling function as for the data points
        for val in size_legend_values:
            # Use the same scaling function defined above
            # If the value is outside the actual data range, clamp it to the range
            clamped_val = min(max(val, min_expr), max_expr)
            size = scale_point_size(clamped_val, min_expr, max_expr, min_size, max_size)
            
            # Log the actual size being used for the legend point
            logging.debug(f"Legend point {val} TPM scaled to size {size}")
            
            # Convert area to diameter for Line2D (sqrt of area * 2)
            marker_diameter = 2 * np.sqrt(size / np.pi)
            
            size_legend_elements.append(
                plt.Line2D([0], [0], marker='o', color='w', 
                          markerfacecolor='gray', markersize=marker_diameter,
                          label=f'{val:.0f} TPM')
            )
        
        # Position legend
        plt.legend(handles=size_legend_elements, 
                  title="Expression Level",
                  loc='center left', 
                  bbox_to_anchor=(1.15, 0.5),
                  frameon=False,
                  title_fontsize=12,
                  fontsize=12)
        
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        
        # Add labels and title with actual condition names
        plt.xlabel('Log2 Fold Change', fontsize=12)
        plt.ylabel('Total expression from novel transcripts (%)', fontsize=12)
        plt.title('Novel Transcript Usage vs Expression Change between High Risk and Low Risk Phenotypes', fontsize=18)
        
        plt.grid(True, alpha=0.3)
        
        # Use tighter layout settings
        plt.tight_layout()
        
        # Save figure to parent directory instead of gene_visualizations_dir
        output_path = os.path.join(parent_dir, "novel_transcript_expression_plot.pdf")
        plt.savefig(output_path, dpi=300, bbox_inches='tight', pad_inches=0.5)
        plt.close()
        
        logging.debug(f"Novel transcript expression plot saved to {output_path}")


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
            label=f"Down-regulated in ({target_label})",
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
            self.output_path / f"volcano_plot_{feature_type}.pdf"  # Changed from .png to .pdf
        )
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
        plot_path = self.output_path / f"ma_plot_{feature_type}.pdf"  # Changed from .png to .pdf
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

    
    def plot_pca(
        self,
        pca_df: pd.DataFrame,
        title: str,
        output_prefix: str,
        explained_variance: np.ndarray,
        loadings: np.ndarray,
        feature_names: List[str]
    ) -> Path:
        """Plot PCA scatter plot, scree plot, and loadings."""
        plt.figure(figsize=(8, 6))

        # Extract variance info from title for axis labels only
        pc1_var = title.split("PC1 (")[1].split("%)")[0]
        pc2_var = title.split("PC2 (")[1].split("%)")[0]

        # Get clean title without PCs and variance
        base_title = title.split(' Level PCA: ')[0]
        comparison = title.split(': ')[1].split('PC1')[0].strip()
        clean_title = f"{base_title} Level PCA: {comparison}"

        # Update group labels in the DataFrame
        condition_mapping = {'Target': title.split(": ")[1].split(" vs ")[0],
                            'Reference': title.split(" vs ")[1].split("PC1")[0].strip()}
        pca_df['group'] = pca_df['group'].map(condition_mapping)

        # Create scatter plot
        sns.scatterplot(x='PC1', y='PC2', hue='group', data=pca_df, s=100)
        plt.xlabel(f'PC1 ({pc1_var}%)')
        plt.ylabel(f'PC2 ({pc2_var}%)')
        plt.title(clean_title)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.tight_layout()

        scatter_plot_path = self.output_path / f"{output_prefix}.pdf" # Changed from .png to .pdf
        plt.savefig(scatter_plot_path)
        plt.close()

        # --- Scree Plot ---
        self._plot_scree(explained_variance, output_prefix)

        # --- Loadings ---
        self._output_loadings(loadings, feature_names, output_prefix)

        return scatter_plot_path # Return path to scatter plot

    def _plot_scree(self, explained_variance: np.ndarray, output_prefix: str) -> Path:
        """Plot scree plot of explained variance."""
        plt.figure(figsize=(8, 6))
        num_components = len(explained_variance)
        component_numbers = range(1, num_components + 1)

        plt.bar(component_numbers, explained_variance * 100)
        plt.xlabel('Principal Component')
        plt.ylabel('Percentage of Explained Variance')
        plt.title('Scree Plot')
        plt.xticks(component_numbers) # Ensure all component numbers are labeled
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.tight_layout()

        scree_plot_path = self.output_path / f"scree_{output_prefix}.pdf" # Changed from .png to .pdf
        plt.savefig(scree_plot_path)
        plt.close()
        return scree_plot_path

    def _output_loadings(self, loadings: np.ndarray, feature_names: List[str], output_prefix: str, top_n: int = 10) -> Path:
        """Output top N loadings for PC1 and PC2."""
        # Generate column names dynamically based on the number of components
        num_components = loadings.shape[0] # Get the number of components from loadings shape
        pc_columns = [f'PC{i+1}' for i in range(num_components)]

        loadings_df = pd.DataFrame(loadings.T, index=feature_names, columns=pc_columns) # Use dynamic column names

        output_path = self.output_path / f"loadings_{output_prefix}.txt"
        with open(output_path, 'w') as f:
            f.write("PCA Loadings (Top {} Features for PC1 and PC2):\n\n".format(top_n))
            for pc_name in ['PC1', 'PC2']:
                f.write(f"\n--- {pc_name} ---\n")
                # Sort by absolute value of loading
                top_loadings = loadings_df.sort_values(by=pc_name, key=lambda x: x.abs(), ascending=False).head(top_n)
                for gene, loading in top_loadings[pc_name].items(): # Iterate over series items
                    f.write(f"{gene}:\t{loading:.4f}\n") # Tab-separated for readability
        return output_path

