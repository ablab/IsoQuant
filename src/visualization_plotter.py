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
        self.display_threshold = filter_transcripts
        
        # Explicitly set reference and target conditions
        self.ref_conditions = ref_conditions if ref_conditions else []
        self.target_conditions = target_conditions if target_conditions else []
        
        # Log conditions for debugging
        if self.ref_conditions or self.target_conditions:
            expected_conditions = set(self.ref_conditions + self.target_conditions)
            actual_conditions = set(self.updated_gene_dict.keys())
            if expected_conditions != actual_conditions:
                logging.warning(f"Mismatch between provided conditions and keys in updated_gene_dict. "
                                f"Expected: {sorted(list(expected_conditions))}, Found: {sorted(list(actual_conditions))}")
            else:
                 logging.info(f"Plotting with ref conditions: {self.ref_conditions} and target conditions: {self.target_conditions}")
        else:
            logging.warning("No ref_conditions or target_conditions set, plots will include all conditions found in updated_gene_dict")
        
        # Log the threshold value if provided (for context)
        if self.display_threshold is not None:
            logging.info(f"Transcript data assumes upstream filtering with TPM >= {self.display_threshold}")
        
        # Ensure output directories exist
        if self.gene_visualizations_dir:
            os.makedirs(self.gene_visualizations_dir, exist_ok=True)
        if self.read_assignments_dir: # Check if read_assignments_dir is not None
            os.makedirs(self.read_assignments_dir, exist_ok=True)

    def plot_transcript_map(self):
        """Plot transcript structure using pre-filtered gene data."""
        if not self.gene_visualizations_dir:
            logging.warning("No gene_visualizations_dir provided. Skipping transcript map plotting.")
            return

    

        for gene_name_or_id in self.gene_names:  # gene_names list contains gene names (symbols)
            gene_data = None  # Initialize gene_data to None
            
            # Find the gene in the pre-filtered dictionary.
            # We only need one instance of the gene structure, as it should be consistent.
            # Iterate through conditions until the gene is found.
            for condition, genes in self.updated_gene_dict.items():
                for gene_id, gene_info in genes.items():
                    # Compare gene names (case-insensitive matching)
                    if "name" in gene_info and gene_info["name"].upper() == gene_name_or_id.upper():
                        gene_data = gene_info
                        # No need to log which condition it came from, as it's pre-filtered.
                        break # Found gene info 
                if gene_data:
                    break  # Found gene, stop searching conditions

            if not gene_data:
                logging.warning(f"Gene '{gene_name_or_id}' not found in the provided gene dictionary.")
                continue  # Skip to the next gene if not found

            # Get chromosome info and calculate buffer
            chromosome = gene_data.get("chromosome", "Unknown")
            start = gene_data.get("start", 0)
            end = gene_data.get("end", 0)
            
            # Find the actual min/max coordinates of all exons
            min_exon_start = min(exon["start"] for transcript in gene_data["transcripts"].values() 
                               for exon in transcript["exons"])
            max_exon_end = max(exon["end"] for transcript in gene_data["transcripts"].values() 
                             for exon in transcript["exons"])
            
            # Calculate buffer (10% of total width)
            width = max(end, max_exon_end) - min(start, min_exon_start)
            buffer = width * 0.10  # Increased from 5% to 10%
            plot_start = min(start, min_exon_start) - buffer
            plot_end = max(end, max_exon_end) + buffer

            # REMOVED FILTERING LOGIC - Directly use transcripts from gene_data
            filtered_transcripts = gene_data["transcripts"]
            
            # Skip plotting if no transcripts are present (this might happen if upstream filtering removed all)
            if not filtered_transcripts:
                logging.warning(f"No transcripts found for gene {gene_name_or_id} in the input data. Skipping plot.")
                continue
                
            # Calculate plot height based on number of filtered transcripts
            num_transcripts = len(filtered_transcripts)
            plot_height = max(10, num_transcripts * 0.6) # Increased base height and multiplier
            # Use INFO level for starting plot creation, DEBUG for saving it.
            logging.debug(f"Creating transcript map for gene '{gene_name_or_id}' with {num_transcripts} transcripts.")

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

                # Sort exons based on strand direction
                exons = sorted(transcript_info["exons"], 
                             key=lambda x: x["start"] if gene_data["strand"] == "+" else -x["start"])

                # Exon blocks with color based on reference status
                for exon_idx, exon in enumerate(exons, 1):
                    exon_length = exon["end"] - exon["start"]
                    if self.ref_only: # Check ref_only flag
                        exon_color = "skyblue" # If ref_only, always treat as reference
                        exon_alpha = 1.0
                    else:
                        is_reference_exon = exon["exon_id"].startswith("E") # Original logic
                        exon_color = "skyblue" if is_reference_exon else "red"
                        exon_alpha = 1.0 if is_reference_exon else 0.6
                    
                    # Add exon rectangle
                    rect = plt.Rectangle(
                        (exon["start"], i - 0.4),
                        exon_length,
                        0.8,
                        color=exon_color,
                        alpha=exon_alpha
                    )
                    ax.add_patch(rect)
                    
                # Store y-axis label information
                y_ticks.append(i)
                # Get transcript name with fallback options
                transcript_name = (transcript_info.get("name") or 
                                  transcript_info.get("transcript_id") or 
                                  transcript_id)

                y_labels.append(f"{transcript_name}")

            # Set up the plot formatting with just chromosome
            gene_display_name = gene_data.get("name", gene_name_or_id) # Fallback to ID if no name
            
            # Update title to include TPM threshold if applied
            if self.display_threshold is not None:
                title = f"Transcript Structure - {gene_display_name} (Chromosome {chromosome}) (Input filtered at TPM >= {self.display_threshold})"
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

            plt.tight_layout(rect=[0.05, 0, 0.9, 1]) # Give more space on left (0.05) and right (1-0.9=0.1)
            plot_path = os.path.join(
                self.gene_visualizations_dir, f"{gene_name_or_id}_splicing.pdf" # Changed from .png to .pdf
            )
            plt.savefig(plot_path, bbox_inches='tight', dpi=300)
            plt.close(fig)


   
    def plot_transcript_usage(self):
        """Visualize transcript usage for each gene across conditions from pre-filtered data."""
        if not self.gene_visualizations_dir:
            logging.warning("No gene_visualizations_dir provided. Skipping transcript usage plotting.")
            return

        # The input updated_gene_dict is assumed to be pre-filtered.
        
        for gene_name_or_id in self.gene_names:  # gene_names list contains gene names (symbols)
            gene_data_per_condition = {}  # Store gene transcript data per condition
            found_gene_any_condition = False  # Flag if gene found in any condition

            # Iterate directly through the pre-filtered dictionary
            for condition, genes in self.updated_gene_dict.items():
                condition_gene_data = None
                for gene_id, gene_info in genes.items():
                    # Compare gene names (case-insensitive matching)
                    if "name" in gene_info and gene_info["name"].upper() == gene_name_or_id.upper():
                        condition_gene_data = gene_info.get("transcripts", {}) # Get transcripts, default to empty dict
                        found_gene_any_condition = True
                        #logging.debug(f"Found gene {gene_name_or_id} data for condition {condition}")
                        break  # Found gene in this condition, break inner loop
                if condition_gene_data is not None: # Store even if empty, to represent the condition
                    gene_data_per_condition[condition] = condition_gene_data

            if not found_gene_any_condition:
                logging.warning(f"Gene {gene_name_or_id} not found in any condition within the pre-filtered updated_gene_dict.")
                continue  # Skip to the next gene if not found

            if not gene_data_per_condition:
                logging.warning(f"No transcript data available for gene {gene_name_or_id} across conditions. Skipping plot.")
                continue

            # --- Reorder conditions: Reference first, then Target ---
            all_conditions = list(gene_data_per_condition.keys())
            ref_conditions_present = sorted([c for c in all_conditions if c in self.ref_conditions])
            target_conditions_present = sorted([c for c in all_conditions if c in self.target_conditions])
            # Include any other conditions found in the data but not specified as ref/target (shouldn't happen with pre-filtering, but safe)
            other_conditions_present = sorted([c for c in all_conditions if c not in self.ref_conditions and c not in self.target_conditions])
            
            conditions = ref_conditions_present + target_conditions_present + other_conditions_present
            # --- End Reordering ---

            n_bars = len(conditions)
            
            if n_bars == 0:
                logging.warning(f"No conditions to plot for gene {gene_name_or_id}.") 
                continue
            
            
            fig, ax = plt.subplots(figsize=(12, 8))
            index = np.arange(n_bars)
            bar_width = 0.35
            opacity = 0.8
            
            # Determine unique transcripts across all plotted conditions for consistent coloring
            all_transcript_ids = set()
            for condition in conditions:
                all_transcript_ids.update(gene_data_per_condition[condition].keys())
            unique_transcripts = sorted(list(all_transcript_ids))
            transcript_to_color_idx = {tid: idx for idx, tid in enumerate(unique_transcripts)}
            colors = plt.cm.plasma(np.linspace(0, 1, num=len(unique_transcripts)))

            bottom_val = np.zeros(n_bars)
            plotted_labels = set() # To avoid duplicate legend entries
            
            for i, condition in enumerate(conditions):
                transcripts = gene_data_per_condition[condition]
                if not transcripts:  # Skip if no transcript data for this condition
                    continue
                    
                # Sort transcripts for consistent stacking order (optional but good practice)
                sorted_transcript_items = sorted(transcripts.items(), key=lambda item: item[0])
                
                for transcript_id, transcript_info in sorted_transcript_items:
                    color_idx = transcript_to_color_idx.get(transcript_id, 0) # Fallback index 0
                    color = colors[color_idx % len(colors)]
                    value = transcript_info["value"]
                    # Get transcript name with fallback options
                    transcript_name = (transcript_info.get("name") or
                                     transcript_info.get("transcript_id") or
                                     transcript_id)
                    
                    label = transcript_name if transcript_name not in plotted_labels else ""
                    if label:
                        plotted_labels.add(label)
                        
                    ax.bar(
                        i,
                        float(value),
                        bar_width,
                        bottom=bottom_val[i],
                        alpha=opacity,
                        color=color,
                        label=label,
                    )
                    bottom_val[i] += float(value)

            ax.set_xlabel("Condition")
            ax.set_ylabel("Transcript Usage (TPM)")
            # Find a representative gene name (assuming transcripts exist in at least one condition)
            first_condition_with_data = next((cond for cond in conditions if gene_data_per_condition[cond]), None)
            gene_display_name = gene_name_or_id # Default to ID
            if first_condition_with_data:
                # Attempt to get gene name from the first transcript entry in the first condition with data
                first_transcript_info = next(iter(gene_data_per_condition[first_condition_with_data].values()), None)
                if first_transcript_info:
                    # Assuming gene name might be stored within transcript info, or fallback
                    # This part might need adjustment based on your actual data structure
                    # If gene name isn't in transcript info, you might need to fetch it differently
                    pass # Placeholder - logic to get gene name needs review based on structure

            # Updated title - Include threshold if available
            if self.display_threshold is not None:
                 ax.set_title(f"Transcript Usage for {gene_display_name} by Condition (Input filtered at TPM >= {self.display_threshold})")
            else:
                ax.set_title(f"Transcript Usage for {gene_display_name} by Condition")
                
            ax.set_xticks(index)
            ax.set_xticklabels(conditions)
            
            # Update legend handling to use plotted_labels
            handles, labels = ax.get_legend_handles_labels()
            if handles: # Only show legend if there are items to show
                ax.legend(
                    handles,
                    labels,
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
        
        # Input data is assumed to be pre-filtered, so no need to check ref/target conditions here.
        # Plot for all sample groups found in the data.
        
        for title, data in zip(titles, self.reads_and_class):
            if isinstance(data, dict):
                # Check if the dictionary values are also dictionaries (indicating separate sample groups)
                if data and isinstance(next(iter(data.values()), None), dict):
                    # Separate sample data case (e.g. {'Mutants': {...}, 'WildType': {...}})
                    logging.debug(f"Creating separate pie charts for samples in '{title}'")
                    for sample_name, sample_data in data.items():
                        # No filtering needed here, plot for every sample found
                        self._create_pie_chart(f"{title} - {sample_name}", sample_data)
                elif data: # Check if data is not empty before proceeding
                    # Combined data case or single sample group provided directly
                    logging.debug(f"Creating combined pie chart for '{title}'")
                    self._create_pie_chart(title, data)
                else:
                     logging.warning(f"Empty data dictionary provided for pie chart '{title}'. Skipping.")
            else:
                logging.warning(f"Skipping unexpected data type for pie chart '{title}': {type(data)}")

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
        plt.savefig(plot_path, bbox_inches='tight', dpi=300)
        plt.close()

    def plot_read_length_effects(self, length_effects):
        """
        Plot how read length relates to (a) assignment uniqueness and (b) FSM/ISM/mono classification.
        Saves two bar charts into read_assignments_dir.
        """
        if not self.read_assignments_dir:
            logging.warning("No read_assignments_dir provided. Skipping length effects plotting.")
            return

        bins = length_effects['bins']
        totals = length_effects['totals']

        # Assignment uniqueness plot
        df_a_rows = []
        for b in bins:
            row = {'bin': b, **length_effects['by_bin_assignment'][b], 'TOTAL': totals[b]}
            df_a_rows.append(row)
        df_a = pd.DataFrame(df_a_rows)
        if df_a.empty:
            logging.warning("No data available for assignment uniqueness plot; skipping.")
            return
        df_a.set_index('bin', inplace=True)

        # Determine assignment categories dynamically and ensure columns exist
        assignment_keys = length_effects.get('assignment_keys', [])
        if not assignment_keys:
            assignment_keys = [c for c in df_a.columns if c != 'TOTAL']
        for key in assignment_keys:
            if key not in df_a.columns:
                df_a[key] = 0

        # Normalize to percentages per bin
        for col in assignment_keys:
            df_a[col] = np.where(df_a['TOTAL'] > 0, df_a[col] / df_a['TOTAL'] * 100.0, 0.0)

        # Preferred column order if present
        preferred_order = ['UNIQUE', 'AMBIGUOUS', 'OTHER', 'INCONSISTENT', 'UNASSIGNED']
        ordered_cols = [c for c in preferred_order if c in assignment_keys] + [c for c in assignment_keys if c not in preferred_order]
        if not ordered_cols:
            logging.warning("No assignment columns to plot after normalization; skipping.")
            return

        ax = df_a[ordered_cols].plot(kind='bar', stacked=True, figsize=(12,6), colormap='tab20')
        ax.set_ylabel('Percentage of reads')
        ax.set_title('Read assignment uniqueness by read length')
        ax.legend(title='Assignment')
        plt.tight_layout()
        out1 = os.path.join(self.read_assignments_dir, 'read_length_vs_assignment_uniqueness.pdf')
        plt.savefig(out1, bbox_inches='tight', dpi=300)
        plt.close()

    def plot_read_length_histogram(self, hist_data):
        """
        Plot a histogram of read lengths using precomputed bin edges/counts.
        """
        if not self.read_assignments_dir:
            logging.warning("No read_assignments_dir provided. Skipping length histogram plot.")
            return

        edges = hist_data.get('edges', [])
        counts = hist_data.get('counts', [])
        total = hist_data.get('total', 0)
        if not edges or not counts:
            logging.warning("Empty histogram data; skipping.")
            return

        # Build midpoints for bar plotting
        mids = [(edges[i] + edges[i+1]) / 2.0 for i in range(len(counts))]
        widths = [edges[i+1] - edges[i] for i in range(len(counts))]

        plt.figure(figsize=(12,6))
        plt.bar(mids, counts, width=widths, align='center', color='steelblue', edgecolor='black')
        plt.xlabel('Read length (bp)')
        plt.ylabel('Read count')
        plt.title(f'Read length histogram (total n={total:,})')
        plt.tight_layout()
        outp = os.path.join(self.read_assignments_dir, 'read_length_histogram.pdf')
        plt.savefig(outp, bbox_inches='tight', dpi=300)
        plt.close()

    def plot_read_length_vs_assignment(self, length_vs_assignment):
        """
        Plot read-length bins vs assignment_type and vs classification as stacked bar charts.
        Saves two PDFs into read_assignments_dir.
        """
        if not self.read_assignments_dir:
            logging.warning("read_assignments_dir not set; skipping length vs assignment plots")
            return

        import pandas as pd
        import matplotlib.pyplot as plt

        bins = length_vs_assignment.get('bins', [])
        a_counts = length_vs_assignment.get('assignment', {})
        c_counts = length_vs_assignment.get('classification', {})

        # Build DataFrames
        def to_df(counts_dict):
            rows = []
            for (b, key), val in counts_dict.items():
                rows.append({'bin': b, 'key': key, 'count': val})
            df = pd.DataFrame(rows)
            if df.empty:
                return df
            pivot = df.pivot_table(index='bin', columns='key', values='count', aggfunc='sum', fill_value=0)
            # Ensure bin order
            pivot = pivot.reindex(bins, axis=0).fillna(0)
            return pivot

        df_a = to_df(a_counts)
        df_c = to_df(c_counts)

        def plot_stacked(pivot_df, title, filename):
            if pivot_df.empty:
                logging.warning(f"No data for plot: {title}")
                return
            ax = pivot_df.plot(kind='bar', stacked=True, figsize=(12, 6))
            ax.set_xlabel('Read length bin')
            ax.set_ylabel('Read count')
            ax.set_title(title)
            plt.tight_layout()
            out = os.path.join(self.read_assignments_dir, filename)
            plt.savefig(out)
            plt.close()
            logging.info(f"Saved plot: {out}")

        plot_stacked(df_a, 'Read length vs assignment_type', 'length_vs_assignment_type.pdf')
        plot_stacked(df_c, 'Read length vs classification', 'length_vs_classification.pdf')

    def plot_novel_transcript_contribution(self):
        """
        Creates a plot showing the percentage of expression from novel transcripts.
        - Y-axis: Percentage of expression from novel transcripts (combined across conditions)
        - X-axis: Expression log2 fold change between conditions
        - Point size: Overall expression level
        - Color: Red (target) to Blue (reference) indicating which condition contributes more to novel transcript expression
        Assumes input updated_gene_dict is already filtered appropriately.
        """
        logging.info("Creating novel transcript contribution plot")
        
        # Skip if we don't have reference vs target conditions defined
        if not (hasattr(self, 'ref_conditions') and self.ref_conditions and 
                hasattr(self, 'target_conditions') and self.target_conditions):
            logging.warning("Cannot create novel transcript plot: missing reference or target conditions")
            return
        
        # Get actual condition labels
        ref_label = "+".join(self.ref_conditions)
        target_label = "+".join(self.target_conditions)
        
        # Track all unique genes across all conditions
        all_genes = {}  # Dictionary to track gene_id -> gene_info mapping across conditions
        
        # Collect all genes present in the (presumably pre-filtered) input dictionary
        for condition, genes in self.updated_gene_dict.items():
            for gene_id, gene_info in genes.items():
                gene_name = gene_info.get('name', gene_id)
                if gene_id not in all_genes:
                    all_genes[gene_id] = {'name': gene_name, 'conditions': {}}
                
                # Store condition-specific data only if it's a ref or target condition
                if condition in self.ref_conditions or condition in self.target_conditions:
                     all_genes[gene_id]['conditions'][condition] = gene_info
        
        logging.info(f"Total unique genes found across relevant conditions: {len(all_genes)}")
        
        # Prepare data storage for the main plot
        plot_data = []
        
        # Process each gene from all_genes
        for gene_id, gene_data in all_genes.items():
            gene_name = gene_data['name']
            conditions_data = gene_data['conditions'] # Contains only ref/target conditions now
            
            # Calculate expression for each condition group
            ref_total_exp = {cond: 0 for cond in self.ref_conditions}
            target_total_exp = {cond: 0 for cond in self.target_conditions}
            ref_novel_exp = {cond: 0 for cond in self.ref_conditions}
            target_novel_exp = {cond: 0 for cond in self.target_conditions}
            
            gene_has_any_transcript = False # Check if the gene has any transcripts in ref/target
            
            # Process each relevant condition for the gene
            for condition, gene_info in conditions_data.items():
                transcripts = gene_info.get('transcripts', {})
                if not transcripts:
                    continue
                 
                gene_has_any_transcript = True # Mark that this gene has data
                    
                # Check if this condition is in our condition groups (redundant check now, but safe)
                is_ref = condition in self.ref_conditions
                is_target = condition in self.target_conditions
                
                for transcript_id, transcript_info in transcripts.items():
                    # Improved novel transcript identification - transcript is novel if not from Ensembl
                    transcript_is_reference = transcript_id.startswith("ENST")
                    is_novel = not transcript_is_reference
                    
                    value = float(transcript_info.get("value", 0))
                    
                    # REMOVED Filtering by TPM threshold - Now process all transcripts present
                    if is_ref:
                        ref_total_exp[condition] += value
                        if is_novel:
                            ref_novel_exp[condition] += value
                    
                    if is_target:
                        target_total_exp[condition] += value
                        if is_novel:
                            target_novel_exp[condition] += value
            
            # Only proceed if the gene had transcripts in the relevant conditions
            if gene_has_any_transcript:
                # Calculate average expression for each condition group
                ref_novel_pct = 0
                target_novel_pct = 0
                ref_expr_total = 0
                target_expr_total = 0
                ref_novel_expr_total = 0
                target_novel_expr_total = 0
                
                # Sum up expression values across conditions
                num_ref_conditions_with_expr = 0
                for cond in self.ref_conditions:
                    cond_total_exp = ref_total_exp.get(cond, 0)
                    cond_novel_exp = ref_novel_exp.get(cond, 0)
                    ref_expr_total += cond_total_exp
                    ref_novel_expr_total += cond_novel_exp
                    if cond_total_exp > 0:
                        ref_novel_pct += (cond_novel_exp / cond_total_exp) * 100
                        num_ref_conditions_with_expr += 1
                
                num_target_conditions_with_expr = 0
                for cond in self.target_conditions:
                    cond_total_exp = target_total_exp.get(cond, 0)
                    cond_novel_exp = target_novel_exp.get(cond, 0)
                    target_expr_total += cond_total_exp
                    target_novel_expr_total += cond_novel_exp
                    if cond_total_exp > 0:
                        target_novel_pct += (cond_novel_exp / cond_total_exp) * 100
                        num_target_conditions_with_expr += 1
                
                # Average the condition-specific percentages (for color coding only)
                ref_novel_pct /= num_ref_conditions_with_expr or 1
                target_novel_pct /= num_target_conditions_with_expr or 1
                
                # Calculate overall novel percentage (for y-axis)
                combined_expr_total = ref_expr_total + target_expr_total
                combined_novel_expr_total = ref_novel_expr_total + target_novel_expr_total
                
                # Check for non-zero total expression before calculating percentages and fold change
                if combined_expr_total > 0: 
                    # Calculate log2 fold change using the total expression values
                    # Add pseudocount to avoid division by zero or log(0)
                    pseudocount = 1e-6 # Small value to add
                    log2fc = np.log2((target_expr_total + pseudocount) / (ref_expr_total + pseudocount))
                    
                    # Calculate novel transcript contribution difference (for color)
                    novel_pct_diff = target_novel_pct - ref_novel_pct
                    
                    # Calculate overall novel percentage (for y-axis)
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
        
        # Create dataframe
        df = pd.DataFrame(plot_data)
        
        if df.empty:
            logging.warning("No data available for novel transcript plot after processing.") # Adjusted warning
            return
            
        # Get the parent directory of gene_visualizations_dir
        parent_dir = os.path.dirname(self.gene_visualizations_dir)
        
        # Save the CSV to parent directory instead of gene_visualizations_dir
        csv_path = os.path.join(parent_dir, "novel_transcript_expression_data.csv")
        df.to_csv(csv_path, index=False)
        logging.info(f"Novel transcript expression data saved to {csv_path}")
        
        # Log the number of genes used in the plot
        logging.info(f"Number of genes included in novel transcript plot: {len(df)}")
            
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
        # Handle case where expression_values might be empty or all zero
        if len(expression_values) == 0 or expression_values.max() == expression_values.min():
            max_expr = 1
            min_expr = 0
            logging.warning("Cannot determine expression range for point scaling; using default [0, 1].")
        else:
            max_expr = expression_values.max()
            min_expr = expression_values.min()
        
        # Log the actual min and max expression values for reference
        logging.debug(f"Expression range in data: min={min_expr}, max={max_expr}")
        
        # Define the scaling function that will be used for both data points and legend
        def scale_point_size(expr_value, min_expr, max_expr, min_size, max_size, power=0.3):
            """Scale expression values to point sizes using the same formula for data and legend"""
            # Normalize the expression value to [0,1] range
            if max_expr == min_expr:  # Avoid division by zero or invalid range
                normalized = 0.5 # Default to middle size if range is zero
            else:
                # Clamp value to range before normalizing to handle potential outliers from pseudocounts
                clamped_value = np.clip(expr_value, min_expr, max_expr)
                normalized = (clamped_value - min_expr) / (max_expr - min_expr)
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
        plt.figtext(0.92, 0.72, f'Blue = higher (%) in {ref_label}', fontsize=12, ha='center')
        plt.figtext(0.92, 0.75, f'Red = higher (%) in {target_label}', fontsize=12, ha='center')
        
        # Add size legend directly to the plot
        # Create legend elements for different sizes with new values: 50, 500, 5000
        size_legend_values = [50, 500, 5000]
        size_legend_elements = []
        
        # Calculate sizes for legend using EXACTLY the same scaling function as for the data points
        for val in size_legend_values:
            # Use the same scaling function defined above
            size = scale_point_size(val, min_expr, max_expr, min_size, max_size)
            
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
        top_n: int = 60,  # Increased from 10 to 20
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
            self.output_path / f"volcano_plot_{feature_type}.pdf"
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

