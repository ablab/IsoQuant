import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

class PlotOutput:
    def __init__(self, updated_gene_dict, gene_names, output_directory, reads_and_class=None, filter_transcripts=None, conditions=False, use_counts=False):
        self.updated_gene_dict = updated_gene_dict
        self.gene_names = gene_names
        self.output_directory = output_directory
        self.reads_and_class = reads_and_class
        self.filter_transcripts = filter_transcripts
        self.conditions = conditions
        self.use_counts = use_counts

        # Create visualization subdirectory if it doesn't exist
        self.visualization_dir = os.path.join(self.output_directory, "visualization")
        os.makedirs(self.visualization_dir, exist_ok=True)

    def plot_transcript_map(self):
        # Get the first condition's gene dictionary
        first_condition = next(iter(self.updated_gene_dict))
        gene_dict = self.updated_gene_dict[first_condition]
        
        for gene_name in self.gene_names:
            if gene_name in gene_dict:
                gene_data = gene_dict[gene_name]
                num_transcripts = len(gene_data['transcripts'])
                plot_height = max(3, num_transcripts * 0.3)  # Adjust the height dynamically

                fig, ax = plt.subplots(figsize=(12, plot_height))  # Adjust height dynamically

                if self.filter_transcripts is not None:
                    ax.set_title(f"Transcripts of Gene: {gene_data['name']} on Chromosome {gene_data['chromosome']} with value over {self.filter_transcripts}")
                else:
                    ax.set_title(f"Transcripts of Gene: {gene_data['name']} on Chromosome {gene_data['chromosome']}")
                    
                ax.set_xlabel("Chromosomal position")
                ax.set_ylabel("Transcripts")
                ax.set_yticks(range(num_transcripts))
                ax.set_yticklabels([f"{transcript_id}" for transcript_id in gene_data['transcripts'].keys()])

                ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))  # Ensure genomic positions are integers
                ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{int(x)}'))  # Format x-axis ticks as integers

                # Plot each transcript
                for i, (transcript_id, transcript_info) in enumerate(gene_data['transcripts'].items()):
                    # Determine the direction based on the gene's strand information
                    direction_marker = '>' if gene_data['strand'] == '+' else '<'

                    # Add a direction marker to indicate the direction of the transcript
                    marker_pos = transcript_info['end'] + 100 if gene_data['strand'] == '+' else transcript_info['start'] - 100
                    ax.plot(marker_pos, i, marker=direction_marker, markersize=5, color="blue")

                    # Draw the line for the whole transcript
                    ax.plot([transcript_info['start'], transcript_info['end']], [i, i], color="grey", linewidth=2)

                    # Exon blocks
                    for exon in transcript_info['exons']:
                        exon_length = exon['end'] - exon['start']
                        ax.add_patch(plt.Rectangle((exon['start'], i - 0.4), exon_length, 0.8, color="skyblue"))

                ax.set_xlim(gene_data['start'], gene_data['end'])
                ax.invert_yaxis()  # First transcript at the top

                plt.tight_layout()
                plot_path = os.path.join(self.visualization_dir, f'{gene_name}_splicing.png')
                plt.savefig(plot_path)  # Saving plot by gene name
                plt.close(fig)

    def plot_transcript_usage(self):
        """
        Visualize transcript usage for each gene in gene_names across different conditions.
        """
        
        for gene_name in self.gene_names:
            gene_data = {}
            for condition, genes in self.updated_gene_dict.items():
                if gene_name in genes:
                    gene_data[condition] = genes[gene_name]['transcripts']

            if not gene_data:
                print(f"Gene {gene_name} not found in the data.")
                continue
            
            conditions = list(gene_data.keys())
            n_bars = len(conditions)
            
            fig, ax = plt.subplots(figsize=(12, 8))
            index = np.arange(n_bars)
            bar_width = 0.35
            opacity = 0.8
            
            for sample_type, transcripts in gene_data.items():
                print(f"Sample Type: {sample_type}")
                for transcript_id, transcript_info in transcripts.items():
                    print(f"  Transcript ID: {transcript_id}, Value: {transcript_info['value']}")
            # Adjusting the colors for better within-bar comparison
            max_transcripts = max(len(gene_data[condition]) for condition in conditions)
            colors = plt.cm.plasma(np.linspace(0, 1, num=max_transcripts))  # Using plasma for better color gradation
            
            bottom_val = np.zeros(n_bars) 
            for i, condition in enumerate(conditions):
                transcripts = gene_data[condition]
                for j, (transcript_id, transcript_info) in enumerate(transcripts.items()):
                    color = colors[j % len(colors)]
                    value = transcript_info['value']
                    plt.bar(i, float(value), bar_width, bottom=bottom_val[i], alpha=opacity, color=color, label=transcript_id if i == 0 else "")
                    bottom_val[i] += float(value)
            
            plt.xlabel('Sample Type')
            plt.ylabel('Transcript Usage (TPM)')
            plt.title(f'Transcript Usage for {gene_name} by Sample Type')
            plt.xticks(index, conditions)
            plt.legend(title="Transcript IDs", bbox_to_anchor=(1.05, 1), loc='upper left')
            
            plt.tight_layout()
            plot_path = os.path.join(self.visualization_dir, f'{gene_name}_transcript_usage_by_sample_type.png')
            plt.savefig(plot_path)
            plt.close(fig)




def make_pie_chart():

    data = {
        "ambiguous": 236646,
        "inconsistent": 1212565,
        "intergenic": 6886,
        "noninformative": 130493,
        "unique": 745194,
        "unique_minor_difference": 79178,
        }
    labels = data.keys()
    sizes = data.values()
    total = sum(sizes)
    print(total)
    plt.pie(sizes, labels=labels, autopct='%1.1f%%')
    plt.axis('equal')
    plt.title(f"Total: {total}")
    plt.show()
    plt.savefig('read_assignment_pie_chart.png')







def visualize_transcript_usage_single_gene(gene_data, gene_name):
    """

    :param gene_data: Dict containing transcript usage data for a given gene across different sample types.
    :param gene_name: The gene to visualize.
    """
    
    if gene_name not in gene_data:
        print(f"Gene {gene_name} not found in the data.")
        return
    
    sample_types = gene_data[gene_name].keys()
    n_bars = len(sample_types)
    
    fig, ax = plt.subplots(figsize=(10, 7))
    index = np.arange(n_bars)
    bar_width = 0.35
    opacity = 0.8
    
    # Adjusting the colors for better within-bar comparison
    max_transcripts = max(len(gene_data[gene_name][sample]) for sample in sample_types)
    colors = plt.cm.plasma(np.linspace(0, 1, num=max_transcripts))  # Using plasma for better color gradation
    
    bottom_val = np.zeros(n_bars) 
    for i, sample_type in enumerate(sample_types):
        transcripts = gene_data[gene_name][sample_type]
        for j, (transcript_id, value) in enumerate(transcripts):
            color = colors[j % len(colors)]
            plt.bar(i, float(value), bar_width, bottom=bottom_val[i], alpha=opacity, color=color, label=transcript_id if i == 0 else "")
            bottom_val[i] += float(value)
    
    plt.xlabel('Sample Type')
    plt.ylabel('Transcript Usage (TPM)')
    plt.title(f'Transcript Usage for {gene_name} by Sample Type')
    plt.xticks(index, sample_types)
    plt.legend(title="Transcript IDs", bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.show()
    plt.savefig(f'{gene_name}_transcript_usage_by_sample_type_ref.png') 