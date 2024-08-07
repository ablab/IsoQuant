import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pprint


class PlotOutput:
    def __init__(
        self,
        updated_gene_dict,
        gene_names,
        output_directory,
        create_visualization_subdir=False,
        reads_and_class=None,
        filter_transcripts=None,
        conditions=False,
        use_counts=False,
    ):
        self.updated_gene_dict = updated_gene_dict
        self.gene_names = gene_names
        self.output_directory = output_directory
        self.reads_and_class = reads_and_class
        self.filter_transcripts = filter_transcripts
        self.conditions = conditions
        self.use_counts = use_counts

        # Create visualization subdirectory if specified
        if create_visualization_subdir:
            self.visualization_dir = os.path.join(
                self.output_directory, "visualization"
            )
            os.makedirs(self.visualization_dir, exist_ok=True)
        else:
            self.visualization_dir = self.output_directory

    def plot_transcript_map(self):
        # Get the first condition's gene dictionary
        first_condition = next(iter(self.updated_gene_dict))
        gene_dict = self.updated_gene_dict[first_condition]

        for gene_name in self.gene_names:
            if gene_name in gene_dict:
                gene_data = gene_dict[gene_name]
                num_transcripts = len(gene_data["transcripts"])
                plot_height = max(
                    3, num_transcripts * 0.3
                )  # Adjust the height dynamically

                fig, ax = plt.subplots(
                    figsize=(12, plot_height)
                )  # Adjust height dynamically

                if self.filter_transcripts is not None:
                    ax.set_title(
                        f"Transcripts of Gene: {gene_data['name']} on Chromosome {gene_data['chromosome']} with value over {self.filter_transcripts}"
                    )
                else:
                    ax.set_title(
                        f"Transcripts of Gene: {gene_data['name']} on Chromosome {gene_data['chromosome']}"
                    )

                ax.set_xlabel("Chromosomal position")
                ax.set_ylabel("Transcripts")
                ax.set_yticks(range(num_transcripts))
                ax.set_yticklabels(
                    [
                        f"{transcript_id}"
                        for transcript_id in gene_data["transcripts"].keys()
                    ]
                )

                ax.xaxis.set_major_locator(
                    ticker.MaxNLocator(integer=True)
                )  # Ensure genomic positions are integers
                ax.xaxis.set_major_formatter(
                    ticker.FuncFormatter(lambda x, pos: f"{int(x)}")
                )  # Format x-axis ticks as integers

                # Plot each transcript
                for i, (transcript_id, transcript_info) in enumerate(
                    gene_data["transcripts"].items()
                ):
                    # Determine the direction based on the gene's strand information
                    direction_marker = ">" if gene_data["strand"] == "+" else "<"
                    marker_pos = (
                        transcript_info["end"] + 100
                        if gene_data["strand"] == "+"
                        else transcript_info["start"] - 100
                    )
                    ax.plot(
                        marker_pos,
                        i,
                        marker=direction_marker,
                        markersize=5,
                        color="blue",
                    )

                    # Draw the line for the whole transcript
                    ax.plot(
                        [transcript_info["start"], transcript_info["end"]],
                        [i, i],
                        color="grey",
                        linewidth=2,
                    )

                    # Exon blocks
                    for exon in transcript_info["exons"]:
                        exon_length = exon["end"] - exon["start"]
                        ax.add_patch(
                            plt.Rectangle(
                                (exon["start"], i - 0.4),
                                exon_length,
                                0.8,
                                color="skyblue",
                            )
                        )

                ax.set_xlim(gene_data["start"], gene_data["end"])
                ax.invert_yaxis()  # First transcript at the top

                plt.tight_layout()
                plot_path = os.path.join(
                    self.visualization_dir, f"{gene_name}_splicing.png"
                )
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
                    gene_data[condition] = genes[gene_name]["transcripts"]

            if not gene_data:
                print(f"Gene {gene_name} not found in the data.")
                continue

            conditions = list(gene_data.keys())
            n_bars = len(conditions)

            fig, ax = plt.subplots(figsize=(12, 8))
            index = np.arange(n_bars)
            bar_width = 0.35
            opacity = 0.8

            # for sample_type, transcripts in gene_data.items():
            # print(f"Sample Type: {sample_type}")
            # for transcript_id, transcript_info in transcripts.items():
            # print(
            # f"  Transcript ID: {transcript_id}, Value: {transcript_info['value']}"
            # )
            # Adjusting the colors for better within-bar comparison
            max_transcripts = max(len(gene_data[condition]) for condition in conditions)
            colors = plt.cm.plasma(
                np.linspace(0, 1, num=max_transcripts)
            )  # Using plasma for better color gradation

            bottom_val = np.zeros(n_bars)
            for i, condition in enumerate(conditions):
                transcripts = gene_data[condition]
                for j, (transcript_id, transcript_info) in enumerate(
                    transcripts.items()
                ):
                    color = colors[j % len(colors)]
                    value = transcript_info["value"]
                    plt.bar(
                        i,
                        float(value),
                        bar_width,
                        bottom=bottom_val[i],
                        alpha=opacity,
                        color=color,
                        label=transcript_id if i == 0 else "",
                    )
                    bottom_val[i] += float(value)

            plt.xlabel("Sample Type")
            plt.ylabel("Transcript Usage (TPM)")
            plt.title(f"Transcript Usage for {gene_name} by Sample Type")
            plt.xticks(index, conditions)
            plt.legend(
                title="Transcript IDs", bbox_to_anchor=(1.05, 1), loc="upper left"
            )

            plt.tight_layout()
            plot_path = os.path.join(
                self.visualization_dir,
                f"{gene_name}_transcript_usage_by_sample_type.png",
            )
            plt.savefig(plot_path)
            plt.close(fig)

    def make_pie_charts(self):
        """
        Create pie charts for transcript alignment classifications and read assignment consistency.
        Handles both combined and separate sample data structures.
        """
        print("self.reads_and_class structure:")
        pprint.pprint(self.reads_and_class)

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
        plot_path = os.path.join(self.visualization_dir, f"{file_title}_pie_chart.png")
        plt.savefig(plot_path, bbox_inches="tight", dpi=300)
        plt.close()
