#!/usr/bin/env python3

from src.post_process import OutputConfig, DictionaryBuilder
from src.plot_output import PlotOutput
import argparse
from src.process_dict import simplify_and_sum_transcripts
from src.gene_model import rank_and_visualize_genes


class FindGenesAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values is None:
            values = 100  # Default value when the flag is used without a value
        setattr(namespace, self.dest, values)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Visualize your IsoQuant output.")
    parser.add_argument(
        "output_directory", type=str, help="Directory containing IsoQuant output files."
    )
    parser.add_argument(
        "--viz_output",
        type=str,
        help="Optional directory to save visualization output files, defaults to the main output directory.",
        default=None,
    )
    parser.add_argument(
        "--gtf",
        type=str,
        help="Optional path to a GTF file if unable to be extracted from IsoQuant log",
        default=None,
    )
    parser.add_argument(
        "--counts", action="store_true", help="Use counts instead of TPM files."
    )
    parser.add_argument(
        "--ref_only",
        action="store_true",
        help="Use only reference transcript quantification instead of transcript model quantification.",
    )
    parser.add_argument(
        "--filter_transcripts",
        type=float,
        help="Filter transcripts by minimum value occurring in at least one condition.",
        default=None,
    )
    parser.add_argument(
        "--gene_list",
        type=str,
        required=True,
        help="Path to a .txt file containing a list of genes, each on its own line.",
    )
    parser.add_argument(
        "--find_genes",
        nargs="?",
        const=100,
        type=int,
        help="Find genes with the highest combined rank and visualize them. Optionally specify the number of top genes to evaluate (default is 100).",
    )
    parser.add_argument(
        "--known_genes_path",
        type=str,
        help="Path to a CSV file containing known target genes.",
        default=None,
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    output = OutputConfig(
        args.output_directory,
        use_counts=args.counts,
        ref_only=args.ref_only,
        gtf=args.gtf,
    )
    dictionary_builder = DictionaryBuilder(output)
    gene_list = dictionary_builder.read_gene_list(args.gene_list)
    update_names = not all(gene.startswith("ENS") for gene in gene_list)
    gene_dict = dictionary_builder.build_gene_transcript_exon_dictionaries()
    reads_and_class = (
        dictionary_builder.build_read_assignment_and_classification_dictionaries()
    )

    if output.conditions:
        gene_file = (
            output.gene_grouped_tpm
            if not output.use_counts
            else output.gene_grouped_counts
        )
    else:
        gene_file = output.gene_tpm if not output.use_counts else output.gene_counts

    updated_gene_dict = dictionary_builder.update_gene_dict(gene_dict, gene_file)
    if update_names:
        print("Updating Ensembl IDs to gene symbols.")
        updated_gene_dict = dictionary_builder.update_gene_names(updated_gene_dict)

    if output.ref_only or not output.extended_annotation:
        print("Using reference-only based quantification.")
        if output.conditions:
            updated_gene_dict = dictionary_builder.update_transcript_values(
                updated_gene_dict,
                (
                    output.transcript_grouped_tpm
                    if not output.use_counts
                    else output.transcript_grouped_counts
                ),
            )
        else:
            updated_gene_dict = dictionary_builder.update_transcript_values(
                updated_gene_dict,
                (
                    output.transcript_tpm
                    if not output.use_counts
                    else output.transcript_counts
                ),
            )
    else:
        print("Using transcript model quantification.")
        if output.conditions:
            updated_gene_dict = dictionary_builder.update_transcript_values(
                updated_gene_dict,
                (
                    output.transcript_model_grouped_tpm
                    if not output.use_counts
                    else output.transcript_model_grouped_counts
                ),
            )
        else:
            updated_gene_dict = dictionary_builder.update_transcript_values(
                updated_gene_dict,
                (
                    output.transcript_model_tpm
                    if not output.use_counts
                    else output.transcript_model_counts
                ),
            )

    if args.filter_transcripts is not None:
        print(
            f"Filtering transcripts with minimum value {args.filter_transcripts} in at least one condition."
        )
        updated_gene_dict = dictionary_builder.filter_transcripts_by_minimum_value(
            updated_gene_dict, min_value=args.filter_transcripts
        )
    else:
        updated_gene_dict = dictionary_builder.filter_transcripts_by_minimum_value(
            updated_gene_dict
        )

    # Visualization output directory decision
    viz_output_directory = args.viz_output if args.viz_output else args.output_directory

    if args.find_genes:
        print("Finding genes.")
        simple_gene_dict = simplify_and_sum_transcripts(updated_gene_dict)
        path = rank_and_visualize_genes(
            simple_gene_dict,
            viz_output_directory,
            args.find_genes,
            known_genes_path=args.known_genes_path,
        )
        gene_list = dictionary_builder.read_gene_list(path)

    # dictionary_builder.save_gene_dict_to_json(updated_gene_dict, viz_output_directory)
    plot_output = PlotOutput(
        updated_gene_dict,
        gene_list,
        viz_output_directory,
        create_visualization_subdir=(viz_output_directory == args.output_directory),
        reads_and_class=reads_and_class,
        filter_transcripts=args.filter_transcripts,
        conditions=output.conditions,
        use_counts=args.counts,
    )
    plot_output.plot_transcript_map()
    plot_output.plot_transcript_usage()
    plot_output.make_pie_charts()


if __name__ == "__main__":
    main()
