from src.post_process import OutputConfig, DictionaryBuilder
from src.plot_output import PlotOutput
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description="Visualize your IsoQuant output.")
    parser.add_argument("output_directory", type=str, help="Directory containing IsoQuant output files.")
    parser.add_argument("--gtf", type=str, help="Optional path to a GTF file if unable to be extracted from IsoQuant log", default=None)
    parser.add_argument("--counts", action="store_true", help="Use counts instead of TPM files.")
    parser.add_argument("--ref_only", action="store_true", help="Use only reference transcript quantification instead of transcript model quantification.")
    parser.add_argument("--filter_transcripts", type=float, help="Filter transcripts by minimum value occuring in at least one condition.", default=None)
    parser.add_argument("--gene_list", type=str, required=True, help="Path to a .txt file containing a list of genes, each on its own line.")
    return parser.parse_args()


def main():
    args = parse_arguments()
    output = OutputConfig(args.output_directory, use_counts=args.counts, ref_only=args.ref_only, gtf=args.gtf)
    dictionary_builder = DictionaryBuilder(output)
    gene_list = dictionary_builder.read_gene_list(args.gene_list)
    update_names = not all(gene.startswith("ENS") for gene in gene_list)
    gene_dict = dictionary_builder.build_gene_transcript_exon_dictionaries()
    reads_and_class = dictionary_builder.build_read_assignment_and_classification_dictionaries()

    if output.conditions:
        gene_file = output.gene_grouped_tpm if not output.use_counts else output.gene_grouped_counts
    else:
        gene_file = output.gene_tpm if not output.use_counts else output.gene_counts

    updated_gene_dict = dictionary_builder.update_gene_dict(gene_dict, gene_file)

    if update_names:
        print("Updating gene names to gene symbols.")
        updated_gene_dict = dictionary_builder.update_gene_names(updated_gene_dict)

    if output.ref_only or not output.extended_annotation:
        print("Using reference-only based quantification.")
        if output.conditions:
            updated_gene_dict = dictionary_builder.update_transcript_values(updated_gene_dict, output.transcript_grouped_tpm if not output.use_counts else output.transcript_grouped_counts)
        else:
            updated_gene_dict = dictionary_builder.update_transcript_values(updated_gene_dict, output.transcript_tpm if not output.use_counts else output.transcript_counts)
    else:
        print("Using transcript model quantification.")
        if output.conditions:
            updated_gene_dict = dictionary_builder.update_transcript_values(updated_gene_dict, output.transcript_model_grouped_tpm if not output.use_counts else output.transcript_model_grouped_counts)
        else:
            updated_gene_dict = dictionary_builder.update_transcript_values(updated_gene_dict, output.transcript_model_tpm if not output.use_counts else output.transcript_model_counts)

    if args.filter_transcripts is not None:
        print(f"Filtering transcripts with minimum value {args.filter_transcripts} in at least one condition.")
        updated_gene_dict = dictionary_builder.filter_transcripts_by_minimum_value(updated_gene_dict, min_value=args.filter_transcripts)
    else:
        updated_gene_dict = dictionary_builder.filter_transcripts_by_minimum_value(updated_gene_dict)

    plot_output = PlotOutput(updated_gene_dict, gene_list, args.output_directory, reads_and_class, filter_transcripts=args.filter_transcripts, conditions=output.conditions, use_counts=args.counts)
    plot_output.plot_transcript_map()
    plot_output.plot_transcript_usage()


if __name__ == "__main__":
    main()
