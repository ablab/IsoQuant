#!/usr/bin/env python3

import argparse
import sys
import logging
from src.visualization_output_config import OutputConfig
from src.visualization_dictionary_builder import DictionaryBuilder
from src.visualization_plotter import PlotOutput
from src.visualization_differential_exp import DifferentialAnalysis
from src.visualization_gsea import GSEAAnalysis
from src.visualization_simple_ranker import SimpleGeneRanker
from pathlib import Path


def setup_logging(viz_output_dir: Path) -> None:
    """Configure centralized logging for all visualization processes."""
    log_file = viz_output_dir / "visualize.log"

    # Create formatters
    file_formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(module)s - %(funcName)s - %(levelname)s - %(message)s'
    )
    console_formatter = logging.Formatter('%(levelname)s: %(message)s')

    # File handler - detailed logging
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(file_formatter)

    # Console handler - less detailed
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO) # Console output at INFO level
    console_handler.setFormatter(console_formatter)

    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG) # Root logger at DEBUG level
    root_logger.handlers = []  # Clear existing handlers
    root_logger.addHandler(file_handler)
    root_logger.addHandler(console_handler)


    logging.info("Initialized centralized logging system")
    logging.debug(f"Log file location: {log_file}")


def setup_viz_output(output_directory: str, viz_output: str = None) -> Path:
    """Set up visualization output directory."""
    if viz_output:
        viz_output_dir = Path(viz_output)
    else:
        viz_output_dir = Path(output_directory) / "visualization"
    viz_output_dir.mkdir(parents=True, exist_ok=True)
    return viz_output_dir


class FindGenesAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values is None:
            values = 100  # Default if flag used without value
        setattr(namespace, self.dest, values)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Visualize your IsoQuant output.")

    # Positional Argument
    parser.add_argument(
        "output_directory", type=str, help="Directory containing IsoQuant output files."
    )

    # Optional Arguments
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
        "--gsea",
        action="store_true",
        help="Perform GSEA analysis on differential expression results",
    )
    parser.add_argument(
        "--technical_replicates",
        type=str,
        help="Technical replicate specification. Can be a file path (.txt/.csv) with 'sample,group' format, or inline format 'sample1:group1,sample2:group1,sample3:group2'",
        default=None,
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--gene_list",
        type=str,
        help="Path to a .txt file containing a list of genes to evaluate.",
    )
    group.add_argument(
        "--find_genes",
        nargs="?",
        const=100,
        type=int,
        help="Find top genes with highest combined rank (default 100).",
    )

    args = parser.parse_args()

    if args.find_genes is not None:
        output = OutputConfig(
            args.output_directory,
            ref_only=args.ref_only,
            gtf=args.gtf,
        )
        if output.conditions:
            gene_file = output.transcript_grouped_tpm
        else:
            gene_file = output.transcript_tpm

        if not gene_file or not Path(gene_file).is_file():
            print(f"Error: Grouped TPM/Counts file not found at {gene_file}.")
            sys.exit(1)

        with open(gene_file, "r") as f:
            header = f.readline().strip().split("\t")

        if len(header) < 2:
            print(
                "Error: The grouped TPM/Counts file does not contain condition information."
            )
            sys.exit(1)

        available_conditions = header[1:]
        if not available_conditions:
            print("Error: No conditions found in the grouped TPM/Counts file.")
            sys.exit(1)

        args.available_conditions = available_conditions

    return args


def select_conditions_interactively(args):
    print("\nAvailable conditions:")
    for idx, condition in enumerate(args.available_conditions, 1):
        print(f"{idx}. {condition}")

    def get_selection(prompt, max_selection, exclude=[]):
        while True:
            try:
                choices = input(prompt)
                choice_indices = [int(x.strip()) for x in choices.split(",")]
                if all(1 <= idx <= max_selection for idx in choice_indices):
                    selected = [
                        args.available_conditions[idx - 1]
                        for idx in choice_indices
                        if args.available_conditions[idx - 1] not in exclude
                    ]
                    if not selected:
                        print("No valid conditions selected. Please try again.")
                        continue
                    return selected
                else:
                    print(f"Please enter numbers between 1 and {max_selection}.")
            except ValueError:
                print("Invalid input. Please enter numbers separated by commas.")

    max_idx = len(args.available_conditions)
    args.reference_conditions = get_selection(
        "\nEnter refs (comma-separated): ", max_idx
    )
    selected_refs = set(args.reference_conditions)
    args.target_conditions = get_selection(
        "\nEnter targets (comma-separated): ", max_idx, exclude=selected_refs
    )

    print("\nSelected Reference Conditions:", ", ".join(args.reference_conditions))
    print("Selected Target Conditions:", ", ".join(args.target_conditions), "\n")





def main():
    # First, parse just the output directory argument to set up logging
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("output_directory", type=str, nargs='?')
    parser.add_argument("--viz_output", type=str, default=None)
    
    # Parse just these arguments first
    first_args, _ = parser.parse_known_args()
    
    # Initialize output directory early
    if not first_args.output_directory:
        print("Error: Output directory is required.")
        sys.exit(1)
        
    # Set up visualization directory early
    viz_output_dir = setup_viz_output(first_args.output_directory, first_args.viz_output)
    
    # Set up logging immediately to capture all operations
    setup_logging(viz_output_dir)
    logging.info("Starting IsoQuant visualization pipeline")
    
    # Now parse the full arguments with the real parser
    args = parse_arguments()
    
    # If find_genes is specified, get conditions interactively
    if args.find_genes is not None:
        select_conditions_interactively(args)

    logging.info("Reading IsoQuant parameters.")
    output = OutputConfig(
        args.output_directory,
        ref_only=args.ref_only,
        gtf=args.gtf,
        technical_replicates=args.technical_replicates,
    )
    dictionary_builder = DictionaryBuilder(output)
    logging.debug("OutputConfig details:")
    logging.debug(f"Output directory: {output.output_directory}")
    logging.debug(f"Reference only: {output.ref_only}")

    # Ask user about read assignments (optional)
    use_read_assignments = (
        input("Do you want to look at read_assignment data? (y/n): ")
        .strip()
        .lower()
        .startswith("y")
    )

    # If gene_list was given, read it; might use later for some optional steps
    if args.gene_list:
        logging.info(f"Reading gene list from {args.gene_list}")
        gene_list = dictionary_builder.read_gene_list(args.gene_list)
        # Decide if you need to rename Genes -> Symbol
        update_names = not all(gene.startswith("ENS") for gene in gene_list)
    else:
        gene_list = None
        update_names = True

    min_val = args.filter_transcripts if args.filter_transcripts is not None else 1.0
    logging.info(f"FLOW_DEBUG: Building updated_gene_dict with:")
    logging.info(f"  min_value: {min_val}")
    logging.info(f"  reference_conditions: {getattr(args, 'reference_conditions', None)}")
    logging.info(f"  target_conditions: {getattr(args, 'target_conditions', None)}")
    
    updated_gene_dict = dictionary_builder.build_gene_dict_with_expression_and_filter(
        min_value=min_val,
        reference_conditions=getattr(args, 'reference_conditions', None),
        target_conditions=getattr(args, 'target_conditions', None)
    )
    
    logging.info(f"FLOW_DEBUG: updated_gene_dict created:")
    logging.info(f"  type: {type(updated_gene_dict)}")
    logging.info(f"  keys (conditions): {list(updated_gene_dict.keys()) if updated_gene_dict else 'None'}")
    if updated_gene_dict:
        for condition, genes in updated_gene_dict.items():
            logging.info(f"  condition '{condition}': {len(genes)} genes")
            sample_genes = list(genes.keys())[:3]
            if sample_genes:
                for gene_id in sample_genes:
                    gene_info = genes[gene_id]
                    logging.info(f"    gene '{gene_id}': name='{gene_info.get('name', 'MISSING')}', keys={list(gene_info.keys())}")
                    if 'transcripts' in gene_info:
                        logging.info(f"      transcripts: {len(gene_info['transcripts'])} items")
            break  # Only show details for first condition

    # Debug: log whether gene_dict keys are Ensembl IDs or gene names
    if updated_gene_dict:
        sample_condition = next(iter(updated_gene_dict))
        sample_keys = list(updated_gene_dict[sample_condition].keys())[:5]
        logging.info(
            "Sample gene_dict keys for condition '%s': %s", sample_condition, sample_keys
        )

    # 2. If read assignments are desired, build those as well (cached)
    if use_read_assignments:
        logging.info("Building read assignment and classification dictionaries.")
        reads_and_class = (
            dictionary_builder.build_read_assignment_and_classification_dictionaries()
        )
    else:
        reads_and_class = None

    # 3. If user wants to find top genes (--find_genes), choose method based on replicate availability
    if args.find_genes is not None:
        ref_str = "_".join(x.upper().replace(" ", "_") for x in args.reference_conditions)
        target_str = "_".join(x.upper().replace(" ", "_") for x in args.target_conditions)
        main_dir_name = f"find_genes_{ref_str}_vs_{target_str}"
        base_dir = viz_output_dir / main_dir_name if not args.viz_output else viz_output_dir
        base_dir.mkdir(exist_ok=True)

        tech_rep_dict = output.technical_replicates_dict
        replicate_ok = output.check_biological_replicates_for_conditions(
            args.reference_conditions, args.target_conditions
        )

        if replicate_ok:
            logging.info("Finding genes via DESeq2 (replicates detected).")
            diff_analysis = DifferentialAnalysis(
                output_dir=output.output_directory,
                viz_output=base_dir,
                ref_conditions=args.reference_conditions,
                target_conditions=args.target_conditions,
                updated_gene_dict=updated_gene_dict,
                ref_only=args.ref_only,
                dictionary_builder=dictionary_builder,
                tech_rep_dict=tech_rep_dict,
            )
            gene_results, transcript_results, _, deseq2_df = diff_analysis.run_complete_analysis()

            if args.gsea:
                gsea = GSEAAnalysis(output_path=base_dir)
                target_label = f"{'+'.join(args.target_conditions)}_vs_{'+'.join(args.reference_conditions)}"
                gsea.run_gsea_analysis(deseq2_df, target_label)

            # Path to DESeq2-derived top genes
            top_n = args.find_genes
            contrast_label = f"{'+'.join(args.target_conditions)}_vs_{'+'.join(args.reference_conditions)}"
            top_genes_filename = f"genes_of_top_{top_n}_DE_transcripts_{contrast_label}.txt"
            find_genes_list_path = gene_results.parent / top_genes_filename
            logging.info(f"Reading gene list generated by differential analysis from: {find_genes_list_path}")
            
            logging.info(f"FLOW_DEBUG: DESeq2 path - reading from file: {find_genes_list_path}")
            if find_genes_list_path.exists():
                with open(find_genes_list_path, 'r') as f:
                    file_contents = f.read().strip().split('\n')
                logging.info(f"FLOW_DEBUG: DESeq2 file has {len(file_contents)} lines, first 5: {file_contents[:5]}")
            else:
                logging.error(f"FLOW_DEBUG: DESeq2 gene list file does not exist: {find_genes_list_path}")
            
            gene_list = dictionary_builder.read_gene_list(find_genes_list_path)
            logging.info(f"FLOW_DEBUG: DESeq2 gene_list after dictionary_builder.read_gene_list:")
            logging.info(f"  type: {type(gene_list)}")
            logging.info(f"  length: {len(gene_list) if gene_list else 'None'}")
            logging.info(f"  content (first 10): {gene_list[:10] if gene_list else 'None'}")
        else:
            logging.info("No biological replicates detected – using SimpleGeneRanker.")
            logging.info(f"FLOW_DEBUG: Creating SimpleGeneRanker with:")
            logging.info(f"  output_dir: {output.output_directory}")
            logging.info(f"  ref_conditions: {args.reference_conditions}")
            logging.info(f"  target_conditions: {args.target_conditions}")
            logging.info(f"  ref_only: {args.ref_only}")
            logging.info(f"  updated_gene_dict keys: {list(updated_gene_dict.keys()) if updated_gene_dict else 'None'}")
            
            simple_ranker = SimpleGeneRanker(
                output_dir=output.output_directory,
                ref_conditions=args.reference_conditions,
                target_conditions=args.target_conditions,
                ref_only=args.ref_only,
                updated_gene_dict=updated_gene_dict,
            )
            
            logging.info(f"FLOW_DEBUG: Calling simple_ranker.rank(top_n={args.find_genes})")
            gene_list = simple_ranker.rank(top_n=args.find_genes)
            logging.info(f"FLOW_DEBUG: SimpleGeneRanker returned gene_list with {len(gene_list)} genes")
            logging.info(f"FLOW_DEBUG: Gene list type: {type(gene_list)}")
            logging.info(f"FLOW_DEBUG: Gene list content (first 10): {gene_list[:10] if gene_list else 'EMPTY'}")
            
            # Write gene list to file for reproducibility
            contrast_label = f"{'+'.join(args.target_conditions)}_vs_{'+'.join(args.reference_conditions)}"
            top_genes_filename = f"genes_of_top_{args.find_genes}_simple_{contrast_label}.txt"
            simple_list_path = base_dir / top_genes_filename
            import pandas as _pd
            _pd.Series(gene_list).to_csv(simple_list_path, index=False, header=False)
            logging.info(f"Simple gene list written to {simple_list_path}")
            logging.info(f"FLOW_DEBUG: File contents verification:")
            try:
                with open(simple_list_path, 'r') as f:
                    file_contents = f.read().strip().split('\n')
                logging.info(f"FLOW_DEBUG: File has {len(file_contents)} lines, first 5: {file_contents[:5]}")
            except Exception as e:
                logging.error(f"FLOW_DEBUG: Error reading written file: {e}")
    else:
        base_dir = viz_output_dir

    # 5. Set up output directories
    gene_visualizations_dir = base_dir / "gene_visualizations"
    gene_visualizations_dir.mkdir(exist_ok=True)

    if use_read_assignments:
        read_assignments_dir = base_dir / "read_assignments"
        read_assignments_dir.mkdir(exist_ok=True)
    else:
        read_assignments_dir = None # Set to None if not used

    # 6. Plotting with PlotOutput
    logging.info(f"FLOW_DEBUG: Creating PlotOutput with:")
    logging.info(f"  gene_names type: {type(gene_list)}")
    logging.info(f"  gene_names length: {len(gene_list) if gene_list else 'None'}")
    logging.info(f"  gene_names content (first 10): {gene_list[:10] if gene_list else 'None'}")
    logging.info(f"  updated_gene_dict keys: {list(updated_gene_dict.keys()) if updated_gene_dict else 'None'}")
    logging.info(f"  conditions: {output.conditions}")
    logging.info(f"  filter_transcripts: {min_val}")
    logging.info(f"  ref_only: {args.ref_only}")
    
    plot_output = PlotOutput(
        updated_gene_dict=updated_gene_dict,
        gene_names=gene_list,
        gene_visualizations_dir=str(gene_visualizations_dir),
        read_assignments_dir=str(read_assignments_dir) if read_assignments_dir else None,
        reads_and_class=reads_and_class,
        filter_transcripts=min_val,
        conditions=output.conditions,
        ref_only=args.ref_only,
        ref_conditions=args.reference_conditions if hasattr(args, "reference_conditions") else None,
        target_conditions=args.target_conditions if hasattr(args, "target_conditions") else None,
    )
    

    plot_output.plot_transcript_map()
    
    plot_output.plot_transcript_usage()
  

    if use_read_assignments:
        plot_output.make_pie_charts()


if __name__ == "__main__":
    main()
