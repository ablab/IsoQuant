# Visualization

IsoQuant provides a visualization tool to help interpret and explore the output data. The goal of this visualization is to create informative plots that represent transcript usage and splicing patterns for genes of interest. Additionally, we provide global transcript and read assignment statistics from the IsoQuant analysis.

## Running the visualization tool

To run the visualization tool, use one of the following commands:

```bash
# Visualize a predefined list of genes
python visualize.py <output_directory> --gene_list <gene_list> [options]

# Automatically find the top N most differentially expressed genes
python visualize.py <output_directory> --find_genes [N] [options]
```

## Command line options

* `output_directory` (required): Directory containing IsoQuant output files.
* `--gene_list`: Path to a .txt file containing a list of genes, each on its own line. Mutually exclusive with `--find_genes`.
* `--find_genes [N]`: Automatically select the top **N** genes with the highest combined differential-expression rank between chosen conditions (default 100 if *N* is omitted).
* `--viz_output`: Optional directory to save visualization output files. Defaults to `<output_directory>/visualization`.
* `--gtf`: Optional path to a GTF file if it cannot be extracted from the IsoQuant log.
* `--ref_only`: Use only reference transcript quantification instead of transcript model quantification.
* `--filter_transcripts <float>`: Minimum expression value a transcript must reach in at least one condition to be included in plots (default 1.0).
* `--gsea`: Perform Gene Set Enrichment Analysis on differential expression results (requires `--find_genes`).
* `--technical_replicates`: Specify technical replicate groupings as a file (`sample,group`) or inline (`sample1:group1,sample2:group1`).


## Output

The visualization tool can generate the following outputs:

1. Transcript usage profiles: For each gene, a plot showing the relative usage of different transcripts across conditions or samples.

2. Gene-specific transcript maps: Visual representation of the different splicing patterns of transcripts for each gene, allowing easy comparison of exon usage and alternative splicing events.

3. Global read assignment consistency: A summary plot showing the overall consistency of read assignments across all genes and transcripts analyzed (enabled interactively).

4. Global transcript alignment classifications: A chart representing the distribution of different transcript alignment categories (e.g., full splice match, incomplete splice match, novel isoforms) across the entire dataset.

5. Differential expression tables and volcano plots when `--find_genes` is used, with optional GSEA pathway visualizations if `--gsea` is supplied.

These visualizations and reports provide valuable insights into transcript diversity, splicing patterns, differential expression, and the overall quality of the IsoQuant analysis.
