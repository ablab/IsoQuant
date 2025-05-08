# Visualization

IsoQuant provides a visualization tool to help interpret and explore the output data. The goal of this visualization is to create informative plots that represent transcript usage and splicing patterns for genes of interest. Additionally, we provide global transcript and read assignment statistics from the IsoQuant analysis.

## Installing supplementary R packages

Visualization now includes various differential expression analysis that relies on the following R packages via `rpy2` library.

- DESeq2
- ggplot2
- ggrepel
- RColorBrewer
- BiocManager
- clusterProfiler
- org.Hs.eg.db

The last two are only available only via Bioconductor.

To install these packages simply run the following Python script:

```
intsall_r_packages.py
```

Or install them directly via R.

## Running the visualization tool

To run the visualization tool, use the following command:

```bash

visualize.py <output_directory> --gene_list <gene_list> [options]

```

## Command line options

* `output_directory` (required): Directory containing IsoQuant output files.
* * `--gene_list` (required): Path to a .txt file containing a list of genes, each on its own line.
* `--viz_output`: Optional directory to save visualization output files. Defaults to the main output directory if not specified.
* `--gtf`: Optional path to a GTF file if it cannot be extracted from the IsoQuant log.
* `--counts`: Use counts instead of TPM files for visualization.
* `--ref_only`: Use only reference transcript quantification instead of transcript model quantification.
* `--filter_transcripts`: Filter transcripts by minimum value occurring in at least one condition.


## Output

The visualization tool generates the following plots based on the IsoQuant output:

1. Transcript usage profiles: For each gene specified in the gene list, a plot showing the relative usage of different transcripts across conditions or samples.

2. Gene-specific transcript maps: Visual representation of the different splicing patterns of transcripts for each gene, allowing easy comparison of exon usage and alternative splicing events.

3. Global read assignment consistency: A summary plot showing the overall consistency of read assignments across all genes and transcripts analyzed.

4. Global transcript alignment classifications: A chart or plot representing the distribution of different transcript alignment categories (e.g., full splice match, incomplete splice match, novel isoforms) across the entire dataset.

These visualizations provide valuable insights into transcript diversity, splicing patterns, and the overall quality of the IsoQuant analysis.
