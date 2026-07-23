# IsoQuant output files

IsoQuant output files will be stored in `<output_dir>`, which is set by the user.
If the output directory was not specified the files are stored in `isoquant_output`.

IsoQuant consists of two stages, which generate its own output:
1. Reference-based analysis. Runs only if reference annotation is provided. Performs read-to-isoform assignment,
splice site correction and abundance quantification for reference genes/transcripts.
2. Transcript discovery. Reconstructs transcript models and performs abundance quantification for discovered isoforms.

## Reference-based analysis output

_Will be produced only if a reference gene annotation is provided._

* `SAMPLE_ID.read_info.tsv.gz` - TSV file with unified per-read information including assignments, exon coordinates, and barcode/UMI data (default output, gzipped by default);
* `SAMPLE_ID.read_assignments.tsv.gz` - legacy TSV file with read to isoform assignments (only with `--large_output read_assignments`);
* `SAMPLE_ID.corrected_reads.bed.gz` - BED file with corrected read alignments (only with `--large_output corrected_bed`);
* `SAMPLE_ID.transcript_counts.tsv` - TSV file with raw read counts for reference transcript;
* `SAMPLE_ID.gene_counts.tsv` - TSV file with raw read counts for reference genes;
* `SAMPLE_ID.transcript_tpm.tsv` - TSV file with reference transcript expression in TPM;
* `SAMPLE_ID.gene_tpm.tsv` - TSV file with reference gene expression in TPM;


If `--sqanti_output` is set, IsoQuant will produce output in [SQANTI](https://github.com/ConesaLab/SQANTI3)-like format:

* `SAMPLE_ID.novel_vs_known.SQANTI-like.tsv` - discovered novel transcripts vs reference transcripts (similar, but not identical to SQANTI `classification.txt`);

If `--analysis exon_quantification` is set (or the deprecated `--count_exons`), exon, splice junction and intron retention counts will be produced:

* `SAMPLE_ID.exon_counts.tsv` - region-based exon counts: overlapping reference exons are grouped into regions, and each region reports per-variant inclusion counts plus one region-level exclusion count;
* `SAMPLE_ID.exon_splice_site_counts.tsv` - exon splice-site counts: per-candidate full / left / right splice-site support and per-region exclusion / ambiguous counts, one row per feature and group;
* `SAMPLE_ID.splice_junction_counts.tsv` - reference splice junction inclusion/exclusion read counts;

The legacy per-exon inclusion/exclusion counts (previous IsoQuant exon format) are no longer produced by default. Use `--old_exon_count_format` to additionally output them as `SAMPLE_ID.old_exon_counts.tsv` (deprecated, will be removed in a future release).

As part of exon quantification (`--analysis exon_quantification`; deprecated aliases `--count_exons` / `--count_intron_retentions`), intron retention counts are also produced:

* `SAMPLE_ID.intron_retention_counts.tsv` - intron retention event counts per reference intron (same format as exon/splice junction counts);

#### PolyA / TSS site prediction

Whenever a gene annotation is supplied via `--genedb`, IsoQuant accumulates
per-transcript histograms of read polyA positions, detects peaks, filters them
with a pre-trained XGBoost classifier, and classifies each retained peak as
`Known` (within 10 bp of the annotated transcript end) or `Novel`:

* `SAMPLE_ID.polyA_prediction.tsv` - predicted polyA cleavage sites per reference transcript.

If `--fl_data` is also supplied (reads represent full-length transcripts), the
same machinery is applied to read start positions:

* `SAMPLE_ID.TSS_prediction.tsv` - predicted transcription start sites per reference transcript.

Columns: `chromosome`, `transcript_id`, `gene_id`, `prediction` (genomic
coordinate of the peak), `counts` (number of reads supporting the peak window),
`flag` (`Known` / `Novel`). When `--read_group` is set, an additional file per
grouping strategy is produced with two extra columns, `counts_byGroup` and
`group_id`:

* `SAMPLE_ID.polyA_prediction_grouped_<strategy>`
* `SAMPLE_ID.TSS_prediction_grouped_<strategy>` (only with `--fl_data`)

If `--read_group` is set or multiple files are provided, the per-group expression values for reference features will be also computed:

#### Default grouped counts in linear format

Grouped count file names contain the grouping strategy as `<strategy>` (e.g. `file_name`, `barcode`, `barcode_spot`, or `file0_col1` for a `file:...:0:1` group). When several `--read_group` strategies are given, one set of files is produced per strategy.

* `SAMPLE_ID.gene_grouped_<strategy>_counts.linear.tsv`
* `SAMPLE_ID.transcript_grouped_<strategy>_counts.linear.tsv`
* `SAMPLE_ID.exon_grouped_<strategy>_counts.linear.tsv`
* `SAMPLE_ID.exon_splice_site_grouped_<strategy>_counts.linear.tsv`
* `SAMPLE_ID.splice_junction_grouped_<strategy>_counts.linear.tsv`
* `SAMPLE_ID.intron_retention_grouped_<strategy>_counts.linear.tsv` (only with exon quantification)
* `SAMPLE_ID.old_exon_grouped_<strategy>_counts.linear.tsv` (only with `--old_exon_count_format`)

The region-based `exon` and `exon_splice_site` grouped counts are produced in linear format only (they are not converted to matrix/MTX format).

The exon splice-site counts carry a `group_id` column (one row per feature and group; `NA` when ungrouped). To reconstruct the per-molecule group-list format (one entry per read, e.g. barcodes for downstream cell/cell-type aggregation) use `isoquant_lib/scripts/exon_splice_site_to_group_lists.py`.

Note that grouped counts can be converted to any format using `{IsoQuant installation folder}/isoquant_lib/quantification/convert_grouped_counts.py`.
The script accepts the following arguments:

`--output` or `-o`
    Output prefix name;

`--input` or `-i`
    Path to counts files in linear IsoQuant format;

`--genedb` or `-g`
    Gene annotation in gffutils `.db` format (can be found in IsoQuant log), feature names will be used instead of IDs if provided; works only for genes and transcripts;

`--feature_type {gene,transcript,exon,intron}`
    Feature type to be converted [gene, transcript, exon, intron]; annotation lookup applies only to genes/transcripts;

`--output_format {mtx,matrix}` or `-f {mtx,matrix}`
    Output format; `matrix` is a simple TSV matrix (not recommended for large matrices), `mtx` is a Seurat-compatible MTX format;
`--tpm`
    Convert counts to TPM (works only for genes and transcripts);

`--gzip`
    Gzip output files.


#### Other formats
By default, IsoQuant converts grouped counts with small number of groups/samples (<=100) to standard matrix format; 
larger matrices (e.g. for single-cell experiments) will be saved to MTX.
See [options](cmd.md#specific-output-options) for details.

* `SAMPLE_ID.gene_grouped_<strategy>_counts.tsv` - grouped gene counts in standard matrix format;
* `SAMPLE_ID.transcript_grouped_<strategy>_counts.tsv` - grouped transcript counts in standard matrix format;
* `SAMPLE_ID.gene_grouped_<strategy>_tpm.tsv` - grouped gene TPM values in standard matrix format;
* `SAMPLE_ID.transcript_grouped_<strategy>_tpm.tsv` - grouped transcript TPM values in standard matrix format;
* `SAMPLE_ID.splice_junction_grouped_<strategy>_counts.tsv` - grouped splice junction counts in standard matrix format; row IDs are `chr:start-end:strand`, each cell holds `include,exclude` (comma-separated); produced with exon quantification;
* `SAMPLE_ID.intron_retention_grouped_<strategy>_counts.tsv` - grouped intron retention counts in standard matrix format; same layout as splice junction counts; produced with exon quantification;
* `SAMPLE_ID.old_exon_grouped_<strategy>_counts.tsv` - grouped legacy per-exon counts in standard matrix format; same layout; only with `--old_exon_count_format`;

* `SAMPLE_ID.gene_grouped_<strategy>_counts.matrix.mtx`, `SAMPLE_ID.gene_grouped_<strategy>_counts.features.tsv`, `SAMPLE_ID.gene_grouped_<strategy>_counts.barcodes.tsv` - grouped gene counts in Seurat-compatible MTX format;
* `SAMPLE_ID.transcript_grouped_<strategy>_counts.matrix.mtx`, `SAMPLE_ID.transcript_grouped_<strategy>_counts.features.tsv`, `SAMPLE_ID.transcript_grouped_<strategy>_counts.barcodes.tsv` - grouped transcript counts in Seurat-compatible MTX format;
* `SAMPLE_ID.gene_grouped_<strategy>_tpm.matrix.mtx`, `SAMPLE_ID.gene_grouped_<strategy>_tpm.features.tsv`, `SAMPLE_ID.gene_grouped_<strategy>_tpm.barcodes.tsv` - grouped gene TPM values in Seurat-compatible MTX format;
* `SAMPLE_ID.transcript_grouped_<strategy>_tpm.matrix.mtx`, `SAMPLE_ID.transcript_grouped_<strategy>_tpm.features.tsv`, `SAMPLE_ID.transcript_grouped_<strategy>_tpm.barcodes.tsv` - grouped transcript TPM values in Seurat-compatible MTX format;
* `SAMPLE_ID.splice_junction_grouped_<strategy>_counts.include.matrix.mtx`, `SAMPLE_ID.splice_junction_grouped_<strategy>_counts.exclude.matrix.mtx`, `SAMPLE_ID.splice_junction_grouped_<strategy>_counts.features.tsv`, `SAMPLE_ID.splice_junction_grouped_<strategy>_counts.barcodes.tsv` - grouped splice junction counts in Seurat-compatible MTX format; one features file and one barcodes file are shared between the include and exclude matrices; produced with exon quantification;
* `SAMPLE_ID.intron_retention_grouped_<strategy>_counts.include.matrix.mtx`, `SAMPLE_ID.intron_retention_grouped_<strategy>_counts.exclude.matrix.mtx`, `SAMPLE_ID.intron_retention_grouped_<strategy>_counts.features.tsv`, `SAMPLE_ID.intron_retention_grouped_<strategy>_counts.barcodes.tsv` - grouped intron retention counts in Seurat-compatible MTX format; produced with exon quantification;
* `SAMPLE_ID.old_exon_grouped_<strategy>_counts.include.matrix.mtx`, `SAMPLE_ID.old_exon_grouped_<strategy>_counts.exclude.matrix.mtx`, `SAMPLE_ID.old_exon_grouped_<strategy>_counts.features.tsv`, `SAMPLE_ID.old_exon_grouped_<strategy>_counts.barcodes.tsv` - grouped legacy per-exon counts in Seurat-compatible MTX format; only with `--old_exon_count_format`;


## Transcript discovery output

_Produced only when `transcript_discovery` is among the requested `--analysis` values (on by default with `--genedb` in bulk mode, and without `--genedb`). Not produced when transcript discovery is disabled (e.g. `--analysis quant`, or the deprecated `--no_model_construction`)._

File names typically contain `transcript_model` in their name.

* `SAMPLE_ID.transcript_models.gtf` - GTF file with discovered expressed transcript (both known and novel transcripts);
* `SAMPLE_ID.transcript_model_reads.tsv.gz` - which reads contributed to which transcript models, in the same unified read_info format as `SAMPLE_ID.read_info.tsv` (columns: read_id, chr, strand, gene_id, isoform_id, assignment types, exons, polyA/CAGE/barcode/UMI/groups, etc.); the `isoform_id` column holds the assigned transcript model, and reads on no model appear as noninformative (gzipped by default);
* `SAMPLE_ID.discovered_transcript_counts.tsv` - raw read counts for discovered transcript models (corresponds to `SAMPLE_ID.transcript_models.gtf`);
* `SAMPLE_ID.discovered_gene_counts.tsv` - raw read counts for discovered genes (corresponds to `SAMPLE_ID.transcript_models.gtf`);
* `SAMPLE_ID.discovered_transcript_tpm.tsv` - expression of discovered transcripts models in TPM (corresponds to `SAMPLE_ID.transcript_models.gtf`);
* `SAMPLE_ID.discovered_gene_tpm.tsv` - expression of discovered genes in TPM (corresponds to `SAMPLE_ID.transcript_models.gtf`);
* `SAMPLE_ID.extended_annotation.gtf` - GTF file with the entire reference annotation plus all discovered novel transcripts;


If `--read_group` is set, the per-group counts for discovered transcripts will be also computed:

* `SAMPLE_ID.discovered_transcript_grouped_<strategy>_counts.linear.tsv`
* `SAMPLE_ID.discovered_gene_grouped_<strategy>_counts.linear.tsv`

Similarly to the reference-based counts, these counts are converted to other formats as described [above](#other-formats).


If multiple experiments are provided, aggregated expression matrices will be placed in `<output_dir>`:

* `combined_gene_counts.tsv`
* `combined_gene_tpm.tsv`
* `combined_transcript_counts.tsv`
* `combined_transcript_tpm.tsv`


## Fusion detection output

_Produced only when `fusion` is among the requested `--analysis` values (deprecated alias: `--fusion`)._

Fusion detection runs after the isoform pipeline and reports candidate gene fusions for each input BAM file:

* `fusion_<BAM_NAME>.tsv` - detected gene fusions for the respective input BAM, placed directly in `<output_dir>` (one file per BAM);

See [output formats](formats.md#fusion-prediction-format) for a description of the columns.


Additionally, an `isoquant.log` log file will be saved to the output directory.  

If raw reads were provided, BAM file(s) will be stored in `<output_dir>/<SAMPLE_ID>/aux/`.  
In case `--keep_tmp` option was specified this directory will also contain temporary files.

