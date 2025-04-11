# IsoQuant output files

IsoQuant output files will be stored in `<output_dir>`, which is set by the user.
If the output directory was not specified the files are stored in `isoquant_output`.

IsoQuant consists of two stages, which generate its own output:
1. Reference-based analysis. Runs only if reference annotation is provided. Performs read-to-isofrom assignment,
splice site correction and abundance quantification for reference genes/transcripts.
2. Transcript discovery. Reconstructs transcript models and performs abundance quantification for discovered isoforms.

## Reference-based analysis output

_Will be produced only if a reference gene annotation is provided._

* `SAMPLE_ID.read_assignments.tsv.gz` - TSV file with read to isoform assignments (gzipped by default);
* `SAMPLE_ID.corrected_reads.bed.gz` - BED file with corrected read alignments (gzipped by default);
* `SAMPLE_ID.transcript_counts.tsv` - TSV file with raw read counts for reference transcript;
* `SAMPLE_ID.gene_counts.tsv` - TSV file with raw read counts for reference genes;
* `SAMPLE_ID.transcript_tpm.tsv` - TSV file with reference transcript expression in TPM;
* `SAMPLE_ID.gene_tpm.tsv` - TSV file with reference gene expression in TPM;


If `--sqanti_output` is set, IsoQuant will produce output in [SQANTI](https://github.com/ConesaLab/SQANTI3)-like format:

* `SAMPLE_ID.novel_vs_known.SQANTI-like.tsv` - discovered novel transcripts vs reference transcripts (similar, but not identical to SQANTI `classification.txt`);

If `--count_exons` is set, exon and intron counts will be produced:

* `SAMPLE_ID.exon_counts.tsv` - reference exon inclusion/exclusion read counts;
* `SAMPLE_ID.intron_counts.tsv` - reference intron inclusion/exclusion read counts;

If `--read_group` is set or multiple files are provided, the per-group expression values for reference features will be also computed:

#### Default grouped counts in linear format
* `SAMPLE_ID.gene_grouped_counts.linear.tsv`
* `SAMPLE_ID.transcript_grouped_counts.linear.tsv`
* `SAMPLE_ID.exon_grouped_counts.tsv`
* `SAMPLE_ID.intron_grouped_counts.tsv`

Note, that grouped counts can be converted to any format using `src/convert_grouped_counts.py`.

#### Other formats 
By default, IsoQuant converts grouped counts with small number of groups/samples (<=100) to standard matrix format; 
larger matrices (e.g. for single-cell experiments) will be saved to MTX.
See [options](cmd.md#specific-output-options) for details.

* `SAMPLE_ID.gene_grouped_counts.tsv` - grouped gene counts in standard matrix format;
* `SAMPLE_ID.transcript_grouped_counts.tsv` - grouped transcript counts in standard matrix format;
* `SAMPLE_ID.gene_grouped_tpm.tsv` - grouped gene TPM values in standard matrix format;
* `SAMPLE_ID.transcript_grouped_tpm.tsv` - grouped TPM values counts in standard matrix format;

* `SAMPLE_ID.gene_grouped_counts.matrix.mtx`, `SAMPLE_ID.gene_grouped_counts.features.tsv`, `SAMPLE_ID.gene_grouped_counts.barcodes.tsv` - grouped gene counts in Seurat-compatible MTX format;
* `SAMPLE_ID.transcript_grouped_counts.matrix.mtx`, `SAMPLE_ID.transcript_grouped_counts.features.tsv`, `SAMPLE_ID.transcript_grouped_counts.barcodes.tsv` - grouped transcript counts in Seurat-compatible MTX format;
* `SAMPLE_ID.gene_grouped_tpm.matrix.mtx`, `SAMPLE_ID.gene_grouped_tpm.features.tsv`, `SAMPLE_ID.gene_grouped_tpm.barcodes.tsv` - grouped gene TPM values in Seurat-compatible MTX format;
* `SAMPLE_ID.transcript_grouped_tpm.matrix.mtx`, `SAMPLE_ID.transcript_grouped_tpm.features.tsv`, `SAMPLE_ID.transcript_grouped_tpm.barcodes.tsv` - grouped transcript TPM values in Seurat-compatible MTX format;


## Transcript discovery output

_Will not be produced if `--no_model_construction` is set._

File names typically contain `transcript_model` in their name.

* `SAMPLE_ID.transcript_models.gtf` - GTF file with discovered expressed transcript (both known and novel transcripts);
* `SAMPLE_ID.transcript_model_reads.tsv.gz` - TSV file indicating which reads contributed to transcript models (gzipped by default);
* `SAMPLE_ID.transcript_model_tpm.tsv` - expression of discovered transcripts models in TPM (corresponds to `SAMPLE_ID.transcript_models.gtf`);
* `SAMPLE_ID.transcript_model_counts.tsv` - raw read counts for discovered transcript models (corresponds to `SAMPLE_ID.transcript_models.gtf`);
* `SAMPLE_ID.extended_annotation.gtf` - GTF file with the entire reference annotation plus all discovered novel transcripts;


If `--read_group` is set, the per-group counts for discovered transcripts will be also computed:

* `SAMPLE_ID.transcript_model_grouped_counts.tsv`
* `SAMPLE_ID.transcript_model_grouped_tpm.tsv`


If multiple experiments are provided, aggregated expression matrices will be placed in `<output_dir>`:

* `combined_gene_counts.tsv`
* `combined_gene_tpm.tsv`
* `combined_transcript_counts.tsv`
* `combined_transcript_tpm.tsv`

Additionally, an `isoquant.log` log file will be saved to the output directory.  

If raw reads were provided, BAM file(s) will be stored in `<output_dir>/<SAMPLE_ID>/aux/`.  
In case `--keep_tmp` option was specified this directory will also contain temporary files.

