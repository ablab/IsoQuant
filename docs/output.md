# IsoQuant output files

IsoQuant output files will be stored in `<output_dir>`, which is set by the user.
If the output directory was not specified the files are stored in `isoquant_output`.

IsoQuant consists of two stages, which generate its own output:
1. Reference-based analysis. Runs only if reference annotation is provided. Performs read-to-isofrom assignment,
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

If `--count_exons` is set, exon and intron counts will be produced:

* `SAMPLE_ID.exon_counts.tsv` - reference exon inclusion/exclusion read counts;
* `SAMPLE_ID.splice_junction_counts.tsv` - reference intron inclusion/exclusion read counts;

If `--read_group` is set or multiple files are provided, the per-group expression values for reference features will be also computed:

#### Default grouped counts in linear format
* `SAMPLE_ID.gene_grouped_counts.linear.tsv`
* `SAMPLE_ID.transcript_grouped_counts.linear.tsv`
* `SAMPLE_ID.exon_grouped_counts.linear.tsv`
* `SAMPLE_ID.splice_junction_grouped_counts.linear.tsv`

Note that grouped counts can be converted to any format using `{IsoQuant intsllation folder}/isoquant_lib/convert_grouped_counts.py`.
The script accepts the following arguments:

`--output` or `-o`
    Output prefix name;

`--input` or `-i`
    Path to counts files in linear IsoQuant format;

`--genedb` or `-g`
    Gene annotation in gffutils `.db` format (can be found in IsoQuant log), feature names will be used instead of IDs if provided; works only for genes and trascncripts;

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

* `SAMPLE_ID.gene_grouped_counts.tsv` - grouped gene counts in standard matrix format;
* `SAMPLE_ID.transcript_grouped_counts.tsv` - grouped transcript counts in standard matrix format;
* `SAMPLE_ID.gene_grouped_tpm.tsv` - grouped gene TPM values in standard matrix format;
* `SAMPLE_ID.transcript_grouped_tpm.tsv` - grouped TPM values counts in standard matrix format;
* `SAMPLE_ID.exon_grouped_counts.tsv` - grouped exon counts in standard matrix format; row IDs are `chr:start-end:strand`, each cell holds `include,exclude` (comma-separated);
* `SAMPLE_ID.intron_grouped_counts.tsv` - grouped intron counts in standard matrix format; same layout as exon counts;

* `SAMPLE_ID.gene_grouped_counts.matrix.mtx`, `SAMPLE_ID.gene_grouped_counts.features.tsv`, `SAMPLE_ID.gene_grouped_counts.barcodes.tsv` - grouped gene counts in Seurat-compatible MTX format;
* `SAMPLE_ID.transcript_grouped_counts.matrix.mtx`, `SAMPLE_ID.transcript_grouped_counts.features.tsv`, `SAMPLE_ID.transcript_grouped_counts.barcodes.tsv` - grouped transcript counts in Seurat-compatible MTX format;
* `SAMPLE_ID.gene_grouped_tpm.matrix.mtx`, `SAMPLE_ID.gene_grouped_tpm.features.tsv`, `SAMPLE_ID.gene_grouped_tpm.barcodes.tsv` - grouped gene TPM values in Seurat-compatible MTX format;
* `SAMPLE_ID.transcript_grouped_tpm.matrix.mtx`, `SAMPLE_ID.transcript_grouped_tpm.features.tsv`, `SAMPLE_ID.transcript_grouped_tpm.barcodes.tsv` - grouped transcript TPM values in Seurat-compatible MTX format;
* `SAMPLE_ID.exon_grouped_counts.include.matrix.mtx`, `SAMPLE_ID.exon_grouped_counts.exclude.matrix.mtx`, `SAMPLE_ID.exon_grouped_counts.features.tsv`, `SAMPLE_ID.exon_grouped_counts.barcodes.tsv` - grouped exon counts in Seurat-compatible MTX format; one features file and one barcodes file are shared between the include and exclude matrices;
* `SAMPLE_ID.intron_grouped_counts.include.matrix.mtx`, `SAMPLE_ID.intron_grouped_counts.exclude.matrix.mtx`, `SAMPLE_ID.intron_grouped_counts.features.tsv`, `SAMPLE_ID.intron_grouped_counts.barcodes.tsv` - grouped intron counts in Seurat-compatible MTX format;


## Transcript discovery output

_Will not be produced if `--no_model_construction` is set._

File names typically contain `transcript_model` in their name.

* `SAMPLE_ID.transcript_models.gtf` - GTF file with discovered expressed transcript (both known and novel transcripts);
* `SAMPLE_ID.transcript_model_reads.tsv.gz` - TSV file indicating which reads contributed to transcript models (gzipped by default);
* `SAMPLE_ID.discovered_transcript_counts.tsv` - raw read counts for discovered transcript models (corresponds to `SAMPLE_ID.transcript_models.gtf`);
* `SAMPLE_ID.discovered_gene_counts.tsv` - raw read counts for discovered genes (corresponds to `SAMPLE_ID.transcript_models.gtf`);
* `SAMPLE_ID.discovered_transcript_tpm.tsv` - expression of discovered transcripts models in TPM (corresponds to `SAMPLE_ID.transcript_models.gtf`);
* `SAMPLE_ID.discovered_gene_tpm.tsv` - expression of discovered genes in TPM (corresponds to `SAMPLE_ID.transcript_models.gtf`);
* `SAMPLE_ID.extended_annotation.gtf` - GTF file with the entire reference annotation plus all discovered novel transcripts;


If `--read_group` is set, the per-group counts for discovered transcripts will be also computed:

* `SAMPLE_ID.discovered_transcript_grouped_counts.linear.tsv`
* `SAMPLE_ID.discovered_gene_grouped_counts.linear.tsv`

Similarly to the reference-based counts, these counts are converted to other formats as described [above](#other-formats).


If multiple experiments are provided, aggregated expression matrices will be placed in `<output_dir>`:

* `combined_gene_counts.tsv`
* `combined_gene_tpm.tsv`
* `combined_transcript_counts.tsv`
* `combined_transcript_tpm.tsv`

Additionally, an `isoquant.log` log file will be saved to the output directory.  

If raw reads were provided, BAM file(s) will be stored in `<output_dir>/<SAMPLE_ID>/aux/`.  
In case `--keep_tmp` option was specified this directory will also contain temporary files.

