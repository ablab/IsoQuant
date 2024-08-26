# IsoQuant command line options


## Basic options
`--output` (or `-o`)
    Output folder, will be created automatically.

Note: if your output folder is located on a shared disk, use `--genedb_output` for storing
reference annotation database.

`--help` (or `-h`)
    Prints help message.

`--full_help`
    Prints all available options (including hidden ones).

`--test`
    Runs IsoQuant on the toy data set.   


## Input options
`--data_type` or `-d`
    Type of data to process, supported values are:  `pacbio_ccs` (same as `pacbio`), `nanopore` (same as `ont`)
and  `assembly` (same as `transcripts`). This option affects the algorithm parameters.

Note, that for novel mono-exonic transcripts are not reported for ONT data by default, use `--report_novel_unspliced true`.

`--reference` or `-r`
    Reference genome in FASTA format (can be gzipped), required even when BAM files are provided.

`--index`
    Reference genome index for the specified aligner (`minimap2` by default),
can be provided only when raw reads are used as an input (constructed automatically if not set).

`--genedb` or `-g`
    Gene database in gffutils database format or GTF/GFF format (can be gzipped).
If you use official gene annotations we recommend to set `--complete_genedb` option.

`--complete_genedb`
    Set this flag if gene annotation contains transcript and gene meta-features.
Use this flag when providing official annotations, e.g. GENCODE.
This option will set `disable_infer_transcripts` and `disable_infer_genes` gffutils options,
which dramatically speeds up gene database conversion (see more [here](https://daler.github.io/gffutils/autodocs/gffutils.create.create_db.html)).

### Providing input reads via command line option:

`--fastq`
    Input FASTQ/FASTA file(s), can be gzipped;  a single GTF will be generated for all files. If multiple files are provided,
expression tables with "per-file" columns will be computed. See more about [input data](input.md).


`--bam`
    Sorted and indexed BAM file(s); a single GTF will be generated for all files. If multiple files are provided,
expression tables with "per-file" columns will be computed. See more about [input data](input.md).


### Providing input reads via YAML configuration file:

`--yaml`
    Path to dataset description file in [YAML](https://www.redhat.com/en/topics/automation/what-is-yaml) format. The file should contain a list with `data format` property,
which can be `fastq` or `bam` and an individual entry for experiment.
Each experiment is represented as set of parameters (e.g. in curly brackets):
- `name` - experiment name, string (optional);
- `long read files` - a list of paths to long read files matching the specified format;
- `lables` - a list labels for long read files for expression table (optional, must be equal to the number of long read files)
- `illumina bam` - a list of paths to short read BAM files for splice site correction (optional).

All paths should be either absolute or relative to the YAML file.
See more in [examples](examples.md).

### Providing input reads via dataset description file (deprecated since 3.4)

`--bam_list` (_deprecated since 3.4_)
    Text file with list of BAM files, one file per line. Each file must be sorted and indexed.
Leave empty line or experiment name starting with # between the experiments.
For each experiment IsoQuant will generate a individual GTF and count tables.
You may also give a label for each file specifying it after a colon (e.g. `/PATH/TO/file.bam:replicate1`).

`--fastq_list` (_deprecated since 3.4_)
    Text file with list of FASTQ/FASTA files (can be gzipped),  one file per line.
Leave empty line or experiment name starting with # between the experiments.
For each experiment IsoQuant will generate a individual GTF and count tables.
You may also give a label for each file specifying it after a colon (e.g. `/PATH/TO/file.fastq:replicate1`).

### Other input options:
`--stranded`
    Reads strandness type, supported values are: `forward`, `reverse`, `none`.

`--fl_data`
    Input sequences represent full-length transcripts; both ends of the sequence are considered to be reliable.

`--prefix` or `-p`
    Prefix for all output files and sub-folder name. `OUT` if not set.

`--labels` or `-l`
    Sets space-separated sample names. Make sure that the number of labels is equal to the number of files.
Input file names are used as labels if not set.

`--read_group`
 Sets a way to group feature counts (e.g. by cell type). Available options are:

 * `file_name`: groups reads by their original file names (or file name labels) within an experiment.
This option makes sense when multiple files are provided.
This option is designed for obtaining expression tables with a separate column for each file.
If multiple BAM/FASTQ files are provided and `--read_group` option is not set, IsoQuant will set `--read_group:file_name`
by default.
 * `tag`: groups reads by BAM file read tag: set `tag:TAG`, where `TAG` is the desired tag name
(e.g. `tag:RG` with use `RG` values as groups, `RG` will be used if unset);
 * `read_id`: groups reads by read name suffix: set `read_id:DELIM` where `DELIM` is the
symbol/string by which the read id will be split
(e.g. if `DELIM` is `_`, for read `m54158_180727_042959_59310706_ccs_NEU` the group will set as `NEU`);
 * `file`: uses additional file with group information for every read: `file:FILE:READ_COL:GROUP_COL:DELIM`,
where `FILE` is the file name, `READ_COL` is column with read ids (0 if not set),
`GROUP_COL` is column with group ids (1 if not set),
`DELIM` is separator symbol (tab if not set). File can be gzipped.


### Output options

`--sqanti_output`
    Produce comparison between novel and known transcripts in SQANTI-like format.
    Will take effect only when reference annotation is provided.

`--check_canonical`
    Report whether read or constructed transcript model contains non-canonical splice junction (requires more time).

`--count_exons`
    Perform exon and intron counting in addition to gene and transcript counting.
    Will take effect only when reference annotation is provided.

`--bam_tags`
    Comma separated list of BAM tags that will be imported into `read_assignments.tsv`.

## Pipeline options

`--resume`
    Resume a previously unfinished run. Output folder with previous run must be specified.
    Allowed options are `--threads` and `--debug`, other options cannot be changed.
    IsoQuant will run from the beginning if the output folder does not contain the previous run.

`--force`
    force to overwrite the folder with previous run.

`--threads` or `-t`
    Number of threads to use, 16 by default.

`--clean_start`
    Do not use previously generated gene database, genome indices or BAM files, run pipeline from the very beginning (will take more time).

`--no_model_construction`
    Do not report transcript models, run read assignment and quantification of reference features only.

`--run_aligner_only`
    Align reads to the reference without running IsoQuant itself.


## Algorithm parameters

### Quantification

`--transcript_quantification` Transcript quantification strategy;
`--gene_quantification` Gene quantification strategy;

Available options for quantification:

* `unique_only` - use only reads that are uniquely assigned and consistent with a transcript/gene
(i.e. flagged as unique/unique_minor_difference), default fot transcript quantification;
* `with_ambiguous` - in addition to unique reads, ambiguously assigned consistent reads are split between features with equal weights 
(e.g. 1/2 when a read is assigned to 2 features simultaneously);
* `unique_splicing_consistent` - uses uniquely assigned reads that do not contradict annotated splice sites
(i.e. flagged as unique/unique_minor_difference or inconsistent_non_intronic), default for gene quantification;
* `unique_inconsistent` - uses uniquely assigned reads allowing any kind of inconsistency;
* `all` - all of the above.


### Read to isoform matching:

`--matching_strategy` A preset of parameters for read-to-isoform matching algorithm, should be one of:

* `exact` - delta = 0, all minor errors are treated as inconsistencies;  
* `precise` - delta = 4, only minor alignment errors are allowed, default for PacBio data;  
* `default` - delta = 6, alignment errors typical for Nanopore reads are allowed, short novel introns are treated as deletions;   
* `loose` - delta = 12, even more serious inconsistencies are ignored, ambiguity is resolved based on nucleotide similarity.

Matching strategy is chosen automatically based on specified data type.
However, the parameters will be overridden if the matching strategy is set manually.

### Read alignment correction:

`--splice_correction_strategy` A preset of parameters for read alignment correction algorithms, should be one of:

* `none` - no correction is applied;  
* `default_pacbio` - optimal settings for PacBio CCS reads;
* `default_ont` - optimal settings for ONT reads;
* `conservative_ont` - conservative settings for ONT reads, only incorrect splice junction and skipped exons are fixed;
* `assembly` - optimal settings for a transcriptome assembly;    
* `all` - correct all discovered minor inconsistencies, may result in overcorrection.

This option is chosen automatically based on specified data type, but will be overridden if set manually.

### Transcript model construction:
`--model_construction_strategy` A preset of parameters for transcript model construction algorithm, should be one of

* `reliable` - only the most abundant and reliable transcripts are reported, precise, but not sensitive;  
* `default_pacbio` - optimal settings for PacBio CCS reads;
* `sensitive_pacbio` - sensitive settings for PacBio CCS reads, more transcripts are reported possibly at a cost of precision;
* `fl_pacbio` - optimal settings for full-length PacBio CCS reads, will be used if `--data_type pacbio_ccs` and `--fl_data` options are set;
* `default_ont` - optimal settings for ONT reads, novel mono-exonic transcripts are not reported (use `--report_novel_unspliced true`);
* `sensitive_ont` - sensitive settings for ONT reads, more transcripts are reported possibly at a cost of precision (including novel mono-exonic isoforms);
* `assembly` - optimal settings for a transcriptome assembly: input sequences are considered to be reliable and each transcript to be represented only once, so abundance is not considered;    
* `all` - reports almost all novel transcripts, loses precision in favor to recall.

This option is chosen automatically based on specified data type, but will be overridden if set manually.


`--report_novel_unspliced` Report novel mono-exonic transcripts (set `true` or `false`).
The default value is `false` for Nanopore data and `true` for other data types.
The main explanation that some aligners report a lot of false unspliced alignments
for ONT reads.


`--report_canonical`
    Strategy for reporting novel transcripts based on canonical splice sites, should be one of:

* `auto` - automatic selection based on the data type and model construction strategy (default); 
* `only_canonical` - report novel transcripts, which contain only canonical splice sites;
* `only_stranded` - report novel transcripts, for which the strand can be unambiguously derived using splice sites and 
presence of a polyA tail, allowing some splice sites to be non-canonical;
* `all` -- report all transcript model regardless of their splice sites.


`--polya_requirement` Strategy for using polyA tails during transcript model construction, should be one of:

* `auto` - default behaviour: polyA tails are required if at least 70% of the reads have polyA tail; 
polyA tails are always required for 1/2-exon transcripts when using ONT data (this is caused by elevated number of false 1/2-exonic alignments reported by minimap2); 
* `never` - polyA tails are never required; use this option **at your own risk** as it may noticeably increase false discovery rate, especially for ONT data;
* `always` - reported transcripts are always required to have polyA support in the reads.

Note, that polyA tails are always required for reporting novel unspliced isoforms. 



## Hidden options

Options below are shown only with `--full_help` option.
We recommend _not_ to modify these options unless you are clearly aware of their effect.

`--no_gzip`
    Do not compress large output files.

`--no_gtf_check`
    Do not perform input GTF checks.

`--no_secondary`
    Ignore secondary alignments.

`--aligner`
    Force to use this alignment method, can be `starlong` or `minimap2`; `minimap2` is currently used as default. Make sure the specified aligner is in the `$PATH` variable.

`--no_junc_bed`
    Do not use gene annotation for read mapping.

`--junc_bed_file`
    Annotation in BED12 format produced by `paftools.js gff2bed` (can be found in `minimap2`), will be created automatically if not given.

`--delta`
    Delta for inexact splice junction comparison, chosen automatically based on data type (e.g. 4bp for PacBio, 6pb for ONT).

`--genedb_output`
    If your output folder is located on a shared storage (e.g. NFS share), use this option to set another path
    for storing the annotation database, because SQLite database cannot be created on a shared disks.
    The folder will be created automatically.

`--high_memory`
    Cache read alignments instead for making several passes over a BAM file, noticeably increases RAM usage, 
but may improve running time when disk I/O is relatively slow.

`--min_mapq`
    Filers out all alignments with MAPQ less than this value (will also filter all secondary alignments, as they typically have MAPQ = 0).

`--inconsistent_mapq_cutoff`
    Filers out inconsistent alignments with MAPQ less than this value (works when the reference annotation is provided, default is 5).

`--simple_alignments_mapq_cutoff`
    Filers out alignments with 1 or 2 exons and MAPQ less than this value (works only in annotation-free mode, default is 1).

`--normalization_method`
    Method for normalizing non-grouped counts into TPMs:

* `simple` - standard method, scale factor equals to 1 million divided by the counts sum (default);
* `usable_reads` - includes all reads assigned to a feature including the ones that were filtered out
during quantification (i.e. inconsistent or ambiguous);
scale factor equals to 1 million divided by the number of all assigned reads.
In this case the sum of all gene/transcript TPMs may not add up to 1 million.
Experiments with simulated data show that this method could give more accurate estimations.
However, normalization method does not affect correlation/relative proportions.

`--counts_format`
    Output format for grouped counts:
* `matrix` - usual format with genes as rows and groups as columns;
* `linear` - linear format, each line contains gene/transcript id, group name and respective count value (no TPM output); 
* `both` - output counts in both formats (default).