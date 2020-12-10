![TeamCity Simple Build Status](http://chihua.cab.spbu.ru:3000/app/rest/builds/buildType:(id:IsoQuant_SimpleTest)/statusIcon)
![Python version](https://img.shields.io/badge/python-3.7-blue)
![License](https://img.shields.io/badge/licence-GPLv2-blue)
# IsoQuant 1.1 manual

1. [About IsoQuant](#sec1) </br>
    1.1. [Supported data types](#sec1.1)</br>
    1.2. [Supported reference data](#sec1.2)</br>
2. [Installation](#sec2)</br>
    2.1. [Installing from conda](#sec2.1)</br>
    2.2. [Manual installation and requirements](#sec2.2)</br>
    2.3. [Verifying your installation](#sec2.3)</br>
3. [Running IsoQuant](#sec3)</br>
    3.1. [IsoQuant input](#sec3.1)</br>
    3.2. [Command line options](#sec3.2)</br>
    3.3. [IsoQuant output](#sec3.3)</br>
4. [Citation](#sec4)</br>
5. [Feedback and bug reports](#sec5)</br>

**Quick start:**  

*   IsoQuant can be downloaded from [https://github.com/ablab/IsoQuant](https://github.com/ablab/IsoQuant) or installed via conda:

        conda install -c bioconda isoquant

    Or:

        conda install -c bioconda -c conda-forge -c isoquant isoquant

*   If installing manually, you will need Python3 (preferably 3.7), [gffutils](https://pythonhosted.org/gffutils/installation.html), [pysam](https://pysam.readthedocs.io/en/latest/index.html), [pybedtools](https://daler.github.io/pybedtools/), [biopython](https://biopython.org/) and some other common Python libraries to be installed. See `requirements.txt` for details. You will also need to have [minimap2](https://github.com/lh3/minimap2) and [samtools](http://www.htslib.org/download/) to be in your `$PATH` variable.
  
*   To run IsoQuant on raw FASTQ/FASTA files use the following command

        python IsoQuant.py 
        --reference /PATH/TO/reference_genome.fasta --genedb /PATH/TO/gene_annotation.gtf 
        --fastq /PATH/TO/sample1.fastq.gz /PATH/TO/sample2.fastq.gz 
        --data_type (pacbio_css|pacbio_raw|nanopore) -o OUTPUT_FOLDER


*   To run IsoQuant on aligned reads (make sure your BAM is sorted and indexed) use the following command:

        python IsoQuant.py --genedb /PATH/TO/gene_annotation.gtf 
        --fastq /PATH/TO/sample1.sorted.bam /PATH/TO/sample2.sorted.bam 
        --data_type (pacbio_css|pacbio_raw|nanopore) -o OUTPUT_FOLDER

<a name="sec1"></a>
# About IsoQuant

IsoQuant is a tool for reference-based analysis of long RNA reads, such as PacBio or Oxford Nanopores. IsoQuant maps reads to the reference genome and assigns them to the annotated isoforms based on their intron and exon structure. IsoQuant is also capable of discovering various modifications, such as intron retention, alternative splice sites, skipped exons etc. IsoQuant further performs gene, isoform, exon and intron quantification. If reads are grouped (e.g. according to cell type), counts are reported according to the provided grouping. In addition, IsoQuant generates discovered transcript models, including novel ones.

IsoQuant version 1.1 was released under GPLv2 on December 11th, 2020 and can be downloaded from [https://github.com/ablab/IsoQuant](https://github.com/ablab/IsoQuant).


<a name="sec1.1"></a>
## Supported data types

IsoQuant support all kinds of long RNA reads:
* PacBio CCS / HiFi
* PacBio IsoSeq
* ONT dRNA
* ONT cDNA

Reads must be provided in FASTQ or FASTA format (can be gzipped). If you have already aligned your reads to the reference genome, simply provide sorted and indexed BAM files.

<a name="sec1.2"></a>
## Supported reference data

Reference genome should be provided in multi-FASTA format. When BAM files are provided, reference genome is not mandatory, but can be used to count some additional statistics.

Gene annotation can be provided in GFF/GTF format. In this case it will be converted to [gffutils](https://pythonhosted.org/gffutils/installation.html) database. Information on converted databases will be stored in your `~/.config/IsoQuant/db_config.json` to increase speed of future runs. You can also provide gffutils database manually. Make sure that chromosome/scaffold names are identical in FASTA file and gene annotation.

Pre-constructed aligner index can also be provided to increase mapping time.

<a name="sec2"></a>
# Installation
IsoQuant requires a 64-bit Linux system or Mac OS and Python (3.7 and higher) to be pre-installed on it. 
You will also need 
* [gffutils](https://pythonhosted.org/gffutils/installation.html) 
* [pysam](https://pysam.readthedocs.io/en/latest/index.html) 
* [biopython](https://biopython.org/)
* [pybedtools](https://daler.github.io/pybedtools/)
* [pandas](https://pandas.pydata.org/)
* [numpy](https://numpy.org/)
* [minimap2](https://github.com/lh3/minimap2) 
* [samtools](http://www.htslib.org/download/) 
* [STAR](https://github.com/alexdobin/STAR) (optional)

<a name="sec2.1"></a>
## Installing from conda
Isoquant can be installed with conda:
```bash
conda install -c bioconda isoquant
```

If this command does not work, it means that bioconda is not updated yet. Try installing via:
```bash
conda create -n isoquant python=3.7
conda activate isoquant
conda install -c bioconda -c conda-forge -c isoquant isoquant
```

<a name="sec2.2"></a>
## Manual installation and requirements
To obtain IsoQuant you can download repository and install requirements.  
Clone IsoQuant repository and switch to latest release:
```bash
git clone https://github.com/ablab/IsoQuant.git
cd IsoQuant
git checkout latest
```
Install requirements:
```bash
pip install -r requirements.txt
```

You also need [minimap2](https://github.com/lh3/minimap2) to be in the `$PATH` variable.
  
<a name="sec2.3"></a>
## Verifying your installation 
To verify IsoQuant installation type
```bash
python3 isoquant.py --test
```
to run on toy dataset.  
If the installation is successful, you will find the following information at the end of the log:
```bash
=== IsoQuant pipeline finished === 
=== TEST PASSED CORRECTLY ===
```

<a name="sec3"></a>
# Running IsoQuant
<a name="sec3.1"></a>
## IsoQuant input
To run IsoQuant, you should provide:
* gene annotation in gffutils database or GTF/GFF format;
* reads in FASTA/FASTQ (can be gzipped) or sorted and indexed BAM;
* reference sequence in FASTA format (optional when BAMs are provided).  

By default, each file with reads is treated as a separate sample. To group multiple files into a single sample, provide a text files with paths to your FASTQ/FASTA/BAM files. Provide each file in a separate line, leave blank lines between samples.

<a name="sec3.2"></a>
## IsoQuant command line options

### Basic options
`--output` (or `-o`) 
    Output folder, will be created automatically.  

`--help` (or `-h`) 
    Prints help message.

`--full_help` 
    Prints all available options.

`--test` 
    Runs IsoQuant on the toy data set.   


### Input options
`--data_type` or `-d`
    Type of data to process, supported types are: `assembly`, `pacbio_ccs`, `pacbio_raw`, `nanopore`. This option affects some of the algorithm parameters.

`--genedb` or `-g`
    Gene database in gffutils database format or GTF/GFF format. If you use official gene annotations we recommend to set `--complete_genedb` option.

`--complete_genedb`
    Set this flag if gene annotation contains transcript and gene metafeatures. Use this flag when providing official annotations, e.g. GENCODE. This option will set `disable_infer_transcripts` and `disable_infer_genes` gffutils options, which dramatically speed up gene database conversion (see more [here](https://pythonhosted.org/gffutils/autodocs/gffutils.create_db.html?highlight=disable_infer_transcripts)).


`--reference` or `-r`
    Reference genome in FASTA format, should be provided  when raw reads are used as an input and to compute some additional stats.

#### Using mapped reads as input:
To provide aligned reads use one of the following options:

`--bam`
    Sorted and indexed BAM file(s); each file will be treated as a separate sample.

`--bam_list` 
    Text file with list of BAM files, one file per line, leave empty line between samples. 

#### Using raw read as an input:  
To provide read sequences use one of the following options:

`--fastq` 
    Input FASTQ/FASTA file(s); each file will be treated as a separate sample.
  
`--fastq_list` 
    Text file with list of FASTQ/FASTA files, one file per line, leave empty line between samples.  

#### Other input options:
`--stranded`
    Reads strandness type, supported values are: `forward`, `reverse`, `none`.

`--has_polya`
    Set this flag if reads were not poly-A trimmed; poly-A tails will be detected and further required for transcript model construction.

`--fl_data`
    Input sequences represent full-length transcripts; both ends of the sequence are considered to be reliable.

`--labels` or `-l`
    Sets space-separated sample names; make sure that the number of labels is equal to the number of samples; input file names are used if not set.

`--read_group`
 Sets a way to group feature counts (e.g. by cell type). Available options are: 
 * by BAM file read tag: set `tag:TAG`, where `TAG` is the desired tag name (e.g. `tag:RG` with use `RG` values as groups);
 * by read name suffix: set `read_id:DELIM` where `DELIM` is the symbol/string by which the read id will be split (e.g. if `DELIM` is `_`, for read `m54158_180727_042959_59310706_ccs_NEU` the group will set as `NEU`);
 * using additional file with group information for every read: `file:FILE:READ_COL:GROUP_COL:DELIM`, where `FILE` is the file name, `READ_COL` is column with read ids (0 if not set), `GROUP_COL` is column with group ids (1 if not set), `DELIM` is separator symbol (tab if not set).


### Pipeline and output options

`--clean_start`
    Do not use previously generated gene database, genome indices or BAM files, run pipeline from the very beginning (will take more time).

`--sqanti_output`
    Produce SQANTI-like TSV output (requires more time).

`--count_exons`
    Perform exon and intron counting in addition to gene and transcript counting.

`--use_secondary`
    Do not ignore secondary alignments (may significantly affect the results).

`--aligner`
    Force to use this alignment method, can be `starlong` or `minimap2`; chosen automatically based on data type if not set.

`--index`
    Reference genome index for the specified aligner, can be provided only when raw reads are used as an input (constructed automatically if not set).

`--run_aligner_only` 
    Align reads to reference without running IsoQuant itself.

`--threads` or `-t`
    Number of threads to use, 16 by default (currently affects the alignment stage only). 

`--keep_tmp` 
    Do not remove temporary files in the end.

### Algorithm parameters

#### Read to isoform matching:

`--matching-strategy` A preset of parameters for read to isoform matching algorithm, should be one of 
* `exact` - delta = 0, all minor errors are treated as inconsistencies;  
* `precise` - delta = 3, only minor alignment errors are allowed;  
* `default` - delta = 6, alignment errors typical for Nanopore reads are allowed, short novel introns are treated as deletions;   
* `loose` - delta = 12, even more serious inconsistencies are ignored, ambiguity is resolved based on nucleotide similarity.

Matching strategy is chosen automatically based on specified data type. However, the parameters will be overridden if the matching strategy is set manually.

You can manually set some of the parameters (will override options in the preset):

`--delta` 
    Delta for inexact splice junction comparison, chosen automatically based on data type.  

`--max-intron-shift` 
    Set maximum length for intron shift that will be treated as misalignment.

`--max-missed-exon-len` 
    Set maximum length for skipped exon that will be treated as misalignment.

#### Transcript model construction:

`--model_construction_strategy` A preset of parameters for transcript model construction algorithm, should be one of 
* `reliable` - only the most abundant and reliable transcripts are reported, precise, but not sensitive; intron retention is not reported;  
* `default` - a just trade-off between precision and recall for usual long-read dataset, intron retention is reported;   
* `all` - report most of detected modification as novel transcripts, loses precision in favor to recall; intron retention is reported;
* `fl` - input reads are considered as full-length transcripts; intron retention is reported;
* `assembly` - input sequences are considered to be reliable and each transcript to be represented only once, so abundance is not requires; intron retention is reported;

Transcript model construction strategy is chosen automatically based on specified data type. However, parameters will be overridden if set manually.

You can manually set some of the parameters (will override options in the preset):

`--report_intron_retention` 
    Report intron retention events as novel transcript models.

`--collapse_subisoform` 
    Collapse isoforms whose intron chain is a subsequence of another intron chain.

`--min_ref_fsm_supporting_reads` 
    Set a minimal number of full splice match reads that support known isoform.

`--min_ref_supporting_reads` 
    Set a minimal number of matching reads that support known isoform.

`--min_novel_fsm_supporting_reads` 
    Set a minimal number of full splice match reads that support novel isoform.

`--min_novel_supporting_reads` 
    Set a minimal number of reads that support a novel isoform.

`--min_reads_supporting_tsts` 
    Set a minimal number of reads that support isoform terminal sites.


### Examples

* Mapped raw PacBio reads in BAM format; not poly-A trimmed; pre-converted gene annotation:

```bash
isoquant.py -d pacbio_raw --has_polya --bam mapped_reads.bam --genedb annotation.db --output output_dir 
```

* Nanopore dRNA reads; not poly-A trimmed; official annotation in GTF format, used sample label instead of file name:
```bash
isoquant.py -d nanopore --has_polya --stranded forward --fastq ONT.raw.fastq.gz --reference reference.fasta --genedb annotation.gtf --complete_genedb --output output_dir --threads 8 --labels My_ONT
```

* PacBio FL reads, poly-A trimmed; custom annotation in GTF format, which contains only exon features:
```bash
python3 isoquant.py -d pacbio_ccs --fl_data --fastq CCS.fastq --reference reference.fasta --genedb genes.gtf --output output_dir --threads 8
```

<a name="sec3.3"></a>
## IsoQuant output

### Output files

IsoQuant output files will be stored in `<output_dir>`, which is set by the user. If the output directory was not specified the files are stored in `isoquant_output`.   
Output directory will contain one folder per sample with the following files:  

* `SAMPLE_ID.read_assignments.tsv` - TSV file with each read to isoform assignments;
* `SAMPLE_ID.transcript_counts.tsv` - TSV file with isoform counts;
* `SAMPLE_ID.gene_counts.tsv` - TSV file with gene counts;
* `SAMPLE_ID.transcript_models.gtf` - GTF file with constructed transcript models;
* `SAMPLE_ID.transcript_models_reads.tsv` - TSV file indicating which reads contributed to transcript models;
* `SAMPLE_ID.transcript_models_counts.tsv` - counts for constructed transcript models;
* `SAMPLE_ID.mapped_reads.bed` - coordinates of mapped reads in BED format.

If `--sqanti_output` is set, IsoQuant will save read assignments in [SQANTI](https://github.com/ConesaLab/SQANTI3#class)-like format:
* `SAMPLE_ID.SQANTI-like.tsv`

If `--count_exons` is set, exon and intron counts will be produced:
* `SAMPLE_ID.exon_counts.tsv`
* `SAMPLE_ID.intron_counts.tsv`

If `--read_group` is set, the per-group counts will be also computed:
* `SAMPLE_ID.exon_grouped_counts.tsv`
* `SAMPLE_ID.gene_grouped_counts.tsv`
* `SAMPLE_ID.intron_grouped_counts.tsv`
* `SAMPLE_ID.transcript_grouped_counts.tsv`

If multiple samples are provided, aggregated expression matrices will be placed in `<output_dir>`:
* `combined_gene_counts.tsv`
* `combined_transcript_counts.tsv`

Additionally, a log file will be saved to the directory.  
* <output_dir>/isoquant.log   

In case `--keep_tmp` option was specified output directory will also contain temporary files  
* <output_dir>/tmp/  

### Output file formats

Although most output files include headers that describe the data, a brief explanation of the output files is provided below.

#### Read to isoform assignment

Tab-separated values, the columns are:

* `read_id` - read id;
* `isoform_id` - isoform ids from the annotation, can contain several ids separated by comma;
* `assignment_type` - assignment type, can be:
    - `unique` - reads was unambiguously assigned to a single known isoform;
    - `unique_minor_difference` - read was assigned uniquely but has alignment artifacts;
    - `inconsistent` - read was matched with inconsistencies, closest match(es) are reported;
    - `ambiguous` - read was assigned to multiple isoforms equally well;
    - `noninfomative` - reads is intronic/intergenic.
* `assignment_events` - list of detected inconsistencies; for each assigned isoform a list of detected inconsistencies relative to the respective isoform is stored; values in each list are separated by `+` symbol, lists are separated by comma, the number of lists equals to the number of assigned isoforms; possible inconsistencies events are:
    - consistent events:
        - `none` / `.` / `undefined` - no special event detected;
        - `mono_exon_match` mono-exonic read matched to mono-exonic transcript;
        - `fsm` - full splice match;
        - `ism_5/3` - incomplete splice match, truncated on 5'/3' side;
        - `ism_internal` - incomplete splice match, truncated on both sides;
        - `mono_exonic` - mono-exonic read matching spliced isoform;
    - alignment artifacts:
        - `intron_shift` - intron that seems to be shifted due to misalignment (typical for Nanopores);
        - `exon_misallignment` - short exon that seems to be missed due to misalignment  (typical for Nanopores);
        - `fake_terminal_exon_5/3` - short terminal exon at 5'/3' end that looks like an alignment artifact (typical for Nanopores);  
        - `exon_elongation_5/3` - minor exon extension at 5'/3' end (not exceeding 30bp);
    - intron retentions:
        - `intron_retention` - intron retention;
        - `unspliced_intron_retention`  - intron retention by mono-exonic read;
        - `incomplete_intron_retention_5/3` - terminal exon at 5'/3' end partially covers adjacent intron;
    - significant inconsistencies:
        - `major_exon_elongation_5/3` - significant exon extension at 5'/3' end (exceeding 30bp);
        - `extra_intron_5/3` - additional intron on the 5'/3' end of the isoform;
        - `extra_intron_known` - read contains additional annotated intron in the middle of exon;
        - `extra_intron` - read contains additional unannotated intron in the middle of exon;
        - `alt_donor_site_known` - read contains alternative annotated donor site;
        - `alt_acceptor_site_known` - read contains alternative annotated acceptor site;
        - `alt_donor_site_novel` - read contains alternative unannotated donor site;
        - `alt_acceptor_site_novel` - read contains alternative unannotated acceptor site;
        - `intron_migration` - read contains alternative annotated intron of approximately the same length as the isoform;
        - `intron_alternation_known` - read contains alternative annotated intron;
        - `intron_alternation_novel` - read contains alternative unannotated intron; 
        - `mutually_exclusive_exons_known` - read contains different exon(s) comparing to the isoform, all resulting introns are annotated;
        - `mutually_exclusive_exons_novel` - read contains different exon(s) comparing to the isoform, at least one of introns is unannotated;
        - `exon_skipping_known_intron` - read skips exon(s) comparing to the isoform, resulting intron is annotated;
        - `exon_skipping_novel_intron` - read skips exon(s) comparing to the isoform, resulting intron is unannotated;
        - `exon_gain_known` - read contains additional exon(s) comparing to the isoform, all resulting introns are annotated;
        - `exon_gain_novel` - read contains additional exon(s) comparing to the isoform, at least one of introns is unannotated;
        - `alternative_structure_known` - reads has different intron chain that does not fall into any of categories above, all read's introns are annotated;
        - `alternative_structure_novel`  - reads has different intron chain that does not fall into any of categories above, at least one of read's introns is unannotated;
    - alternative transcription start / end (reported when CAGE data / poly-A tails are present):
        - `alternative_polya_site` - read has alternative polyadenylation site;
        - `fake_polya_site` - poly-A tail detected but seems to be originated from A-rich intronic region;
        - `alternative_tss` - alternative transcription start site.
* `polyA_found` - whether poly-A/T was detected at 5' end (True or False);
* `aligned_blocks` - list of coordinates for normalized aligned blocks (1-based, indels excluded);
* `intron_profile`, `split_exon_profile` - supplementary information intended for internal use; so called read profile, i.e. which of the gene's known introns and exons match the read.

#### Gene and transcript count format

Tab-separated values, the columns are:

* `feature_id` - genomic feature ID;
* `group_id` - read group if provided (NA by default);
* `count` - number of reads that were assigned to this feature (float value);

#### Exon and intron count format

Tab-separated values, the columns are:

* `chr` - chromosome ID;
* `start` - feature leftmost 1-based positions;
* `end` - feature rightmost 1-based positions;
* `strand` - feature strand;
* `flags` - symbolic feature flags, can contain the following characters:
    - `X` - terminal feature;
    - `I` - internal feature;
    - `T` - feature appears as both terminal and internal in different isoforms;
    - `S` - feature has similar positions to some other feature;
    - `C` - feature is contained in another feature;
    - `U` - unique feature, appears only in a single known isoform;
    - `M` - feature appears in multiple different genes.
* `gene_ids` - list if gene ids feature belong to;
* `group_id` - read group if provided (NA by default);
* `include_counts` - number of reads that include this feature;
* `exclude_counts` - number of reads that span, but do not include this feature; 

#### Transcript models format

Constructed transcript models are stored in usual [GTF format](https://www.ensembl.org/info/website/upload/gff.html). Currently, only exon features are listed. Metafeatures, such as genes and transcripts will be added later.

Transcript ids given in the `attribute` field have the following format: `transcript_###.TYPE`, where ### is the unique number (not necessarily consecutive) and TYPE can be one of the following:
* known - previously annotated transcripts;
* nic - novel in catalog, new transcript that contains only annotated introns;
* nnic - novel not in catalog, new transcript that contains unannotated introns.

The `attribute` field also contains `gene_id` (matches reference gene id), `reference_gene_id` (same value) and `reference_transcript_id` - the most similar known isoform. 

<a name="sec4"></a>
## Citation
Manuscript is in preparation.

<a name="sec5"></a>
## Feedback and bug reports
Your comments, bug reports, and suggestions are very welcome. They will help us to further improve IsoQuant. If you have any troubles running IsoQuant, please send us isoquant.log from the <output_dir> directory. 

You can leave your comments and bug reports at our [GitHub repository tracker](https://github.com/ablab/IsoQuant/issues) or send them via email: isoquant.rna@gmail.com.

