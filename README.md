![TeamCity Simple Build Status](http://chihua.cab.spbu.ru:3000/app/rest/builds/buildType:(id:IsoQuant_SimpleTest)/statusIcon)
![Python version](https://img.shields.io/badge/python-3.7-blue)
![License](https://img.shields.io/badge/licence-GPLv2-blue)
# IsoQuant 1.0 manual

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

*   If installing manually, you will need Python3, [gffutils](https://pythonhosted.org/gffutils/installation.html), [pysam](https://pysam.readthedocs.io/en/latest/index.html), [pyfaidx](https://pypi.org/project/pyfaidx/), [biopython](https://biopython.org/) and some other common Python libraries to be installed. See `requirements.txt` for details.
  
*   To run IsoQuant on raw FASTQ files use the following command

        python IsoQuant.py \
        --reference /PATH/TO/reference_genome.ta --genedb /PATH/TO/gene_annotation.gtf \
        --fastq /PATH/TO/sample1.fastq /PATH/TO/sample2.fastq \
        --data_type (pacbio_css|pacbio_raw|nanopore) -o OUTPUT_FOLDER


*   To run IsoQuant on aligned reads (make sure your BAM is sorted and indexed) use the following command:

        python IsoQuant.py --genedb /PATH/TO/gene_annotation.gtf \
        --fastq /PATH/TO/sample1.sorted.bam /PATH/TO/sample2.sorted.bam \
        --data_type (pacbio_css|pacbio_raw|nanopore) -o OUTPUT_FOLDER

<a name="sec1"></a>
# About IsoQuant

IsoQuant is a tool for reference-based analysis of long RNA reads, such as PacBio or Oxford Nanopores. IsoQuant maps reads to the reference genome and assigns them to the annotated isoforms based on their intron and exon structure. IsoQuant is also capable of discovering various modifications, such as intron retention, alternative splice sites, skipped exons etc. IsoQuant further performs gene, isoform, exon and intron quantification. If reads are grouped (e.g. according to cell type), counts are reported according to the provided grouping. In addition, IsoQuant generates discovered transcript models, including novel ones.

IsoQuant version 1.0 was released under GPLv2 on July 11th, 2020 and can be downloaded from [https://github.com/ablab/IsoQuant](https://github.com/ablab/IsoQuant).


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

Pre-constructed aligner index can be also provided to increase mapping time.

<a name="sec2"></a>
# Installation
IsoQuant requires a 64-bit Linux system or Mac OS and Python (3.6 and higher) to be pre-installed on it. 
You will also need 
* [gffutils](https://pythonhosted.org/gffutils/installation.html) 
* [pysam](https://pysam.readthedocs.io/en/latest/index.html) 
* [biopython](https://biopython.org/)
* [pyfaidx](https://pypi.org/project/pyfaidx/)
* [pandas](https://pandas.pydata.org/)
* [numpy](https://numpy.org/)

<a name="sec2.1"></a>
## Installing from conda
Isoquant can be installed with conda:
```bash
conda install -c bioconda isoquant
```
<a name="sec2.2"></a>
## Manual installation and requirements
To obtain IsoQuant you can download repository and install requirements.  
Clone IsoQuant repository:
```bash
git clone https://github.com/ablab/IsoQuant.git
```
Enter IsoQuant directory and install requirements:
```bash
cd IsoQuant && pip install -r requirements.txt
```

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
To run IsoQuant from the command line, enter IsoQuant directory
```bash
cd path/to/IsoQuant
```
and type
```bash
python3 isoquant.py -d assembly --bam <alignment.bam> --genedb <genes.gtf.db> --output <output_dir> --threads 4 
```

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
    Gene database in gffutils database format or GTF/GFF format.

`--reference` or `-r`
    Reference genome in FASTA format, should be provided  when raw reads are used as an input and to compute some additional stats.

#### Using alignment as input
To provide aligned reads use one of the following options:

`--bam`
    Sorted and indexed BAM file(s); each file will be treated as a separate sample.

`--bam_list` 
    Text file with list of BAM files, one file per line, leave empty line between samples. 

#### Using row read as an input:  
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
    Sets sample names; make sure that the number of labels is equal to the number of samples; input file names are used if not set.

`--read_group`
     Sets a way to group feature counts (e.g. by cell type). Available options are: 
     - by BAM file read tag: set `tag:TAG`, where `TAG` is the desired tag name (e.g. `tag:RG` with use `RG` values as groups);
     - by read name suffix: set `read_id:DELIM` where `DELIM` is the symbol/string by which the read id will be split (e.g. if `DELIM` is `_`, for read `m54158_180727_042959_59310706_ccs_NEU` the group will set as `NEU`);
     - using additional file with group information for every read: `file:FILE:READ_COL:GROUP_COL:DELIM`, where `FILE` is the file name, `READ_COL` is column with read ids (0 if not set), `GROUP_COL` is column with group ids (1 if not set), `DELIM` is separator symbol (tab if not set).


### Pipeline and output options

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

`--path_to_aligner`
    Directory with the aligner binary, `$PATH` is used by default.

`--run_aligner_only` 
    Align reads to reference without running IsoQuant itself.

`--threads` or `-t`
    Number of threads to use, 16 by default (currently affects the alignment stage only). 

`--keep_tmp` 
    Do not remove temporary files in the end.

### Algorithm parameters

#### Read to isoform matching:

`--matching-strategy` A preset of parameters for read to isoform matching algorithm, should be one of 
    * `exact` - delta = 0, all minor errors and exon elongations are treated as inconsistencies;  
    * `precise` - delta = 3, only minor alignment errors and exon elongations are allowed, cutoff = 30;  
    * `default` - delta = 6, alignment errors and exon elongations are allowed, cutoff = 100;   
    * `loose` - delta = 12, resolve ambiguous matches based on nucleotide similarity, allow extra introns/alternative TSS/polyA sites, minor errors and exon elongation allowed, cutoff = 300.
    Matching strategy is chosen automatically based on specified data type. However, parameters will be overridden if matching strategy is set manually.

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
    * `all` - report most of detected modification as novel transcripts, looses precition in favor of recall; intron retention is reported;
    * `fl` - input reads are considered as full-length transcripts; intron retention is reported;
    * `assembly` - input sequences are considered to be reliable and each transcript to be represented only once, so abundance is not requires; intron retention is reported;
    Transcript model construction strategy is chosen automatically based on specified data type. However, parameters will be overridden if set manually.

You can manually set some of the parameters (will override options in the preset):
`--report_intron_retention` 
    Report intron retention events as novel transcript models.

`--collapse_subisoform` 
    Collapse isoforms whose intron chain is a subsequence of other intron chain.

`--min_ref_fsm_supporting_reads` 
    Set minimal number of full splice match reads that support known isoform.

`--min_ref_supporting_reads` 
    Set minimal number of matching reads that support known isoform.

`--min_novel_fsm_supporting_reads` 
    Set minimal number of full splice match reads that support novel isoform.

`--min_novel_supporting_reads` 
    Set minimal number of reads that support novel isoform.

`--min_reads_supporting_tsts` 
    Set minimal number of reads that support isoform terminal sites.


### Examples

In case you have 

* Aligment in BAM format

```bash
python3 isoquant.py -d assembly --bam alignment.bam --genedb genes.gtf.db --output output_dir --threads 4 
```

* Unaligned row reads
```bash
python3 isoquant.py -d raw_long_reads --fastq row_reads.fastq --reference reference.fasta --genedb genes.gtf.db --output output_dir --threads 4 
```

* High-quality reads
```bash
python3 isoquant.py -d hq_long_reads --fastq hq_reads.fastq --reference reference.fasta --genedb genes.gtf.db --output output_dir --threads 4 
```

<a name="sec3.3"></a>
## IsoQuant output
IsoQuant output files will be stored in in <output_dir>, which is set by the user. If output directory was not specified the files are stored in `isoquant_output` directory.   
Output directory will contain one folder per sample with three tsv files:  
* <output_dir>/sample_dir/sample.altered_reads.tsv  
* <output_dir>/sample_dir/sample.unmatched_reads.tsv  
* <output_dir>/sample_dir/sample.assigned_reads.tsv  

Additionally log file will be saved to the directory.  
* <output_dir>/isoquant.log   

In case `--keep_tmp` option was specified output directory will also contain temporary files  
* <output_dir>/tmp/  


<a name="sec4"></a>
## Citation
BioArxiv paper is in progress.

<a name="sec5"></a>
## Feedback and bug reports
Your comments, bug reports, and suggestions are very welcome. They will help us to further improve IsoQuant. If you have any troubles running IsoQuant, please send us isoquant.log from the <output_dir> directory. 

You can leave your comments and bug reports at our [GitHub repository tracker](https://github.com/ablab/IsoQuant/issues) or send them via e-mail: isoquant.rna@gmail.com.
