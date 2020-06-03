# IsoQuant 0.1.0 manual

1. [About IsoQuant](#sec1) </br>
    1.1. [Supported data types](#sec1.1)</br>
    1.2. [IsoQuant pipeline](#sec1.2)</br>
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

**For impatient people:**  

*   IsoQuant can be downloaded from [https://github.com/ablab/IsoQuant](https://github.com/ablab/IsoQuant) or installed via conda:
        conda install -c bioconda isoquant (in progress)

*   You will need Python3, [gffutils](https://pythonhosted.org/gffutils/installation.html), [pysam](https://pysam.readthedocs.io/en/latest/index.html) and [biopython](https://biopython.org/) to be installed.
  
*   To run IsoQuant on raw FASTQ files use the following command

        python IsoQuant.py \
        --reference /PATH/TO/reference_genome.ta --gtf /PATH/TO/gene_annotation.gtf \
        --fastq /PATH/TO/sample1.fastq /PATH/TO/sample2.fastq \
        -o OUTPUT_FOLDER


*   To run IsoQuant on aligned reads use the following command

        python IsoQuant.py --gtf /PATH/TO/gene_annotation.gtf \
        --fastq /PATH/TO/sample1.sorted.bam /PATH/TO/sample2.sorted.bam \
        -o OUTPUT_FOLDER

<a name="sec1"></a>
# About IsoQuant

IsoQuant is a tool for reference-based analysis of long RNA reads, such as IsoSeq or Oxford Nanopore. IsoQuant maps reads to the reference genome and assigns them to annotated isoforms based on their intron and exon structure. IsoQuant is also capable of discovering various modifications, such as intron retention, alternative splice sites, skipped exons etc. Beside read-to-isofrom assignments, IsoQuant performs gene and isoform quantification, and computes various informative statics.

IsoQuant version 0.1.0 was released under GPLv2 on May 1st, 2020 and can be downloaded from [https://github.com/ablab/IsoQuant](https://github.com/ablab/IsoQuant).


<a name="sec1.1"></a>
## Supported data types

IsoQuant support all kinds of long RNA reads:
* IsoSeq CCS / HiFi
* Raw IsoSeq
* ONT dRNA
* ONT cDNA

Reads must be provided in FASTQ format (can be gzipped). If you have your data already aligned to  or in sorted and

<a name="sec1.2"></a>
## IsoQuant pipeline

<a name="sec2"></a>
# Installation
IsoQuant requires a 64-bit Linux system or Mac OS and Python (3.5 and higher) to be pre-installed on it. 
You will also need 
* [gffutils](https://pythonhosted.org/gffutils/installation.html) 
* [pysam](https://pysam.readthedocs.io/en/latest/index.html) 
* [biopython](https://biopython.org/)

<a name="sec2.1"></a>
## Installing from conda (in progress)
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
## Verifying your installation (in progress)
To verify IsoQuant installation type
```bash
python3 path/to/isoquant.py --test
```
to run on toy dataset.  
If the installation is successful, you will find the following information at the end of the log:
```bash
=== IsoQuant pipeline finished === 

========= TEST PASSED CORRECTLY.
```

<a name="sec3"></a>
# Running IsoQuant
<a name="sec3.1"></a>
## IsoQuant input
IsoQuant takes following files as an input:
* one or several alignments in indexed and sorted BAM file(s);
* file with DNA features in gffutils .db format or GTF/GFF format;   

Optionally, you can provide reads and reference sequence instead of BAM alignment: 
* reads in FASTA or FASTQ format;
* reference sequence in FASTA format.  

In this case they will be aligned and converted to BAM before the main pipeline.

<a name="sec3.2"></a>
## IsoQuant command line options
To run SPAdes from the command line, enter IsoQuant directory
```bash
cd path/to/IsoQuant
```
and type
```bash
python3 isoquant.py -d assembly --bam <alignment.bam> --genedb <genes.gtf.db> --output <output_dir> --threads 4 
```

### Basic options
`--output` (or `-o`) Output folder, will be created automatically.  
`--data_type` (or `-d`) Type of data to process, supported types are: `assembly`, `raw_long_reads`, `hq_long_reads`.  
`--threads` (or `-t`) Number of threads to use, 16 by default.   

`--help` (or `-h`) Prints help message.  
`--test` Runs IsoQuant on the toy data set (in progress).   


### Input options

`--genedb` (or `-g`) Gene database in gffutils .db format or GTF/GFF format.

#### Using alignment as input
`--bam` Sorted and indexed BAM file(s)  
or `--bam_list` Text file with list of BAM files, one file per line, leave empty line between samples. 




#### Using row read as an input:  
`--fastq` Input FASTQ file(s), each file will be treated as a separate sample  
or `--fastq_list` Text file with list of FASTQ files, one file per line, leave empty line between samples.  
`--reference` Reference genome in FASTA format.  

### Matching options
The easiest way to set matching parameters is to specify strategy.  
`--matching-strategy` Strategy that will be used, should be one of 
* `exact` - delta = 0, all minor errors and exon elongation ignored;  
* `precise` - delta = 3, only small intron shift and exon elongation allowed, cutoff = 30;  
* `default` - delta = 6, minor errors and exon elongation allowed, cutoff = 100;   
* `loose` - delta = 12, resolve ambiguous, allow extra introns/alternative TSS/polyA sites, minor errors and exon elongation allowed, cutoff = 300.  

Or you can manually set some of the parameters:  
`--delta` Delta for inexact splice junction comparison, chosen automatically based on data type.  
`--max-exon-extension` Set maximum length for exon elongation.  
`--max-intron-shift` Set maximum length for intron shift.  
`--max-missed-exon-len` Set maximum length for skipped exon.  


### Advanced options
`--labels` Sample names to be used.  
`--keep_tmp` Do not remove temporary files in the end.  
`--prefix` Prefix for output files.  
`--read_info` Text file with tab-separated information about input reads, according to which counts are groupped, e.g. cell type, barcode, etc.

#### Alignment options when reads in FASTA/FASTQ format are used
`--run_aligner_only` Align reads to reference without isoform assignment.  
`--aligner` Force to use this alignment method, can be `starlong` or `minimap2`, chosen based on data type if not set.  
`--index` Genome index for specified aligner.  
`--path_to_aligner` Folder with the aligner, $PATH is used by default.

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

You can leave your comments and bug reports at our GitHub repository tracker or send them via e-mail: isoquant.support@cab.spbu.ru.