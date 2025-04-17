[![BioConda Install](https://img.shields.io/conda/dn/bioconda/isoquant.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/isoquant)
[![Python version](https://img.shields.io/badge/python-3.8-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/licence-GPLv2-blue)](https://www.gnu.org/licenses/old-licenses/gpl-2.0)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/ablab/IsoQuant)](https://github.com/ablab/IsoQuant/releases/)
[![GitHub Downloads](https://img.shields.io/github/downloads/ablab/IsoQuant/total.svg?style=social&logo=github&label=Download)](https://github.com/ablab/IsoQuant/releases)
[![UnitTests](https://github.com/ablab/IsoQuant/actions/workflows/Unit_tests.yml/badge.svg)](https://github.com/ablab/IsoQuant/actions/workflows/Unit_tests.yml)
[![User manual](https://github.com/ablab/IsoQuant/actions/workflows/docs.yml/badge.svg)](https://ablab.github.io/IsoQuant/)


# IsoQuant 3.6

[Full IsoQuant documentation can found here](https://ablab.github.io/IsoQuant/).
Information in this README is given only for convenience and is not a full user manual.

* [Citation information](#citation)
* [Feedback and bug reports](#feedback-and-bug-reports)
* [Quick start examples](#quick-start)


## About IsoQuant

IsoQuant is a tool for the genome-based analysis of long RNA reads, such as PacBio or
Oxford Nanopores. IsoQuant allows to reconstruct and quantify transcript models with
high precision and decent recall. If the reference annotation is given, IsoQuant also
assigns reads to the annotated isoforms based on their intron and exon structure.
IsoQuant further performs annotated gene, isoform, exon and intron quantification.
If reads are grouped (e.g. according to cell type), counts are reported according to the provided grouping.

Latest IsoQuant version can be downloaded from [github.com/ablab/IsoQuant/releases/latest](https://github.com/ablab/IsoQuant/releases/latest).

Full IsoQuant documentation is available at [ablab.github.io/IsoQuant](https://ablab.github.io/IsoQuant/).

## Supported sequencing data

IsoQuant support all kinds of long RNA data:
* PacBio CCS
* ONT dRNA / ONT cDNA
* Assembled / corrected transcript sequences

Reads must be provided in FASTQ or FASTA format (can be gzipped). If you have already aligned your reads to the reference genome, simply provide sorted and indexed BAM files.
IsoQuant expect reads to contain polyA tails. For more reliable transcript model construction do not trim polyA tails.

IsoQuant can also take aligned Illumina reads to correct long-read spliced alignments. However, short reads are _not_
used to discover transcript models or compute abundances.


## Supported reference data

Reference genome is mandatory and should be provided in multi-FASTA format (can be gzipped).

Reference gene annotation is not mandatory, but is likely to increase precision and recall.
It can be provided in GFF/GTF format (can be gzipped).

Pre-constructed `minimap2` index can also be provided to increase mapping time.


## Citation
The paper describing IsoQuant algorithms and benchmarking is available at [10.1038/s41587-022-01565-y](https://doi.org/10.1038/s41587-022-01565-y).

To try IsoQuant you can use the data that was used in the publication [zenodo.org/record/7611877](https://zenodo.org/record/7611877).


## Feedback and bug reports
Your comments, bug reports, and suggestions are very welcome. They will help us to further improve IsoQuant. If you have any troubles running IsoQuant, please send us `isoquant.log` from the `<output_dir>` directory.

You can leave your comments and bug reports at our [GitHub repository tracker](https://github.com/ablab/IsoQuant/issues) or send them via email: isoquant.rna@gmail.com.



## Quick start

*   Full IsoQuant documentation is available at [ablab.github.io/IsoQuant](https://ablab.github.io/IsoQuant/).

*   IsoQuant can be downloaded from [github.com/ablab/IsoQuant](https://github.com/ablab/IsoQuant) or installed via conda:

        conda create -c conda-forge -c bioconda -n isoquant python=3.8 isoquant

*   If installing manually, you will need Python3 (3.8 or higher), [gffutils](https://pythonhosted.org/gffutils/installation.html), [pysam](https://pysam.readthedocs.io/en/latest/index.html), [pybedtools](https://daler.github.io/pybedtools/), [biopython](https://biopython.org/) and some other common Python libraries to be installed. See `requirements.txt` for details. You will also need to have [minimap2](https://github.com/lh3/minimap2) and [samtools](http://www.htslib.org/download/) to be in your `$PATH` variable.

*   Verify your installation by running:

        isoquant.py --test

*   To run IsoQuant on raw FASTQ/FASTA files use the following command

        isoquant.py --reference /PATH/TO/reference_genome.fasta \
        --genedb /PATH/TO/gene_annotation.gtf \
        --fastq /PATH/TO/sample1.fastq.gz /PATH/TO/sample2.fastq.gz \
        --data_type (assembly|pacbio_ccs|nanopore) -o OUTPUT_FOLDER

    For example, using the toy data provided within this repository,

        ./isoquant.py --reference tests/toy_data/MAPT.Mouse.reference.fasta \
        --genedb tests/toy_data/MAPT.Mouse.genedb.gtf \
        --fastq tests/toy_data/MAPT.Mouse.ONT.simulated.fastq \
        --data_type nanopore -o toy_data_out


* To run IsoQuant on aligned reads (make sure your BAM is sorted and indexed) use the following command:

        isoquant.py --reference /PATH/TO/reference_genome.fasta \
        --genedb /PATH/TO/gene_annotation.gtf \
        --bam /PATH/TO/sample1.sorted.bam /PATH/TO/sample2.sorted.bam \
        --data_type (assembly|pacbio_ccs|nanopore) -o OUTPUT_FOLDER

    For example, using the toy data provided within this repository,

        ./isoquant.py --reference tests/toy_data/MAPT.Mouse.reference.fasta \
        --genedb tests/toy_data/MAPT.Mouse.genedb.gtf \
        --fastq tests/toy_data/MAPT.Mouse.ONT.simulated.fastq \
        --data_type nanopore -o toy_data_out

* If using official annotations containing `gene` and `transcript` features use `--complete_genedb` to save time.

* Using reference annotation is optional since version 3.0, you may preform de novo transcript discovery without providing `--genedb` option':

        isoquant.py --reference /PATH/TO/reference_genome.fasta \
        --fastq /PATH/TO/sample1.fastq.gz /PATH/TO/sample2.fastq.gz \
        --data_type (assembly|pacbio|nanopore) -o OUTPUT_FOLDER

* If multiple files are provided, IsoQuant will create a single output annotation and a single set of gene/transcript expression tables.
