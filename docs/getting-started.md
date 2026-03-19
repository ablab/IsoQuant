# Quick start

*   IsoQuant can installed via pip:

        pip install isoquant

*   Via conda (bioconda channel):

        conda create -c conda-forge -c bioconda -n isoquant python=3.12 isoquant

*   Or from GitHub:

        git clone https://github.com/ablab/IsoQuant.git 
        cd IsoQuant
        git checkout latest
        pip install -e .

Installation typically takes no more than a few minutes.

*   If running simply from [the source archive](https://github.com/ablab/IsoQuant/releases/), 
you will need Python3 (3.8 or higher), [gffutils](https://pythonhosted.org/gffutils/installation.html), [pysam](https://pysam.readthedocs.io/en/latest/index.html), [biopython](https://biopython.org/), [pyfaidx](https://pypi.org/project/pyfaidx/),
 [ssw-py](https://pypi.org/project/ssw-py/), [editdistance](https://pypi.org/project/editdistance/) and some other common Python libraries to be installed. See `requirements.txt` for details. 
You will also need to have [minimap2](https://github.com/lh3/minimap2) and [samtools](http://www.htslib.org/download/) to be in your `$PATH` variable.
All required Python libraries can be installed via: 

        pip install -r requirements.txt

*   Verify your installation by running (typically takes less than 1 minute):

        isoquant --test

*   To run IsoQuant on raw FASTQ/FASTA files, use the following command

        isoquant --reference /PATH/TO/reference_genome.fasta \
        --genedb /PATH/TO/gene_annotation.gtf \
        --fastq /PATH/TO/sample1.fastq.gz /PATH/TO/sample2.fastq.gz \
        --data_type (assembly|pacbio_ccs|nanopore) -o OUTPUT_FOLDER

    For example, using the toy data provided within this repository,

        ./isoquant --fastq /home/andreyp/ablab/IsoQuant/isoquant_tests/simple_data/chr9.4M.ont.sim.fq.gz \
        --reference /home/andreyp/ablab/IsoQuant/isoquant_tests/simple_data/chr9.4M.fa.gz \
        --genedb /home/andreyp/ablab/IsoQuant/isoquant_tests/simple_data/chr9.4M.gtf.gz \
        --data_type nanopore --complete_genedb -p TEST_DATA --output isoquant_test 


* To run IsoQuant on aligned reads (make sure your BAM is sorted and indexed) use the following command:

        isoquant --reference /PATH/TO/reference_genome.fasta \
        --genedb /PATH/TO/gene_annotation.gtf \
        --bam /PATH/TO/sample1.sorted.bam /PATH/TO/sample2.sorted.bam \
        --data_type (assembly|pacbio_ccs|nanopore) -o OUTPUT_FOLDER

* If using official annotations containing `gene` and `transcript` features use `--complete_genedb` to save time.

* Using reference annotation is optional since version 3.0, you may preform de novo transcript discovery without providing `--genedb` option':

        isoquant --reference /PATH/TO/reference_genome.fasta \
        --fastq /PATH/TO/sample1.fastq.gz /PATH/TO/sample2.fastq.gz \
        --data_type (assembly|pacbio|nanopore) -o OUTPUT_FOLDER

* If multiple files are provided, IsoQuant will create a single output annotation and a single set of gene/transcript expression tables.