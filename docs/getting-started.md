# Quick start

*   IsoQuant can be downloaded from [https://github.com/ablab/IsoQuant](https://github.com/ablab/IsoQuant) or installed via conda:

        conda create -c conda-forge -c bioconda -n isoquant python=3.8 isoquant

*   If installing manually, you will need Python3 (3.8 or higher), [gffutils](https://pythonhosted.org/gffutils/installation.html), [pysam](https://pysam.readthedocs.io/en/latest/index.html), [pyfaidx](https://pypi.org/project/pyfaidx/) and some other common Python libraries to be installed. See `requirements.txt` for details. You will also need to have [minimap2](https://github.com/lh3/minimap2) and [samtools](http://www.htslib.org/download/) to be in your `$PATH` variable.

*   Verify your installation by running:

        isoquant.py --test

*   To run IsoQuant on raw FASTQ/FASTA files use the following command

        isoquant.py --reference /PATH/TO/reference_genome.fasta \
        --genedb /PATH/TO/gene_annotation.gtf \
        --fastq /PATH/TO/sample1.fastq.gz /PATH/TO/sample2.fastq.gz \
        --data_type (assembly|pacbio_ccs|nanopore) -o OUTPUT_FOLDER

    For example, using the toy data provided within this repository,

        ./isoquant.py --fastq /home/andreyp/ablab/IsoQuant/tests/simple_data/chr9.4M.ont.sim.fq.gz \
        --reference /home/andreyp/ablab/IsoQuant/tests/simple_data/chr9.4M.fa.gz \
        --genedb /home/andreyp/ablab/IsoQuant/tests/simple_data/chr9.4M.gtf.gz \
        --data_type nanopore --complete_genedb -p TEST_DATA --output isoquant_test 


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