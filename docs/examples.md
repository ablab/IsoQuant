# Examples

* Mapped PacBio CCS reads in BAM format; pre-converted gene annotation:

```bash
isoquant.py -d pacbio_ccs --bam mapped_reads.bam \
 --genedb annotation.db --output output_dir
```

* Nanopore dRNA stranded reads; official annotation in GTF format, use custon prefix for output:
```bash
isoquant.py -d nanopore --stranded forward --fastq ONT.raw.fastq.gz \
 --reference reference.fasta --genedb annotation.gtf --complete_genedb \
 --output output_dir --prefix My_ONT
```

* Nanopore cDNA reads; no reference annotation:
```bash
isoquant.py -d nanopore --fastq ONT.cDNA.raw.fastq.gz \
 --reference reference.fasta --output output_dir --prefix My_ONT_cDNA
```

* PacBio FL reads; custom annotation in GTF format, which contains only exon features:
```bash
isoquant.py -d pacbio_ccs --fl_data --fastq CCS.fastq \
 --reference reference.fasta --genedb genes.gtf --output output_dir
```

* Nanopore cDNA reads, multiple samples/replicas within a single experiment; official annotation in GTF format:
```bash
isoquant.py -d nanopore --bam ONT.cDNA_1.bam ONT.cDNA_2.bam ONT.cDNA_3.bam \
 --reference reference.fasta --genedb annotation.gtf --complete_genedb --output output_dir
 --predix ONT_3samples --labels A1 A2 A3
```

* ONT cDNA reads; 2 experiments with 3 replicates; official annotation in GTF format:
```bash
isoquant.py -d nanopore --yaml dataset.yaml  \
 --complete_genedb --genedb genes.gtf \
 --reference reference.fasta --output output_dir
```

dataset.yaml file :

```
[
  data format: "fastq",
  {
    name: "Experiment1",
    long read files: [
      "/PATH/TO/SAMPLE1/file1.fastq",
      "/PATH/TO/SAMPLE1/file2.fastq",
      "/PATH/TO/SAMPLE1/file3.fastq"
    ],
    labels: [
      "Replicate1",
      "Replicate2",
      "Replicate3"
    ]
  },
  {
    name: "Experiment1",
    long read files: [
      "/PATH/TO/SAMPLE2/file1.fastq",
      "/PATH/TO/SAMPLE2/file2.fastq",
      "/PATH/TO/SAMPLE2/file3.fastq"
    ],
    labels: [
      "Replicate1",
      "Replicate2",
      "Replicate3"
    ]
  }
]

```


IsoQuant will produce 2 sets of resulting files (including annotations and expression tables), one for each experiment.
Output sub-folder will be named `Experiment1` and `Experiment2`.
Expression tables will have columns "Replicate1", "Replicate2" and "Replicate3".


* ONT cDNA reads; 1 experiment with 2 replicates, each replicate has 2 files; official annotation in GTF format:
```bash
isoquant.py -d nanopore --yaml dataset.yaml  \
  --complete_genedb --genedb genes.gtf \
 --reference reference.fasta --prefix MY_SAMPLE \
 --output output_dir  
```

`dataset.yaml` file :


```
[
  data format: "fastq",
  {
    name: "Experiment1",
    long read files: [
      "/PATH/TO/SAMPLE1/file1.fastq",
      "/PATH/TO/SAMPLE1/file2.fastq",
      "/PATH/TO/SAMPLE1/file3.fastq",
      "/PATH/TO/SAMPLE1/file3.fastq"
    ],
    labels: [
      "Replicate1",
      "Replicate1",
      "Replicate2",
      "Replicate2"
    ]
  }
]

```


IsoQuant will produce one output sub-folder `Experiment1`.
Expression tables will have columns "Replicate1" and "Replicate2".
Files having identical labels will be treated as a single replica (and thus the counts will be combined).
