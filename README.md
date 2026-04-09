
[Full IsoQuant documentation can found here](https://ablab.github.io/IsoQuant/).

This is a feature branch developed for complementing IsoQuant with fusion transcript detection.
It is a separate pipeline triggered by `--fusion` flag in the IsoQuant input. Example command:
```
python3 isoquant.py -d nanopore --fusion --fastq ~/Thesis/real_data/SRR12048357_Ont.fastq \
 --genedb ~/Thesis/refs/gencode.v49.chr_patch_hapl_scaff.annotation.db \
 --reference ~/Thesis/refs/GRCh38.primary_assembly.genome.fa --output ~/Thesis/Isoquant/AML
```
