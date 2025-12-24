[![BioConda Install](https://img.shields.io/conda/dn/bioconda/isoquant.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/isoquant)
[![Python version](https://img.shields.io/badge/python-3.8-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/licence-GPLv2-blue)](https://www.gnu.org/licenses/old-licenses/gpl-2.0)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/ablab/IsoQuant)](https://github.com/ablab/IsoQuant/releases/)
[![GitHub Downloads](https://img.shields.io/github/downloads/ablab/IsoQuant/total.svg?style=social&logo=github&label=Download)](https://github.com/ablab/IsoQuant/releases)
[![UnitTests](https://github.com/ablab/IsoQuant/actions/workflows/Unit_tests.yml/badge.svg)](https://github.com/ablab/IsoQuant/actions/workflows/Unit_tests.yml)
[![User manual](https://github.com/ablab/IsoQuant/actions/workflows/docs.yml/badge.svg)](https://ablab.github.io/IsoQuant/)


# IsoQuant

[Full IsoQuant documentation can found here](https://ablab.github.io/IsoQuant/).

This is a feature branch developed for complementing IsoQuant with fusion transcript detection.
It is a separate pipeline triggered by `--fusion` flag in the IsoQuant input. Example command:
```
python3 isoquant.py -d nanopore --fusion --fastq ~/Thesis/real_data/SRR12048357_Ont.fastq \
 --genedb ~/Thesis/refs/gencode.v49.chr_patch_hapl_scaff.annotation.db \
 --reference ~/Thesis/refs/GRCh38.primary_assembly.genome.fa --output ~/Thesis/Isoquant/AML
```
