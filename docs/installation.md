# Installation

IsoQuant requires a 64-bit Linux system or Mac OS and Python (3.8 and higher) to be pre-installed on it.
You will also need
* [gffutils](https://pythonhosted.org/gffutils/installation.html)
* [pysam](https://pysam.readthedocs.io/en/latest/index.html)
* [biopython](https://biopython.org/)
* [pybedtools](https://daler.github.io/pybedtools/)
* [pyfaidx](https://pypi.org/project/pyfaidx/)
* [pandas](https://pandas.pydata.org/)
* [pyyaml](https://pypi.org/project/PyYAML/)
* [minimap2](https://github.com/lh3/minimap2)
* [samtools](http://www.htslib.org/download/)
* [STAR](https://github.com/alexdobin/STAR) (optional)

## Installing from conda
IsoQuant can be installed with conda:
```bash
conda install -c bioconda isoquant
```

## Manual installation and requirements
To obtain IsoQuant you can download repository and install requirements.  
Clone IsoQuant repository and switch to the latest release:
```bash
git clone https://github.com/ablab/IsoQuant.git
cd IsoQuant
git checkout latest
```
Install requirements:
```bash
pip install -r requirements.txt
```

You also need [samtools](http://www.htslib.org/download/) and [minimap2](https://github.com/lh3/minimap2) to be in the `$PATH` variable.

## Verifying your installation
To verify IsoQuant installation type
```bash
isoquant.py --test
```
to run on toy dataset.  
If the installation is successful, you will find the following information at the end of the log:
```bash
=== IsoQuant pipeline finished ===
=== TEST PASSED CORRECTLY ===
```