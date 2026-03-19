# Installation

IsoQuant requires a 64-bit Linux system or Mac OS and Python (3.8 and higher) to be pre-installed on it.
Note that IsoQuant 3.11.x and earlier versions are not well compatible with Python 3.14.
No special hardware is required. 

You will also need

-   [gffutils](https://pythonhosted.org/gffutils/installation.html)
-   [pysam](https://pysam.readthedocs.io/en/latest/index.html)
-   [pyfaidx](https://pypi.org/project/pyfaidx/)
-   [ssw-py](https://pypi.org/project/ssw-py/)
-   [editdistance](https://pypi.org/project/editdistance/)
-   [biopython](https://biopython.org/)
-   [numba](https://numba.pydata.org/)
-   [minimap2](https://github.com/lh3/minimap2)
-   [samtools](http://www.htslib.org/download/)
-   [STAR](https://github.com/alexdobin/STAR) (optional)

and some other common Python libraries. 
The full requirements list can be found in [requirements.txt](https://github.com/ablab/IsoQuant/blob/master/requirements.txt).

## Installing via pip

Simply run:

```bash
pip install isoquant
```

## Installing via conda
IsoQuant can be installed with conda using the bioconda channel:

```bash
conda create -c conda-forge -c bioconda -n isoquant python=3.12 isoquant
```

## Installing from GitHub

To get IsoQuant from GitHub, you can clone the repository:
```bash
git clone https://github.com/ablab/IsoQuant.git 
cd IsoQuant
git checkout latest
```

Or download the latest release [here](https://github.com/ablab/IsoQuant/releases).

IsoQuant can then be installed using pip:
```bash
pip install -e .
```

Alternatively, you can install the requirements manually:
```bash
pip install -r requirements.txt
```

Add [samtools](http://www.htslib.org/download/) and [minimap2](https://github.com/lh3/minimap2) to the `$PATH` variable and run IsoQuant without intalling it.

Typically, the whole installation takes a few minutes regardless of the method.

## Verifying your installation
To verify IsoQuant installation run
```bash
isoquant --test
```
to run on toy dataset. This should typically take less than 1 minute.  


If the installation is successful, you will find the following information at the end of the log:
```bash
=== IsoQuant pipeline finished ===
=== TEST PASSED CORRECTLY ===
```
