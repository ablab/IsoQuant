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


## Installing from conda
IsoQuant can be installed with conda:
```bash
conda create -c conda-forge -c bioconda -n isoquant python=3.12 isoquant
```
Typically, conda installation takes a few minutes.

## Installing from GitHub
To obtain IsoQuant you can download the repository and install the requirements.  
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
Typically, package installation takes a few minutes.

You also need [samtools](http://www.htslib.org/download/) and [minimap2](https://github.com/lh3/minimap2) to be in the `$PATH` variable.

## Verifying your installation
To verify IsoQuant installation run
```bash
isoquant.py --test
```
to run on toy dataset. This should typically take less than 1 minute.  


If the installation is successful, you will find the following information at the end of the log:
```bash
=== IsoQuant pipeline finished ===
=== TEST PASSED CORRECTLY ===
```
