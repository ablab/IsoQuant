#!/usr/bin/env python3
import rpy2
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector

print("rpy2 version is:", rpy2.__version__)

# List of R packages to install
r_package_names = ('DESeq2', 'ggplot2', 'ggrepel', 'RColorBrewer', 'BiocManager')

# Get R's utility package
utils = rpackages.importr('utils')

# Select CRAN mirror (optional, but recommended for reproducibility)
utils.chooseCRANmirror(ind=1) # Select the first mirror in the list

# Install R packages if not already installed
packages_to_install = [pkg for pkg in r_package_names if not rpackages.isinstalled(pkg)]

if packages_to_install:
    print(f"Installing R packages: {', '.join(packages_to_install)}")
    utils.install_packages(StrVector(packages_to_install))
else:
    print("All required R packages are already installed.")

biocondunctor = rpackages.importr('BiocManager')
print("Checking Bioconductor packages.")
for bio_pack in ['clusterProfiler', 'org.Hs.eg.db']:
    if not rpackages.isinstalled(bio_pack):
        biocondunctor.install()

print("Done")
