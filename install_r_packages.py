import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector

# List of R packages to install
r_package_names = ('DESeq2', 'ggplot2', 'ggrepel', 'RColorBrewer', 'clusterProfiler', 'org.Hs.eg.db')

# Get R's utility package
utils = rpackages.importr('utils')

# Select CRAN mirror (optional, but recommended for reproducibility)
utils.chooseCRANmirror(ind=1) # Select the first mirror in the list

# Function to check if R package is installed
def is_installed(package_name):
    return package_name in rpackages.packages()

# Install R packages if not already installed
packages_to_install = [pkg for pkg in r_package_names if not is_installed(pkg)]

if packages_to_install:
    print(f"Installing R packages: {', '.join(packages_to_install)}")
    utils.install_packages(StrVector(packages_to_install))
else:
    print("All required R packages are already installed.") 