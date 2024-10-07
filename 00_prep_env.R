# metabopathia
# This code is written by Kinza Rian (rian.kinza@gmail.com)
# This code was developed under R version 4.3.1
# Step 1: Environment Setup

# Check if BiocManager is available; if not, install it
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install hipathia package from Bioconductor
if (!requireNamespace("hipathia", quietly = TRUE)) {
  BiocManager::install("hipathia")
}

# Install igraph package
if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}

# Install SummarizedExperiment package from Bioconductor
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  BiocManager::install("SummarizedExperiment")
}

# Install Matrix package
if (!requireNamespace("Matrix", quietly = TRUE)) {
  install.packages("Matrix")
}

# Install preprocessCore package
if (!requireNamespace("preprocessCore", quietly = TRUE)) {
  install.packages("preprocessCore")
}

# Install KEGGREST package from Bioconductor
if (!requireNamespace("KEGGREST", quietly = TRUE)) {
  BiocManager::install("KEGGREST")
}

# Install optparse package
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}

# Install R.utils package
if (!requireNamespace("R.utils", quietly = TRUE)) {
  install.packages("R.utils")
}

# Install gplots package
if (!requireNamespace("gplots", quietly = TRUE)) {
  install.packages("gplots")
}

# Check if all packages are installed successfully
if(all(c("BiocManager", "hipathia", "igraph", "SummarizedExperiment", "Matrix", "preprocessCore", "KEGGREST", "optparse","R.utils","gplots") %in% .packages())) {
  cat("All required dependencies and packages have been successfully installed.\n")
} else {
  stop("Installation failed for one or more packages. Please check the error messages and try installing the packages manually.")
}
# External bugs with compatibility bioc release (- 3.17) 
# install.packages("devtools")
# devtools::install_version("dbplyr", version = "2.3.4")