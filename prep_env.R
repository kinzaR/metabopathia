# metabopathia
# step1: env setup
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hipathia")
# annotation for Metabolites

BiocManager::install("KEGGREST")
