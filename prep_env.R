# metabopathia
# step1: env setup
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hipathia")
install.packages("igraph")
BiocManager::install("SummarizedExperiment")
install.packages("Matrix")
install.packages("preprocessCore")
# annotation for Metabolites

BiocManager::install("KEGGREST")
