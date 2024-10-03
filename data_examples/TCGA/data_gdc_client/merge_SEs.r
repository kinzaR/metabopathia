suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(SEtools))
source("utils.R")
cancer_list <- c("BLCA","BRCA","COAD","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","THCA","UCEC")
cancer_list_id <- paste("TCGA", cancer_list, sep = "-")
my_dir <- "processed_data/"
data_list <- lapply(cancer_list_id, load_data, my_dir)
merged_data <- mergeSEs(ll = data_list)
saveRDS(object = merged_data, file = file.path(my_dir,"counts_merged.rds"))