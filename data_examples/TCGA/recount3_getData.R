cancer_list <- c("BLCA","BRCA","COAD","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","THCA","UCEC")

# BiocManager::install("recount3")
library("recount3")
## Find all available human projects
human_projects <- available_projects()
## Find the project you are interested in,
proj_info <- subset(
  human_projects,
  file_source == "tcga" & project %in% cancer_list
)
## Create a RangedSummarizedExperiment (RSE) object at the gene level
rse_gene_tcga <- lapply(rownames(proj_info), function(p) create_rse(proj_info[p,], type = "gene"))
names(rse_gene_tcga) <- proj_info$project
# this was for debuging : memory problem !
# rse_gene_tcga <- list()
# for (p in rownames(proj_info)) {
#   print(proj_info[p,])
#   rse_gene_tcga[[proj_info[p,"project"]]] <-  create_rse(proj_info[p,], type = "gene")
# }
saveRDS(object=rse_gene_tcga, file= "TCGA/rse_gene_tcga.rds")