# these has to be received args!
# # This are my cancers 
# cancer_list <- c("BLCA","BRCA","COAD","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","THCA","UCEC")
# cancer_list_id <- paste("TCGA", cancer_list, sep = "-")
###extract info from arguments
#### INPUT DATA
args <- commandArgs()
cancer <- args[grep("-cancer",args)+1] # cancer <- "BLCA"
cancer_id <- paste("TCGA", cancer, sep = "-")
my_cache_dir <- args[grep("-output_folder",args)+1] #"/media/kinza/DATA/gdcdata_cache_dir"
# Sys.getenv("HOME")
disk.usage <- function(path = my_cache_dir) {
  if(length(system("which df", intern = TRUE, ignore.stderr =  TRUE))) {
    cmd <- sprintf("df -h %s", path)
    exec <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
    exec <- strsplit(exec[length(exec)], "[ ]+")[[1]]
    exec <- as.character(exec[3:4])
    structure(exec, names = c("used", "available"))
  } else {
    stop("'df' command not found")
  }
}

du <- disk.usage()
message("Disk usage :\n\t")
print(du)
## Load Libraries
if (!require("BiocManager"))
  install.packages("BiocManager")
# if (!require("GenomicDataCommons"))
#   BiocManager::install('GenomicDataCommons')
if (!require("TCGAbiolinks"))
  BiocManager::install("TCGAbiolinks")

library(dplyr)
library(SummarizedExperiment)
library(TCGAbiolinks)
# Step 01: get data from TCGA
query <- GDCquery(project = cancer_id,
                 data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification",
                 workflow.type = "STAR - Counts",
                 access = "open")
GDCdownload(query, directory = my_cache_dir)
exp_data <- GDCprepare(query,
                       directory = my_cache_dir,
                       save = T,
                       save.filename = file.path(my_cache_dir,cancer_id, paste0("counts_",cancer_id,".RData")))#cancer_id folderis created priviosly while downloading data by GDCdownload() function.
# Step 02: Here check if there is a batch effect in the data 
# To be disscused : removing low-count for differential expression and normalize data is not relevant for hipathia, alcontrario is wrong !
library(edgeR)
# SEE utils.R
# Draft of using GenomicDataCommons package 
# # Check connectivity and status
# # GenomicDataCommons::status()
# # stopifnot(GenomicDataCommons::status()$status=="OK")
# 
# # Get data using filters
# file_cases <- files() |> filter(~ cases.project.project_id==cancer_id &
#                             type == 'gene_expression' &
#                             data_category == "Transcriptome Profiling" &
#                             analysis.workflow_type == 'STAR - Counts' &
#                             access == 'open') |>
#   GenomicDataCommons::select(c('file_name','cases.case_id')) |>
#   GenomicDataCommons::expand('cases.tissue_source_site')|>
#   response_all()
# # Find data
# ge_manifest <- file_cases %>% manifest()
# file_cases$results$cases_id<-unlist(file_cases$results$cases)
# # Here I checked if is the same AS WEB MANEFIST
# # from_web_manifest <- read.table("data_gdc_client/gdc_manifest.2024-06-20.txt", sep = "\t", header = T)
# # all(ge_manifest$id == from_web_manifest$id) # TRUE
# # all(ge_manifest$md5sum == from_web_manifest$md5) # TRUE
# # For my local: changing cache dir
# # original_cache_dir <- rappdirs::app_dir(appname = "GenomicDataCommons")$cache()
# gdc_set_cache(
#   directory = paste0(my_cache_dir,"_tmp"),
#   verbose = TRUE,
#   create_without_asking = F
# )
# #check if cache dir changed 
# message("Cache dir, where data is saved, is : ",gdc_cache()) # "/media/kinza/DATA/gdcdata_cache_dir" here I have more space!
# # Download Data
# fnames <- lapply(ge_manifest$id, gdcdata)
# 
# # Metadata queries
# ## Clinical data
# case_ids <- file_cases$results$cases
# clindat <- gdc_clinical(case_ids)
# names(clindat)## [1] "demographic" "diagnoses"   "exposures"   "main"
# head(clindat[["main"]])
# head(clindat[["diagnoses"]])
# 
# # other filter 
# resp = cases() |> filter(~ project.project_id==cancer_id &
#                            samples.sample_type=='Solid Tissue Normal') |>
#   GenomicDataCommons::select(c(default_fields(cases()),'samples.sample_type')) |>
#   response_all()
# count(resp)
# 
# ## General metadata queries
# expands = c("demographic","diagnoses","exposures",
#             "annotations")
# clinResults = cases() %>%
#   GenomicDataCommons::filter(project.project_id == cancer_id)  %>%
#   GenomicDataCommons::expand(expands) %>%
#   results_all()
# str(clinResults[[1]],list.len=6)
# exp <- read.table("/media/kinza/DATA/gdcdata_cache_dir_tmp/48c72fa3-43b5-4a3f-a712-afce30a963b4/867e266d-0161-4011-820b-e7ae874d7887.rna_seq.augmented_star_gene_counts.tsv", sep = "\t")
# ##  chr [1:50] "69eced5b-1e76-45c9-bc9c-2aa71a921c57" ...
# # or listviewer::jsonedit(clinResults)
