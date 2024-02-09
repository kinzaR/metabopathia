# Please dont forget to setwd to the principal WD metabopathia folder

################## config
# Load pathways
pathways <- readRDS("pathways/mgi_metabopathia_hsa_1_20240208.rds")
species <-  "hsa"
# data file path
metData_fname <- ""
# genData_fname <- "data_examples/GSE253171/GSE253171_Processed_data.txt"
# genData_fname <- "data_examples/GSE242284/processed_gene_fpkm.txt"
genData_fname <- "data_examples/GSE215221/GSE215221_process_data.txt"

## Lib
library(data.table)
########## data coverage for integrative data  
# read data
metData <- fread(
  file = metData_fname,
  header = T, sep = "\t"
) %>% as.data.frame(.)
genData <- fread(
  file = genData_fname,
  header = T, sep = "\t"
) %>% as.data.frame(.)
### annotation to EntrezIds and Metabolite Kegg id
metabolite_names <-
# gene_names <- genData$GeneID  %>%  .[!is.na(.)]
# gene_names <- genData$`GeneCard Gene ID` %>%  .[!is.na(.)]
gene_names <- genData$affy_snp_id %>%  .[!is.na(.)]
  
new_ids <- gsub("\\.[0123456789]+$", "", gene_names)
xref <- hipathia:::load_xref(species)
gene_names_trans <- hipathia:::translate_ids(ids = new_ids, xref = xref)
sum(!gene_names_trans$is_na)
# gene_names_trans2 <- gene_names_trans[!gene_names_trans$is_na, drop = FALSE]

### validator 
# 1- check the number of mutual samples

# 2- check the coverge of genes in Signaling pathways (SP)
all_genes<-pathways$all.genes
sum(all_genes %in%  gene_names)/length(all_genes)
# 3- check the coverage of metabolite in SP
all_metabolites<-pathways$all.metabolite
