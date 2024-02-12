# Please dont forget to setwd to the principal WD metabopathia folder

################## config
# Load pathways
pathways <- readRDS("pathways/mgi_metabopathia_hsa_1_20240208.rds")
species <-  "hsa"
# data file path
# metData_fname <- "data_examples/GSE207088/non_polar.csv"
metData_fname <- "data_examples/GSE207088/polar_name_metaboanalyst.csv"
# genData_fname <- "data_examples/GSE253171/GSE253171_Processed_data.txt"
# genData_fname <- "data_examples/GSE242284/processed_gene_fpkm.txt"
# genData_fname <- "data_examples/GSE215221/GSE215221_process_data.txt"
# genData_fname <- "data_examples/GSE159857/GSE159857_REZ.July2017.coding_genes.qn.submission.csv"
genData_fname <- "data_examples/GSE207088/RNA_hNESCs_TPM.csv"

## Lib
library(data.table)
########## data coverage for integrative data  
# read data
metData <- fread(
  file = metData_fname,
  header = T, sep = "\t"
) %>% as.data.frame(.)
# metdata <- read.csv(metData_fname, header=T, row.names=1, sep=",")
genData <- fread(
  file = genData_fname,
  header = T, sep = "\t"
) %>% as.data.frame(.)
### annotation to EntrezIds and Metabolite Kegg id
metabolite_names <- metdata$KEGG %>%  .[!is.na(.)]

# gene_names <- genData$GeneID  %>%  .[!is.na(.)]
# gene_names <- genData$`GeneCard Gene ID` %>%  .[!is.na(.)]
# gene_names <- genData$affy_snp_id %>%  .[!is.na(.)]
# gene_names <- genData$GeneSymbol %>%  .[!is.na(.)]
gene_names <- genData$Gene %>%  .[!is.na(.)]
  
new_ids <- gsub("\\.[0123456789]+$", "", gene_names)
xref <- hipathia:::load_xref(species)
gene_names_trans <- hipathia:::translate_ids(ids = new_ids, xref = xref)
sum(!gene_names_trans$is_na)
gene_names <- gene_names_trans$translation
# gene_names_trans2 <- gene_names_trans[!gene_names_trans$is_na, drop = FALSE]

### validator 
# 1- check the number of mutual samples

# 2- check the coverge of genes in Signaling pathways (SP)
all_genes<-pathways$all.genes
sum(all_genes %in%  gene_names)/length(all_genes)
# 3- check the coverage of metabolite in SP
all_metabolites<-pathways$all.metabolite
sum(all_metabolites %in%  metabolite_names)/length(all_metabolites)
all_metabolites[all_metabolites %in%  metabolite_names]
