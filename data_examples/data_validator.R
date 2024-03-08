# Please dont forget to setwd to the principal WD metabopathia folder

################## config
# Load pathways
pathways <- readRDS("pathways/metabo_mgi_v1.0.0.RDS")
species <-  "hsa"
# data file path
metData_fname <- "data_examples/CCLE19/from_source/CCLE_metabolomics_20190502.csv"
# metData_fname <- "data_examples/GSE207088/non_polar.csv"
# metData_fname <- "data_examples/GSE207088/polar_name_metaboanalyst.csv"
# genData_fname <- "data_examples/GSE253171/GSE253171_Processed_data.txt"
# genData_fname <- "data_examples/GSE242284/processed_gene_fpkm.txt"
# genData_fname <- "data_examples/GSE215221/GSE215221_process_data.txt"
# genData_fname <- "data_examples/GSE159857/GSE159857_REZ.July2017.coding_genes.qn.submission.csv"
# genData_fname <- "data_examples/GSE207088/RNA_hNESCs_TPM.csv"
genData_fname <- "data_examples/CCLE19/from_source/CCLE_RNAseq_genes_counts_20180929.gct.gz"
## meta data
meta_data_fname <- "data_examples/CCLE19/from_source/Cell_lines_annotations_20181226.txt"
## Lib
library(dplyr)
library(data.table)
library(tidyr)

########## data coverage for integrative data  
# read data
metData <- fread(
  file = metData_fname,
  header = T
) %>% as_tibble() %>% column_to_rownames("CCLE_ID")
metData_info <- metData %>% select(1)
metData<-metData %>% select(-"DepMap_ID") %>% t()
# metdata <- read.csv(metData_fname, header=T, row.names=1, sep=",")
genData <- fread(
  file = genData_fname,
  header = T, sep = "\t"
) %>% as.data.frame(.)


### annotation to EntrezIds and Metabolite Kegg id
metabolite_names <- rownames(metData) %>% as_tibble() %>% drop_na() %>% set_colnames("Query")
write_delim(metabolite_names,"data_examples/CCLE19/metabolite_original_names.tsv")
# Then using https://www.metaboanalyst.ca/MetaboAnalyst/Secure/utils/NameMapView.xhtml I did the quick /manual annotatio
name_maped <- read_csv("data_examples/CCLE19/name_map.csv")
annotated <- metabolite_names %>% left_join(name_maped, "Query")
metData_annot <- metData %>% as_tibble() %>%
  mutate(KEGG = annotated$KEGG) %>% drop_na() %>% 
# annotator<- read_csv("pathways/all_metabolites_annotations.csKEGG# annotator<- read_csv("pathways/all_metabolites_annotations.csv")
# intersect(annotator$Match , metabolite_names)
# library(KEGGREST)
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
