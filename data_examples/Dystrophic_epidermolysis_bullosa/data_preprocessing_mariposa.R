# data preprocessing : Dystrophic_epidermolysis_bullosa
metab_data <- read.csv2("data_examples/Dystrophic_epidermolysis_bullosa/Orina.csv")
metab_data <- t(metab_data)
colnames(metab_data) <- metab_data[1,]
metab_data <- metab_data[-c(1,2),]
# check how much metabolites are in metbMGI
metabo_pathways <- readRDS("metabo_pathways_hsa_120.RDS")
intersect(metabo_pathways$all.metabolite, rownames(metab_data))
intersect( rownames(metab_data), metabo_pathways$all.metabolite)
unique(metab_data[,48]) %in% rownames(metab_data)
setdiff(unique(metab_data[,48]), metabo_pathways$all.metabolite)
intersect(unique(metab_data[,48]), metabo_pathways$all.metabolite)
