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

# check if there is these metabolite in other no-simple pathways 
library(hipathia)
species<- "hsa"
metabo_pathways <- readRDS("metabo_pathways_hsa_120.RDS")
pathways <- load_pathways(species = species)
complex_pathways <- load_pathways(species = species, pathways_list = setdiff(names(pathways$pathigraphs),names(metabo_pathways$pathigraphs)))
cp <- complex_pathways$pathigraphs$hsa04024
cp <- complex_pathways$pathigraphs$hsa04742
allmetabo <- c()
library("stringr")
for(cp in complex_pathways$pathigraphs){
  cat("\n *",cp$path.id, "\n")
  graph_att <-do.call(what = cbind, args = vertex_attr(cp$graph)) %>% as.data.frame(.)
  
  detected_metabolites <- graph_att[graph_att$shape=="circle","tooltip"] %>% 
    str_extract("(?<=href=http://www.kegg.jp/dbget-bin/www_bget\\?)[^>]+") %>% 
    .[!is.na(.)]
  detected_metabolites <- unique(detected_metabolites)
  detected_metabolites <- unlist(strsplit(x = detected_metabolites, split = "\\+"))
  if(any(rownames(metab_data) %in% detected_metabolites)){
    cat("These metabolite: ",intersect(rownames(metab_data) , detected_metabolites),"\n were found in ",cp$path.id)
  }
  # for the moment I will not see the complexity between gene/metabolite in the same node 
  # geneMetabo <- unlist(graph_att$genesList) %>% .[!is.na(.) & . != "/"] %>% .[is.na(as.numeric(.))]
  # cat(unique(geneMetabo))
  # if(length(geneMetabo)>0) detected_metabolites <- c(detected_metabolites, geneMetabo)
  allmetabo <- c(allmetabo, detected_metabolites)
}
allmetabo <- unique(allmetabo)
intersect(allmetabo,
          rownames(metab_data))
intersect(rownames(metab_data),
          allmetabo)
View(metab_data[c("C00042", "C00186", "C00158"),])
