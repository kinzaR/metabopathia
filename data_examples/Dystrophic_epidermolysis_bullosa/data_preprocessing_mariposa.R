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
rownames(metab_data) %>% .[startsWith(., "C")] 
######################## datos de suero
metab_data_suero <- read.csv2("data_examples/Dystrophic_epidermolysis_bullosa/Suero_KEGGID-Conc_signif_final.csv") %>% 
  .[!is.na(.$KEGG_ID),]
rownames(metab_data_suero) <- metab_data_suero$KEGG_ID
metab_data_suero <- metab_data_suero[,-c(1,2)]
metabo_pathways <- readRDS("metabo_pathways_hsa_120.RDS")
intersect(metabo_pathways$all.metabolite,
          rownames(metab_data_suero))
intersect(rownames(metab_data_suero), metabo_pathways$all.metabolite
          )
allmetabo <- c()
library("stringr")
for(cp in metabo_pathways$pathigraphs){
  cat("\n *",cp$path.id, "\n")
  graph_att <-do.call(what = cbind, args = vertex_attr(cp$graph)) %>% as.data.frame(.)
  
  detected_metabolites <- graph_att[graph_att$shape=="circle","tooltip"] %>% 
    str_extract("(?<=href=http://www.kegg.jp/dbget-bin/www_bget\\?)[^>]+") %>% 
    .[!is.na(.)]
  detected_metabolites <- unique(detected_metabolites)
  detected_metabolites <- unlist(strsplit(x = detected_metabolites, split = "\\+"))
  if(any(rownames(metab_data_suero) %in% detected_metabolites)){
    cat("These metabolite: ",intersect(rownames(metab_data_suero) , detected_metabolites),"\n were found in ",cp$path.id)
  }
  # for the moment I will not see the complexity between gene/metabolite in the same node 
  # geneMetabo <- unlist(graph_att$genesList) %>% .[!is.na(.) & . != "/"] %>% .[is.na(as.numeric(.))]
  # cat(unique(geneMetabo))
  # if(length(geneMetabo)>0) detected_metabolites <- c(detected_metabolites, geneMetabo)
  allmetabo <- c(allmetabo, detected_metabolites)
}
allmetabo <- unique(allmetabo)
intersect(rownames(metab_data_suero), rownames(metab_data))
intersect(metab_data_ampolla$KEGG, rownames(metab_data))
##################### ampolla data
metab_data_ampolla <- read.csv2("data_examples/Dystrophic_epidermolysis_bullosa/Ampolla_KEGGID-Conc_signif_final.csv") %>% 
  .[!is.na(.$KEGG),]
# rownames(metab_data_ampolla) <- metab_data_ampolla$KEGG
metab_data_ampolla <- metab_data_ampolla[,-c(1)]
intersect(metab_data_ampolla$KEGG, metabo_pathways$all.metabolite)
intersect(intersect(rownames(metab_data_suero), rownames(metab_data)), metab_data_ampolla$KEGG)

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
  if(any(metab_data_ampolla$KEGG %in% detected_metabolites)){
    cat("These metabolite: ",intersect(metab_data_ampolla$KEGG , detected_metabolites),"\n were found in ",cp$path.id)
  }
  # for the moment I will not see the complexity between gene/metabolite in the same node 
  # geneMetabo <- unlist(graph_att$genesList) %>% .[!is.na(.) & . != "/"] %>% .[is.na(as.numeric(.))]
  # cat(unique(geneMetabo))
  # if(length(geneMetabo)>0) detected_metabolites <- c(detected_metabolites, geneMetabo)
  allmetabo <- c(allmetabo, detected_metabolites)
}
allmetabo <- unique(allmetabo)


exp_file <- "data_examples/Dystrophic_epidermolysis_bullosa/counts_TMM_normalization.txt"
met_file <-"data_examples/Dystrophic_epidermolysis_bullosa/Suero_KEGGID-Conc_signif_final.csv"
##read dta 
exp <- read.table(exp_file,header=T,sep="\t",stringsAsFactors=F,row.names = 1,comment.char="",check.names=F)
metabo_data<- read.csv2(met_file) %>% .[!is.na(.$KEGG_ID),-1]
rownames(metabo_data) <- metabo_data$KEGG_ID
metabo_data<- metabo_data[,-1] 
rn<-rownames(metabo_data)
metabo_data<-sapply(metabo_data, as.numeric)
rownames(metabo_data)<-rn
#fake_ to be removed
rownames(metabo_data)[1] <- "C00025"
# boxplot(metabo_data)
metabo_data <- metabo_data[,c("Suero_C01_08", "Suero_C02","Suero_C03","Suero_C04","Suero_C05","Suero_C06","Suero_C07" ,
                        "Suero_P1_V1", "Suero_P2_V1", "Suero_P3_V1","Suero_P4_V1","Suero_P5_V1","Suero_P6_V1","Suero_P7_V1")]
colnames(metabo_data) <- sapply(colnames(metabo_data), function(x){
  (unlist(strsplit(x, "_"))[2])
})%>% unname()
exp <- (exp[,c("M1.V1", "M2.V1", "M3.V1", "M4.V1", "M5.V1", "M6.V1", "M7.V1", "M8.V1", "C02", "C03", "C04","C06", "C07", "C08")])
colnames(exp)[1:8]<- paste0("P",seq_along(1:8))
colnames(exp)[colnames(exp)=="C08"] <- "C01"
to_keep <-intersect(colnames(metabo_data),colnames(exp))
to_keep<-to_keep[to_keep!="P5"]
design_data <- data.frame(sample=to_keep, group=c(rep("control",6 ),rep("visit1",6 )))
#out put
design_file <- "data_examples/Dystrophic_epidermolysis_bullosa/integration_design.tsv"
write.table(x = design_data, file = design_file, append = F, quote = F, sep = "\t", row.names = F, col.names = F)
write.table(x = exp[,to_keep], file ="data_examples/Dystrophic_epidermolysis_bullosa/counts_TMM_normalization.tsv", append = F,quote = F, sep = "\t", row.names = T,col.names = T)
write.table(x = metabo_data[,to_keep], file ="data_examples/Dystrophic_epidermolysis_bullosa/metabolite_suero.tsv", append = F,quote = F, sep = "\t", row.names = T,col.names = T)
