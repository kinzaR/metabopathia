library(data.table)
metabo<-fread("data_examples/CCLE19/CCLE_metabolomics_20190502(1).csv")
View(metabo)
metID<- colnames(metabo) %>%.[-c(1,2)] %>% unique
write.table(metID, "data_examples/CCLE19/metabIDs.tsv", quote = F, append = F, row.names = F,col.names = F)
#####
metabo_pathways$all.metabolite
ids<-fread("data_examples/CCLE19/metaboIdsAnnotated.csv")
keggID<-unique(ids$KEGG)
length(intersect(metabo_pathways$all.metabolite,keggID))/length(unique(metabo_pathways$all.metabolite))
View(metabo)
expdata<-fread("data_examples/CCLE19/normalized_expression_CCLE.tsv")
intersect(colnames(expdata),metabo$CCLE_ID)
