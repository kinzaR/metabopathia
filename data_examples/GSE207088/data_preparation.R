# annotation was performed using https://www.metaboanalyst.ca/MetaboAnalyst/home.xhtml
metabo_data <-read.table(file = "data_examples/GSE207088/polar_metab.tsv", sep = "\t", row.names = 1, header =  T)
metab_annotation <- read.csv("data_examples/GSE207088/polar_name_metaboanalyst.csv", header=T, row.names=1, sep=",")
metabo_data$ID<-metab_annotation[rownames(metabo_data),"KEGG"]
metabo_data<-metabo_data[!is.na(metabo_data$ID),]
rownames(metabo_data) <- metabo_data$ID
metabo_data$ID <- NULL
write.table(metabo_data, file = "data_examples/GSE207088/toBeUsed/polar_metab_normalized_annotated.tsv",sep = "\t", append = F, quote = F,row.names = T,col.names = T)
#######
# design
des <- data.frame(sample=colnames(metabo_data), group=c(rep("CTRL",3),rep("IPD",3)))
write.table(x = des, file = "data_examples/GSE207088/toBeUsed/design.tsv", sep = "\t", quote = F, row.names = F, col.names = F, append = F)
####
#RNA
rna<-read.table("data_examples/GSE207088/RNA_hNESCs_TPM.csv",header = T, sep = ",",row.names = 1)
colnames(rna) <- colnames(metabo_data)
write.table(rna, file = "data_examples/GSE207088/toBeUsed/RNA_hNESCs_TPM.tsv", quote = F, append = F,sep = "\t", row.names = T,col.names = T)
