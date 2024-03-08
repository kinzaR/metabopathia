library(data.table)
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
#### re-normalize
expData_row<- fread("data_examples/GSE207088/GSE207088_feature_counts_rnaseq_hg19.txt", data.table = F, header = T)
# Set the first column as row names
rownames(expData_row) <- expData_row$Gene
expData_row$Gene<-NULL
library(edgeR)
#TMM normalization
exp_data_normalized<- DGEList(counts=as.matrix(expData_row))
exp_data_normalized<- calcNormFactors(exp_data_normalized,method = "TMM")
##get the TMM normalized matrix (counts*1000000/(lib.size*norm.factors)), 
exp_data_normalized<-cpm(exp_data_normalized,log = T,prior.count = 3)# we have in other pipline the pcount fixed to 3
write.table(exp_data_normalized, file = "data_examples/GSE207088/toBeUsed/RNA_hNESCs_TMM.tsv", quote = F, append = F,sep = "\t", row.names = T,col.names = T)
################
# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE207088", "file=GSE207088_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)
exp_data_normalized<- DGEList(counts=as.matrix(tbl))
exp_data_normalized<- calcNormFactors(exp_data_normalized,method = "TMM")
exp_data_normalized<-cpm(exp_data_normalized,log = T,prior.count = 3)
colnames(exp_data_normalized)
colnames(exp_data_normalized)<-c("IPD_1","IPD_2","IPD_3","CTRL_1","CTRL_2","CTRL_3")
# colnames(exp_data_normalized[,c("GSM6278356","GSM6278357","GSM6278358", "GSM6278359","GSM6278360","GSM6278361")])
write.table(exp_data_normalized, file = "data_examples/GSE207088/toBeUsed/RNA_hNESCs_TMM.tsv", quote = F, append = F,sep = "\t", row.names = T,col.names = T)
########33

# box-and-whisker plot
dat<-exp_data_normalized
par(mar=c(7,4,2,1))
boxplot(dat, boxwex=0.7, notch=T, main="GSE207088", ylab="tmm", outline=F, las=2)

# UMAP plot (dimensionality reduction)
require(umap)
dat <- dat[!duplicated(dat), ] # first remove duplicates
ump <- umap(t(dat), n_neighbors = 3, random_state = 123)
plot(ump$layout, main="GSE207088 UMAP plot, nbrs =3", xlab="", ylab="", pch=20, cex=1.5)
library(car)
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)



