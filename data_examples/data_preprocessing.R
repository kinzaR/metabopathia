# metabolomics data preprocessing from scaled and imputed metabolomics datasets of paired BRCA (65 tumor and 65 normal
# breast tissue) samples were downloaded from the supplementary materials of a publication [https://www.jci.org/articles/view/71180 ].
install.packages("readxl")
library(readxl)
excel_file <- "data_examples/JCI71180sd2(1).xlsx"
sheet_names <- excel_sheets(excel_file)
# we are interested by the scaled&Imputed data  second sheet.
metabo_data <- read_excel(excel_file, sheet = sheet_names[2])
# Design
metabo_des <- metabo_data[c(1:3),]
colnames(metabo_des) <- NULL
metabo_des <- metabo_des[, !colSums(is.na(metabo_des)) > 0]
metabo_des <- t(metabo_des)
table(metabo_des[,2])
dim(metabo_des)
write.table(metabo_des, file = "data_examples/metabo_brca_design.tsv", append = F,quote = F, sep = "\t", row.names = F, col.names = F)
# metabo available info & annotation 
metabo_info <- metabo_data[,c(1:9)]
colnames(metabo_info)<-metabo_info[4,]
metabo_info <- metabo_info[-c(1:4),]
write.table(metabo_info, file = "data_examples/metabo_brca_info.tsv", append = F,quote = F, sep = "\t", row.names = F, col.names = T)
# Dataset
metabo_data<-metabo_data[,-c(2:10)]
colnames(metabo_data)<- c("BIOCHEMICAL", metabo_data[1,-1])
metabo_data<-metabo_data[-c(1:4),]
# rownames(metabo_data)<- metabo_data$BIOCHEMICAL
write.table(metabo_data, file = "data_examples/metabo_brca_data.tsv", append = F,quote = F, sep = "\t", row.names = F, col.names = T)

