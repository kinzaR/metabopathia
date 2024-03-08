metData_from_metabolizer <- read.table("data_examples/metabolizer_as_DS/module_vals.txt", header = T, row.names = 1)

metData_from_metabolizer$new_rownames<-gsub(".*_","",rownames(metData_from_metabolizer))

metData_from_metabolizer_new<-metData_from_metabolizer %>% distinct(new_rownames, .keep_all = T)
rownames(metData_from_metabolizer_new)<- metData_from_metabolizer_new$new_rownames
metData_from_metabolizer_new$new_rownames<-NULL
write.table(metData_from_metabolizer_new, "data_examples/metabolizer_as_DS/inferedmetabolic_data.tsv", quote = F, append = F, sep = "\t")

intersect(metData_from_metabolizer_new$new_rownames , pathways$all.metabolite)
