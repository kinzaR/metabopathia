library(dplyr) 
species <- "hsa"
load(paste0(species,"_module_data_Dec2016.RData"))
m <- hsa_module_data$M00001_C00022
metabolizer_metabolites<-lapply(hsa_module_data, function(m){
  m$KEGG_met_path_node$KEGG_met_path_node_info[m$KEGG_met_path_node$KEGG_met_path_node_info$type== "compound","graphics_name"]
}) %>% unlist %>% unname() %>% unique()

suero<-read.csv2("../../data_examples/Dystrophic_epidermolysis_bullosa/Suero_KEGGID-Conc_signif_final.csv")
intersect(suero$KEGG_ID, metabolizer_metabolites) %>% unique %>% length()
length(unique(suero$KEGG_ID))

ampolla<-read.csv2("../../data_examples/Dystrophic_epidermolysis_bullosa/Ampolla_KEGGID-Conc_signif_final.csv")
intersect(ampolla$KEGG, metabolizer_metabolites)  %>% unique %>% length()
length(unique(ampolla$KEGG))

orina<-read.csv("../../data_examples/Dystrophic_epidermolysis_bullosa/Orina.csv", sep = "\t") %>% t() %>% as.data.frame()
colnames(orina)<- orina[1,]
orina <- orina[-c(1,2),]
intersect(rownames(orina), metabolizer_metabolites)  %>% unique %>% length()
length(rownames(orina))
