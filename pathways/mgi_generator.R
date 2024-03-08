codebase <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(codebase)
species <- "hsa"
# here I will load only easy pathways 
source("filter_path.R")
pathways_list <- get_easy_pathways(hipathia::load_pathways(species))
pathways <- hipathia::load_pathways(species = species, pathways_list = pathways_list)
source("utils.R")
metabo_pathways <- add_metabolite_to_mgi(pathways)
saveRDS(object = metabo_pathways, file = paste0("metabo_pathways_hsa_",length(pathways_list),"2.RDS"))
