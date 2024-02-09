# NB: please setwd to the main WD metabopathia
# generating MGIs
source("src/utils.R")
source("src/metabopathia.R")
# vars
spe <- "hsa"
pathways_list <- NULL
v <- "v.1"
current_date <- format(Sys.Date(), "%Y%m%d")
### load pathways with pathway list:  preprocessed KEGG pathway
message("Loading pathways...") #I have to load from prepared MGI already ! time consuming 
# Load the pre-processed pathway object
pathways <- hipathia::load_pathways(species = spe, pathways_list = pathways_list)
## adaptation of the MGI: to be removed, because I have to load it already prepared
metabo_pathways <- add_metabolite_to_mgi(pathways)
## have to filter easy pathwys
saveRDS(object = metabo_pathways, file = paste0("pathways/mgi_metabopathia_",spe,"_",v,"_",current_date,".rds"))
