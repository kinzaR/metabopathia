## Selecting and Identifying MGI pathways without any complication. 
# Simple pathways are :
#   - Have only entrezID or NA in the genesList
# - If there is metabolites: only one per metaboliteID; no + & no hsa: 
library(hipathia)
species <- "hsa"
p <- hipathia::load_pathways(species = species)
pgs <- p$pathigraphs
i <-pgs$
for(i in pgs){
  # condition 1: check the genesList is it is only entrezId or / or NA
  print(i$path.name)
  
  # condition 2: check metabolite
}
