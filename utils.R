# ### to be removed 
# pathigraphs <- pathways$pathigraphs
# x <- pathways$pathigraphs$hsa04071
# mgi <- pathways
# ##############################################

library(igraph)
library(stringr)

all_needed_metabolites <- function (pathigraphs) 
{
  metabolite <- unique(unlist(sapply(pathigraphs, function(x) {
    tooltip<-V(x$graph)$tooltip[which(V(x$graph)$shape == 
                                      "circle")]
    # Extract the kegg id
    str_extract(tooltip, "(?<=href=http://www.kegg.jp/dbget-bin/www_bget\\?)[^>]+")
  })))
  return(metabolite[!is.na(metabolite)])
}
add_param_metaboID <- function(ig){
  # ig <-pathigraphs[[pathway]]$graph
  V(ig)$metaboID <- str_extract(V(ig)$tooltip, "(?<=href=http://www.kegg.jp/dbget-bin/www_bget\\?)[^>]+")
  return(ig)
}
metabo_graphs <-function(pathigraphs){
  newpathigraphs <- list()
  for (pathway in names(pathigraphs)) {
    newpathigraphs[[pathway]]$graph <- add_param_metaboID(pathigraphs[[pathway]]$graph)
    # no, because thse has no fnction nodes !!!!!!!!!!!!!!!!
    subs <- hipathia:::create_subgraphs(newpathigraphs[[pathway]]$graph)
    newpathigraphs[[pathway]]$subgraphs <- subs[[1]]
    newpathigraphs[[pathway]]$subgraphs.mean.length <- subs[[2]]
  }
}
add_metabolite_to_mgi <- function(mgi, verbose = FALSE){
  newmgi <- list()
  newmgi$species <- mgi$species
  newmgi$all.genes <- mgi$all.genes
  newmgi$all.metabolite <- all_needed_metabolites(mgi$pathigraphs)
  newmgi$path.norm  <- mgi$path.norm ## to change?
  newmgi$eff.norm  <- mgi$eff.norm ## to change?
  newmgi$pathigraphs <- metabo_graphs(mgi$pathigraphs)
  return(newmgi)
}