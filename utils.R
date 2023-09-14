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
add_param_metaboID_subgraphs <- function(subgraphs){
  # subgraphs <-pathigraphs[[pathway]]$subgraphs
  for (sg in names(subgraphs)){
    subgraphs[[sg]] <-add_param_metaboID(subgraphs[[sg]])
  }
  return(subgraphs)
}

metabo_graphs <-function(pathigraphs){
  newpathigraphs <- list()
  for (pathway in names(pathigraphs)) {
    newpathigraphs[[pathway]]$graph <- add_param_metaboID(pathigraphs[[pathway]]$graph)
    newpathigraphs[[pathway]]$subgraphs <- add_param_metaboID_subgraphs(pathigraphs[[pathway]]$subgraphs)
    newpathigraphs[[pathway]]$subgraphs.mean.length <- pathigraphs[[pathway]]$subgraphs.mean.length
    newpathigraphs[[pathway]]$effector.subgraphs <- add_param_metaboID_subgraphs(pathigraphs[[pathway]]$effector.subgraphs)
    newpathigraphs[[pathway]]$path.name <- pathigraphs[[pathway]]$path.name
    newpathigraphs[[pathway]]$path.id <- pathigraphs[[pathway]]$path.id
    newpathigraphs[[pathway]]$label.id <- pathigraphs[[pathway]]$label.id
    newpathigraphs[[pathway]]$subgraphs_funs <- add_param_metaboID_subgraphs(pathigraphs[[pathway]]$subgraphs_funs)
    newpathigraphs[[pathway]]$effector.subgraphs_funs <- add_param_metaboID_subgraphs(pathigraphs[[pathway]]$effector.subgraphs_funs)
    newpathigraphs[[pathway]]$rl <- pathigraphs[[pathway]]$rl
    newpathigraphs[[pathway]]$fixed <- pathigraphs[[pathway]]$fixed
  }
  return(newpathigraphs)
}
add_metabolite_to_mgi <- function(mgi, verbose = FALSE){
  newmgi <- list()
  newmgi$species <- mgi$species
  newmgi$all.genes <- mgi$all.genes
  newmgi$all.metabolite <- all_needed_metabolites(mgi$pathigraphs)
  newmgi$path.norm  <- mgi$path.norm ## to change?
  newmgi$eff.norm  <- mgi$eff.norm ## to change?
  newmgi$pathigraphs <- metabo_graphs(mgi$pathigraphs)
  newmgi$all.labelids <- mgi$all.labelids
  newmgi$group.by <- mgi$group.by
  return(newmgi)
}