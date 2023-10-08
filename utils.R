library(igraph)
library(stringr)
library(hipathia)
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
add_metabolite_to_mgi <- function(mgi, verbose = FALSE, basal.value = 0.5){
  newmgi <- list()
  newmgi$species <- mgi$species
  newmgi$all.genes <- mgi$all.genes
  newmgi$all.metabolite <- all_needed_metabolites(mgi$pathigraphs)
  newmgi$pathigraphs <- metabo_graphs(mgi$pathigraphs)
  newmgi$all.labelids <- mgi$all.labelids
  
  genes.vals.05 <- matrix(basal.value, ncol = 2, nrow = length(newmgi$all.genes), 
                          dimnames = list(newmgi$all.genes, c("1", "2")))
  metabolites.vals.05 <- matrix(basal.value, ncol = 2, nrow = length(newmgi$all.metabolite), 
                          dimnames = list(newmgi$all.metabolite, c("1", "2")))
  meta.05 <- NULL
  meta.05$pathigraphs <- newmgi$pathigraphs
  meta.05$all.labelids <- newmgi$all.labelids
  source("metabopathia.R")
  results.05 <- metabopathia(genes.vals.05, metabolites.vals.05, meta.05, test = FALSE, 
                         verbose = FALSE)
  results.dec.05 <- metabopathia(genes.vals.05, metabolites.vals.05, meta.05, decompose = TRUE, 
                             test = FALSE, verbose = FALSE)
  
  newmgi$path.norm <- assay(results.dec.05, "paths")[, 1]
  newmgi$eff.norm <- assay(results.05, "paths")[, 1]
  
  #check here if there is any changes? Yes, specially in pathways that has metabolites
  # newmgi$path.norm  <- mgi$path.norm ## to change? -> yes
  # newmgi$eff.norm  <- mgi$eff.norm ## to change? -> yes
  newmgi$group.by <- mgi$group.by
  return(newmgi)
}

add_missing_metabolites <- function(metabo_vals, metabolites, default = NULL){
  if (is.null(default)) 
    default <- stats::median(metabo_vals) # here I have to add more sophisticated method
  missing_metabolites <- setdiff(metabolites, rownames(metabo_vals))
  if (length(missing_metabolites > 0)) {
    fakemat <- default + matrix(0, nrow = length(missing_metabolites), 
                                ncol = ncol(metabo_vals))
    rownames(fakemat) <- missing_metabolites
    colnames(fakemat) <- colnames(metabo_vals)
    metabo_vals <- rbind(metabo_vals, fakemat)
    message("Added missing metabolites: ", length(missing_metabolites), 
            " (", round(length(missing_metabolites)/nrow(metabo_vals) * 
                          100, digits = 2), "%)")
  }
  return(metabo_vals)
}
translate_metab_matrix <- function (metabo_vals, species, verbose = TRUE) 
{
  
  
  # Print the mapping
  print(kegg_ids)
  ###
  xref <- hipathia:::load_xref(species)
  
  new_ids <- gsub("^\\s+|\\s+", " ", trimws(rownames(metabo_vals)))
  tt <- hipathia:::translate_ids(new_ids, xref)
  exp2 <- exp[!tt$is_na, , drop = FALSE]
  valid_translation <- tt$translation[!tt$is_na]
  raw_exp3 <- by(exp2, valid_translation, colMeans, na.rm = TRUE)
  if (ncol(exp2) > 1) {
    exp3 <- do.call("rbind", raw_exp3)
  }
  else {
    exp3 <- matrix(raw_exp3, ncol = 1)
    rownames(exp3) <- names(raw_exp3)
    colnames(exp3) <- colnames(exp2)
  }
  if (verbose == TRUE) {
    cat("translated ids = ", tt$translated_ids_count, " (", 
        format(digits = 2, tt$translated_ids_ratio), ") \n", 
        sep = "")
    cat("untranslated ids = ", tt$untranslated_ids_count, 
        " (", format(digits = 2, tt$untranslated_ids_ratio), 
        ") \n", sep = "")
    cat("multihit ids = ", sum(tt$duplicated_ids_count), 
        " (", format(digits = 2, tt$duplicated_ids_ratio), 
        ") \n", sep = "")
  }
  attr(exp3, "translation") <- tt
  return(exp3)
}

test_metabo_pathways_object <- function (pathways) 
{
  hasall <- length(pathways) == 8 | length(pathways) == 9
  if (length(pathways) == 8) {
    byuser <- FALSE
  }
  else if (length(pathways) == 9) {
    byuser <- pathways$by.user
  }
  else if(is.null(pathways$by.user)){
    byuser <- FALSE
  }
  spec <- hipathia:::is_accepted_species(pathways$species) || byuser == 
    TRUE
  isigraph <- is(pathways$pathigraphs[[1]]$graph, "igraph")
  if (!hasall == TRUE | !spec == TRUE | !isigraph == TRUE) 
    stop("Pathways object not allowed")
}
