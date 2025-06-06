suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(stringr))

all_needed_metabolites <- function (pathigraphs)
{
  metabolite <- unique(unlist(sapply(pathigraphs, function(x) {
    tooltip<-V(x$graph)$tooltip[which(V(x$graph)$shape == 
                                      "circle")]
    # Extract the kegg id
    met_ids<- str_extract(tooltip, "(?<=href=https?://www.kegg.jp/dbget-bin/www_bget\\?)[^>]+")
    if(length(met_ids)>0){
      met_ids<- met_ids %>% strsplit("\\+") %>% unlist %>% 
        .[!grepl("^hsa:", .)] %>% 
        .[!grepl("^D[0-9]{5}$", .)]
    }
    
  })))
  return(metabolite[!is.na(metabolite)])
}
add_param_metaboID <- function(ig){
  # ig <-pathigraphs[[pathway]]$graph
 list_metabs <- str_extract(V(ig)$tooltip, "(?<=href=https?://www.kegg.jp/dbget-bin/www_bget\\?)[^>]+") %>%
    strsplit("\\+") 
 V(ig)$metaboID  <- sapply(list_metabs, function(m_l){
   m_l %>% gsub("^hsa:", "", .) %>%  .[!grepl("^D[0-9]{5}$", .)] %>% paste(collapse = ",/,")
 })
 V(ig)$metaboID[V(ig)$metaboID == "NA" | V(ig)$metaboID == ""] <- NA
  return(ig)
}
add_param_metaboID_subgraphs <- function(subgraphs){
  subgraphs<-lapply(subgraphs, add_param_metaboID)
  return(subgraphs)
}

metabo_graphs <-function(pathigraphs){
  return(lapply(pathigraphs, function(pathigraph){
    cat("|")
    newpathigraph <- list()
    newpathigraph$graph <- add_param_metaboID(pathigraph$graph)
    newpathigraph$subgraphs <- add_param_metaboID_subgraphs(pathigraph$subgraphs)
    newpathigraph$subgraphs.mean.length <- pathigraph$subgraphs.mean.length
    newpathigraph$effector.subgraphs <- add_param_metaboID_subgraphs(pathigraph$effector.subgraphs)
    newpathigraph$path.name <- pathigraph$path.name
    newpathigraph$path.id <- pathigraph$path.id
    newpathigraph$label.id <- pathigraph$label.id
    newpathigraph$subgraphs_funs <- add_param_metaboID_subgraphs(pathigraph$subgraphs_funs)
    newpathigraph$effector.subgraphs_funs <- add_param_metaboID_subgraphs(pathigraph$effector.subgraphs_funs)
    newpathigraph$rl <- pathigraph$rl
    newpathigraph$fixed <- pathigraph$fixed
    newpathigraph
  }))
}
add_metabolite_to_mgi <- function(mgi, verbose = FALSE, basal.value = 0.5, basal.met.value = 1){
  newmgi <- list()
  newmgi$species <- mgi$species
  newmgi$all.genes <- mgi$all.genes

  newmgi$all.metabolite <- all_needed_metabolites(mgi$pathigraphs)
  if(verbose) message(length(newmgi$all.metabolite), " detected metabolites")
  newmgi$pathigraphs <- metabo_graphs(mgi$pathigraphs)
  if(verbose) message("metabo_graphs done!")
  newmgi$all.labelids <- mgi$all.labelids
  
  genes.vals.05 <- matrix(basal.value, ncol = 2, nrow = length(newmgi$all.genes), 
                          dimnames = list(newmgi$all.genes, c("1", "2")))
  metabolites.vals.05 <- matrix(basal.met.value, ncol = 2, nrow = length(newmgi$all.metabolite), 
                          dimnames = list(newmgi$all.metabolite, c("1", "2")))
  meta.05 <- NULL
  meta.05$pathigraphs <- newmgi$pathigraphs
  meta.05$all.labelids <- newmgi$all.labelids
  results.05 <- metabopathia(genes.vals.05, metabolites.vals.05, metaginfo = meta.05, test = FALSE, 
                         verbose = FALSE)
  if(verbose) message("results.05 done!")
  results.dec.05 <- metabopathia(genes.vals.05, metabolites.vals.05, meta.05, decompose = TRUE, 
                             test = FALSE, verbose = FALSE)
  if(verbose) message("results.dec.05 done!")
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

# Get the script directory
getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}
get_main_dir <- function(){
  # Check if running in RStudio
  is_rstudio <- Sys.getenv("RSTUDIO") == "1"
  if (is_rstudio) {
    codebase <- dirname(rstudioapi::getActiveDocumentContext()$path)
  } else {
    codebase <-  getScriptPath()
  }
  return(codebase)
}
# These are functions in the new version of hipathia that are not avail for old version ! that why I add here manually 
local_get_node_type <- function(metaginfo){
  types <- c("compound", "gene")
  names(types) <- c("circle", "rectangle")
  df <- lapply(metaginfo$pathigraphs, function(pg) {
    d <- data.frame(name = V(pg$graph)$name, label = V(pg$graph)$label, 
                    type = types[V(pg$graph)$shape])
    d <- d[!grepl("_func", d$name), ]
  })
  alltypes <- do.call(rbind, df)
  return(alltypes)
}

local_get_path_nodes <- function (metaginfo, path.names, decompose){
  pathnodes <- sapply(path.names, function(name) {
    pathway <- unlist(strsplit(name, "-"))[2]
    if (decompose == TRUE) {
      nodes <- V(metaginfo$pathigraphs[[pathway]]$subgraphs[[name]])$name
    }
    else {
      nodes <- V(metaginfo$pathigraphs[[pathway]]$effector.subgraphs[[name]])$name
    }
    allnodes <- paste(nodes, collapse = ", ")
  })
}
