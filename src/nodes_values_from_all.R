# ig <- pathigraph$graph
# Pathways are structured as graphs, consisting of interconnected nodes and edges. 
# Certain nodes may encompass multiple genes, representing distinct isoforms of the protein or members within the same gene family, among other possibilities.
# As each gene has its unique level of expression, the initial step involves summarizing this information into a score.
# Which will characterizes the overall value of the node.
# In the case of metabolite nodes, two methods are proposed to align the values of these nodes with the methodology:
# - Binarization: Conversion of values into binary form.
# - Metabolite/Protein Ratio of Activation Association: Establishing a ratio between metabolite and protein activation values.
# Check if the library is loaded, if not, load it
if (!require("hipathia", quietly = TRUE)) {
  library("hipathia")
}
nodes_values_from_all <- function(genes_vals, metabo_vals, ig, summ = "per90"){
  # Assuming that colnames(genes_vals) == colnames(metabo_vals)
  nodes <- igraph::as_data_frame(ig, what = "vertices") %>% 
    dplyr::select(c("name", "shape", "genesList","metaboID")) %>% 
    filter(!grepl("_func", name))
  nodes_vals <- do.call(rbind,
                        lapply(nodes$name, function(n) get_node_value(nodes[n,], genes_vals, metabo_vals, summ = "per90")))
  rownames(nodes_vals) <- nodes$name 
  return(nodes_vals)
}
get_node_value <- function(node, genes_vals, metabo_vals, summ = "per90"){
  node_type <- get_node_type(node)
  if(is.na(node_type)) return(rep(1, ncol(genes_vals))) #Drug node or Unkown type ! 
  if (node_type == "complex") {
    complex_list <- get_complex_lists(node$genesList[[1]])
    return(do.call(rbind,lapply(complex_list, function(component){
      if(length(component)>1){
        gv_l1 <- genes_vals[component, , drop = FALSE]
        hipathia:::summarize_probabilities(gv_l1, summ)
      }else{
        genes_vals[component, , drop = FALSE]
      }
    })) %>% colMins(., na.rm = TRUE))
  }
  if(node_type == "protein_family"){
    gv <- genes_vals[node$genesList[[1]], , drop = FALSE]
    return(hipathia:::summarize_probabilities(gv,summ))
  }
  if(node_type == "protein"){
    return(data.matrix(genes_vals[node$genesList[[1]], , drop = FALSE]))
  }
  if(node_type == "metabolite"){
    #Here normal metabolites
    return(data.matrix(metabo_vals[node$metaboID, , drop = FALSE]))
  }
  if(node_type =="complex_metabolite"){
    complex_list <- strsplit(node$metaboID, ",/,")[[1]]
    return(data.matrix(metabo_vals[complex_list, , drop = FALSE]) %>%
      colMins(., na.rm = TRUE))
    
  }
  #if no type found then is a drug or whatever:/
  return(rep(1, ncol(genes_vals)))
}
get_node_type <- function(node){
  if ("/" %in% node$genesList[[1]]) return("complex") # Here some complex has metabolites as well need to be detected
  if(length(node$genesList[[1]]) > 1) return("protein_family")
  if(length(node$genesList[[1]]) == 1 && !is.na(node$genesList[[1]])) return("protein")
  if(!is.na(node$metaboID)){
    if(!grepl("/", node$metaboID)) return("metabolite")
    if(is.na(node$genesList[[1]])) return("complex_metabolite")
  }
  return(NA)
}
get_complex_lists <- function(complex) {
  #it does what hipathia:::get_genes_lists(complex) does, but in one line!
  return(split(complex[complex != "/"], cumsum(complex == "/")[complex != "/"]))
}

