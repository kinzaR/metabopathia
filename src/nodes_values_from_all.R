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
  genes_list <- V(ig)$genesList
  metabolites_list <- V(ig)$metaboID
  names(genes_list) <- V(ig)$name
  names(metabolites_list) <- V(ig)$name
  genes_list <- genes_list[!grepl("_func", names(genes_list))]
  
  # Assuming that colnames(genes_vals) == colnames(metabo_vals)
  nodes_vals <- matrix(NA, nrow = length(names(genes_list)), 
                       ncol = ncol(genes_vals), dimnames = list(names(genes_list), 
                                                                colnames(genes_vals)))
  for (node_name in names(genes_list)) {
    genes <- genes_list[[node_name]]
    if ("/" %in% genes) {
      lists <- hipathia:::get_genes_lists(genes)
      probabilities_mat <- matrix(NA, nrow = 0, ncol = ncol(genes_vals))
      for (list1 in lists) {
        if (length(list1) > 1) {
          gv_l1 <- genes_vals[list1, , drop = FALSE]
          prob <- hipathia:::summarize_probabilities(gv_l1, summ)
        }
        else {
          prob <- genes_vals[list1, , drop = FALSE]
        }
        probabilities_mat <- rbind(probabilities_mat, 
                                   prob)
      }
      nodes_vals[node_name, ] <- colMins(probabilities_mat, 
                                         na.rm = TRUE)
    }
    else {
      glist <- genes_list[[node_name]]
      if (length(glist) > 1) {
        gv <- genes_vals[glist, , drop = FALSE]
        nodes_vals[node_name, ] <- hipathia:::summarize_probabilities(gv, 
                                                           summ)
      }
      else if (length(glist) == 1 && !is.na(glist)) {
        dm <- data.matrix(genes_vals[glist, , drop = FALSE])
        nodes_vals[node_name, ] <- dm
      }
      else {
        # if metabolite ID exist 
        if(!is.na(metabolites_list[[node_name]])){
          metabolite <- metabolites_list[[node_name]]
          dm <- data.matrix(metabo_vals[metabolite, , drop = FALSE])
          nodes_vals[node_name, ] <- dm
        }else{
          nodes_vals[node_name, ] <- rep(1, ncol(nodes_vals))
        }
      }
    }
  }
  return(nodes_vals)
}
