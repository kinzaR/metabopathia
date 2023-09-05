library(hipathia)
pathways <- hipathia::load_pathways(species = "hsa", pathways_list = c("hsa04071"))
ig <- pathways$pathigraphs$hsa04071$graph
genes_vals <- hipathia::brca_data
##
# here I will try  to add metabolite ID into geneList

###
nodes_values_from_all <- function(genes_vals, ig, summ = "per90"){
  genes_list <- V(ig)$genesList
  names(genes_list) <- V(ig)$name
  genes_list <- genes_list[!grepl("_func", names(genes_list))]
  
  nodes_vals <- matrix(NA, nrow = length(names(genes_list)), 
                       ncol = ncol(genes_vals), dimnames = list(names(genes_list), 
                                                                colnames(genes_vals)))
  for (node_name in names(genes_list)) {
    genes <- genes_list[[node_name]]
    if ("/" %in% genes) {
      lists <- get_genes_lists(genes)
      probabilities_mat <- matrix(NA, nrow = 0, ncol = ncol(genes_vals))
      for (list1 in lists) {
        if (length(list1) > 1) {
          gv_l1 <- genes_vals[list1, , drop = FALSE]
          prob <- summarize_probabilities(gv_l1, summ)
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
        nodes_vals[node_name, ] <- rep(1, ncol(nodes_vals))
      }
    }
  }
  return(nodes_vals)
}
