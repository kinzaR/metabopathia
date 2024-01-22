# nodes_vals<-res$nodes.vals
# subgraph_list<-pathigraph$effector.subgraphs

all_path_values <- function (nodes_vals, subgraph_list, method = "maxmin", maxnum = 100, 
          tol = 1e-06, divide = FALSE, response_tol = 0) 
{
  path_vals <- matrix(0, ncol = ncol(nodes_vals), nrow = length(subgraph_list), 
                      dimnames = list(names(subgraph_list), colnames(nodes_vals)))
  signal_dif <- list()
  for (path in names(subgraph_list)) {
    dec_name <- unlist(strsplit(path, "\\-"))
    if (length(dec_name) == 4) {
      ininodes <- paste("N", dec_name[2], dec_name[3], 
                        sep = "-")
      endnode <- paste("N", dec_name[2], dec_name[4], 
                       sep = "-")
    }
    else if (length(dec_name) == 3) {
      endnode <- paste("N", dec_name[2], dec_name[3], 
                       sep = "-")
      sl <- subgraph_list[[path]]
      ininodes <- V(sl)$name[!V(sl)$name %in% get.edgelist(sl)[, 
                                                               2]]
    }
    else {
      stop("Error: Unknown path ID")
    }
    res <- hipathia:::path_value(nodes_vals, subgraph_list[[path]], 
                      ininodes, endnode, method, maxnum = maxnum, tol = tol, 
                      divide = divide, response_tol = response_tol)
    # For the moment I wil not change the way how Hipathia propagate the signal !
    path_vals[path, ] <- res[[1]]
    signal_dif[[path]] <- res[[2]]
  }
  return(list(path_vals, signal_dif))
}