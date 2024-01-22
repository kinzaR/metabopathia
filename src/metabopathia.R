# This function is designed to calculate the activation level of each circuit (subpathway) within a pathway for each sample in the experiment. 
# As inputs, in addition to the pathways object, it takes node values. 
#   - For gene nodes, this involves the matrix of gene expression, 
#   - and for metabolite nodes, it requires the metabolic concentration matrix. 
# The output consists of the computed activation values for all loaded circuits.
### Rules for signal computation
# A protein passes the signal through the following factors:
#   1-The protein must be present.
#   2-Another protein must activate it.
#   3-If it is an enzyme and the previous node is a metabolite, a threshold of concentration has to be reached.
# This method computes signal transduction (alternative terminology for metabolic pathways needs to be proposed) based on these steps:
#   1-Quantification of the presence of a particular gene as a normalized value between 0 and 1 (as a proxy for protein presence).
#   2-Check if there is an enzyme that interacts with the metabolite captured in the proposed data.
#   3-Compute the signal value passing through a node, taking into account:
#     a. The level of expression of each gene inside the node.
#     b. The metabolite concentration of each metabolite that interacts with this node.
#     c. The intensity of the signal arriving at it, either activations or inhibitions.
# The signal value of the circuit is the signal value through the last node of the circuit
## Imputation approch:
# The missing values for genes and metabolites which are needed by the method to compute the signal are added by the function,
# assigning to each sample the median of the matrix. 
# Please note: A high ratio of missingness may lead to unrepresentative results.
source("src/utils.R")
source("src/nodes_values_from_all.R")
#The iterative algorithm that calculate the signal
metabopathia <- function (genes_vals, metabo_vals, metaginfo, uni.terms = FALSE, GO.terms = FALSE, 
          custom.terms = NA, sel_assay = 1, decompose = FALSE, maxnum = 100, 
          verbose = TRUE, tol = 1e-06, test = TRUE) {

  if (is(genes_vals, "SummarizedExperiment")) {
    coldata <- colData(genes_vals)
    genes_vals <- assay(genes_vals, sel_assay)
  }
  else {
    cols <- colnames(genes_vals)
    coldata <- data.frame(cols = cols, stringsAsFactors = FALSE)
  }
  if (test == TRUE) {
    if (is.null(genes_vals)) 
      stop("Missing input genes matrix")
    if (is.null(metabo_vals)) 
      stop("Missing input metabolites matrix")
    if (is.null(metaginfo)) 
      stop("Missing pathways object")
    if(ncol(genes_vals) != ncol(metabo_vals)){
      stop("Input genes matrix and meatbolites matrix has different column length")
    }
    hipathia:::test_matrix(genes_vals)
    hipathia:::test_matrix(metabo_vals)
    test_metabo_pathways_object(metaginfo)
    hipathia:::test_tolerance(tol)
  }
  pathigraphs <- metaginfo$pathigraphs
  genes_vals <- hipathia:::add_missing_genes(genes_vals, genes = metaginfo$all.genes)
  metabo_vals <- add_missing_metabolites(metabo_vals, metabolites = metaginfo$all.metabolite)
  results <- list()
  if (verbose == TRUE) 
    cat("Computing pathways...\n")
  
  # pathigraph <- pathigraphs$hsa04010
  results$by.path <- lapply(pathigraphs, function(pathigraph) {
    res <- list()
    # suppressWarnings(res$nodes.vals <- hipathia:::nodes_values_from_genes(genes_vals, 
    #                                                            pathigraph$graph)) # here I have to change: metabolite vals are 1 instead of 0.5 !
    suppressWarnings(res$nodes.vals <- nodes_values_from_all(genes_vals, metabo_vals,
                                                               pathigraph$graph)) # if not found metab ID is now 1 , exceptions ready to detect drugs!
    if (decompose == FALSE) {
      respaths <- hipathia:::all_path_values(res$nodes.vals, pathigraph$effector.subgraphs, 
                                  maxnum = maxnum, tol = tol)
    }
    else {
      respaths <- hipathia:::all_path_values(res$nodes.vals, pathigraph$subgraphs, 
                                  maxnum = maxnum, tol = tol)
    }
    res$path.vals <- respaths[[1]]
    res$convergence <- respaths[[2]]
    return(res)
  })
  nodes <- do.call("rbind", lapply(results$by.path, function(x) x$nodes.vals))
  nodes_rd <- DataFrame(metaginfo$all.labelids[rownames(nodes), 
  ], node.name = hipathia:::get_node_names(metaginfo, rownames(nodes)), 
  node.type = hipathia:::get_node_type(metaginfo)$type, node.var = apply(nodes, 
                                                              1, var))
  nodes_se <- SummarizedExperiment(list(nodes = nodes), rowData = nodes_rd, 
                                   colData = coldata)
  paths <- do.call("rbind", lapply(results$by.path, function(x) x$path.vals))
  paths_rd <- DataFrame(path.ID = rownames(paths), path.name = hipathia::get_path_names(metaginfo, 
                                                                              rownames(paths)), path.nodes = hipathia:::get_path_nodes(metaginfo, 
                                                                                                                            rownames(paths), decompose = decompose), decomposed = decompose)
  paths_se <- SummarizedExperiment(list(paths = paths), rowData = paths_rd, 
                                   colData = coldata)
  se_list <- list(nodes = nodes_se, paths = paths_se)
  
  ## This has to be checked again
  if (uni.terms == TRUE) {
    if (verbose == TRUE) 
      cat("\nComputing Uniprot terms...\n")
    unis_se <- hipathia:::quantify_funs(paths_se, metaginfo, "uniprot")
    se_list$uni.terms <- unis_se
  }
  if (GO.terms == TRUE) {
    if (verbose == TRUE) 
      cat("\nComputing GO terms...\n")
    gos_se <- hipathia:::quantify_funs(paths_se, metaginfo, "GO")
    se_list$GO.terms <- gos_se
  }
  if (!is.na(custom.terms)) {
    if (verbose == TRUE) 
      cat("\nComputing custom terms...\n")
    custom_se <- hipathia:::quantify_funs(paths_se, metaginfo, dbannot)
    se_list$custom.terms <- custom_se
  }
  resmae <- MultiAssayExperiment(se_list)
  if (verbose == TRUE) 
    message("DONE")
  return(resmae)
}
