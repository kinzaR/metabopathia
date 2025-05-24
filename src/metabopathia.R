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
# Computation of Signal:
# The intensity of the propagated signal is calculated using an iterative algorithm, starting from the first node(s) in a sub-pathway until reaching the last node. (the concept of a sub-pathway/circiut has to be introduced clearly before this!)
# Biologically speaking, the first node represents the receptors within a cell, and the last nodes represent effector protein. 
# In subsequent steps, we will annotate them with cellular functions and infer the functional activity of the studied condition within the cell.
# The initial input signal arriving at the first node is set to 1, assuming that the signal reaching the receptors is at its maximum value.
# As mentioned before, all values must be scaled between 0 and 1. Here, 0 signifies inactivity, and 1 denotes full activation. However, these values may not be meaningful until placed in a comparative context.
# Then, for each node 'n' in the subpathway, the signal value will be the product of the normalized/scaled value of that node calculated previously, the set of arrived inhibition signal values, and activation signal values.
## parameters:
# decompose: indicates whether to use effector subpathways or decomposed subpathways.
# Some comments :
#   very good as you mention early intro for the sub pathway
# clarification of node types -> more description as you mention is good idea. you mentioned that first node represents receptors within a cell and the last represent a TF , consider adding adescription to the other types that might be present their roles
# Biological speaking:  why are we doing this? imagine its like playing detective in the cell, by crunching those numbers, we can unveil cellular mysteries and solve biological puzzles, and may beeven find out why your cells have been mis behaving- real applications of this computation in undersytanding rhe cellular functions or disease mechanism
# scaling points : elaborate why values must be scaled between 0 and 1? what scales represet biologicaly  discuss how 0 signifying in activity and 1 full activation
# you would say about the inhibition and activation. elabort how this can be determined? are those experimental data? predicted/calculated values or both? the idea how these signal influence a signal values
# the iterative process : describe the iterative process in more details here... what are the steps in the propagating the signal from point x to point z ? any math equation govern the process
# explain how the calculated signal values at each node are used for tits functional annotation ? what the downstream analysis?
#   feedback loops, refer if the method consider that bcoz it can play a crucial role in signal stability and propagation
# Clear confusing in cocept : effector protein vs TF
source("src/utils.R")
source("src/nodes_values_from_all.R")
#The iterative algorithm that calculate the signal
metabopathia <- function (genes_vals, metabo_vals, metaginfo, uni.terms = FALSE, GO.terms = FALSE, 
          custom.terms = NA, sel_assay = 1, decompose = FALSE, maxnum = 100, 
          verbose = TRUE, tol = 1e-06, test = TRUE) {
  if (is(genes_vals, "SummarizedExperiment")) {
    coldata <- colData(genes_vals)
    genes_vals <- assay(genes_vals, sel_assay)
  }else {
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
  metabo_vals <- add_missing_metabolites(metabo_vals, 
                                         metabolites = metaginfo$all.metabolite,
                                         default = 1)
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
  nodes_rd <- DataFrame(row.names = rownames(nodes),metaginfo$all.labelids[rownames(nodes), 
  ], node.name = hipathia:::get_node_names(metaginfo, rownames(nodes)), 
  # node.type = hipathia:::get_node_type(metaginfo)$type, node.var = apply(nodes, 1, var))
  node.type = local_get_node_type(metaginfo)$type, node.var = apply(nodes, 1, var))
  nodes_se <- SummarizedExperiment(list(nodes = nodes), rowData = nodes_rd, 
                                   colData = coldata)
  paths <- do.call("rbind", lapply(results$by.path, function(x) x$path.vals))
  paths_rd <- DataFrame(row.names =rownames(paths), path.ID = rownames(paths), path.name = hipathia::get_path_names(metaginfo, 
                                                                              #rownames(paths)), path.nodes = hipathia:::get_path_nodes(metaginfo, rownames(paths), decompose = decompose), decomposed = decompose)
                                                                              rownames(paths)), path.nodes = local_get_path_nodes(metaginfo, rownames(paths), decompose = decompose), decomposed = decompose)
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
