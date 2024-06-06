## Needed lib fro metabolizer pipeline
library(preprocessCore)
library(igraph)
library(limma)
library(R.utils,quietly=T,verbose=F,warn.conflicts=F)
library("gplots")
## Metabolizer package 
metabolizer_files <- "metabolic_pathways/from_metabolizer"
source(file.path(metabolizer_files,"pretty_functions_Sep2016.R"))
# ## Here I have to recod metabolizer
infer_met_from_metabolizer <- function(data_set, all_sig_metabolites, species="hsa", output_folder, verbose=F){
  exp <- data_set$genes_vals
  design <- data_set$des
  output_folder <- filePath(output_folder,"metabolizer_res")
  # Load hsa_module_data 
  load(file.path(metabolizer_files, species, paste0(species,"_module_data_Dec2016.RData")))
  # HERE: if check.names=T then TCGA-6D-AA2E-01A changes into TCGA.6D.AA2E.01A.
  # HERE: EXP data has to be already xrefed here! pre-data has to do this task before
  ################# All module genes ################# 
  # HERE: this need optimization
  all_metabolites <- sort(unique(unlist(sapply(hsa_module_data, function(x) {
    nodes <- c(x$graphobject$SIF[,1],x$graphobject$SIF[,3])
    metabolites <- nodes[c(grep("C",nodes),grep("G",nodes))]
    return(metabolites)}))))
  
  all_module_genes_vec <- get.All.module.genes(hsa_module_data)
  
  ######### get all module rxns ###########
  
  all_module_rxn_vec <- get.All.module.rxns(hsa_module_data)
  
  ########### calculate rxn vals from gene exp ########### 
  
  rxn_gene_mat <- get.rxn.gene.matrix(hsa_module_data)
  
  gene_exp_mat <- get.gene.exp.of.module.genes(exp, all_module_genes_vec, min.exp=0.0001)
  
  rxn_vals_mat <- get.RXNvals.from.genes(all_module_rxn_vec, gene_exp_mat, rxn_gene_mat)
  metabolite_matrix <- mat.or.vec(nr = length(all_metabolites), nc = ncol(rxn_vals_mat))
  rownames(metabolite_matrix) <- all_metabolites
  metabolite_matrix[,] <- 1
  
  # probably "discard_compounds" is not any more important because I will use metabolites with default value 1.# Cankut descesion ! need to be reconcidred ?
  #hsa_module_data <- discard_compounds(hsa_module_data)
  #HERE : compute.node.signal2 is declered twice
  results_module_activities <- methyways(hsa_module_data,rxn_vals_mat=rxn_vals_mat,expbased=T,fluxbased=F,verbose = F, default_value=0.5,metabolitematrix = metabolite_matrix)
  
  #### do it as function
  results_module_activities_list <- sapply(results_module_activities, function(x){ x[[1]][1] })
  results_module_activities_matrix <- do.call("rbind",results_module_activities_list)
  
  # write.table(x=results_module_activities_matrix,file=paste0(output_folder,"/module_vals.txt"),quote=F,sep="\t")
  
  #******************************************************************************************************#
  # rownames(results_module_activities_matrix) <- gsub(".*_","",rownames(results_module_activities_matrix))
  # comp_names<-gsub(".*_","",rownames(results_module_activities_matrix))
  # dup_comp_names <- duplicated(comp_names) %>% comp_names[.]
  # dup_original_names<-rowAnys(sapply(dup_comp_names, function(d) {(grepl(d, rownames(results_module_activities_matrix)))})) %>% rownames(results_module_activities_matrix)[.]
  # results_module_activities_matrix[dup_original_names,] %>% View
  # comp_dup_inSigPaths<-dup_comp_names %in% metabo_pathways$all.metabolite %>% dup_comp_names[.]
  # ss<-sapply(metabo_pathways$pathigraphs, function(ig){
  #   who_is_here <- comp_dup_inSigPaths %in% vertex_attr(ig$graph, name = "metaboID") %>% .[!is.na(.)]
  #   if(any(who_is_here)){
  #     return(paste0(ig$path.name,":", comp_dup_inSigPaths[who_is_here]))
  #   }
  # })%>% unlist
  # hsa04930                           hsa05200                           hsa05211 
  # "Type II diabetes mellitus:C00022"        "Pathways in cancer:C00122"      "Renal cell carcinoma:C00122" 
  ###HERE: I have to decide about duplicated metabolite activities!
  #******************************************************************************************************#
  results_module_activities_matrix_curated <- results_module_activities_matrix %>% as.data.frame() %>% mutate(keggID= gsub(".*_","",rownames(results_module_activities_matrix))) %>% group_by(keggID) %>%
    summarise(across(everything(), mean)) %>% tibble::column_to_rownames(var = "keggID")
  metabo_vals <- results_module_activities_matrix_curated %>% filter(rownames(.) %in% all_sig_metabolites)# keep only relevant metabolites
  return(metabo_vals)
}
