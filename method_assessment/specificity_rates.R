#!/usr/bin/env /opt/R/4.3.1/bin/Rscript
# Script for assessing the Specificity of Metabopathia Vs Hipathia: False Positive Rate 
# Author: Kinza Rian

# Load necessary libraries
message("You are using :",R.version.string)
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(hipathia))
suppressPackageStartupMessages(library(dplyr))

# Get the script directory
getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}
# Check if running in RStudio
is_rstudio <- Sys.getenv("RSTUDIO") == "1"
if (is_rstudio) {
  codebase <- dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  codebase <-  getScriptPath()
}
# Set the working directory to the script directory
setwd(codebase)

# Source necessry files
source("method_assessment/utils.R")
source("src/opt_validator.R")
source("src/report.R")
source("src/data_pre.R")
source("src/utils.R")
source("src/metabopathia.R")

# Define command-line options
option_list <- list(
  make_option(c("-s","--spe"), type = "character", default = "hsa",
              help = "Species variable. Allowed choices: 'hsa', 'mmu', 'rno'. (default: hsa)"),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE, 
              help = "Enable verbose mode."),
  make_option(c("-e", "--exp_file"), type = "character", default = NULL,
              help = "Path to the expression file."),
  make_option(c("-m", "--met_file"), type = "character", default = NULL,
              help = "Path to the metabolomics concentration file."),
  make_option(c("-t", "--met_type"), type = "character", default = "inferred",
              help = "Allowed values: 
              inferred: Infer production of metabolites from Metabolizer (doi: 10.1038/s41540-019-0087-2).
              perturbations: Study perturbance in metabolite concentration.
              concentration_matrix: Using metabolomics concentration matrix."),
  make_option(c("-d", "--design_file"), type = "character", default = NULL,
              help = "Path to the design file."),
  make_option("--group1", type = "character", default = NULL,
              help = "Label of the first group."),
  make_option("--group2", type = "character", default = NULL,
              help = "Label of the second group to be compared (Reference condition)."),
  make_option(c("-p", "--pathways_list"), type = "character", default = NULL,
              help = "Vector of the IDs of the pathways to load. By default, all available pathways are loaded. Example: '04014,04015'."),
  make_option("--paired", action = "store_true", default = FALSE,
              help = "Boolean, whether the samples to be compared are paired. If TRUE, function wilcoxsign_test from package coin is used. If FALSE, function wilcox.test from package stats is used."),
  make_option("--decompose", action = "store_true", default = FALSE,
              help = "Boolean, whether to compute the values for the decomposed subpathways. By default, effector subpathways are computed."),
  make_option("--design_type", type = "character", default = "categorical",
              help = "Type of design. Allowed values: 'categorical' or 'continuous'. Default is 'categorical'."),
  make_option("--adjust", action = "store_true", default = TRUE,
              help = "Boolean, whether to adjust the p.value with Benjamini-Hochberg FDR method. Default is TRUE."),
  make_option("--conf.level", type = "double", default = 0.05,
              help = "Level of significance. By default 0.05."),
  make_option("--difexp", action = "store_true", default = TRUE,
              help = "Boolean, whether to perform differential expression analysis."),
  make_option("--GO.terms", action = "store_true", default = FALSE,
              help = "Boolean, whether to compute functional analysis with Gene Ontology terms."),
  make_option("--uni.terms", action = "store_true", default = FALSE,
              help = "Boolean, whether to compute functional analysis with Uniprot keywords."),
  make_option("--custom.terms", type = "character", default = NA,
              help = "Path to a file containing a data.frame with the custom annotation of the genes to the functions. First column are gene symbols, second column the functions."),
  make_option("--analysis", type = "character", default = "overlay",
              help = "Type of analysis. Allowed values:'overlay','ORA', 'compare', 'predictor_test', 'predictor_train', 'variant_interpreter', 'drug_repurposing'. Default is 'compare'.\n\n\
                      Differential Signaling 'compare': Provides an estimation of significant cell signaling activity changes across different conditions.\n\
                      Predictor 'predictor': Allows you to train a prediction-test and test it with different data.\n\
                      Variant Functional Interpretation 'variant_interpreter': Provides an estimation of the potential impact of genomic variation on cell signaling and cell functionality.\n\
                      Drug Repurposing 'drug_repurposing': drug repurposing."
  ),
  make_option("--output_folder", type = "character", default = "tmp",
              help = "Output folder path. Default is 'tmp'."),
  make_option("--hipathia", action = "store_true", default = FALSE,
              help="Enable calculation using HiPathia for result comparison between MetaboPathia and HiPathia."),
  make_option("--example", action = "store_true", default = FALSE,
              help = "Load variables from the example config file (src/example1.R)."),
  make_option("--ready_mgi", action = "store_true", default = TRUE,
              help = "Load ready mgi pre-processed."),
  make_option("--senario", type = "character", default = "i",
              help = "senario of data selection"),
  make_option("--N", type = "integer", default = 20,
              help = "number of selected samples")
)
# Parse command-line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# be carefull this is a forced hardcoded working-example to be removed !!!
if(opt$example || is_rstudio) source("src/example1.R")
if(opt$verbose){
  message("The recieved options are :")
  str(opt)
}
if(opt$met_type == "inferred") source("src/local_metabolizer.R")
### Piplines ###
#NOTE: Check if changing this to a concatination string is less expensive and easier! 
# source(paste0("src/",opt$analysis,"_pipeline.R")) # after all integrations!
switch (opt$analysis,
        "overlay" = source("src/overlay_pipeline.R"),
        "compare" = source("src/compare_pipeline.R") #NOTE: also  needed for hipathia =T
)
# Validate options
validate_options(opt)

# Extract values from options
attach(opt) # some recomendations are to avoid the use of attach but still not convincing for me in this situation!

# Set default values if not specified by the user
if (is.null(opt$pathways_list)) { #if (!opt$pathways_list) { the other one is more acurate!
  pathways_list <- NULL
}

output_folder <- create_output_folder(output_folder, verbose = verbose)
status("  0", "Analysis is started", output_folder)
###### Workflow #####
## Step 1: Data pre-processing
if(verbose) message("Loading data all data...") 
#dbplyr V=2.3.4 because leatest has a bug
data_set <- load_data(exp_file, met_file = ifelse(exists("met_file"), met_file, NA), met_type, design_file, group1, group2, output_folder, design_type, analysis, spe, verbose) # analysis here is not used!
status(" 20", "Data loaded successfully", output_folder)

## Step 2: MGI preparation & filtering
### from data filter pathways_list ( think about filtering by metabolimocs type in order to avoid high prensentage of missingness)
# ...
### load pathways with pathway list:  preprocessed KEGG pathway
if(verbose) message("Loading pathways...") #I have to load from prepared MGI already ! time consuming 
# Load the pre-processed pathway object
if(!ready_mgi | hipathia | analysis=="overlay")
  pathways <- hipathia::load_pathways(species = spe, pathways_list = pathways_list)
# modules <- 
## adaptation of the MGI: to be removed, because I have to load it already prepared
if(analysis!="overlay") { # new mgi metabo will not be needed for overlay!
  if(ready_mgi) {
    #metabo_pathways <-readRDS("pathways/metabo_mgi_v1.0.0.RDS")# this had bug
    metabo_pathways <-readRDS("pathways/metabo_mgi_v1.1.0.RDS")
  }else metabo_pathways <-add_metabolite_to_mgi(pathways, verbose = verbose, basal.met.value = 1) #saveRDS(object = metabo_pathways , file = "pathways/metabo_mgi_v1.1.0.RDS")
}
status(" 40", "Pathways loaded successfully", output_folder)

## Here I have to do met_vals
if(is.null(data_set$metabo_vals)){
  if(exists("metabo_pathways")){
    all_metabolites <- metabo_pathways$all.metabolite
  }else{
    if(exists("pathways")){
      all_metabolites <- all_needed_metabolites(pathways$pathigraphs)
    }else{
      stop("Pathway meta-graph info was not loaded successfully!")
    }
  }
  # HERE: be carfll with cond1 and 2 because the reference for metabolizer maybe is not the same as hipathia !
  data_set$metabo_vals <- switch (met_type,
                                  "inferred" = {
                                    source("src/local_metabolizer.R")
                                    infer_met_from_metabolizer(data_set, all_sig_metabolites=all_metabolites, species=spe, output_folder, verbose)},
                                  "perturbations" = {
                                    source("src/get_perturbed_met_vals.R")
                                    get_perturbed_met_vals(data_set$metabo_vals, data_set$des, group1, all_metabolites, verbose)
                                  }
  )
}
## Generate sample from all data
################################################################################
############################## specificity_rates ###############################
################################################################################
# for senario in ("i", "ii", "iii") 
# specificity_rates <- lapply(c("i", "ii", "iii") , function(senario){
specificity_rates <- lapply(setNames(c("i"),c("i")) , function(senario){
  # lapply(c(20, 50, 70, 100) , function(N){
  lapply(setNames(c(20), c(20)) , function(N){
    # here iterations c(1:100)
    # lapply(setNames(c(1:100), c(1:100)), function(it){
    lapply(setNames(c(1:2), c(1:2)), function(it){
      if(verbose) message("***Senario:",senario," **N: ",N," *it:",it)
      new_data_set <- get_N_data_perSenario(data_set, N=N, senario=senario) # genes scales per N samples
      ## Step 3: Signal propagation : Pathway activation computation
      if(verbose) message("Propagating signaling...")
      ## Functiona annotation 
      if(!is.na(custom.terms)){
        custom.terms <- read.table(custom.terms, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
      }
      if(analysis!="overlay")
        ## Here I have to add data_set$metabo_vals
        metdata <- metabopathia(genes_vals = new_data_set$genes_vals, metabo_vals = as.matrix(new_data_set$metabo_vals),
                                metaginfo = metabo_pathways,
                                uni.terms = uni.terms, GO.terms = GO.terms, custom.terms = custom.terms,
                                decompose = decompose, verbose=verbose)
      if(hipathia | analysis=="overlay"){
        # This has to be adapted to the web version 4.1.2
        hdata <- hipathia(genes_vals = new_data_set$genes_vals,
                          pathways,
                          uni.terms = uni.terms, GO.terms = GO.terms, custom.terms = custom.terms,
                          decompose = decompose, verbose=verbose)
      }
      
      ### path vals extraction
      ####metabopathia
      if(analysis!="overlay"){
        met_path_vals <- get_paths_data(metdata, matrix = T)
        met_path_vals <- normalize_paths(met_path_vals, metabo_pathways)
        met_node_vals <- hipathia::get_nodes_data(metdata, matrix = T) %>% as.data.frame()
      }
      ###hipathia
      if(hipathia | analysis=="overlay"){
        h_path_vals <- get_paths_data(hdata, matrix = T)
        h_path_vals <- normalize_paths(h_path_vals, pathways)
        h_node_vals <- hipathia::get_nodes_data(hdata, matrix = T) %>% as.data.frame()
      }
      #Here to check that only infered metabolites nodes are changes: (met_node_vals == h_node_vals) %>% rowSums() %>% .[.!=20] %>% names %>% metabo_pathways$all.labelids[.,"label"] %>% .[!grepl("\\*",.)] %>% unique %>% length()
      status(" 60", "Signal propagation computed successfully", output_folder)
      
      ## Step 4: Scenarios & pipelines : 
      # Here : maybe switch will be cleaner
      ## Step 4.0: Differential metabolic overlay
      if(analysis=="overlay"){
        stop("NOT implemeted yet!")
      }
      ## Step 4.1: Differential Activity Analysis
      # I have here to discuss why I have used the statiticakl test different for each!
      if(analysis=="compare"){
        met_results <- compare_pipeline(metdata, groups=new_data_set$des$group, 
                                        expdes=unique(new_data_set$des$group)[1], g2 = unique(new_data_set$des$group)[2],
                                        path.method = "wilcoxon", node.method = "limma", fun.method = "wilcoxon",
                                        order = FALSE, paired = FALSE, adjust = adjust, conf.level = 0.05, sel_assay = 1)
        if(hipathia){
          hi_results <- compare_pipeline(hdata, groups=new_data_set$des$group, expdes=unique(new_data_set$des$group)[1], g2 = unique(new_data_set$des$group)[2],
                                         path.method = "wilcoxon", node.method = "limma", fun.method = "wilcoxon",
                                         order = FALSE, paired = paired, adjust = adjust, conf.level = 0.05, sel_assay = 1)
          
        }
        status(" 80", "Differential Activity Analysis completed successfully", output_folder)
      }
      
      ################ Start : To be removed after testing the hipathia function with 3_terms
      # ## Step 5: Functional analysis: hipathia::quantify_terms()
      # uniprot_vals <- quantify_terms(results, pathways, dbannot = "uniprot") is hipathia:::quantify_funs() and it encapsulated inside hipathia/metabopathia function 
      ################ End : To be removed after testing the hipathia function with 3_terms
      
      ## Step 6: Results 
      ## visualization functions
      # Results overview: Summary of UP & DOWN nodes, paths and functions
      met_results$DAoverview <- DAoverview_plotless(DAdata = met_results, conf.level = opt$conf.level,adjust = opt$adjust)
      if(hipathia) {
        hi_results$DAoverview <- DAoverview_plotless(DAdata = hi_results, conf.level = opt$conf.level,adjust = opt$adjust)
      }
      
      # Results summary by pathway: summary of the number of paths altered in the n most altered pathways
      met_results$DAsummary <- hipathia::DAsummary(DAdata = met_results, n = 10) # Top altered pathways
      if(hipathia) {
        hi_results$DAsummary <- hipathia::DAsummary(DAdata = hi_results, n = 10) # Top altered pathways
      }
      return(met_results)
    })
  })
})
