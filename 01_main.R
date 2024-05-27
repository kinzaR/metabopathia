#!/usr/bin/env /opt/R/4.3.1/bin/Rscript
### !/usr/bin/env Rscript
# Main Script for metabopathia project: 01_main.R
# Author: Kinza Rian
# Description: This script will be the core of Metabopathia web-server with the provided omics.
# Usage: Rscript 01_main.R --species hsa or dirretly ./01_main.R --species hsa

# Load necessary libraries
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
# Load the opt_validator.R file
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
              help = "Load ready mgi pre-processed.")
  )
# Parse command-line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# be carefull this is a forced hardcoded working-example to be removed !!!
if(opt$example || is_rstudio) source("src/example1.R")
if(opt$verbose){
  message("The recieved options are :")
  str(opt)
}
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
if(verbose) message("Loading data...") 
#dbplyr V=2.3.4 because leatest has a bug
data_set <- data_pre(exp_file, met_file, design_file, group1, group2, output_folder, design_type, analysis, spe, verbose) # analysis here is not used!
status(" 20", "Data preprocessed successfully", output_folder)

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
    metabo_pathways <-readRDS("pathways/metabo_mgi_v1.0.0.RDS")
  }else metabo_pathways <-add_metabolite_to_mgi(pathways)
}

status(" 40", "Pathways loaded successfully", output_folder)

## Step 3: Signal propagation : Pathway activation computation
if(verbose) message("Propagating signaling...")
## Functiona annotation 
if(!is.na(custom.terms)){
  custom.terms <- read.table(custom.terms, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
}
if(analysis!="overlay")
  metdata <- metabopathia(genes_vals = data_set$genes_vals, metabo_vals = data_set$metabo_vals,
                        metaginfo = metabo_pathways,
                        uni.terms = uni.terms, GO.terms = GO.terms, custom.terms = custom.terms,
                        decompose = decompose, verbose=verbose)
if(hipathia | analysis=="overlay"){
  hdata <- hipathia(genes_vals = data_set$genes_vals,
                  pathways,
                  uni.terms = uni.terms, GO.terms = GO.terms, custom.terms = custom.terms,
                  decompose = decompose, verbose=verbose)
}

### path vals extraction
####metabopathia
if(analysis!="overlay"){
  met_path_vals <- get_paths_data(metdata, matrix = T)
  met_path_vals <- normalize_paths(met_path_vals, metabo_pathways)
}
###hipathia
if(hipathia | analysis=="overlay"){
  h_path_vals <- get_paths_data(hdata, matrix = T)
  h_path_vals <- normalize_paths(h_path_vals, pathways)
}

status(" 60", "Signal propagation computed successfully", output_folder)

## Step 4: Scenarios & pipelines : 
# Here : maybe switch will be cleaner
## Step 4.0: Differential metabolic overlay
if(analysis=="overlay"){
  # 23 may : will caluculate normal comparaison then in the visualization metabolite will join the party
  # # ## this will calculate only hipathia and add the colors of metabolites
  # met_results <- overlay_pipeline(metdata, groups=data_set$des$group, expdes=group1, g2 = group2,
  #                                 path.method = "wilcoxon", node.method = "limma", fun.method = "wilcoxon",
  #                                 order = FALSE, paired = paired, adjust = adjust, conf.level = 0.05, sel_assay = 1)
  # 
  # if(hipathia){
  # met_results <- hi_results <- compare_pipeline_overlaying(hdata, metabo_vals = data_set$metabo_vals, groups=data_set$des$group, expdes=group1, g2 = group2,
  #Here met_resultre are not same as hipathia_results
  met_results <- compare_pipeline_overlaying(hdata,metabo_vals=data_set$metabo_vals, groups=data_set$des$group, expdes=group1, g2 = group2,
                                   path.method = "wilcoxon", node.method = "limma", fun.method = "wilcoxon",
                                   order = FALSE, paired = paired, adjust = adjust, conf.level = 0.05, sel_assay = 1)
  # }
  if(hipathia){
    hi_results <- hipathia::DAcomp(hidata = hdata, groups = data_set$des$group , expdes = group1, g2 = group2,
                                   path.method = "wilcoxon", node.method = "limma", fun.method = "wilcoxon",
                               order = FALSE, paired = paired, adjust = adjust, conf.level = 0.05, sel_assay = 1)
    
  }
  status(" 80", "metabolomics overlay with a Differential Activity Analysis completed successfully", output_folder)

}
## Step 4.1: Differential Activity Analysis
# I have here to discuss why I have used the statiticakl test different for each!
if(analysis=="compare"){
  met_results <- compare_pipeline(metdata, groups=data_set$des$group, expdes=group1, g2 = group2,
                                  # path.method = "wilcoxon", node.method = "limma", fun.method = "wilcoxon",
                                  path.method = "limma", node.method = "limma", fun.method = "limma",
                                  order = FALSE, paired = paired, adjust = adjust, conf.level = 0.05, sel_assay = 1)
  if(hipathia){
    hi_results <- compare_pipeline(hdata, groups=data_set$des$group, expdes=group1, g2 = group2,
                                 path.method = "wilcoxon", node.method = "limma", fun.method = "wilcoxon",
                                 order = FALSE, paired = paired, adjust = adjust, conf.level = 0.05, sel_assay = 1)
  }
  status(" 80", "Differential Activity Analysis completed successfully", output_folder)
}
## (Upcoming Features)Step 4.2: Drug repurposing | MAchine learning | variante interpreter ....

################ Start : To be removed after testing the hipathia function with 3_terms
# ## Step 5: Functional analysis: hipathia::quantify_terms()
# uniprot_vals <- quantify_terms(results, pathways, dbannot = "uniprot") is hipathia:::quantify_funs() and it encapsulated inside hipathia/metabopathia function 
################ End : To be removed after testing the hipathia function with 3_terms

# Save the entire workspace
save.image(file = file.path(output_folder,"workspace.RData"))
## Step 6: Results visualization

## visualization functions
# Results overview: Summary of UP & DOWN nodes, paths and functions
hipathia::DAoverview(DAdata = met_results)
ggplot2::ggsave(file.path(output_folder, "DAoverview_metabopathia.png")) # Results overview
if(hipathia) {
  hipathia::DAoverview(DAdata = hi_results) # Results overview
  ggplot2::ggsave(file.path(output_folder, "DAoverview_hipathia.png")) 
}

# Results summary by pathway: summary of the number of paths altered in the n most altered pathways
hipathia::DAsummary(DAdata = met_results, n = 10) # Top altered pathways
ggplot2::ggsave(file.path(output_folder, "DAsummary_metabopathia.png"))
if(hipathia) {
  hipathia::DAsummary(DAdata = hi_results, n = 10) # Top altered pathways
  ggplot2::ggsave(file.path(output_folder, "DAsummary_hipathia.png")) 
}
##NOTE: if overlay maybe an Enrichment analysis will be relevant here !

# Top results per feature: Top 10 altered features per class (nodes, paths, functions)
# I have to check if there is some altered or not before , othways it will give an error conf.level =0.5 to fore results (has to be removed offcorse)
# hipathia::DAtop(DAdata = met_results, n = 10, conf.level = 05) # top n differtially activated nodes, paths and functions, and plots a dot plot with that info.
# ggplot2::ggsave(file.path(output_folder, "DAtop_metabopathia.png")) # Results overview
# if(hipathia) {
#   hipathia::DAtop(DAdata = hi_results, n = 10, conf.level = 0.5)
#   ggplot2::ggsave(file.path(output_folder, "DAtop_hipathia.png"))
# }

# Save the entire workspace
save.image(file = file.path(output_folder,"workspace.RData"))
servr::daemon_stop() # kill & close
# DAreport() to easily create a report
# Save and serve all results to browser
## Overlay: The report here is different
if(analysis=="overlay"){
  MTreport <- DAreport_overlay(met_results, pathways, path = output_folder, output_folder = "metabopathia_report", verbose = verbose, adjust = adjust, conf.level = conf.level)
}else{
  MTreport <- DAreport(met_results, metabo_pathways, path = output_folder, output_folder = "metabopathia_report", verbose = verbose, adjust = adjust, conf.level = conf.level)
}
status("100", paste0("HTML report created successfully in the ",file.path(codebase,MTreport)), output_folder)
message("Press Ctrl + C to stop serving the report...\n")
serve_report(file.path(codebase,MTreport), port = servr::random_port(), browser = T, daemon = T)
if(hipathia){
  HIreport <- DAreport(hi_results, pathways, path = output_folder, output_folder = "hipathia4comp_report",verbose = verbose)
  message("Press Ctrl + C to stop serving the report...\n")
  serve_report(file.path(codebase,HIreport),port = servr::random_port(), browser = T, daemon = F)
}

# Pause and wait for user input
if(!is_rstudio){
  cat("Press Enter to finish...")
  invisible(readLines("stdin", n=1))
}
servr::daemon_stop() # kill & close


if(is_rstudio){
  ### Pathway differential activation plot: pathway viewer and other figures (interactivity! visNetwork? we compatibility)
  hipathia::DApathway(name = "hsa04015", pathways = metabo_pathways, DAdata = met_results) # Pathway viewer plot 
  hipathia::DApathway(name = "hsa04015", pathways = pathways, DAdata = hi_results) # Pathway viewer plot 
  # DApathway(metabo_pathways$pathigraphs$hsa04720$path.id, metabo_pathways, DAdata)
  # 100*(intersect(metabo_pathways$all.genes , rownames(genes_vals)) %>% length(.) ) / length(metabo_pathways$all.genes)
  ## other version !
  # metabopathia
  colors_de_m <- node_color_per_de(results = metdata, metaginfo = metabo_pathways, 
                                 group = data_set$des$group, expdes = opt$group1, g2 = opt$group2)
  newComp_m <- as.data.frame(met_results$paths)
  rownames(newComp_m)<- newComp_m$ID
  hipathia::pathway_comparison_plot(comp = newComp_m , metaginfo = metabo_pathways, pathway = "hsa04015", 
                                    node_colors = colors_de_m, conf = 0.5)
  
  # hipathia
  colors_de <- node_color_per_de(results = hdata, metaginfo = pathways, 
                                 group = data_set$des$group, expdes = group1, g2 = group2)
  newComp <- as.data.frame(hi_results$paths)
  rownames(newComp)<- newComp$ID
  hipathia::pathway_comparison_plot(comp = newComp , metaginfo = pathways, pathway = "hsa04015", 
                                    node_colors = colors_de, conf = 0.5)
}



