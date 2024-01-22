#!/usr/bin/env /opt/R/4.3.1/bin/Rscript
### !/usr/bin/env Rscript
# Main Script for metabopathia project: 01_main.R
# Author: Kinza Rian
# Description: This script will be the core of Metabopathia web-server with the provided omics.
# Usage: Rscript 01_main.R --species hsa or dirretly ./01_main.R --species hsa

# Load necessary libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(hipathia))

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
  make_option("--difexp", action = "store_true", default = TRUE,
              help = "Boolean, whether to perform differential expression analysis."),
  make_option("--GO.terms", action = "store_true", default = FALSE,
              help = "Boolean, whether to compute functional analysis with Gene Ontology terms."),
  make_option("--uni.terms", action = "store_true", default = FALSE,
              help = "Boolean, whether to compute functional analysis with Uniprot keywords."),
  make_option("--custom.terms", type = "character", default = NULL,
              help = "Path to a file containing a data.frame with the custom annotation of the genes to the functions. First column are gene symbols, second column the functions."),
  make_option("--analysis", type = "character", default = "compare",
              help = "Type of analysis. Allowed values: 'compare', 'predictor_test', 'predictor_train', 'variant_interpreter', 'drug_repurposing'. Default is 'compare'.\n\n\
                      Differential Signaling 'compare': Provides an estimation of significant cell signaling activity changes across different conditions.\n\
                      Predictor 'predictor': Allows you to train a prediction-test and test it with different data.\n\
                      Variant Functional Interpretation 'variant_interpreter': Provides an estimation of the potential impact of genomic variation on cell signaling and cell functionality.\n\
                      Drug Repurposing 'drug_repurposing': drug repurposing."
  ),
  make_option("--output_folder", type = "character", default = "tmp",
              help = "Output folder path. Default is 'tmp'.")
  )
# Parse command-line arguments
opt <- parse_args(OptionParser(option_list = option_list))
# be carefull this is a forced hardcoded contant to be removed !!!
source("src/example1.R")
# Extract values from options
attach(opt) # some recomendations are to avoid the use of attach but still not convincing for me in this situation!
if(verbose){
  message("The recieved options are :")
  str(opt)
}

# Validate options
validate_options(opt)

# Set default values if not specified by the user
if (is.null(opt$pathways_list)) { #if (!opt$pathways_list) { the other one is more acurate!
  pathways_list <- NULL
}

output_folder <- create_output_folder(output_folder, verbose = verbose)
status("  0", "Analysis is started", output_folder)
###### Workflow #####
## Step 1: Data pre-processing
if(verbose) message("Loading data...")
data_set <- data_pre(exp_file, met_file, design_file, group1, group2, output_folder, design_type, analysis, spe, verbose)
status(" 20", "Data preprocessed successfully", output_folder)

## Step 2: MGI preparation & filtering
### from data filter pathways_list ( think about filtering by metabolimocs type in order to avoid high prensentage of missingness)
# ...
### load pathways with pathway list:
if(verbose) message("Loading pathways...") #I have to load from prepared MGI already ! time consuming 
pathways <- hipathia::load_pathways(species = spe, pathways_list = pathways_list)
# modules <- 
## adaptation of the MGI: to be removed, because I have to load it already prepared
metabo_pathways <- add_metabolite_to_mgi(pathways)
status(" 40", "Pathways loaded successfully", output_folder)

## Step 3: Signal propagation : Pathway activation computation
if(verbose) message("Propagating signaling...")
metdata <- metabopathia(genes_vals, metabo_vals, metabo_pathways, uni.terms = uni.terms, GO.terms = GO.terms,
                        decompose = decompose, verbose=verbose)
hdata <- hipathia(genes_vals = genes_vals,pathways, uni.terms = uni.terms, GO.terms = GO.terms,
                        decompose = decompose, verbose=verbose)
### path vals
####metabopathia
met_path_vals <- get_paths_data(metdata, matrix = T)
met_path_vals <- normalize_paths(met_path_vals, metabo_pathways)
###hipathia
h_path_vals <- get_paths_data(hdata, matrix = T)
h_path_vals <- normalize_paths(h_path_vals, pathways)

status(" 60", "Signal propagation computed successfully", output_folder)



## Step 4: Scenarios & pipelines :
## Step 4.1: Differential Activity Analysis
if(analysis=="compare") compare_pipeline()
# hipathia::DAcomp(hidata = , groups = , expdes = , g2 = , )
## (Upcoming Features)Step 4.2: Drug repurposing | MAchine learning | variante interpreter ....
## Step 5: Functional analysis: hipathia::quantify_terms()
## Step 6: Results
### pathway viewer and other figures (interactivity!)




