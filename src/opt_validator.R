# an option validator
# Function to validate command-line options
############### Config ###############
# Define constants
allowed_species <- c("hsa", "mmu", "rno")
required_files <- c("exp_file", "met_file", "design_file")
required_options <- c("group1", "group2")
allowed_design_types <- c("categorical", "continuous")
allowed_analysis <- c("overlay", "compare", "predictor_test", "predictor_train", "variant_interpreter", "drug_repurposing")
ready_analisis <- c("overlay", "compare")
# to be tested in the future:
## pathways_list from available pathways
## custom.terms if it is checked as a file
############### Validator Functions ############### 
# Function to check if required options are provided
check_required <- function(required_options, opt) {
  for (opt_name in required_options) {
    if (is.null(opt[[opt_name]])) {
      stop(paste("Error: '", opt_name, "' is missing.", sep = ""))
    }
  }
}
# Function to validate species: # Check if the provided species is in the allowed choices
validate_species <- function(spe, allowed_species) {
  if (!spe %in% allowed_species) {
    stop("Invalid species specified. Please choose from 'hsa', 'mmu', or 'rno'.")
  }
}

# Function to validate file existence and non-emptiness
validate_file <- function(file_name, opt) {
  if (is.null(opt[[file_name]]) || !file.exists(opt[[file_name]]) || file.size(opt[[file_name]]) == 0) {
    stop(paste("Error: '", file_name, "' is missing, does not exist, or is empty.", sep = ""))
  }
}

# Function to validate design type
validate_design_type <- function(design_type, allowed_design_types) {
  if (!design_type %in% allowed_design_types) {
    stop("Invalid design_type specified. Allowed values: ",paste(allowed_design_types, collapse = ", "),".")
  }
}

# Function to validate chosen analysis pipeline
validate_analysis <- function(analysis, allowed_analysis, ready_analisis) {
  if (!analysis %in% allowed_analysis) {
    stop("Invalid analysis specified. Allowed values: ",paste(allowed_analysis, collapse = ", "),".")
  } else if (!analysis %in% ready_analisis) {
    stop("Analysis not ready for execution. Ready implemented piplines are: ",paste(ready_analisis, collapse = ", "),".")
  }
}

############### Main validation function ############### 
validate_options <- function(opt, verbose = FALSE) {
  # Check if the required options are provided
  check_required(required_options, opt)
  # Validate species
  validate_species(opt$spe, allowed_species)
  # Check if the required files exist and are non-empty
  for (file_name in required_files) {
    validate_file(file_name, opt)
  }
  # Validate design type
  validate_design_type(opt$design_type, allowed_design_types)
  
  # Validate custom annotation file
  if(!is.na(opt$custom.terms)){
    validate_file("custom.terms", opt)
  }
  
  # Validate chosen analysis pipeline
  validate_analysis(opt$analysis, allowed_analysis, ready_analisis)
  
  if (verbose) message("All options are ok")
}
