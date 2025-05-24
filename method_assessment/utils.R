library(dplyr)
########################################### LOAD DATA ###########################################
load_data <- function(exp_file, met_file=NULL,met_type, design_file, group1, group2, output_folder, design_type="categorical", analysis="compare", spe="hsa", verbose=FALSE){
    #TODO: remove output_folder from here ! no need
  data_set <- list()
  ### Load data from files
  # Load expression
  exp <- read.table(exp_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1, comment.char = "", check.names = FALSE)
  # Load metabolite concentrations
  if(!is.na(met_file))
    metabo_data <- read.table(met_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1, comment.char = "", check.names = FALSE)
  # Load design
  des <- read.table(design_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(des) <- c("sample", c("group", "value")[(design_type == "continuous") + 1])
  rownames(des) <- des$sample
  
  # Filter only the specified groups
  if(design_type == "categorical"){
    des <- des %>% dplyr::filter(group %in% c(group1, group2))
  }
  des <- des %>% filter(sample %in% colnames(exp))
  # Subset expression and metabolite data based on the filtered design
  exp <- exp %>% select(des$sample)
  if(met_type == "concentration_matrix") metabo_data <- metabo_data[, des$sample]
  #metabo_data <- metabo_data %>% dplyr::select(des$sample) %>%
  #  dplyr::relocate(des$sample)
  
  # Check for intersection between design and data
  if(ncol(exp) == 0 || (met_type == "concentration_matrix" && ncol(metabo_data) == 0)){
    stop("ERROR: No intersection between sample names in the Expression matrix ",ifelse(!is.na(met_file),"and/or metabolomics","")," data and the Design matrix; please check your input data ")
  }
  # Check for a sufficient number of samples
  if(ncol(exp) < 3 || (met_type == "concentration_matrix" && ncol(metabo_data) < 3)){
    stop("ERROR: Not enough samples in the Expression/Concentration matrix (at least 3 complete pairs); please check your input data ")
  }
  # Check for consistent data intersection
  if(met_type == "concentration_matrix"){
    if(ncol(exp) != nrow(des) && ncol(metabo_data) != nrow(des)){
      stop("ERROR: Data intersection between design and metabolite data and expression data are not the same")
    }
  }
  ### Additional validation
  if(anyNA(exp) || (met_type == "concentration_matrix" && anyNA(metabo_data))){
    stop("ERROR: Expression data and metabolite data cannot have NAs")
  }
  
  ### Translate names to EntrezGene
  # devtools::install_version("dbplyr", version = "2.3.4")
  trans_exp <- hipathia::translate_data(as.matrix(exp), species = spe, verbose = verbose)
  # Assume metabolite names are KEGG IDs
  
  ### Data scaling & normalization
  # Normalize each dataset independently
  if (max(trans_exp) == min(trans_exp) || (met_type == "concentration_matrix" && max(metabo_data) == min(metabo_data))) {
    stop("ERROR: The min is equal to max in either Expression or Metabolite data !")
  }
  # If all is okay
  data_set$genes_vals <- trans_exp
  if(met_type == "concentration_matrix")
    stop("NOT Implemented yet!")
  if(met_type == "perturbations")
    data_set$metabo_vals <- metabo_data
  data_set$des <- des 
  return(data_set)
}
########################################### get_N_data_perSenario ###########################################
get_N_data_perSenario <- function(data_set, N, senario="i"){
    data_set <- get_N_data(data_set, N)
    if(senario=="i") {
        # Scenario i: N samples from original data with identical groups
        data_set$genes_vals <- hipathia::normalize_data(as.matrix(data_set$genes_vals), by_quantiles = FALSE, 
                                 by_gene = FALSE, percentil = FALSE)
    }
    if(senario=="ii") {
        # Scenario ii- genes: Simulate gene expression from empirical distribution
        data_set$genes_vals <- apply(data_set$genes_vals, 1, function(gene_values) {
          rnorm(N, mean = mean(gene_values), sd = sd(gene_values))
        }) %>% t()
        colnames(data_set$genes_vals) <- data_set$des$sample
        # Values need to be between 0 and 1
        data_set$genes_vals <- hipathia::normalize_data(as.matrix(data_set$genes_vals), by_quantiles = FALSE, by_gene = FALSE, percentil = FALSE)
        # Scenario ii - metabolites: Simulate Metabolites activities from empirical distribution
    	data_set$metabo_vals <- t(apply(data_set$metabo_vals, 1, function(metabolite_vals) {
    	  # Step 1: Simulate new values using rnorm
    	  simulated_vals <- rnorm(N, mean = mean(metabolite_vals), sd = sd(metabolite_vals))
    	  # Step 2: Perform min-max scaling
            if(max(simulated_vals) > min(simulated_vals)){
        	  scaled_vals <- (simulated_vals - min(simulated_vals)) / (max(simulated_vals) - min(simulated_vals))  # Scale to 0-1
        	  return(scaled_vals * (max(metabolite_vals) - min(metabolite_vals)) + min(metabolite_vals))  # Scale to original min-max range
            }else return(simulated_vals)
    	}))
    	# Assign column names
    	colnames(data_set$metabo_vals) <- data_set$des$sample
    }
    if(senario=="iii") {
        # Scenario iii: Simulate from beta distribution with mean 0.5 and sd ~0.05
        data_set$genes_vals <- apply(data_set$genes_vals, 1, function(gene_values) {
          rnorm(N, mean = 0.5, sd = 0.05)
        }) %>% t()
        colnames(data_set$genes_vals) <- data_set$des$sample
        # Values need to be between 0 and 1
        data_set$genes_vals <- hipathia::normalize_data(as.matrix(data_set$genes_vals), by_quantiles = FALSE, by_gene = FALSE, percentil = FALSE)
        # Metabolite
        data_set$metabo_vals <- apply(data_set$metabo_vals, 1, function(metabolite_vals) {
          rnorm(N, mean = 0.5, sd = 0.05)
        }) %>% t()
        colnames(data_set$metabo_vals) <- data_set$des$sample         
    }

    return(data_set)
}
########################################### DAsummary ###########################################
DAoverview_plotless <- function (DAdata, conf.level = 0.05, adjust = TRUE) 
{
  summ <- lapply(names(DAdata), function(feat) {
    data <- DAdata[[feat]]
    if (adjust == TRUE) {
      summdf <- data.frame(feature = feat, total = nrow(data), 
                           sigs = sum(data$FDRp.value < conf.level), UPs = sum(data$FDRp.value < 
                                                                                 conf.level & data$statistic > 0), DOWNs = sum(data$FDRp.value < 
                                                                                                                                 conf.level & data$statistic < 0))
    }
    else {
      summdf <- data.frame(feature = feat, total = nrow(data), 
                           sigs = sum(data$p.value < conf.level), UPs = sum(data$p.value < 
                                                                              conf.level & data$statistic > 0), DOWNs = sum(data$p.value < 
                                                                                                                              conf.level & data$statistic < 0))
    }
  })
  summ <- tibble(do.call(rbind, summ))
  return(summ)
}

########################################### Scenario i ###########################################
# Scenario i 
# Example usage:
# Assuming 'samples' is a vector of BRCA sample names
# samples <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")
# get_N_from_data(samples, 3)  # Randomly select 3 samples from 'samples'
get_N_from_data <- function(samples, N) {
  # Check if N is less than or equal to the number of available samples
  if (N > length(samples)) {
    stop("N is larger than the number of available samples.")
  }
  # Randomly sample N individuals without replacement
  N_random_sample <- sample(samples, N, replace = FALSE) # no samples are repeated
  
  return(N_random_sample)
}
# return same list of data with selected N randomly and identicaly
get_N_data <- function(data_set, N){
  samples <- data_set$des %>% filter(group == sample(unique(group),1))
  new_data_set <- purrr::map(data_set, as.data.frame)
  new_data_set$des <- new_data_set$des %>% filter(sample %in%get_N_from_data(samples$sample,N)) %>%
                                           mutate(group = as.factor(c(rep("group1",length(group)/2), rep("group2",length(group)/2))))
  new_data_set$genes_vals <- new_data_set$genes_vals %>% select(c(new_data_set$des$sample))
  new_data_set$metabo_vals <- new_data_set$metabo_vals %>% select(c(new_data_set$des$sample))
  return(new_data_set)
}
# Scenario ii
# Function to simulate gene expression data for N individuals in Scenario ii
simulate_gene_expression_empirical <- function(N, mu_g, sigma2_g) {
  mapply(function(mu, sigma2) {
    rnorm(N, mean = mu, sd = sqrt(sigma2))  # Normal distribution N(μg, σ2g)
  }, mu_g, sigma2_g)
}