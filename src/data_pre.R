## Data preprocessment
# Metabopathia accepts as input a gene expression matrix and metabolite concentration data. 
# Expression may have been measured using any available sequencing technique, while metabolomics must specify the type of technique used.
# However, Metabopathia assumes that the data has already been normalized to correct for any potential sequencing bias, including batch effect correction.
# The method does not accept missing data or NA values, which is why imputation methods are strongly recommended if needed.
# For the expression matrix, the row names must be the Entrez IDs of the genes in the rows. 
# That's why we use the translate_data() function to transform the gene names.
# For metabolites, their names must be in KEGG IDs. Note that there is no function for this transformation in the current version.
## Data scaling & normalization
# Both data matrices (gene expression and metabolite concentration) have to be scaled between 0 and 1 before computing the subpaths activation values. 
# The normalization has to be performed independently for each data type. 
# In the actual analysis we have used the the percentil to compute the normalized value between 0 and 1, it was chosen in order to remove observed outliers in metabolomics data (heavy-tailed distributions of the metabolomics).
# In order to use same technique for both data , normalization by percentil were used as well for transcriptomic dataset.
# by_gene vs by_metabolite ? by row: yes/no?
data_pre <- function(exp_file, met_file=NULL,met_type, design_file, group1, group2, output_folder, design_type="categorical", analysis="compare", spe="hsa", verbose=FALSE){
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
  
  ### Distribution of received data
  # Save the combined plot to a PNG file
  png(file.path(output_folder, "data_distribution.png"), type="cairo")
  # Set up the layout for the plots
  par(mfrow = c(2, 2))
  boxplot(trans_exp, las = 2, main = "Expression (Before Normalization)",
          xlab = "Genes", ylab = "Expression", col = "lightblue")
  if(met_type == "concentration_matrix")
    boxplot(metabo_data, las = 2, main = "Metabolites (Before Normalization)",
            xlab = "Metabolites", ylab = "Concentration", col = "lightcoral")
  # Metabolite data with outliers
  genes_vals <- normalize_data(as.matrix(trans_exp), percentil = TRUE)
  if(met_type == "concentration_matrix")
    metabo_vals <- normalize_data(as.matrix(metabo_data), percentil = TRUE)
  
  boxplot(genes_vals, las = 2, main = "Expression (After Normalization)",
          xlab = "Genes", ylab = "Expression", col = "blue")
  if(met_type == "concentration_matrix")
    boxplot(metabo_vals, las = 2, main = "Metabolites (After Normalization)",
            xlab = "Metabolites", ylab = "Concentration", col = "red")
  # Close the PNG device
  dev.off()
  par(mfrow = c(1, 1))
  # If all is okay
  data_set$genes_vals <- genes_vals
  if(met_type == "concentration_matrix")
    data_set$metabo_vals <- metabo_vals
  if(met_type == "perturbations")
    data_set$metabo_vals <- metabo_data
  data_set$des <- des 
  return(data_set)
}
