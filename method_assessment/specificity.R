# Load required libraries
library(dplyr)

# Function to simulate individuals for each gene based on normal distribution
met_results_paths <- read.csv(file = "supplementary_files/brca_caseStudy/results/met_results_paths.tsv", sep = "\t")

# Function to assess differential activation (this is a placeholder for actual analysis)
assess_differential_activation <- function(dataset1, dataset2) {
  # Example: here we would perform the actual statistical test (e.g., t-test) to compare datasets
  p_values <- apply(dataset1 - dataset2, 2, function(x) t.test(x)$p.value)
  # Consider p-value < 0.05 as "differentially activated"
  return(sum(p_values < 0.05))
}

# Set simulation parameters
N_values <- c(20, 50, 100, 200, 400)
num_iterations <- 2000
num_genes <- 100  # Example: Number of genes to simulate

# Example empirical distribution parameters (mean and variance for each gene)
set.seed(42)  # For reproducibility
mu <- runif(num_genes, 5, 15)  # Mean gene expression
sigma2 <- runif(num_genes, 0.1, 1)  # Variance in gene expression

# Initialize a list to store results
results <- data.frame(N = rep(N_values, each = num_iterations), false_positives = 0)

# Run simulations
for (N in N_values) {
  for (i in 1:num_iterations) {
    # Generate two datasets of N individuals
    dataset1 <- simulate_individuals(N, mu, sigma2)
    dataset2 <- simulate_individuals(N, mu, sigma2)
    
    # Assess differential activation
    false_positives <- assess_differential_activation(dataset1, dataset2)
    
    # Store the result (if any false positives found)
    results$false_positives[results$N == N & 1:num_iterations == i] <- false_positives
  }
}

# Calculate false positive rate for each N
false_positive_rate <- results %>%
  group_by(N) %>%
  summarise(rate = mean(false_positives > 0))

# Show results
print(false_positive_rate)
