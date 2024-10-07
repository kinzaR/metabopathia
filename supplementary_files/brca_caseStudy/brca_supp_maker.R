### Supp material generator :
## plot metabo_vals.tsv
library(reshape2)
activity_matrix <- read.csv( file = "supplementary_files/brca_caseStudy/metabo_vals.tsv", sep = "\t")
activity_matrix$met <- rownames(activity_matrix)
activity_df <- melt(activity_matrix, id.vars = "met")
colnames(activity_df) <- c("metabolite", "sample", "infActivity")
# Merge with sample information
activity_df$group <- data_set$des$group # reded from here opt$design_file <-"data_examples/TCGA/processed_data/brca_109NX109T_paired/BRCA_109Nx109T_paired_des_v2.tsv"
library(ggplot2)
library(gridExtra)

# Create a function to plot distributions for each metabolite
plot_metabolite_distribution <- function(data, metabolite) {
  ggplot(data[data$metabolite == metabolite, ], aes(x = infActivity, fill = sample)) +
    geom_density(alpha = 0.5) +
    labs(title = metabolite, x = "Activity Level", y = "Density") +
    theme_minimal() +
    theme(legend.position = "top")
}

# Generate a list of plots for each metabolite
metabolites <- unique(activity_df$metabolite)
plots <- lapply(metabolites, plot_metabolite_distribution, data = activity_df)

# Arrange plots in a grid (5 per row)
grid.arrange(grobs = plots, ncol = 5)
