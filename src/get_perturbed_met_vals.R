# group2 is the ref
# is only KO for the moment 
get_perturbed_met_vals <- function(metabo_vals, des, group1, all_metabolites, verbose){
  all_met_df <- as.data.frame(matrix(1, nrow = length(all_metabolites), ncol = length(des$sample),
                          dimnames = list(all_metabolites, des$sample)))
  all_met_df[metabo_vals,]<- ifelse(des$group == group1, 0,1)
  return(all_met_df)
}