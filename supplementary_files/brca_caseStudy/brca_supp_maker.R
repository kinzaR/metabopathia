#TODO: Rewrite this script to generat dirctly supplimentaryu matyerials
### Supp material generator :
######################################### Boxplot of inferred metabolite (15) - Tumor Vs Normal ######################################### 
## plot metabo_vals.tsv
library(reshape2)
activity_matrix <- read.csv( file = "supplementary_files/brca_caseStudy/metabo_vals.tsv", sep = "\t")
metabo_annotation <- read.csv(file = "supplementary_files/brca_caseStudy/inferred_metabolite_brca_metabolizer_v2.csv")
activity_matrix$met <- rownames(activity_matrix)
activity_df <- melt(activity_matrix, id.vars = "met")
colnames(activity_df) <- c("metabolite", "sample", "infActivity")
activity_df <-activity_df %>% left_join(y = metabo_annotation, by = join_by("metabolite"=="Query")) %>% 
  mutate(metabolite= paste0(Match,"(",metabolite,")")) %>% select(c("metabolite", "sample", "infActivity"))
# Merge with sample information
activity_df$group <- data_set$des$group # reded from here opt$design_file <-"data_examples/TCGA/processed_data/brca_109NX109T_paired/BRCA_109Nx109T_paired_des_v2.tsv"
library(ggplot2)
library(gridExtra)


activity_df<- activity_df %>% group_by(group)

# Create a boxplot for the same selected genes divided by sample type after truncation
ggplot(activity_df ,
       aes(x = group, y = infActivity, fill = group)) +
  geom_boxplot(outlier.shape = NA) +  # No outliers displayed in the plot
  labs(title = "Boxplot of inferred metabolite (15) - Tumor Vs Normal",
       x = "Sample Group",
       y = "Inferred activity") +
  theme_minimal() +
  scale_fill_brewer(palette = "Dark2") + 
  facet_wrap(~ metabolite, nrow = 3)  # Separate plots for each gene
#######################################################################################################################################
# Difirentia activity motabolites 
######################################################################################################################################
library(tibble)
library(KEGGREST)
get_pathway_info <- function(path_id) {
  # Retrieve pathway information
  path_info <- KEGGREST::keggGet(path_id)[[1]]
  
  # Return tibble with relevant information
  tibble(
    path_id = path_id,
    name = path_info$NAME ,
    class = if(!is.null(path_info$CLASS)) path_info$CLASS else NA,
    description = if (!is.null(path_info$DESCRIPTION)) paste(path_info$DESCRIPTION, collapse = "- ") else NA,
    compounds = if (!is.null(path_info$COMPOUND)) paste(path_info$COMPOUND, collapse = ", ") else NA
  )
}
pathways_info_table <- lapply(names(metabo_pathways$pathigraphs),get_pathway_info) %>% do.call(what = rbind, args =.)
pathways_info_table <- pathways_info_table %>% rowwise %>% mutate(shortName= metabo_pathways$pathigraphs[[path_id]]$path.name,
                                           numberOfNode= igraph::V(graph = metabo_pathways$pathigraphs[[path_id]]$graph)  %>% length(),
                                           numberOfMetabolites = igraph::vertex_attr(graph = metabo_pathways$pathigraphs[[path_id]]$graph, "shape") %>% .[.=="circle"] %>%length(),
                                           annotatedMetabolite = igraph::vertex_attr(graph = metabo_pathways$pathigraphs[[path_id]]$graph, "metaboID") %>% .[!is.na(.)]  %>%length(),
                                           metaboliteList= igraph::vertex_attr(graph = metabo_pathways$pathigraphs[[path_id]]$graph, "metaboID") %>% .[!is.na(.)]  %>% paste(., collapse = ","))
write.table(pathways_info_table, file = "supplementary_files/pathways_information.tsv", append = F, quote = F, sep = "\t")
#### Metabo
module_paths <- 
  lapply(hsa_module_data, function(m){
    tibble(
      module = m$graphobject$module_name,
      path_id = m$graphobject$modulePATH$V1,
      path_name = m$graphobject$modulePATH$V2
    )
}) %>% do.call(args = ., what = rbind) %>% unique
pathways_module_info_table <- lapply(module_paths$path_id,get_pathway_info)%>% do.call(what = rbind, args =.)
pathways_module_info_table<- pathways_module_info_table %>% unique()
module_paths_extended_info <- module_paths %>% left_join(., pathways_module_info_table, by = join_by("path_id"=="path_id"))
write.table(module_paths_extended_info, file = "supplementary_files/module_paths_extended_info.tsv", append = F, quote = F, sep = "\t")

summary_tibble <-module_paths_extended_info %>%
  group_by(class) %>%
  summarise(
    number_of_pathways = n_distinct(path_id),
    number_of_modules = n_distinct(module)
  )
## number of classes per module
module_paths_extended_info %>%  select(c(module, class)) %>% unique %>%
  group_by(module) %>% summarise(number_of_classes = n_distinct(class)) %>% arrange(desc(number_of_classes)) %>%
  filter(number_of_classes>=2)
