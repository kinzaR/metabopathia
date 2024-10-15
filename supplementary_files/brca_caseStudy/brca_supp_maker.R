#TODO: Rewrite this script to generat dirctly supplimentaryu matyerials
### Supp material generator :
######################################### Boxplot of inferred metabolite (15) - Tumor Vs Normal ######################################### 
## plot metabo_vals.tsv
library(reshape2)
write.table(x = data_set$metabo_vals, file = "supplementary_files/brca_caseStudy/metabo_vals.tsv", append = F, quote = F, sep =  "\t", col.names = T, row.names = T)
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
  labs(title = "Boxplot of inferred metabolite (16) - Tumor Vs Normal",
       x = "Sample Group",
       y = "Inferred activity") +
  theme_minimal() +
  scale_fill_brewer(palette = "Dark2") + 
  facet_wrap(~ metabolite, nrow = 4)  # Separate plots for each gene
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
############# sub-pathways
circuits <- do.call(what = rbind, lapply( pathways$pathigraphs, function(p){
  tibble(path_id = p$path.id,
         numberOfCircuits = length(p$effector.subgraphs),
         circuit_id = names(p$effector.subgraphs)) %>%
    rowwise() %>% 
    mutate(circuit_name = hipathia::get_path_names(metaginfo = pathways,  circuit_id) )
}))
write.table(x = circuits, file = "supplementary_files/circuitsOf146_pathway.tsv", 
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)


# Load ggplot2
library(ggplot2)

# Example data
data <- circuits %>% select(c(path_id,numberOfCircuits)) %>% unique %>% rowwise() %>% mutate(pathway = pathways$pathigraphs[[path_id]]$path.name)

# Load ggplot2
library(ggplot2)

# Create a ggplot with adjustments for readability
ggplot(data, aes(x = numberOfCircuits, y = reorder(pathway, numberOfCircuits))) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Number of Circuits", y = "Pathways", title = "Number of Circuits per Pathway") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 4),  # Reduce the font size for y-axis text
    plot.margin = margin(1, 1, 1, 3, "cm"),  # Increase left margin for better readability
    axis.title.y = element_text(size = 5),  # Y-axis label size adjustment
    axis.title.x = element_text(size = 10)   # X-axis label size adjustment
  ) +
  coord_fixed(ratio = 1.5)  # Adjusting the aspect ratio
#"supplementary_files/circuits_per_pathways_barplot.png"
####################################################################
######draft for diff node vals 
bol <- met_node_vals == h_node_vals
diff_node_values<-which(rowSums(bol) !=218)
diff_nodes <- met_node_vals[diff_node_values,] %>% rownames()

all_atts_xl <- do.call(what = rbind, lapply(metabo_pathways$pathigraphs, function(g){
  igraph::as_data_frame(x = g$graph,"vertices") %>%  
    select(c(name, shape, label, genesList, metaboID, tooltip))
}))
all_atts <- all_atts_xl %>%  select(c(name, shape, label, genesList, metaboID))
metabo_atts <- all_atts %>% filter(shape =="circle")
all_atts_xl %>% 
  filter(name %in% (metabo_atts %>% filter(is.na(metaboID )) %>% .$name)) %>% # these are without metaboID
  View
diff_nodes_att <- all_atts %>% filter(name %in% diff_nodes)
################# save tables 
write.table(x = met_results$nodes , file = "supplementary_files/brca_caseStudy/results/met_results_nodes.tsv", append = F,quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = met_results$paths , file = "supplementary_files/brca_caseStudy/results/met_results_paths.tsv", append = F,quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = hi_results$nodes , file = "supplementary_files/brca_caseStudy/results/hi_results_nodes.tsv", append = F,quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = hi_results$paths , file = "supplementary_files/brca_caseStudy/results/hi_results_paths.tsv", append = F,quote = F, sep = "\t", row.names = F, col.names = T)
################## clinical data
clinical_data <- read.csv("/home/krian/Downloads/clinical.tsv", sep = "\t")
case_ids <- data_set$des %>% rowwise() %>% mutate(case_id= strsplit(sample, "-") %>% .[[1]] %>% .[1:3] %>% paste(collapse = "-")) %>% .$case_id
clinical_data <- clinical_data %>% filter(case_submitter_id %in% case_ids) %>% select(c(case_submitter_id,age_at_index,gender)) %>% unique()
clinical_data %>% dim()
# Load required libraries
library(ggplot2)
library(dplyr)

# Define menopause threshold (e.g., 50 years)
menopause_threshold <- 50

# Add a column to classify menopausal status
clinical_data <- clinical_data %>%
  mutate(menopausal_status = ifelse(gender =="female" ,ifelse(age_at_index >= menopause_threshold, "2-post-menopausal", "1-pre-menopausal"), NA))
  # Create the ggplot graph
  ggplot(clinical_data %>% filter(gender == "female") %>% arrange(age_at_index) %>%
           mutate(age_group = cut(as.numeric(age_at_index), breaks = seq(0, 100, by = 5), right = FALSE)), 
         aes(y = age_group, fill = menopausal_status)) +
    geom_bar(position = "dodge") +  # geom_bar is appropriate for counting frequencies
    labs(
      title = "Frequency of ages by menopausal status (in 5-year bins)",
      x = "Frequency",
      y = "Age group (5-year bins)",
      fill = "Menopausal status"
    ) + 
    facet_wrap(~ menopausal_status) +  # Facet by menopausal status
    theme_minimal()
write.table(x = clinical_data %>% select(-menopausal_status), file = "supplementary_files/brca_caseStudy/clinical_info.tsv", append = F, quote = F, sep = "\t", row.names = F, col.names = T)
