# utils.R
# 
# Description: This file contains utility functions for the project 'Metabopathia' 
# to streamline data analysis, package management, and other common tasks.
#
# Author: Kinza Rian
# Email: rian.kinza@gmail.com
# 
# Last Modified: 29/09/2024
#
# License:  GPL-3.0
#
# Usage: Source this file in your main scripts to access helper functions.
######################################### Functions #########################################
######################################### check_and_install #################################
# Function to load needed packages. If not installed, it installs them for you. 
# If you are working in a conda environment, renv will be used to ensure that (the installed packages will be available for conda envirnmt the next time).
# Example of how to call the function:
# check_and_install(c("ggplot2", "dplyr", "tibble"), c("SummarizedExperiment", "SEtools", "edgeR"))
#############################################################################################
check_and_install <- function(pkgs, bioc_pkgs = c(), condaEnv = FALSE) { # a void function!
  # Install CRAN packages if not already installed
  for (pkg in pkgs) {
    if (!suppressPackageStartupMessages(require(pkg, character.only = TRUE, quietly = TRUE))) {
        if(condaEnv){
            renv::install(pkg)# For use in a conda environment, use renv::install("pkg") instead.
        }else{
            install.packages(pkg, dependencies = TRUE)
        }        
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }
  # Install Bioconductor packages if not already installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    ifelse(condaEnv, renv::install("BiocManager"), install.packages("BiocManager"))
  }
  for (pkg in bioc_pkgs) {
    if (!suppressPackageStartupMessages(require(pkg, character.only = TRUE, quietly = TRUE))) {
      if(condaEnv){
          renv::install(bioc::pkg)
        }else{
          BiocManager::install(pkg) # For use in a conda environment, use renv::install("bioc::pkg") instead.
          }
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }
}
#############################################################################################
######################################### load_data #########################################
# Function to load RData files from a specified directory for a given cancer type.
# The function constructs the file path using the provided 'cancer_id' and 'path2dir',
# checks if the file exists, and if it does, loads the RData file and returns the 'data' object.
# If the file does not exist, it throws an error message indicating the missing file.
#
# Args:
#   cancer_id: A string representing the cancer type (used to construct the file name).
#   path2dir: A string representing the path to the main directory where RData files are stored.
#
# Returns:
#   The 'data' object from the loaded RData file, if it exists.
#
# Raises:
#   An error if the file is not found in the specified location.
#
# Example:
#   data <- load_data("TCGA-BRCA", "processed_data/")
#############################################################################################
load_data <- function(cancer_id, path2dir) {
  file_path <- file.path(path2dir, cancer_id, paste0("counts_", cancer_id, ".RData"))
  if (file.exists(file_path)) {
   load(file_path)
    return(data)
  } else {
    stop(paste("File not found for cancer:", cancer_id))
  }
}
#############################################################################################
############################## summarize_preservation_methods ###############################
# Function to summarize the preservation methods used across different cancer datasets.
# For each dataset in the provided data_list, the function retrieves the 'sample_type' and 
# 'preservation_method', then groups and summarizes the data, including the 'cancer_id'.
#
# Args:
#   data_list: A list of SummarizedExperiment objects loaded from multiple cancer datasets.
#
# Returns:
#   A tibble summarizing the number of samples ('n') grouped by 'sample_type' and 
#   'preservation_method' for each 'cancer_id'.
#
# Example:
#   preservation_summary <- summarize_preservation_methods(data_list)
#############################################################################################
summarize_preservation_methods <- function(data_list) {
  library(dplyr)
  library(tibble)
  summary_table <- lapply(names(data_list), function(i) {
    data <- data_list[[i]]
    cancer_id <- i
    
    # Extract and summarize preservation method and sample type
    colData(data)[, c('sample_type', 'preservation_method')] %>% 
      as_tibble() %>%
      group_by(sample_type, preservation_method) %>%
      summarise(n = n()) %>%
      mutate(cancer_id = cancer_id) %>% 
      arrange(preservation_method)
  }) %>% bind_rows()
  
  return(summary_table)
}
#############################################################################################
####################### plot_preservation_methods_per_sample_types ##########################
change_label <- function(plot) {
  
  plot$plot_env$Pielabel <- 
    plot$plot_env$data2$label %>% 
    stringr::str_replace_all("<br>", "\n") %>% 
    stringr::str_replace("\\(", " \\(")
    plot$plot_env$label2 <- 
    plot$plot_env$dat1$label %>% 
    stringr::str_replace_all("<br>", "\n") %>% 
    stringr::str_replace("\\(", " \\(") 
  plot
}
# Identify unique preservation methods per sample type
plot_Preservation_Method_Distribution<- function(preservation_summary){
    unique_methods <- preservation_summary %>%
        group_by(preservation_method) %>%
        filter(n_distinct(sample_type) == 1) %>%
        ungroup()
    ggplot(preservation_summary, aes(x = sample_type, y = n, fill = preservation_method)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~ cancer_id) +
      labs(title = "Sample Type and Preservation Method Distribution Per Cancer",
           x = "Sample Type", 
           y = "Number of Samples",
           fill = "Preservation Method") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      # Annotate the unique methods
      geom_text(data = unique_methods, aes(label = preservation_method), 
                position = position_dodge(width = 0.9), vjust = -0.5, color="red")
}
boxPlots_preservation_summary<- function(preservation_summary){
    ggplot(preservation_summary, aes(x = cancer_id, y = n, fill = preservation_method)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(title = "Total Sample Counts by Cancer and Preservation Method",
           x = "Cancer ID", 
           y = "Number of Samples",
           fill = "Preservation Method") +
      theme_minimal()
}
plot_preservation_methods_per_sample_types <- function(preservation_summary){
    plots <- purrr::map(unique(preservation_summary$cancer_id),
      ~ preservation_summary %>% 
        dplyr::filter( cancer_id == .x) %>% ggiraphExtra::ggPieDonut(
            ggplot2::aes(pies = sample_type , donuts = preservation_method), 
          interactive = FALSE
        ) + 
        ggplot2::labs(title = .x)
    ) %>% 
    purrr::map(change_label)
    return(plots)
}
plots_preservation_methods <- function(preservation_summary){
    library(ggpubr)
    p1 <-plot_Preservation_Method_Distribution(preservation_summary)
    p2 <-boxPlots_preservation_summary(preservation_summary)
    ggpubr::ggarrange(p1, p2, ncol = 2)
}
arrager <- function(plots2,i,j,k,l){
    ggpubr::ggarrange(plots2[[i]],plots2[[j]],plots2[[k]],plots2[[l]], ncol = 4)
}
piPlots_preservation_methods <- function(preservation_summary,p){
    library(ggpubr)
    plots2 <- plot_preservation_methods_per_sample_types(preservation_summary)
    #ggpubr::ggarrange(plotlist = plots2 ,ncol=4)
    arrager(plots2,1+p,2+p,3+p,4+p)
    }
#############################################################################################
################################### check_and_save_RDSs #####################################
check_and_save_RDSs <- function(data, metadata, my_dir, version="v1", 
                                prefix_countdata_file= "unstranded_counts_data_list_", 
                                prefix_metadata_file ="counts_metadata_list_"){
    # Define version and file paths
    countdata_file <- file.path(my_dir, paste0(prefix_countdata_file, version, ".rds"))
    metadata_file <- file.path(my_dir, paste0(prefix_metadata_file, version, ".rds"))
    
    # Check if the RDS files already exist before saving
    if (!file.exists(countdata_file)) {
      saveRDS(object = data, file = countdata_file)
      message("Saved: ", countdata_file)
    } else {
      message("File already exists: ", countdata_file)
    }
    if (!file.exists(metadata_file)) {
      saveRDS(object = metadata, file = metadata_file)
      message("Saved: ", metadata_file)
    } else {
      message("File already exists: ", metadata_file)
    }
}
#############################################################################################
####################################### load_RDS_data #######################################
# Example usage:
#countdata_list <- load_RDS_data(name="countdata_list",my_dir, version= "v1", prefix_file= "unstranded_counts_data_list_" )
#metadata_list <- load_RDS_data(name="metadata_list",my_dir, version= "v1", prefix_file= "counts_metadata_list_" )
load_RDS_data <- function(name, dir_path, version, prefix_file) {
  # Construct the file path using the directory, prefix, and version
  file_path <- file.path(dir_path, paste0(prefix_file, version, ".rds"))
  
  # Check if the variable is already loaded in the environment
  if (!exists(name, envir = .GlobalEnv)) {
    # If not loaded, read the RDS file and return the value
    data <- readRDS(file_path)
    message(paste(name, "loaded from", file_path))
    return(data)
  } else {
    message(paste(name, "is already loaded in the environment."))
    # Return the existing variable
    return(get(name, envir = .GlobalEnv))
  }
}
#############################################################################################
#############################################################################################
plot_cases_vs_samples <- function(metadata_list){
    samples_per_cancer<-sapply(metadata_list,dim) %>% .[1,]
    cases_per_cancer<-sapply(metadata_list,function(md){
        md %>% rowwise() %>% mutate(participant = unlist(strsplit(barcode, '-'))[3])%>% dplyr::select(participant) %>% unique %>% dim %>% .[1] })
    cancer_S_C<- rbind(cbind(cases_per_cancer,"cases"),cbind(samples_per_cancer,"samples")) %>% as.data.table(keep.rownames = T)
    colnames(cancer_S_C) <- c("cancer_type", "count","Count_type")
    ggplot(cancer_S_C, aes(cancer_type,count, fill=Count_type)) +
      geom_bar(stat = "identity", position = "dodge")  +
    geom_text(aes(label = count), 
                position = position_dodge(width = 0.9), 
                vjust = -0.3, size = 3) + 
      labs(title = "Number of Samples and Cases per Cancer Type", 
           x = "Cancer Type", 
           y = "Count") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
plot_per_cancer_cases_vs_samples <- function(meta_data){
    meta_data <- meta_data %>% mutate(participant_id = sapply(strsplit(as.character(barcode), "-"), function(x) x[3]))
    # Group by sample_type and summarize the number of samples and unique participants
    summary_data <- meta_data %>%
      group_by(sample_type) %>%
      summarise(
        num_samples = n(),
        num_participants = n_distinct(participant_id)
      )
    # Print summary data
    print(summary_data)
    # Reshape data for ggplot (long format)
    summary_data_long <- summary_data %>%
      pivot_longer(cols = c(num_samples, num_participants), 
                   names_to = "type", 
                   values_to = "count")
    ggplot(summary_data_long, aes(x = sample_type, y = count, fill = type)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.6) +
      scale_fill_manual(values = c("num_samples" = "skyblue", "num_participants" = "orange"),
                        labels = c("Samples", "Participants")) +
      labs(x = "Sample Type", y = "Count", title = "Number of Samples and Participants by Sample Type", fill = "Legend") +
      theme_minimal()
}
#############################################################################################
###################################### get_info_counts ######################################
get_info_counts <- function(metadata_list){
    sapply(metadata_list,function(md){
        md %>%rowwise() %>% mutate(part = unlist(strsplit(barcode, '-'))[3])   %>% 
        mutate(group=ifelse(sample_type=="Solid Tissue Normal","N","T")) %>% dplyr::select(part,group) %>% unique %>%
        .[duplicated(.$part),]%>% arrange(part,group) %>% dim %>%.[1] %>%
        as.data.frame %>% 
        mutate(perc_Paired_samples_in_Normals = ./dim(md %>% filter(sample_type=="Solid Tissue Normal"))[1]*100)%>% 
        mutate(perc_Paired_samples_in_Tumors = ./dim(md %>% filter(sample_type!="Solid Tissue Normal")  %>% rowwise() %>% mutate(participant = unlist(strsplit(barcode, '-'))[3])%>% dplyr::select(participant) %>% unique )[1]*100) %>% # some paqired patient have more than tumor sample
        mutate(Tot_part_Tumors = dim(md %>% filter(sample_type!="Solid Tissue Normal")%>%
                         rowwise() %>% mutate(participant = unlist(strsplit(barcode, '-'))[3])%>%
                         dplyr::select(participant) %>% unique)[1]) %>%
        mutate(Tot_Samples_Tumors = dim(md %>% filter(sample_type!="Solid Tissue Normal")%>%
                         dplyr::select(barcode))[1]) %>%
        mutate(Tot_part_Normals = dim(md %>% filter(sample_type=="Solid Tissue Normal")%>%
                         rowwise() %>% mutate(participant = unlist(strsplit(barcode, '-'))[3])%>%
                         dplyr::select(participant) %>% unique)[1]) %>%
        mutate(Tot_Samples_Normals = dim(md %>% filter(sample_type=="Solid Tissue Normal")%>%
                         dplyr::select(barcode))[1]) %>%
        mutate(dup_participat_in_Normals = dim(md %>% filter(sample_type=="Solid Tissue Normal")%>%
                         rowwise() %>% mutate(participant = unlist(strsplit(barcode, '-'))[3])%>%
                         dplyr::select(participant) %>% .[duplicated(.$participant),] )[1]) %>%
        mutate(dup_participat_in_Tumors = dim(md %>% filter(sample_type!="Solid Tissue Normal")%>%
                         rowwise() %>% mutate(participant = unlist(strsplit(barcode, '-'))[3])%>%
                         dplyr::select(participant) %>% .[duplicated(.$participant),] )[1]) %>%
        
        dplyr::select(-'.')
    }) %>% t %>% as.data.frame
}
#############################################################################################
############################ plot_case_vs_samples_per_tissueType ############################
plot_case_vs_samples_per_tissueType <- function(metadata_info){
    # Reshape data for plotting
    data_long <- metadata_info %>% dplyr::select(c(Tot_part_Tumors,Tot_Samples_Tumors,Tot_part_Normals,Tot_Samples_Normals))%>%
        rownames_to_column("cancer_id") %>% pivot_longer(cols = c(Tot_part_Tumors,Tot_Samples_Tumors,Tot_part_Normals,Tot_Samples_Normals),
                              names_to = "Sample_Type", values_to = "Count") %>% mutate(Count=as.double(Count)) # %>% arrange(cancer_id, Sample_Type)
    # Ensure the Sample_Type is a factor with the specified order
    data_long$Sample_Type <- factor(data_long$Sample_Type, 
                                    levels = c("Tot_part_Tumors", "Tot_Samples_Tumors", 
                                               "Tot_part_Normals", "Tot_Samples_Normals"))
        # Plot the data
    ggplot(data_long, aes(x = Sample_Type, y = Count, fill = Sample_Type)) +
      geom_bar(stat = "identity", position = "dodge") +
      # Add text labels showing the count on top of each bar
      geom_text(aes(label = Count), 
                position = position_dodge(width = 0.9), 
                vjust = -0.3, size = 3) +  # Adjust vjust to control the vertical position of labels
      labs(title = "Total Participants vs. Total Samples in Tumor and Normal Data",
           x = "", y = "Count") +
      facet_wrap(~ cancer_id, ncol = 4) +  # This creates a separate plot for each cancer type
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels for better readability
    scale_fill_manual(
        values = c("Tot_part_Tumors" = "red", "Tot_Samples_Tumors" = "darkred",
                   "Tot_part_Normals" = "blue", "Tot_Samples_Normals" = "darkblue"),
        labels = c("Total Participants (Tumors)", "Total Samples (Tumors)", 
                   "Total Participants (Normals)", "Total Samples (Normals)")
      )
}

#############################################################################################
###################################### get_object_size ######################################
get_object_size <- function(){
    sapply(ls(envir = .GlobalEnv), function(x) object.size(get(x))%>%  format( units = "Mb")) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Object") %>%
      dplyr::rename(Size = ".") %>% filter(as.numeric(gsub(" Mb", "", Size))>0)  %>%
      arrange(desc(as.numeric(gsub(" Mb", "", Size)))) 
}
#############################################################################################
######################################## get_batchs #########################################
get_batchs <- function(metadata_list, cancer_id = NULL){
    #Define batches from barcode
    if (is.data.table(metadata_list)){
        if(!is.null(cancer_id)) metadata_list$project_id <- cancer_id
        batchs <- metadata_list
    }else{
        cancer_list_id <- names(metadata_list)
        batchs<-do.call(rbind, metadata_list)
            ####
        # COLORS: some palettes:
        cancer_safe_colorblind_palette12 <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                                     "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888") %>% setNames(.,cancer_list_id)
    }
    batchs[, c("project", "tss", "participant", "sampleVial", "portionAnalyte", "plate", "center") := tstrsplit(barcode, "-", fixed = TRUE)]

    # by plate and by tss: The number of plates is large, making it challenging to create color-blind friendly visuals. A solution needs to be found.
    batchs[, (names(batchs)[-1]) := lapply(.SD, as.factor), .SDcols = names(batchs)[-1]]
    cols_plate<-setNames(rainbow(length(levels(batchs$plate))),levels(batchs$plate)) 
    cols_tss<-setNames(rainbow(length(levels(batchs$tss))),levels(batchs$tss))
    # Colors for sample types:
    ##by center -> HAs to be skipped because theres only one center !
    batchs <- batchs  %>% mutate(groups=case_when(sample_type == "Solid Tissue Normal" ~ "Normal",.default =  "Tumor")) %>%
        mutate(c_by_type= case_when(sample_type == "Solid Tissue Normal" ~ "#009E73", 
                            sample_type == "Primary Solid Tumor" | sample_type == "Primary Tumor" ~ "#E69F00",
                            sample_type == "Recurrent Solid Tumor" | sample_type == "Recurrent Tumor" ~ "#D55E00",
                            sample_type == "Additional - New Primary" ~ "#0072B2",
                            sample_type == "Metastatic" ~ "#CC79A7")) %>%
        mutate(c_by_group = case_when(sample_type == "Solid Tissue Normal" ~ "#0000FF" ,
                                                .default = "#FF0000")) %>%#Then  colors by cancer type or project id:
        mutate(c_by_plate = cols_plate[batchs$plate]) %>%
        mutate(c_by_tss = cols_tss[batchs$tss])
    if(exists("cancer_safe_colorblind_palette12")) 
        batchs %>%  mutate(c_by_cancer = cancer_safe_colorblind_palette12[project_id])
    # Convert all columns except the first one to factors
    batchs[, (names(batchs)[-1]) := lapply(.SD, as.factor), .SDcols = names(batchs)[-1]]
    #scales::show_col(safe_colorblind_palette)
    rownames(batchs) <- batchs$barcode
    return(batchs)
}
#############################################################################################
#################################### plot_Samples_Cancer ####################################
plot_Samples_Cancer <- function(batchs, cancer_names = NULL){
    by_type_n <-batchs %>% group_by(project_id, sample_type, c_by_type) %>% summarise(count = n()) %>% 
        mutate(cancer_full_name = cancer_names[project_id])
    p<-ggplot(by_type_n, aes(y=cancer_full_name, x=count, fill=sample_type)) +
    geom_bar(stat='identity', position='dodge') + 
    geom_text(aes(label=count), 
              position=position_dodge(width=1), 
              hjust=-0.25, vjust=0.5, size=4) +  # Adjust vjust and size as needed
    scale_fill_manual(values = setNames(by_type_n$c_by_type, by_type_n$sample_type))+
    theme(axis.text.y = element_text(size=14),  # Larger y-axis labels
          axis.text.x = element_text(hjust=1, vjust=0.5, size=14),  # Larger x-axis labels
          axis.title.x = element_text(size=18, margin = margin(t = 10)),  # Larger x-axis title
          axis.title.y = element_text(size=18, margin = margin(r = 10)),  # Larger y-axis title
          plot.title = element_text(size=20, face="bold", hjust = 0.5),  # Larger, centered plot title
          legend.title = element_text(size=16),  # Larger legend title
          legend.text = element_text(size=14),  # Larger legend text
          plot.margin = margin(10, 10, 10, 10)) +
    ggtitle("Number and Type of Samples by Cancer Type") +
    labs(fill = "Sample Types", y="Cancer Types", x="Number of Samples")
    return(p)
}
#############################################################################################
################################### plot_by_tumorVsNormal ###################################
plot_by_tumorVsNormal <- function(batchs){
    by_tumorVsNormal_n <- batchs %>% group_by(project_id, groups, c_by_group) %>% summarise(count = n()) %>% 
        mutate(cancer_full_name = cancer_list_names[project_id])
    p<-ggplot(by_tumorVsNormal_n, aes(y=cancer_full_name, x=count, fill=groups)) +
    geom_bar(stat='identity', position='dodge') + 
    geom_text(aes(label=count), 
              position=position_dodge(width=1), 
              hjust=-0.25, vjust=0.5, size=4) +  # Adjust vjust and size as needed 
    #scale_fill_manual(values = setNames(by_tumorVsNormal_n$c_by_group, by_tumorVsNormal_n$groups))+
    scale_fill_manual(values = c("Tumor"='#FF0000',"Normal"='#0000FF')) +
    theme(axis.text.y = element_text(size=14),  # Larger y-axis labels
          axis.text.x = element_text(hjust=1, vjust=0.5, size=14),  # Larger x-axis labels
          axis.title.x = element_text(size=18, margin = margin(t = 10)),  # Larger x-axis title
          axis.title.y = element_text(size=18, margin = margin(r = 10)),  # Larger y-axis title
          plot.title = element_text(size=20, face="bold", hjust = 0.5),  # Larger, centered plot title
          legend.title = element_text(size=16),  # Larger legend title
          legend.text = element_text(size=14),  # Larger legend text
          plot.margin = margin(10, 10, 10, 10)) +
    ggtitle("Number and Type of Samples by Cancer Type: Tumor Vs Normal") +
    labs(fill = "Sample Types", y="Cancer Types", x="Number of Samples")
    return(p)
}
############################################################################################# 
#################################### plot_distPerCancer #####################################
           # Plot the distribution of counts for each project ID
plot_distPerCancer <- function(metadata_list, t_countdata_list){
    data <- lapply(names(metadata_list), function(id){
        t_countdata_list[[id]] %>%
        mutate(sample_type = metadata_list[[id]]$sample_type) %>%
        mutate(project_id = metadata_list[[id]]$project_id) %>%
        gather(key = "gene", value = "expression_level", -c(barcode,sample_type,project_id))
    }) %>% do.call( what = rbind, args = .)     
    ggplot(data, aes(x = expression_level)) +
        geom_histogram() +
        facet_wrap(~ project_id, scales = "free") +
    labs(title = "Distribution of Gene Expression Counts by Project ID by condition",
       x = "Expression Count", y = "Frequency")
}
#############################################################################################
######################################## transposedt ########################################
# Function to clean up output of data.table transpose:
transposedt <- function(dt, varlabel) {
  require(data.table)
  dtrows<-names(dt)
  dtcols <- as.list(c(dt[,1]))
  dtt <- data.table::transpose(dt)
  dtt[, eval(varlabel) := dtrows]
  setnames(dtt, old = names(dtt), new = c(dtcols[[1]], eval(varlabel)))
  dtt = dtt[-1,]
  setcolorder(dtt, c(eval(varlabel), names(dtt)[1:(ncol(dtt) - 1)]))
  return(dtt)
}
#############################################################################################
######################################## do_PCA_4list ########################################
prcomp_list <- function(data_list,batchs){
    lapply(data_list, function(d){
        pr <- prcomp(t(d[,-1]), center = F, scale. = F)
        pr <- as.data.table(pr$x, keep.rownames = T) %>% 
            dplyr::select(rn, PC1,PC2) %>% 
            mutate(sample_type = (batchs %>% filter(barcode %in% rn) %>% dplyr::select("sample_type") %>% pull(sample_type)),
                   plate = (batchs %>% filter(barcode %in% rn) %>% dplyr::select("plate") %>% pull(plate)),
                   tss = (batchs %>% filter(barcode %in% rn) %>% dplyr::select("tss") %>% pull(tss)),
                   project_id = (batchs %>% filter(barcode %in% rn) %>% dplyr::select("project_id") %>% pull(project_id)))
            return(pr)
            }) 
}
plot_pca_per_cancer <- function(pca_data, title = "PCA of Gene Expression Data colored by sample type", color_by = pca_data$sample_type){
    # Plot PCA
    ggplot(pca_data, aes(x = PC1, y = PC2, color = color_by)) +
        geom_point(alpha = 0.6) +
        facet_wrap(~ project_id, scales = "free") +
        labs(title = title, 
             x = "Principal Component 1", 
             y = "Principal Component 2") +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Title centered and bigger
            legend.position = "right" # Keeps the legend on the right or you can remove it if needed
        )
}

do_PCA_4list<-function(data_list,batchs){
    # data_list: a named list of data frames, each df has rownames and colnames
    pca_data_list <- prcomp_list(data_list,batchs)
    pca_data <- do.call(what = rbind, args= pca_data_list)
    # Plot PCA
    p<- plot_pca_per_cancer(pca_data)
    return(p)
}
#############################################################################################
################################# Normalization #############################################
do_normalization <- function(d,batchs,method ="TMM", group_col = "sample_type"){
    grouped_by<-batchs %>% filter(barcode %in% colnames(d[,-1])) %>% dplyr::select(all_of(group_col)) %>% pull() %>% factor
        my_DGEList <- DGEList(counts=d[,-1], group = grouped_by, samples = batchs %>% filter(barcode %in% colnames(d[,-1])))
        dge <- edgeR::normLibSizes(my_DGEList, method=method)
    dge$rn <- d$rn
    return(dge)
}
do_normalization_byList <-function(countdata_list,batchs,method ="TMM", group_col = "sample_type"){
    lapply(countdata_list, do_normalization, batchs, method = "TMM", group_col = "sample_type")
}
get_logCPM <- function(dge_list,log=TRUE){
    lapply(dge_list,function(dge){
        logCPM <- cpm(dge, log=log)
        rownames(logCPM)<-dge$rn
        return(as.data.table(logCPM, keep.rownames = T))
        })
}
#############################################################################################
####################################### cbind_dt_list #######################################
cbind_dt_list<-function(dt_list){
    lapply(dt_list, function(d){
        d<-as.data.frame(d)
        rownames(d) <- d$rn
        d[,-1] }) %>% do.call(what = cbind,args = .) %>% as.data.table(keep.rownames=T)
    }
#############################################################################################
######################################### plot_PCAs #########################################
## Function to plot the 5 plot by batchs
plot_pca <- function(pca_data, title, by, colors){
    p1 <- ggplot(pca_data, aes(x = PC2, y = PC1, color = by)) +
  geom_point(alpha = 0.7) +
  labs(title = title,
       x = "PC 2",
       y = "PC 1") +
    scale_color_manual(values = colors)+
  theme_minimal()+ theme(legend.position="top",
                        plot.margin = margin(1,0.5,0,0.5, "cm"))
    return(p1)
}
plot_PCAs <- function(pca_data, batchs){
    #PCAs
    if ("c_by_cancer" %in% names(batchs)) 
        p1<-plot_pca(pca_data, "PCA of TCGA RNA-seq Data colored by cancer type", batchs$project_id, 
                 setNames(batchs$c_by_cancer , batchs$project_id))
    #p2<-plot_pca(pca_data, "PCA of TCGA RNA-seq Data colored by GCC", pca_data$by_center, setNames(colors$by_center, colors$by_center))
    p3<-plot_pca(pca_data, "PCA of TCGA RNA-seq Data colored by plate", batchs$plate,
                 setNames(batchs$c_by_plate, batchs$plate))+ theme(legend.position="none")
    p4<-plot_pca(pca_data, "PCA of TCGA RNA-seq Data colored by TSS", batchs$tss,
                 setNames(batchs$c_by_tss, batchs$tss))+ theme(legend.position="none")
    p5<-plot_pca(pca_data, "PCA of TCGA RNA-seq Data colored by Types", batchs$sample_type,
                 setNames(batchs$c_by_type, batchs$sample_type))
    p6<-plot_pca(pca_data, "PCA of TCGA RNA-seq Data colored by Tumor-Normal", batchs$groups,
                 setNames(batchs$c_by_group, batchs$groups)) + aes(shape=batchs$sample_type)
    ## all in one
    options(repr.plot.width = 30, repr.plot.height = 10)  # Adjust width and height as needed
    library("ggpubr")
    #p2<-p2 + theme(plot.margin = margin(1,0.5,0,0.5, "cm"))
    #p2_unscaled <- p2_unscaled + theme(plot.margin = margin(1,0.5,0,0.5, "cm"))
    if ("c_by_cancer" %in% names(batchs)) 
        ggarrange(p1, p3, p4, p6,
                    labels = c("A", "B","C", "D"),
                    ncol = 4, nrow = 1)
    else
        ggarrange(p3, p4, p6,
                    labels = c("A", "B","C", "D"),
                    ncol = 3, nrow = 1)
}
#############################################################################################
################################## get_hipathia_ens_genes ###################################
get_hipathia_ens_genes <- function(){
    suppressPackageStartupMessages(library("hipathia"))
    paths<- suppressMessages(hipathia::load_pathways(species = "hsa"))
    paths$all.genes %>% length
    ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    # Retrieve Ensembl IDs using Entrez Gene IDs
    gene_mapping <- biomaRt::getBM(
      attributes = c("entrezgene_id", "ensembl_gene_id_version"),
      filters = "entrezgene_id",
      values = paths$all.genes,
      mart = ensembl
    )
    return(gene_mapping)
}
#############################################################################################         
# Plot the Venn diagram showing only the number of genes
do_venn_diagram <- function(hipathia_genes, genes0){
    library(VennDiagram)
    venn.diagram(
      x = list(Hipathia = hipathia_genes, Genes0 = genes0),
      category.names = c("Hipathia genes", "Lowly Expressed Genes"),
      filename = NULL,  # To prevent saving to file, use NULL for plotting in RStudio
      output = TRUE,
      col = "transparent",
      fill = c("lightblue", "lightgreen"),
      alpha = 0.5,
      cex = 2,  # Text size
      cat.cex = 2,  # Category label size
      cat.pos = c(-20, 20),  # Position of category names
      cat.dist = 0.05,  # Distance of category names from the diagram
      lwd = 2,  # Line width for the circles
      print.mode = "raw"  # Only display the count of genes
)}
#############################################################################################
################################## truncate_at_percentile ###################################
# Hre I was trying  several methods: TODO: rewirte only one eficient one :D
# Function to truncate values at the 99th percentile
truncate_at_percentile_long <- function(df, percentile = 99) {
  # Calculate the 99th percentile for each gene
  pct_99 <- df %>%
    group_by(rn) %>%
    summarise(max_val = quantile(Expression, percentile / 100, na.rm = TRUE))
  
  # Merge the 99th percentile back to the data and truncate values
  df_truncated <- df %>%
    left_join(pct_99, by = "rn") %>%
    mutate(Expression = ifelse(Expression > max_val, max_val, Expression)) %>%
    dplyr::select(-max_val)  # Remove the max_val column
  
  return(df_truncated)
}
# This function is from Hipathia package , when a trancation percentil is sitted to a number
truncate_at_percentile <- function(dt , truncation_percentil = 0.99){
      norm_data <- t(apply(dt, 1, function(x) {
        quan_inf <- stats::quantile(x, 1 - truncation_percentil, na.rm = TRUE)
        x[x < quan_inf] <- quan_inf
        quan_sup <- stats::quantile(x, truncation_percentil,na.rm = TRUE)
        x[x > quan_sup] <- quan_sup
        return(x)
      }))
}
#############################################################################################
do_hc <- function(brca_combat4plate, brca_batchs_paired){
    library(ggplot2)
    library(ggdendro)
    library(plotly)
    
    hc <- hclust(dist(t(brca_combat4plate[,-1])), "ward.D")
    dend <- as.dendrogram(hc)  
    
    dend_data <- dendro_data(dend)
    
    edge_colors <- dend_data$segments %>% dplyr::mutate(colors = case_when(xend<=218 ~ orders[xend],
                                                                                  .default ="black"),
                                                               colors = case_match(colors, "1" ~ "blue", "2" ~ "red"),
                                                      colors =as.factor(colors)) %>%  .$colors 
    p <- ggplot() +
      geom_segment(data = dend_data$segments, aes(x = x, y = y, xend = xend, yend = yend), color =edge_colors)  +
      geom_text(data = dend_data$labels, size = 1,  
                aes(x = x, y = y-200, label = label), 
                colour = brca_batchs_paired$c_by_group[match(dend_data$labels$label, brca_batchs_paired$barcode)]) +
        coord_flip()  
    ggplotly(p) %>%
    layout(
      yaxis = list(
        title = list(text = "Samples",
                     font = list(size = 14),  # Larger font size for better readability
                     standoff = 25)
      ),
      xaxis = list(
        title = list(text = "Labels are colored by real Tumor Vs Control\t Tree is colored by Classification Group",
                     font = list(size = 14),  # Larger font size for better readability
                     standoff = 1)
      )
    )
}
#############################################################################################
pca_plotly <- function(pca_unscaled_uncenterd_brca_combat4plate, brca_batchs_paired){
           library(ggplot2)
    library(plotly)
    p <- ggplot() +
          suppressWarnings(geom_point(data = pca_unscaled_uncenterd_brca_combat4plate$x, show.legend = T,colour = brca_batchs_paired$c_by_group, 
                     aes(x= pca_unscaled_uncenterd_brca_combat4plate$x[,1],
                         y=pca_unscaled_uncenterd_brca_combat4plate$x[,2], text=paste0(brca_batchs_paired$barcode, "(",brca_batchs_paired$groups,")"))))+
        labs(title = paste0("PCA for preprocessed data: "),
                   x = "") 
    embed_notebook(ggplotly(p, tooltip = list("text"))%>%
        layout( 
          yaxis = list(
            title = list(text = "PC2",
                         font = list(size = 14),  # Larger font size for better readability
                         standoff = 25)
          ),
          xaxis = list(
            title = list(text = "PC1",
                         font = list(size = 14),  # Larger font size for better readability
                         standoff = 1)
          )
        ))
    }
           
###################################### trash ################################################           
do_combat <- function(exp_data){
  #Combat 
  data4combat<-assay(exp_data) #dim fom BLCA: 60660 X 431
  plotMDS(data4combat, col=colors$color)
  var.data <- apply(data4combat, 1, var)
  data <- data4combat[which(var.data != 0 ),]  # Maria and Mrin said: eliminamos entras con 0 varianza, porque si no peta ComBat
  plotMDS(data, col=colors$color)
  group <- factor(exp_data$sample_type)
  my_DGEList <- DGEList(counts=data4combat, group=group)
  my_mod = model.matrix(~group, data=colData(exp_data))
  All_cancer_norm <- calcNormFactors(my_DGEList) #, method="upperquartile")
  my_data = cpm(All_cancer_norm, log=TRUE, prior.count=2) # prior.count =1 to have log 1 = 0 then the meaning will not change maybe ?
  # unique(exp_data$preservation_method)
  # unique(exp_data$icd_10_code)
  blocks<- do.call(rbind,lapply(exp_data$barcode, function(bc){
    center <- strsplit(bc,"-")[[1]][7]
    plate <-  strsplit(bc,"-")[[1]][6]
    tss <-  strsplit(bc,"-")[[1]][2]
    return(c(center, plate, tss))
    }
  ))%>% as.data.frame(stringsAsFactors = T)
  colnames(blocks) <- c("center", "plate", "tss") 
  rownames(blocks)<- exp_data$barcode
  blocks$group <-factor (exp_data$sample_type)
  # other colors 
  rownames(colors)<- colnames(exp_data)
  # by center 
  cols_center<-setNames(rainbow(length(levels(blocks$center))),levels(blocks$center))
  colors <- colors %>%  mutate(by_center = cols_center[strsplit(id,"-")[[1]][7]])
  plotMDS(my_data, col= colors$by_center)
  # by plate
  cols_plate<-setNames(rainbow(length(levels(blocks$plate))),levels(blocks$plate))
  colors <- colors %>% rowwise()%>% mutate(by_plate = cols_plate[strsplit(id,"-")[[1]][6]])
  plotMDS(my_data, col= colors$by_plate)
  # by tss
  cols_tss<-setNames(rainbow(length(levels(blocks$tss))),levels(blocks$tss))
  colors <- colors %>%  mutate(by_tss = cols_tss[strsplit(id,"-")[[1]][2]])
  plotMDS(my_data, col= colors$by_tss)
  
  # combat <- ComBat(dat=my_data, batch=blocks$plate, mod=blocks$group, par.prior=TRUE, prior.plots=TRUE)
  combat <- ComBat(dat=my_data, batch=blocks$plate, mod=blocks$group)
  # recheck here !
  plotMDS(combat, col=colors$by_plate)
  plotMDS(combat, col=colors$color)
  
  # combat <- ComBat(dat=my_data, batch=blocks$plate, mod=blocks$group, par.prior=TRUE, prior.plots=TRUE)
  combat <- ComBat(dat=my_data, batch=blocks$tss, mod=blocks$group)
  # recheck here !
  plotMDS(combat, col=colors$by_tss)
  plotMDS(combat, col=colors$color)
  
  # combat <- ComBat(dat=my_data, batch=blocks$plate, mod=blocks$group, par.prior=TRUE, prior.plots=TRUE)
  combat <- ComBat(dat=my_data, batch=blocks$center, mod=blocks$group) # in the curent example BLCA cancer :  only one center no sene here !
  # recheck here !
  plotMDS(combat, col=colors$by_center)
  plotMDS(combat, col=colors$color)
  
  r<-as.data.frame(combat)
  my_names<-getBM(filters = "entrezgene", attributes= c("mgi_symbol","entrezgene"),values=rownames(r),mart= mart)
  r$id <- rownames(r)
  r<-merge(my_names,r, by.y="id",by.x="entrezgene")
  write.table(r, "data_combat", sep = '\t', quote = FALSE, row.names = F) # write RData!
}
## Function to load Rdata from main folder:
load_data <- function(cancer_id, path2dir) {
  file_path <- file.path(path2dir, cancer_id, paste0("counts_", cancer_id, ".RData"))
  if (file.exists(file_path)) {
   load(file_path)
    return(data)
  } else {
    stop(paste("File not found for cancer:", cancer))
  }
}

plotPerCancer <- function (data){
    data %>%
      ggplot(aes(x = gene, y = expression_level)) +
      geom_boxplot() +
      # visualizes the distribution of expression level by gene by tissue type
      # i.e. one set of boxplots for nomal and tumor
      facet_wrap(facets = vars(sample_type)) +
      ylab("Expression level") +
      labs(title = "Gene expression data by tissue type"
           , caption = paste0("Source: ",unique(data$project_id)))
    }
