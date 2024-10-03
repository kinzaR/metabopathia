# TCGA RNA-seq Data Analysis Report

## Introduction

In this report, we describe the data pre-processing steps used to compare the Metabopathia approach with Hipathia. First, a Principal Component Analysis (PCA) and quality control analysis were performed on the RNA-seq data. RNA-seq data were obtained from 12 different cancer types from The Cancer Genome Atlas (TCGA) data portal (https://tcga-data.nci.nih.gov/tcga/). For these cancer types, RNA-seq counts were available for both healthy control and cancer samples:
- Bladder Urothelial Carcinoma (BLCA) [https://doi.org/10.1038/nature12965, https://doi.org/10.1016/j.cell.2017.09.007],
- Breast invasive carcinoma (BRCA) [https://doi.org/10.1038/nature11412, https://doi.org/10.1016/j.cell.2015.09.033],
- Colorectal Adenocarcinoma (COAD) [https://doi.org/10.1038/nature11252],
- Head and Neck squamous cell carcinoma (HNSC) [https://doi.org/10.1038/nature14129],
- (Kidney) Clear Cell Renal Cell Carcinoma (KIRC) [https://doi.org/10.1038/nature12222], 
- (Kidney) Papillary Renal Cell Carcinomaa (KIRP) [10.1056/NEJMoa1505917],
- Liver hepatocellular carcinoma (LIHC) [https://doi.org/10.1016/j.cell.2017.05.046],
- Lung adenocarcinoma (LUAD) [https://doi.org/10.1038/nature13385, https://doi.org/10.1038/ng.3564],
- Lung squamous cell carcinoma (LUSC) [https://doi.org/10.1038/nature11404, https://doi.org/10.1038/ng.3564],
- Prostate adenocarcinoma (PRAD) [https://doi.org/10.1016/j.cell.2015.10.025],
- Thyroid carcinoma (THCA) [https://doi.org/10.1016/j.cell.2014.09.050],
- Uterine Corpus Endometrioid Carcinoma (UCEC) [https://doi.org/10.1038/nature12113]

## Loading Libraries and sourcing needed files 


```R
# Source utils.R file: where functions are pre-definded:
source("utils.R")

# Define packages
cran_packages <- c("ggplot2", "dplyr", "tibble", "data.table", "FactoMineR", "factoextra", "mixOmics", "tidyverse")# Here I have to pay att with libraries conflicts
bioc_packages <- c("SummarizedExperiment", "SEtools", "edgeR") #EDASeq

# Call the function to install or load packages
check_and_install(cran_packages, bioc_packages, condaEnv = T)
```

## Data Acquisition and Preparation

### Define Cancer Types and Paths


```R
#This cancers was mentioned as used cancers to remove batch effect but are not the same
cancer_list <- c("BLCA","BRCA","COAD","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","THCA","UCEC")
cancer_list_names <- c("Bladder Urothelial Carcinoma",
                       "Breast invasive carcinoma",
                       "Colorectal Adenocarcinoma",
                       "Head and Neck squamous cell carcinoma",
                       "Kidney Clear Cell Renal Cell Carcinoma",
                       "Kidney Papillary Renal Cell Carcinomaa",
                       "Liver hepatocellular carcinoma",
                       "Lung adenocarcinoma",
                       "Lung squamous cell carcinoma",
                       "Prostate adenocarcinoma",
                       "Thyroid carcinoma",
                       "Uterine Corpus Endometrioid Carcinoma")
```


```R
cancer_list_id <- paste("TCGA", cancer_list, sep = "-")
names(cancer_list_names) <- cancer_list_id
```

As mentioned earlier, this study focuses on RNA-seq data from 12 cancer types. The data were previously downloaded from The Cancer Genome Atlas (TCGA) using a series of shell scripts to automate the process and handle the large dataset efficiently.

### Data Download Process

To facilitate the download, we used a shell script (get_data_in_parallel.sh) that retrieves data for each cancer type in parallel. This was executed using the Slurm workload manager with the following command (retrieved on June 27, 2024):

```bash
sbatch --job-name=getData --mem=200000 --error='.err_parallel.job' --output='.out_parallel.job' get_data_in_parallel.sh {BLCA,BRCA,COAD,HNSC,KIRC,KIRP,LIHC,LUAD,LUSC,PRAD,THCA,UCEC}
```

```sbatch``` command to submit the job to the Slurm scheduler.
```--job-name``` Option specifies the name of the job as getData for example.
```--mem=200000``` Allocates 200 GB of memory for the job, depends on available resources in the cluster!
```--error``` and ```--output```: To write error and output logs to specific files (.err_parallel.job and .out_parallel.job).

#### Script Details

- **`get_data_in_parallel.sh`**: This script initiates the download process in parallel for the specified cancer types. You can find the script at the following link: [get_data_in_parallel.sh](https://github.com/kinzaR/metabopathia/blob/main/data_examples/TCGA/get_data_in_parallel.sh).

> **Note**: This command may require adaptation depending on the resources available in your computing environment. Adjustments to this command will not affect the reproducibility of this report.


#### Shell Scripts Used

- **`get_data_in_parallel.sh`**: This script manages the parallel execution of the data retrieval for each cancer type by calling the next script for each specified type.  
  [View on GitHub](https://github.com/kinzaR/metabopathia/blob/main/data_examples/TCGA/get_data_in_parallel.sh)

- **`get_data_per_cancer.sh`**: This script is called by the parallel script to handle the data download for a single cancer type.  
  [View on GitHub](https://github.com/kinzaR/metabopathia/blob/main/data_examples/TCGA/get_data_per_cancer.sh)

- **`gdc_getData.R`**: An R script that interacts with the TCGA data portal to fetch the RNA-seq data.  
  [View on GitHub](https://github.com/kinzaR/metabopathia/blob/main/data_examples/TCGA/gdc_getData.R)


#### Data Storage

The downloaded data is organized in a folder named `processed_data`. Within this folder, there are subfolders for each cancer type containing the following:

- Temporary files generated during the download process.
- The key output file: `counts?<cancerID>.RData`, which contains the RNA-seq counts for each cancer type and will be used in subsequent analysis steps.

By organizing the download process in this manner, we ensure that the data acquisition is efficient, reproducible (data retrieved from TCGA portal on June 27, 2024), and easy to manage. This structure allows us to quickly access and process the data for downstream analysis.



```R
my_dir <- "processed_data/"
```


```R
data_list <- lapply(cancer_list_id, load_data, my_dir)
#names(data_list) <- cancer_list_id
names(data_list) <- lapply(data_list, function(d){colData(d)$project_id  %>% unique}) %>% unlist # to be more sure :S
```


```R
countdata_list <- lapply(data_list, assay, i ='unstranded') %>% lapply(data.table::as.data.table,keep.rownames = T)
# rownames are strred in the first column called rn
```


```R
metadata_list <- lapply(data_list, function(d){
    colData(d) %>% data.table::as.data.table() %>% dplyr::select(c(barcode, project_id, sample_type)) 
})
```


```R
#all_metaType in the data
lapply(data_list, function(d){
    colData(d) %>% data.table::as.data.table() %>% dplyr::select(c('sample_type_id','sample_type','specimen_type','tissue_type')) %>% unique
}) %>%do.call(what = rbind, args = .) %>% unique
```


<table class="dataframe">
<caption>A data.table: 6 × 4</caption>
<thead>
	<tr><th scope=col>sample_type_id</th><th scope=col>sample_type</th><th scope=col>specimen_type</th><th scope=col>tissue_type</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>01</td><td>Primary Tumor           </td><td>Solid Tissue</td><td>Tumor </td></tr>
	<tr><td>11</td><td>Solid Tissue Normal     </td><td>Solid Tissue</td><td>Normal</td></tr>
	<tr><td>06</td><td>Metastatic              </td><td>Unknown     </td><td>Tumor </td></tr>
	<tr><td>02</td><td>Recurrent Tumor         </td><td>Solid Tissue</td><td>Tumor </td></tr>
	<tr><td>05</td><td>Additional - New Primary</td><td>Unknown     </td><td>Tumor </td></tr>
	<tr><td>01</td><td>Primary Tumor           </td><td>Unknown     </td><td>Tumor </td></tr>
</tbody>
</table>




```R
# Checking if the .rds files already exist before saving them:
check_and_save_RDSs(data=countdata_list, metadata=metadata_list,
                    my_dir=my_dir, version="v1", 
                    prefix_countdata_file= "unstranded_counts_data_list_", 
                    prefix_metadata_file ="counts_metadata_list_")
```

    File already exists: processed_data//unstranded_counts_data_list_v1.rds
    
    File already exists: processed_data//counts_metadata_list_v1.rds
    



```R
get_object_size()
```


<table class="dataframe">
<caption>A data.frame: 4 × 2</caption>
<thead>
	<tr><th scope=col>Object</th><th scope=col>Size</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>data_list        </td><td>14877 Mb</td></tr>
	<tr><td>countdata_list   </td><td>1678 Mb </td></tr>
	<tr><td>metadata_list    </td><td>0.7 Mb  </td></tr>
	<tr><td>check_and_install</td><td>0.1 Mb  </td></tr>
</tbody>
</table>



### Load previously saved Data (RDS files)

This step is for situations where the counts and metadata have already been saved as RDS files, allowing us to load the existing data and start from there, bypassing the initial steps.


```R
version_of_data <- "v1"
countdata_list <- load_RDS_data(name="countdata_list",my_dir, version= version_of_data, prefix_file= "unstranded_counts_data_list_" )
metadata_list <- load_RDS_data(name="metadata_list",my_dir, version= version_of_data, prefix_file= "counts_metadata_list_" )
```

    countdata_list is already loaded in the environment.
    
    metadata_list is already loaded in the environment.
    


## Metadata Overview and Initial Data Exploration

### Examination of Participant and Sample Counts in Tumor and Normal Data


```R
sapply(countdata_list, dim) %>% .[2,] %>% sum - 12 # 12 rownames column ;)
```


6981


This dataset contains a total of 6,981 samples (first column is for gene names) from 6,251 cases (The number of samples is higher than the number of cases because some patients have more than one sample).
Next, we will load the metadata, which contains information about each sample—barcode, project ID, and the type of tissue sample

In the folowing plot , we can check the number of samples and the number of participants/cases per cancer type:


```R
options(repr.plot.width = 30, repr.plot.height = 10)  
plot_cases_vs_samples(metadata_list)
```


    
![png](output_28_0.png)
    



```R
# Here we will check the percentage of paired samples in normals and in tumor as well
metadata_info<-get_info_counts(metadata_list)

```


```R
metadata_info
```


<table class="dataframe">
<caption>A data.frame: 12 × 8</caption>
<thead>
	<tr><th></th><th scope=col>perc_Paired_samples_in_Normals</th><th scope=col>perc_Paired_samples_in_Tumors</th><th scope=col>Tot_part_Tumors</th><th scope=col>Tot_Samples_Tumors</th><th scope=col>Tot_part_Normals</th><th scope=col>Tot_Samples_Normals</th><th scope=col>dup_participat_in_Normals</th><th scope=col>dup_participat_in_Tumors</th></tr>
	<tr><th></th><th scope=col>&lt;named list&gt;</th><th scope=col>&lt;named list&gt;</th><th scope=col>&lt;named list&gt;</th><th scope=col>&lt;named list&gt;</th><th scope=col>&lt;named list&gt;</th><th scope=col>&lt;named list&gt;</th><th scope=col>&lt;named list&gt;</th><th scope=col>&lt;named list&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>TCGA-BLCA</th><td>         100</td><td>4.679802....</td><td> 406</td><td> 412</td><td> 19</td><td> 19</td><td>0</td><td> 6</td></tr>
	<tr><th scope=row>TCGA-BRCA</th><td>         100</td><td>10.31963....</td><td>1095</td><td>1118</td><td>113</td><td>113</td><td>0</td><td>23</td></tr>
	<tr><th scope=row>TCGA-COAD</th><td>         100</td><td>8.951965....</td><td> 458</td><td> 483</td><td> 41</td><td> 41</td><td>0</td><td>25</td></tr>
	<tr><th scope=row>TCGA-HNSC</th><td>97.72727....</td><td>8.269230....</td><td> 520</td><td> 522</td><td> 44</td><td> 44</td><td>0</td><td> 2</td></tr>
	<tr><th scope=row>TCGA-KIRC</th><td>         100</td><td>13.50844....</td><td> 533</td><td> 542</td><td> 72</td><td> 72</td><td>0</td><td> 9</td></tr>
	<tr><th scope=row>TCGA-KIRP</th><td>         100</td><td>11.03448....</td><td> 290</td><td> 291</td><td> 32</td><td> 32</td><td>0</td><td> 1</td></tr>
	<tr><th scope=row>TCGA-LIHC</th><td>         100</td><td>13.47708....</td><td> 371</td><td> 374</td><td> 50</td><td> 50</td><td>0</td><td> 3</td></tr>
	<tr><th scope=row>TCGA-LUAD</th><td>98.30508....</td><td>11.24031....</td><td> 516</td><td> 541</td><td> 59</td><td> 59</td><td>0</td><td>25</td></tr>
	<tr><th scope=row>TCGA-LUSC</th><td>         100</td><td>10.17964....</td><td> 501</td><td> 502</td><td> 51</td><td> 51</td><td>0</td><td> 1</td></tr>
	<tr><th scope=row>TCGA-PRAD</th><td>         100</td><td>10.46277....</td><td> 497</td><td> 502</td><td> 52</td><td> 52</td><td>0</td><td> 5</td></tr>
	<tr><th scope=row>TCGA-THCA</th><td>         100</td><td>11.68316....</td><td> 505</td><td> 513</td><td> 59</td><td> 59</td><td>0</td><td> 8</td></tr>
	<tr><th scope=row>TCGA-UCEC</th><td>65.71428....</td><td>4.220183....</td><td> 545</td><td> 554</td><td> 35</td><td> 35</td><td>0</td><td> 9</td></tr>
</tbody>
</table>




```R
plot_case_vs_samples_per_tissueType(metadata_info)
```


    
![png](output_31_0.png)
    


### Paired data across cancer types


```R
metadata_info %>% dplyr::select(perc_Paired_samples_in_Normals,perc_Paired_samples_in_Tumors)
```


<table class="dataframe">
<caption>A data.frame: 12 × 2</caption>
<thead>
	<tr><th></th><th scope=col>perc_Paired_samples_in_Normals</th><th scope=col>perc_Paired_samples_in_Tumors</th></tr>
	<tr><th></th><th scope=col>&lt;named list&gt;</th><th scope=col>&lt;named list&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>TCGA-BLCA</th><td>         100</td><td>4.679802....</td></tr>
	<tr><th scope=row>TCGA-BRCA</th><td>         100</td><td>10.31963....</td></tr>
	<tr><th scope=row>TCGA-COAD</th><td>         100</td><td>8.951965....</td></tr>
	<tr><th scope=row>TCGA-HNSC</th><td>97.72727....</td><td>8.269230....</td></tr>
	<tr><th scope=row>TCGA-KIRC</th><td>         100</td><td>13.50844....</td></tr>
	<tr><th scope=row>TCGA-KIRP</th><td>         100</td><td>11.03448....</td></tr>
	<tr><th scope=row>TCGA-LIHC</th><td>         100</td><td>13.47708....</td></tr>
	<tr><th scope=row>TCGA-LUAD</th><td>98.30508....</td><td>11.24031....</td></tr>
	<tr><th scope=row>TCGA-LUSC</th><td>         100</td><td>10.17964....</td></tr>
	<tr><th scope=row>TCGA-PRAD</th><td>         100</td><td>10.46277....</td></tr>
	<tr><th scope=row>TCGA-THCA</th><td>         100</td><td>11.68316....</td></tr>
	<tr><th scope=row>TCGA-UCEC</th><td>65.71428....</td><td>4.220183....</td></tr>
</tbody>
</table>



This section presents a table showing the percentage of paired samples between **Normal** and **Tumor** tissues across various cancer types from the **TCGA** (The Cancer Genome Atlas) dataset.

- **perc_Paired_samples_in_Normals**: The percentage of normal samples that are paired with tumor samples.
- **perc_Paired_samples_in_Tumors**: The percentage of tumor samples that are paired with normal samples.

> **⚠️ Important Notice:**
> - Most cancer types have **100% of their normal samples paired** with corresponding tumor samples (e.g., TCGA-BLCA, TCGA-BRCA, TCGA-COAD).
> - The percentage of tumor samples paired with normal samples is generally lower, with values ranging from **4.22% to 13.51%**, depending on the cancer type.
> - **TCGA-UCEC** has a notably lower percentage of paired normal samples (**65.71%**), which is an outlier compared to the other cancer types.
> - The variability in tumor-normal pairing highlights that while most normal samples are paired, a much smaller proportion of tumor samples have corresponding normal pairs.



### Summary of Preservation Methods Across Sample Types and Cancer Types 


```R
preservation_summary <- suppressMessages(summarize_preservation_methods(data_list))
```


```R
plots_preservation_methods(preservation_summary)
piPlots_preservation_methods(preservation_summary,0) # this need enhancment : code and colors need to be the same coherents
piPlots_preservation_methods(preservation_summary,4)
piPlots_preservation_methods(preservation_summary,8)
```


    
![png](output_37_0.png)
    



    
![png](output_37_1.png)
    



    
![png](output_37_2.png)
    



    
![png](output_37_3.png)
    



```R
preservation_summary %>% filter(cancer_id=="TCGA-BLCA")
```


<table class="dataframe">
<caption>A grouped_df: 5 × 4</caption>
<thead>
	<tr><th scope=col>sample_type</th><th scope=col>preservation_method</th><th scope=col>n</th><th scope=col>cancer_id</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Primary Tumor      </td><td>FFPE   </td><td>  3</td><td>TCGA-BLCA</td></tr>
	<tr><td>Primary Tumor      </td><td>OCT    </td><td>345</td><td>TCGA-BLCA</td></tr>
	<tr><td>Solid Tissue Normal</td><td>OCT    </td><td> 10</td><td>TCGA-BLCA</td></tr>
	<tr><td>Primary Tumor      </td><td>Unknown</td><td> 64</td><td>TCGA-BLCA</td></tr>
	<tr><td>Solid Tissue Normal</td><td>Unknown</td><td>  9</td><td>TCGA-BLCA</td></tr>
</tbody>
</table>



> These graphs need enhancements in visualization for better interpretation. The aim here is to show that the FFPE preservation method was used exclusively for Tumor samples and not for Normal samples in certain cancers.
>
> **e.g. BLCA cancer :** 
> 
> - **3 Primary Tumor samples** have used the **FFPE preservation method**. This method was used exclusively for **Tumor samples** and not for **Normal samples** in some cancers.
> - in BLCA There are **73 samples** with an **Unknown preservation method**, of which **64 are Primary Tumor samples** and **9 are Solid Tissue Normal samples**.



```R
get_object_size()
```


<table class="dataframe">
<caption>A data.frame: 10 × 2</caption>
<thead>
	<tr><th scope=col>Object</th><th scope=col>Size</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>data_list                            </td><td>14877 Mb</td></tr>
	<tr><td>countdata_list                       </td><td>1678 Mb </td></tr>
	<tr><td>get_info_counts                      </td><td>1 Mb    </td></tr>
	<tr><td>metadata_list                        </td><td>0.7 Mb  </td></tr>
	<tr><td>change_label                         </td><td>0.1 Mb  </td></tr>
	<tr><td>check_and_install                    </td><td>0.1 Mb  </td></tr>
	<tr><td>get_object_size                      </td><td>0.1 Mb  </td></tr>
	<tr><td>plot_case_vs_samples_per_tissueType  </td><td>0.1 Mb  </td></tr>
	<tr><td>plot_cases_vs_samples                </td><td>0.1 Mb  </td></tr>
	<tr><td>plot_Preservation_Method_Distribution</td><td>0.1 Mb  </td></tr>
</tbody>
</table>




```R
# No need for data_list, will be removed
rm(data_list)
get_object_size()
```


<table class="dataframe">
<caption>A data.frame: 9 × 2</caption>
<thead>
	<tr><th scope=col>Object</th><th scope=col>Size</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>countdata_list                       </td><td>1678 Mb</td></tr>
	<tr><td>get_info_counts                      </td><td>1 Mb   </td></tr>
	<tr><td>metadata_list                        </td><td>0.7 Mb </td></tr>
	<tr><td>change_label                         </td><td>0.1 Mb </td></tr>
	<tr><td>check_and_install                    </td><td>0.1 Mb </td></tr>
	<tr><td>get_object_size                      </td><td>0.1 Mb </td></tr>
	<tr><td>plot_case_vs_samples_per_tissueType  </td><td>0.1 Mb </td></tr>
	<tr><td>plot_cases_vs_samples                </td><td>0.1 Mb </td></tr>
	<tr><td>plot_Preservation_Method_Distribution</td><td>0.1 Mb </td></tr>
</tbody>
</table>



## Quality Control and Normalization

### Filter Lowly Expressed Genes (skipped):
In this analysis, filtering out lowly expressed genes has been intentionally skipped. Unlike differential expression analysis pipelines, where removing such genes is common, the mechanistic modeling approach utilized here requires the retention of all available data. If lowly expressed genes were excluded, the method would artificially impute missing values (e.g., using 0.5 or another value), which would not accurately represent the biological reality of low expression. Retaining these genes ensures that the model reflects true biological conditions.


```R
# TODO : Get  metabolized genes ! as well 
```


```R
#Get hipathia Ensemble genes
hipathia_genes_table <- get_hipathia_ens_genes()
hipathia_genes<-hipathia_genes_table$ensembl_gene_id_version %>% unique
```


```R
#identify low Expressed genes 
genes0<-sapply(countdata_list, function(d){
    rowSums(d[,-1])
    }) %>% rowSums 
genes0<-which( genes0 ==0) %>% countdata_list[[1]][.,1] %>% .$rn
```


```R
# tradutiona rowSums Vs sofstiucated filterByExpr: TODO: compare the effect of using filterByExpr
#agreedOnRemoving <- sapply(names(countdata_list),function(id){
#    filterByExpr(DGEList(countdata_list[[id]][,-1]), design = (batchs %>% filter(project_id==id) %>% .$groups), 
#                 min.total.count = 1, min.count = 1)
#    }) %>% rowSums 
#agreedOnRemoving <-countdata_list[[1]][which(agreedOnRemoving ==0),rn]
```


```R
# Find common and unique genes
common_genes <- intersect(hipathia_genes, genes0)
only_in_hipathia <- setdiff(hipathia_genes, genes0) ## Here I have to save it somewhere to recheck again !
only_in_genes0 <- setdiff(genes0, hipathia_genes)
do_venn_diagram(hipathia_genes, genes0) %>% grid.draw(.)
```

    Loading required package: futile.logger
    
    
    Attaching package: ‘VennDiagram’
    
    
    The following object is masked from ‘package:ggpubr’:
    
        rotate
    
    



    
![png](output_48_1.png)
    


Great! The intersection between the lowly expressed genes across all cancers and Hipathia is **6 genes**, so we can safely remove the other genes (2331 out of 2337 genes) that have a count of 0 in all samples.


```R
hipathia_null_genes_table<-hipathia_genes_table %>% filter(ensembl_gene_id_version %in% common_genes)
```

These genes will be saved to be tracked in the downstream analysis, allowing for more accurate results.


```R
write.table(x = hipathia_null_genes_table, file = file.path(my_dir, paste0("hipathia_null_genes_",version_of_data,".tsv")), sep = "\t", quote = F,row.names = F, col.names = T)
```


```R
countdata_list <- lapply(countdata_list, function(d){
    d %>% filter(!rn %in% only_in_genes0)
})
```


```R
get_object_size()
```


<table class="dataframe">
<caption>A data.frame: 14 × 2</caption>
<thead>
	<tr><th scope=col>Object</th><th scope=col>Size</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>countdata_list                       </td><td>1613.7 Mb</td></tr>
	<tr><td>get_info_counts                      </td><td>1 Mb     </td></tr>
	<tr><td>metadata_list                        </td><td>0.7 Mb   </td></tr>
	<tr><td>hipathia_genes                       </td><td>0.4 Mb   </td></tr>
	<tr><td>hipathia_genes_table                 </td><td>0.4 Mb   </td></tr>
	<tr><td>only_in_hipathia                     </td><td>0.4 Mb   </td></tr>
	<tr><td>genes0                               </td><td>0.2 Mb   </td></tr>
	<tr><td>only_in_genes0                       </td><td>0.2 Mb   </td></tr>
	<tr><td>change_label                         </td><td>0.1 Mb   </td></tr>
	<tr><td>check_and_install                    </td><td>0.1 Mb   </td></tr>
	<tr><td>get_object_size                      </td><td>0.1 Mb   </td></tr>
	<tr><td>plot_case_vs_samples_per_tissueType  </td><td>0.1 Mb   </td></tr>
	<tr><td>plot_cases_vs_samples                </td><td>0.1 Mb   </td></tr>
	<tr><td>plot_Preservation_Method_Distribution</td><td>0.1 Mb   </td></tr>
</tbody>
</table>




```R
rm(hipathia_null_genes_table, hipathia_genes_table,only_in_hipathia, only_in_genes0, hipathia_genes, do_venn_diagram, genes0,get_info_counts)
```


```R
get_object_size()%>%head
```


<table class="dataframe">
<caption>A data.frame: 6 × 2</caption>
<thead>
	<tr><th></th><th scope=col>Object</th><th scope=col>Size</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>countdata_list                     </td><td>1613.7 Mb</td></tr>
	<tr><th scope=row>2</th><td>metadata_list                      </td><td>0.7 Mb   </td></tr>
	<tr><th scope=row>3</th><td>change_label                       </td><td>0.1 Mb   </td></tr>
	<tr><th scope=row>4</th><td>check_and_install                  </td><td>0.1 Mb   </td></tr>
	<tr><th scope=row>5</th><td>get_object_size                    </td><td>0.1 Mb   </td></tr>
	<tr><th scope=row>6</th><td>plot_case_vs_samples_per_tissueType</td><td>0.1 Mb   </td></tr>
</tbody>
</table>



### Exploratory analysis

#### Defining batchs and colors

Before proceeding with the analysis, we will define a set of colors for each batch to ensure clarity in visualizations. To accommodate colorblindness, we will use the Okabe-Ito palette, recommended by Okabe & Ito (2008). This palette is the default in base R (version 4.0.0 and later) and is widely recognized for its accessibility and effectiveness in presenting results to a broad audience.

In addition to the Okabe-Ito palette, other colors were carefully selected and tested using the "Let's get color blind" plugin. This tool simulates color deficiency based on the model by Gustavo M. Machado et al. ("A Physiologically-based Model for Simulation of Color Vision Deficiency," IEEE Transactions on Visualization and Computer Graphics, 2009), ensuring accessibility for individuals with color vision deficiencies.

First, let's assign colors for each sample type:


```R
#get batchs, with colors by batch
batchs<- get_batchs(metadata_list)
```

Check the [TCGA Barcode web page](https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/) for mor information about Barcode Reading.
Here is the schema of the Creating Barcodes process:

<img src="https://docs.gdc.cancer.gov/Encyclopedia/pages/images/creating_barcodes.png" alt="creating_barcodes" style="display: block; margin-left: auto; margin-right: auto;"/>

The following figure illustrates how metadata identifiers comprise a barcode and shows how to read barcodes. 
A TCGA barcode is composed of a series of identifiers, each uniquely identifying a TCGA data element.

<img src="https://docs.gdc.cancer.gov/Encyclopedia/pages/images/barcode.png" alt="barcode_reading" style="display: block; margin-left: auto; margin-right: auto;"/>


```R
levels(batchs$center)
```


'07'


According to the [GDC TCGA Code Tables](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/center-codes), the code `07` corresponds to the following center:

- **Center Name**: University of North Carolina
- **Center Type**: CGCC
- **Display Name**: UNC


```R
levels(batchs$plate) %>% length # by plate: The number of plates is large, making it challenging to create color-blind friendly visuals. A solution needs to be found.
```


178



```R
levels(batchs$tss) %>% length
```


373



```R
options(repr.plot.width = 20, repr.plot.height = 10)  # Adjust width and height as needed
plot_Samples_Cancer(batchs)
```

    [1m[22m`summarise()` has grouped output by 'project_id', 'sample_type'. You can
    override using the `.groups` argument.



    
![png](output_71_1.png)
    



```R
# Save the plot as a PNG file
#ggsave("figures/samplesTypes_per_cancer.png", plot = p, width = 10, height = 8, dpi = 300)
```


```R
plot_by_tumorVsNormal(batchs)
```

    [1m[22m`summarise()` has grouped output by 'project_id', 'groups'. You can override
    using the `.groups` argument.



    
![png](output_73_1.png)
    


> **⚠️ Important Notice:**
> 
> The data is **imbalanced**—the number of **Tumor samples** is significantly higher than the **Normal (or control) samples**. This imbalance may affect the analysis and should be carefully considered when interpreting the results.


#### Some plots for row data 


```R
get_object_size() %>% head
```


<table class="dataframe">
<caption>A data.frame: 6 × 2</caption>
<thead>
	<tr><th></th><th scope=col>Object</th><th scope=col>Size</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>countdata_list   </td><td>1613.7 Mb</td></tr>
	<tr><th scope=row>2</th><td>batchs           </td><td>2 Mb     </td></tr>
	<tr><th scope=row>3</th><td>metadata_list    </td><td>0.7 Mb   </td></tr>
	<tr><th scope=row>4</th><td>get_batchs       </td><td>0.2 Mb   </td></tr>
	<tr><th scope=row>5</th><td>change_label     </td><td>0.1 Mb   </td></tr>
	<tr><th scope=row>6</th><td>check_and_install</td><td>0.1 Mb   </td></tr>
</tbody>
</table>




```R
#t_countdata_list <- lapply(countdata_list,transposedt, varlabel="barcode")
#plot_distPerCancer(metadata_list,t_countdata_list) # tobe fixed
```

#### Principal component analysis for row data cross cancer type


```R
pca_data <- prcomp_list(countdata_list,batchs) %>% do.call(what = rbind, args= .)
```


```R
plot_pca_per_cancer(pca_data, color_by= pca_data$sample_type)+ 
    scale_color_manual(values = setNames(batchs$c_by_type, batchs$sample_type))
```


    
![png](output_80_0.png)
    



```R
plot_pca_per_cancer(pca_data, color_by= batchs$groups, title = "PCA of Gene Expression Data colored by Normal Vs Tumor") + 
    scale_color_manual(values = setNames(batchs$c_by_group, batchs$groups))
```


    
![png](output_81_0.png)
    



```R
plot_pca_per_cancer(pca_data, color_by= pca_data$plate, title = "PCA of Gene Expression Data colored by plates")+ theme(legend.position = "none") + 
    scale_color_manual(values = setNames(batchs$c_by_plate, batchs$plate))
```


    
![png](output_82_0.png)
    



```R
plot_pca_per_cancer(pca_data, color_by= pca_data$tss, title = "PCA of Gene Expression Data colored by Tissue Source Site (TSS)")+ theme(legend.position = "none")+ 
    scale_color_manual(values = setNames(batchs$c_by_tss, batchs$tss))
```


    
![png](output_83_0.png)
    


#### TMM per subtypes of tissues


```R
#plotMDS(data2plot_per_project_list$`TCGA-KIRC`[,-1])
```


```R
dge_list <- do_normalization_byList(countdata_list,batchs,method ="TMM", group_col = "sample_type") 
```


```R
options(repr.plot.width = 10, repr.plot.height = 10)  # Adjust width and height as needed
par(mfrow = c(3, 4))
v_list<-lapply(dge_list,voom,plot = T) # http://www.genomebiology.com/2014/15/2/R29
options(repr.plot.width = 30, repr.plot.height = 10)  # Adjust width and height as needed
```


    
![png](output_87_0.png)
    



```R
logCPM_list <- get_logCPM(dge_list,log=TRUE)
```

#### PCA for TMM as a QC


```R
pca_logCPM <- prcomp_list(logCPM_list,batchs) %>% do.call(what = rbind, args= .)
plot_pca_per_cancer(pca_logCPM, color_by= pca_logCPM$sample_type)+ 
    scale_color_manual(values = setNames(batchs$c_by_type, batchs$sample_type))
```


    
![png](output_90_0.png)
    



```R
plot_pca_per_cancer(pca_logCPM, color_by= batchs$groups, title = "PCA of Gene Expression Data colored by Normal Vs Tumor") + 
    scale_color_manual(values = setNames(batchs$c_by_group, batchs$groups))
```


    
![png](output_91_0.png)
    



```R
plot_pca_per_cancer(pca_logCPM, color_by= pca_logCPM$plate, title = "PCA of Gene Expression Data colored by plates")+ theme(legend.position = "none") + 
    scale_color_manual(values = setNames(batchs$c_by_plate, batchs$plate))
```


    
![png](output_92_0.png)
    



```R
plot_pca_per_cancer(pca_logCPM, color_by= pca_logCPM$tss, title = "PCA of Gene Expression Data colored by Tissue Source Site (TSS)")+ theme(legend.position = "none")+ 
    scale_color_manual(values = setNames(batchs$c_by_tss, batchs$tss))
```


    
![png](output_93_0.png)
    


> **Note:** While performing the TMM (Trimmed Mean of M-values) normalization, two scenarios were considered:
> 
> 1. The groups were based on tissue subtypes (e.g., metastatic, primary tumor) being considered separately.
> 2. The samples were grouped into broader categories, where all subtypes were grouped under the "tumor" category, and only a distinction between "tumor" and "normal" samples was made.
> 
> as a result of these two approaches, the expression matrix after normalization remains the same in both cases, as no effect of subtype stratification was dectect on the normalization step. in either scenario during the normalization process.


#### Principal Component Analysis for normalized data


```R
logCPM_list_file<-file.path(my_dir,paste0("logCPM_list_file_",version_of_data,".rds"))
saveRDS(object = logCPM_list, file = logCPM_list_file)
```


```R
# chekc the data here
```


```R
merged_count_data <-cbind_dt_list(countdata_list)
```


```R
merged_normalized_data <-cbind_dt_list(logCPM_list)
```


```R
class(merged_count_data)
head(merged_count_data)
#colnames(merged_count_data)
```


```R
pca_count_data <- stats::prcomp(t(merged_count_data[,-1]), scale = F,center = F)
```


```R
pca_normalized <- stats::prcomp(t(merged_normalized_data[,-1]), scale = F,center = F)
```


```R
pca_countX <- as.data.table(pca_count_data$x)
```


```R
pca_normalizedX <- as.data.table(pca_normalized$x)
```

**PCA Plots of Samples to Discover Batch Effects**



```R
plot_PCAs(pca_countX, batchs)
```


    
![png](output_106_0.png)
    



```R
plot_PCAs(pca_normalizedX, batchs)
```


    
![png](output_107_0.png)
    


Principal Component Analysis (PCA) was performed to detect possible batch effects in the RNA-seq counts for 12 cancer types. The PCA plots are colored by various variables to reveal potential sources of variation:

- **Panel A**: Samples are colored by cancer type. This plot reveals the distribution of samples across different cancer types and highlights clustering based on cancer type.

- **Panel B**: Samples are colored by plate. This panel shows a 'clear' batch effect, with distinct clustering of samples from different plates, indicating potential technical variation related to the processing plates.

- **Panel C**: Samples are colored by Tissue Source Site (TSS). This plot provides insight into how TSS affects sample distribution and helps identify any associated batch effects.

- **Panel D**: Samples are colored by normal versus tumor status. This panel illustrates the separation between normal and tumor samples and helps assess whether batch effects might obscure or mimic biological differences.

A clear batch effect was identified in the plots for cancer type and plate. To address and correct these batch effects, the COMBAT algorithm will be applied. COMBAT is a statistical technique used to adjust for batch effects in high-dimensional data, ensuring that the observed biological variations are not confounded by technical artifacts. This adjustment corrects for batch effects related to both the plate and cancer types, facilitating a more accurate comparison of overall tumor data rather than specific tumor types.


### Remove batch effect using COMBAT:

After normalization we need to check for bach effect and then removing it if needed, because normalization does not remove batch effects, which affect specific subsets of genes and may affect different genes in different ways [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3880143/]


```R
my_mod <-model.matrix(~as.factor(groups), data=batchs)
```


```R
var.data <- apply(merged_normalized_data[,-1], 1, var)
data4combat_filtered <- merged_normalized_data[which(var.data != 0 ),]  # Maria P and Mrin mentioned: We removed entries with zero variance because including them would cause COMBAT to fail.
```


```R
var.data %>% summary
```


        Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
     0.08676  0.21392  0.56164  1.22636  1.32804 32.96942 



```R
dim(data4combat_filtered[,-1])
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>58329</li><li>6981</li></ol>




```R
dim(merged_normalized_data[,-1])
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>58329</li><li>6981</li></ol>



> **Important Note**: The removed information needs to be reconsidered. Removing genes with zero variance before applying COMBAT and then imputing them for Pathway Activity analysis introduces an additional layer of potential inaccuracy. These genes do not vary, and their imputation could misslead the analysis. (**NOTHING WERE REMOVED !**)



```R
# Apply ComBat for plate BE
combat4plate <- sva::ComBat(dat = merged_normalized_data[,-1],
                 mod = my_mod,
                  batch = batchs$plate, 
                  par.prior = TRUE, 
                  prior.plots = FALSE)  # Suppress plotting
```

    Found 23922 genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.


    Using the 'mean only' version of ComBat
    
    Found178batches
    
    Note: one batch has only one sample, setting mean.only=TRUE
    
    Adjusting for1covariate(s) or covariate level(s)
    
    Standardizing Data across genes
    
    Fitting L/S model and finding priors
    
    Finding parametric adjustments
    
    Adjusting the Data
    
    



```R
dim(merged_normalized_data[,-1])
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>58329</li><li>6981</li></ol>




```R
dim(combat4plate)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>58329</li><li>6981</li></ol>




```R
head(combat4plate)
```


<table class="dataframe">
<caption>A matrix: 6 × 6981 of type dbl</caption>
<thead>
	<tr><th scope=col>TCGA-BLCA.TCGA-CU-A3KJ-01A-11R-A21D-07</th><th scope=col>TCGA-BLCA.TCGA-K4-A3WU-01B-11R-A23N-07</th><th scope=col>TCGA-BLCA.TCGA-DK-A3IU-01A-11R-A20F-07</th><th scope=col>TCGA-BLCA.TCGA-GV-A40G-01A-11R-A23N-07</th><th scope=col>TCGA-BLCA.TCGA-DK-A3IN-01A-11R-A20F-07</th><th scope=col>TCGA-BLCA.TCGA-SY-A9G0-01A-12R-A38B-07</th><th scope=col>TCGA-BLCA.TCGA-XF-A9SU-01A-31R-A39I-07</th><th scope=col>TCGA-BLCA.TCGA-UY-A78O-01A-12R-A33J-07</th><th scope=col>TCGA-BLCA.TCGA-GD-A2C5-01A-12R-A180-07</th><th scope=col>TCGA-BLCA.TCGA-G2-A2EK-01A-22R-A18C-07</th><th scope=col>⋯</th><th scope=col>TCGA-UCEC.TCGA-E6-A1LX-01A-11R-A14D-07</th><th scope=col>TCGA-UCEC.TCGA-FI-A3PV-01A-11R-A22K-07</th><th scope=col>TCGA-UCEC.TCGA-BG-A0MT-01A-11R-A104-07</th><th scope=col>TCGA-UCEC.TCGA-EY-A549-01A-11R-A27V-07</th><th scope=col>TCGA-UCEC.TCGA-D1-A1O7-01A-11R-A14D-07</th><th scope=col>TCGA-UCEC.TCGA-AJ-A3NF-01A-11R-A22K-07</th><th scope=col>TCGA-UCEC.TCGA-BK-A0C9-01A-11R-A00V-07</th><th scope=col>TCGA-UCEC.TCGA-AJ-A5DW-01A-11R-A27V-07</th><th scope=col>TCGA-UCEC.TCGA-AP-A05J-01A-11R-A00V-07</th><th scope=col>TCGA-UCEC.TCGA-D1-A103-01A-11R-A10J-07</th></tr>
</thead>
<tbody>
	<tr><td> 5.667448</td><td> 8.491573</td><td> 4.062396</td><td> 7.648249</td><td> 6.659637</td><td> 5.579963</td><td> 6.965262</td><td> 8.114996</td><td> 7.569854</td><td> 8.233744</td><td>⋯</td><td> 5.272579</td><td> 6.6572916</td><td> 6.554263</td><td> 6.310736</td><td> 6.629957</td><td> 6.839850</td><td> 5.628924</td><td> 7.109832</td><td> 6.361204</td><td> 6.080712</td></tr>
	<tr><td>-4.147907</td><td>-1.372953</td><td>-3.713166</td><td>-3.099971</td><td>-4.204261</td><td>-2.179601</td><td>-2.999680</td><td>-3.799766</td><td>-3.672981</td><td>-4.094614</td><td>⋯</td><td>-3.509750</td><td>-4.5633311</td><td>-2.915424</td><td>-3.424957</td><td>-3.624429</td><td>-3.529995</td><td>-3.654775</td><td>-2.154001</td><td>-0.768419</td><td>-3.631352</td></tr>
	<tr><td> 5.879755</td><td> 5.127999</td><td> 5.701277</td><td> 5.764180</td><td> 4.657800</td><td> 5.489051</td><td> 5.650628</td><td> 6.182716</td><td> 5.570483</td><td> 5.694698</td><td>⋯</td><td> 4.944719</td><td> 5.3346419</td><td> 4.969617</td><td> 5.109477</td><td> 4.357961</td><td> 5.024670</td><td> 3.403225</td><td> 5.519378</td><td> 4.148964</td><td> 4.255326</td></tr>
	<tr><td> 3.578146</td><td> 4.562290</td><td> 3.643068</td><td> 4.046311</td><td> 3.725897</td><td> 3.726396</td><td> 2.465567</td><td> 3.539136</td><td> 4.516733</td><td> 4.193857</td><td>⋯</td><td> 3.447282</td><td> 3.4959761</td><td> 4.199439</td><td> 4.668052</td><td> 2.777284</td><td> 3.687770</td><td> 3.755970</td><td> 4.095718</td><td> 2.698562</td><td> 3.537804</td></tr>
	<tr><td> 4.191013</td><td> 4.391384</td><td> 3.574144</td><td> 3.005857</td><td> 4.636430</td><td> 3.003469</td><td> 2.456102</td><td> 3.043248</td><td> 4.224812</td><td> 2.863929</td><td>⋯</td><td> 2.871428</td><td> 2.6770370</td><td> 3.589796</td><td> 2.043183</td><td> 2.307135</td><td> 3.050835</td><td> 2.972309</td><td> 3.733101</td><td> 2.492515</td><td> 2.508736</td></tr>
	<tr><td> 2.438952</td><td> 3.440620</td><td> 5.400516</td><td> 1.902097</td><td> 2.939534</td><td> 3.023844</td><td> 3.083837</td><td> 2.237102</td><td> 2.466669</td><td> 2.589599</td><td>⋯</td><td> 1.242578</td><td> 0.4407497</td><td> 2.173327</td><td> 2.965949</td><td> 2.748985</td><td> 3.414149</td><td> 3.132598</td><td> 2.067048</td><td> 1.073637</td><td> 2.648755</td></tr>
</tbody>
</table>




```R
#chack data after plate batch effect removal
pca_unscaled_uncenterd_combat4plate <- prcomp(t(combat4plate), scale = F, center = F)
```


```R
plot_PCAs(as.data.table(pca_unscaled_uncenterd_combat4plate$x), batchs)
```


    
![png](output_122_0.png)
    


> **Important Remark**: Before and after removing batch effects by plate, there is no clear difference observed between the datasets. This step needs to be revisited and redone to ensure proper batch effect correction. Further investigation is required to confirm the effectiveness of the batch effect removal process.



```R
# Apply ComBat for cancer types BE
combat4_plate_and_cancerType <- sva::ComBat(dat=combat4plate, 
                  mod = my_mod,
                  batch=batchs$project_id,
                  par.prior = TRUE, 
                  prior.plots = FALSE)  # Suppress plotting
```

    Found12batches
    
    Adjusting for1covariate(s) or covariate level(s)
    
    Standardizing Data across genes
    
    Fitting L/S model and finding priors
    
    Finding parametric adjustments
    
    Adjusting the Data
    
    



```R
#chack finaly data after plate batch and cancer types  effect removal
pca_unscaled_uncenterd_combat4_plate_and_cancerType <- prcomp(t(combat4_plate_and_cancerType), scale = F, center = F)
plot_PCAs(as.data.table(pca_unscaled_uncenterd_combat4_plate_and_cancerType$x), batchs)
```


    
![png](output_125_0.png)
    


> **Important Remark**: The differences between tumor and control classes are not clearly visible after batch effect correction. This raises concerns that we may be losing important biological effects during the batch effect cleaning process. Further review and re-evaluation are needed to ensure that the biological distinctions between tumor and control samples are accurately preserved and not inadvertently masked by the correction process.



```R
rownames(combat4_plate_and_cancerType) <- merged_normalized_data$rn
```


```R
head(combat4_plate_and_cancerType)
```


<table class="dataframe">
<caption>A matrix: 6 × 6981 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>TCGA-BLCA.TCGA-CU-A3KJ-01A-11R-A21D-07</th><th scope=col>TCGA-BLCA.TCGA-K4-A3WU-01B-11R-A23N-07</th><th scope=col>TCGA-BLCA.TCGA-DK-A3IU-01A-11R-A20F-07</th><th scope=col>TCGA-BLCA.TCGA-GV-A40G-01A-11R-A23N-07</th><th scope=col>TCGA-BLCA.TCGA-DK-A3IN-01A-11R-A20F-07</th><th scope=col>TCGA-BLCA.TCGA-SY-A9G0-01A-12R-A38B-07</th><th scope=col>TCGA-BLCA.TCGA-XF-A9SU-01A-31R-A39I-07</th><th scope=col>TCGA-BLCA.TCGA-UY-A78O-01A-12R-A33J-07</th><th scope=col>TCGA-BLCA.TCGA-GD-A2C5-01A-12R-A180-07</th><th scope=col>TCGA-BLCA.TCGA-G2-A2EK-01A-22R-A18C-07</th><th scope=col>⋯</th><th scope=col>TCGA-UCEC.TCGA-E6-A1LX-01A-11R-A14D-07</th><th scope=col>TCGA-UCEC.TCGA-FI-A3PV-01A-11R-A22K-07</th><th scope=col>TCGA-UCEC.TCGA-BG-A0MT-01A-11R-A104-07</th><th scope=col>TCGA-UCEC.TCGA-EY-A549-01A-11R-A27V-07</th><th scope=col>TCGA-UCEC.TCGA-D1-A1O7-01A-11R-A14D-07</th><th scope=col>TCGA-UCEC.TCGA-AJ-A3NF-01A-11R-A22K-07</th><th scope=col>TCGA-UCEC.TCGA-BK-A0C9-01A-11R-A00V-07</th><th scope=col>TCGA-UCEC.TCGA-AJ-A5DW-01A-11R-A27V-07</th><th scope=col>TCGA-UCEC.TCGA-AP-A05J-01A-11R-A00V-07</th><th scope=col>TCGA-UCEC.TCGA-D1-A103-01A-11R-A10J-07</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ENSG00000000003.15</th><td> 5.240923</td><td> 7.2829385</td><td> 4.080372</td><td> 6.673164</td><td> 5.958337</td><td> 5.177666</td><td> 6.179322</td><td> 7.010650</td><td> 6.616479</td><td> 7.096513</td><td>⋯</td><td> 4.841792</td><td> 6.408050</td><td> 6.291513</td><td> 6.016058</td><td> 6.377131</td><td> 6.614543</td><td> 5.244857</td><td> 6.919921</td><td> 6.07314351</td><td> 5.755877</td></tr>
	<tr><th scope=row>ENSG00000000005.6</th><td>-3.824657</td><td>-0.3473206</td><td>-3.279875</td><td>-2.511472</td><td>-3.895275</td><td>-1.358143</td><td>-2.385796</td><td>-3.388396</td><td>-3.229519</td><td>-3.757873</td><td>⋯</td><td>-3.267034</td><td>-4.490457</td><td>-2.576900</td><td>-3.168572</td><td>-3.400200</td><td>-3.290543</td><td>-3.435438</td><td>-1.692732</td><td>-0.08378953</td><td>-3.408239</td></tr>
	<tr><th scope=row>ENSG00000000419.13</th><td> 5.424674</td><td> 4.7260632</td><td> 5.258813</td><td> 5.317270</td><td> 4.289103</td><td> 5.061591</td><td> 5.211746</td><td> 5.706218</td><td> 5.137266</td><td> 5.252700</td><td>⋯</td><td> 5.205149</td><td> 5.504302</td><td> 5.224252</td><td> 5.331553</td><td> 4.754983</td><td> 5.266489</td><td> 4.022501</td><td> 5.646033</td><td> 4.59463864</td><td> 4.676241</td></tr>
	<tr><th scope=row>ENSG00000000457.14</th><td> 3.633643</td><td> 4.4720291</td><td> 3.688951</td><td> 4.032470</td><td> 3.759512</td><td> 3.759937</td><td> 2.685845</td><td> 3.600412</td><td> 4.433219</td><td> 4.158163</td><td>⋯</td><td> 3.675643</td><td> 3.714771</td><td> 4.280032</td><td> 4.656581</td><td> 3.137272</td><td> 3.868885</td><td> 3.923686</td><td> 4.196688</td><td> 3.07401623</td><td> 3.748381</td></tr>
	<tr><th scope=row>ENSG00000000460.17</th><td> 3.306594</td><td> 3.4591858</td><td> 2.836821</td><td> 2.404045</td><td> 3.645799</td><td> 2.402227</td><td> 1.985382</td><td> 2.432520</td><td> 3.332333</td><td> 2.295960</td><td>⋯</td><td> 2.886295</td><td> 2.708516</td><td> 3.543275</td><td> 2.128827</td><td> 2.370223</td><td> 3.050370</td><td> 2.978555</td><td> 3.674334</td><td> 2.53976206</td><td> 2.554597</td></tr>
	<tr><th scope=row>ENSG00000000938.13</th><td> 2.401667</td><td> 3.3609718</td><td> 5.237980</td><td> 1.887517</td><td> 2.881078</td><td> 2.961823</td><td> 3.019278</td><td> 2.208354</td><td> 2.428211</td><td> 2.545943</td><td>⋯</td><td> 1.987255</td><td> 1.250138</td><td> 2.842889</td><td> 3.571543</td><td> 3.372089</td><td> 3.983571</td><td> 3.724744</td><td> 2.745188</td><td> 1.83194884</td><td> 3.279948</td></tr>
</tbody>
</table>




```R
################To detleted: just to check how data if we start with removing cancer type batch effect:
suppressPackageStartupMessages(library(sva))
#combat CANCER
combat_pre2 <- ComBat(dat=merged_normalized_data[,-1], batch=batchs$project_id,
                  mod = my_mod,
                  par.prior = TRUE, 
                  prior.plots = FALSE)  # Suppress plotting
#combat PLATE
combat_pre3 <- ComBat(dat=combat_pre2, batch=batchs$plate,
                  mod = my_mod,
                  par.prior = TRUE, 
                  prior.plots = FALSE)  # Suppress plotting
#PCA
pca_unscaled_uncenterd_combat_pre2 <- prcomp(t(combat_pre3), scale = F, center = F)
plot_PCAs(as.data.table(pca_unscaled_uncenterd_combat_pre2$x), batchs)
```

    Found12batches
    
    Adjusting for1covariate(s) or covariate level(s)
    
    Standardizing Data across genes
    
    Fitting L/S model and finding priors
    
    Finding parametric adjustments
    
    Adjusting the Data
    
    


    Found 23894 genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.


    Using the 'mean only' version of ComBat
    
    Found178batches
    
    Note: one batch has only one sample, setting mean.only=TRUE
    
    Adjusting for1covariate(s) or covariate level(s)
    
    Standardizing Data across genes
    
    Fitting L/S model and finding priors
    
    Finding parametric adjustments
    
    Adjusting the Data
    
    



    
![png](output_129_3.png)
    



```R
head(combat4_plate_and_cancerType)
```


<table class="dataframe">
<caption>A matrix: 6 × 6981 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>TCGA-BLCA.TCGA-CU-A3KJ-01A-11R-A21D-07</th><th scope=col>TCGA-BLCA.TCGA-K4-A3WU-01B-11R-A23N-07</th><th scope=col>TCGA-BLCA.TCGA-DK-A3IU-01A-11R-A20F-07</th><th scope=col>TCGA-BLCA.TCGA-GV-A40G-01A-11R-A23N-07</th><th scope=col>TCGA-BLCA.TCGA-DK-A3IN-01A-11R-A20F-07</th><th scope=col>TCGA-BLCA.TCGA-SY-A9G0-01A-12R-A38B-07</th><th scope=col>TCGA-BLCA.TCGA-XF-A9SU-01A-31R-A39I-07</th><th scope=col>TCGA-BLCA.TCGA-UY-A78O-01A-12R-A33J-07</th><th scope=col>TCGA-BLCA.TCGA-GD-A2C5-01A-12R-A180-07</th><th scope=col>TCGA-BLCA.TCGA-G2-A2EK-01A-22R-A18C-07</th><th scope=col>⋯</th><th scope=col>TCGA-UCEC.TCGA-E6-A1LX-01A-11R-A14D-07</th><th scope=col>TCGA-UCEC.TCGA-FI-A3PV-01A-11R-A22K-07</th><th scope=col>TCGA-UCEC.TCGA-BG-A0MT-01A-11R-A104-07</th><th scope=col>TCGA-UCEC.TCGA-EY-A549-01A-11R-A27V-07</th><th scope=col>TCGA-UCEC.TCGA-D1-A1O7-01A-11R-A14D-07</th><th scope=col>TCGA-UCEC.TCGA-AJ-A3NF-01A-11R-A22K-07</th><th scope=col>TCGA-UCEC.TCGA-BK-A0C9-01A-11R-A00V-07</th><th scope=col>TCGA-UCEC.TCGA-AJ-A5DW-01A-11R-A27V-07</th><th scope=col>TCGA-UCEC.TCGA-AP-A05J-01A-11R-A00V-07</th><th scope=col>TCGA-UCEC.TCGA-D1-A103-01A-11R-A10J-07</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ENSG00000000003.15</th><td> 5.240923</td><td> 7.2829385</td><td> 4.080372</td><td> 6.673164</td><td> 5.958337</td><td> 5.177666</td><td> 6.179322</td><td> 7.010650</td><td> 6.616479</td><td> 7.096513</td><td>⋯</td><td> 4.841792</td><td> 6.408050</td><td> 6.291513</td><td> 6.016058</td><td> 6.377131</td><td> 6.614543</td><td> 5.244857</td><td> 6.919921</td><td> 6.07314351</td><td> 5.755877</td></tr>
	<tr><th scope=row>ENSG00000000005.6</th><td>-3.824657</td><td>-0.3473206</td><td>-3.279875</td><td>-2.511472</td><td>-3.895275</td><td>-1.358143</td><td>-2.385796</td><td>-3.388396</td><td>-3.229519</td><td>-3.757873</td><td>⋯</td><td>-3.267034</td><td>-4.490457</td><td>-2.576900</td><td>-3.168572</td><td>-3.400200</td><td>-3.290543</td><td>-3.435438</td><td>-1.692732</td><td>-0.08378953</td><td>-3.408239</td></tr>
	<tr><th scope=row>ENSG00000000419.13</th><td> 5.424674</td><td> 4.7260632</td><td> 5.258813</td><td> 5.317270</td><td> 4.289103</td><td> 5.061591</td><td> 5.211746</td><td> 5.706218</td><td> 5.137266</td><td> 5.252700</td><td>⋯</td><td> 5.205149</td><td> 5.504302</td><td> 5.224252</td><td> 5.331553</td><td> 4.754983</td><td> 5.266489</td><td> 4.022501</td><td> 5.646033</td><td> 4.59463864</td><td> 4.676241</td></tr>
	<tr><th scope=row>ENSG00000000457.14</th><td> 3.633643</td><td> 4.4720291</td><td> 3.688951</td><td> 4.032470</td><td> 3.759512</td><td> 3.759937</td><td> 2.685845</td><td> 3.600412</td><td> 4.433219</td><td> 4.158163</td><td>⋯</td><td> 3.675643</td><td> 3.714771</td><td> 4.280032</td><td> 4.656581</td><td> 3.137272</td><td> 3.868885</td><td> 3.923686</td><td> 4.196688</td><td> 3.07401623</td><td> 3.748381</td></tr>
	<tr><th scope=row>ENSG00000000460.17</th><td> 3.306594</td><td> 3.4591858</td><td> 2.836821</td><td> 2.404045</td><td> 3.645799</td><td> 2.402227</td><td> 1.985382</td><td> 2.432520</td><td> 3.332333</td><td> 2.295960</td><td>⋯</td><td> 2.886295</td><td> 2.708516</td><td> 3.543275</td><td> 2.128827</td><td> 2.370223</td><td> 3.050370</td><td> 2.978555</td><td> 3.674334</td><td> 2.53976206</td><td> 2.554597</td></tr>
	<tr><th scope=row>ENSG00000000938.13</th><td> 2.401667</td><td> 3.3609718</td><td> 5.237980</td><td> 1.887517</td><td> 2.881078</td><td> 2.961823</td><td> 3.019278</td><td> 2.208354</td><td> 2.428211</td><td> 2.545943</td><td>⋯</td><td> 1.987255</td><td> 1.250138</td><td> 2.842889</td><td> 3.571543</td><td> 3.372089</td><td> 3.983571</td><td> 3.724744</td><td> 2.745188</td><td> 1.83194884</td><td> 3.279948</td></tr>
</tbody>
</table>




```R
head(combat_pre3)
```


<table class="dataframe">
<caption>A matrix: 6 × 6981 of type dbl</caption>
<thead>
	<tr><th scope=col>TCGA-BLCA.TCGA-CU-A3KJ-01A-11R-A21D-07</th><th scope=col>TCGA-BLCA.TCGA-K4-A3WU-01B-11R-A23N-07</th><th scope=col>TCGA-BLCA.TCGA-DK-A3IU-01A-11R-A20F-07</th><th scope=col>TCGA-BLCA.TCGA-GV-A40G-01A-11R-A23N-07</th><th scope=col>TCGA-BLCA.TCGA-DK-A3IN-01A-11R-A20F-07</th><th scope=col>TCGA-BLCA.TCGA-SY-A9G0-01A-12R-A38B-07</th><th scope=col>TCGA-BLCA.TCGA-XF-A9SU-01A-31R-A39I-07</th><th scope=col>TCGA-BLCA.TCGA-UY-A78O-01A-12R-A33J-07</th><th scope=col>TCGA-BLCA.TCGA-GD-A2C5-01A-12R-A180-07</th><th scope=col>TCGA-BLCA.TCGA-G2-A2EK-01A-22R-A18C-07</th><th scope=col>⋯</th><th scope=col>TCGA-UCEC.TCGA-E6-A1LX-01A-11R-A14D-07</th><th scope=col>TCGA-UCEC.TCGA-FI-A3PV-01A-11R-A22K-07</th><th scope=col>TCGA-UCEC.TCGA-BG-A0MT-01A-11R-A104-07</th><th scope=col>TCGA-UCEC.TCGA-EY-A549-01A-11R-A27V-07</th><th scope=col>TCGA-UCEC.TCGA-D1-A1O7-01A-11R-A14D-07</th><th scope=col>TCGA-UCEC.TCGA-AJ-A3NF-01A-11R-A22K-07</th><th scope=col>TCGA-UCEC.TCGA-BK-A0C9-01A-11R-A00V-07</th><th scope=col>TCGA-UCEC.TCGA-AJ-A5DW-01A-11R-A27V-07</th><th scope=col>TCGA-UCEC.TCGA-AP-A05J-01A-11R-A00V-07</th><th scope=col>TCGA-UCEC.TCGA-D1-A103-01A-11R-A10J-07</th></tr>
</thead>
<tbody>
	<tr><td> 5.238161</td><td> 7.2387378</td><td> 4.069827</td><td> 6.623498</td><td> 5.964624</td><td> 5.149559</td><td> 6.208320</td><td> 6.973849</td><td> 6.667808</td><td> 7.118220</td><td>⋯</td><td> 4.648308</td><td> 6.538783</td><td> 6.522158</td><td> 5.847490</td><td> 6.249317</td><td> 6.754108</td><td> 5.516100</td><td> 6.790012</td><td>6.3798144</td><td> 5.660345</td></tr>
	<tr><td>-3.970997</td><td>-0.6310536</td><td>-3.346015</td><td>-2.764622</td><td>-3.952718</td><td>-1.741980</td><td>-2.339974</td><td>-3.254568</td><td>-2.916656</td><td>-3.946919</td><td>⋯</td><td>-3.211319</td><td>-4.237089</td><td>-2.608840</td><td>-3.811578</td><td>-3.345573</td><td>-3.027368</td><td>-2.940664</td><td>-2.323675</td><td>0.4383793</td><td>-3.364267</td></tr>
	<tr><td> 5.404524</td><td> 4.6698790</td><td> 5.251353</td><td> 5.271146</td><td> 4.265141</td><td> 4.944996</td><td> 5.252229</td><td> 5.669376</td><td> 5.147714</td><td> 5.239577</td><td>⋯</td><td> 5.268608</td><td> 5.600456</td><td> 5.128976</td><td> 5.261239</td><td> 4.845627</td><td> 5.377004</td><td> 4.053155</td><td> 5.556729</td><td>4.5907421</td><td> 4.786749</td></tr>
	<tr><td> 3.645334</td><td> 4.4369364</td><td> 3.698872</td><td> 3.998578</td><td> 3.769240</td><td> 3.671195</td><td> 2.676304</td><td> 3.591198</td><td> 4.477039</td><td> 4.137900</td><td>⋯</td><td> 3.704297</td><td> 3.755691</td><td> 4.312535</td><td> 4.504518</td><td> 3.157500</td><td> 3.912217</td><td> 3.909974</td><td> 4.037427</td><td>3.0470054</td><td> 3.786487</td></tr>
	<tr><td> 3.219324</td><td> 3.2630730</td><td> 2.806984</td><td> 2.200830</td><td> 3.621407</td><td> 2.287604</td><td> 2.051200</td><td> 2.399760</td><td> 3.353775</td><td> 2.153414</td><td>⋯</td><td> 2.884372</td><td> 2.798618</td><td> 3.518298</td><td> 1.898284</td><td> 2.371175</td><td> 3.138568</td><td> 3.107253</td><td> 3.435181</td><td>2.6709042</td><td> 2.561152</td></tr>
	<tr><td> 2.361188</td><td> 3.3817664</td><td> 5.259571</td><td> 1.897793</td><td> 2.885844</td><td> 2.715382</td><td> 3.058456</td><td> 2.160464</td><td> 2.505857</td><td> 2.533101</td><td>⋯</td><td> 2.035387</td><td> 1.330315</td><td> 2.886993</td><td> 3.236559</td><td> 3.435660</td><td> 4.094222</td><td> 3.632679</td><td> 2.400991</td><td>1.7187827</td><td> 3.324470</td></tr>
</tbody>
</table>




```R
errs<-(combat_pre3-combat4_plate_and_cancerType)
```


```R
max(combat4_plate_and_cancerType)
```


42.5814234634876



```R
max(errs)
```


10.3702441344789



```R
min(errs)
```


-9.96790634120659



```R
summary(combat4_plate_and_cancerType[,1])
```


       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    -5.5967 -4.6297 -3.7541 -1.7503  0.2126 12.5609 



```R
summary(combat_pre3[,1])
```


       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    -5.9419 -4.6336 -3.7858 -1.7612  0.2064 12.6339 



```R
summary(errs[,1])
```


          Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
    -0.8779105 -0.0399806 -0.0000189 -0.0108851  0.0002207  1.5615910 



```R
# Here a pca again :P for data after only camner type batch effect removal :
pca_unscaled_uncenterd_combat_CancerType <- prcomp(t(combat_pre2), scale = F, center = F)
plot_PCAs(as.data.table(pca_unscaled_uncenterd_combat_CancerType$x), batchs)
```


    
![png](output_139_0.png)
    


#### PCAs after normalization & Batcheffect

###### PCA of row data:


```R
plot_PCAs(pca_countX, batchs)
```


    
![png](output_142_0.png)
    


###### PCA of normalized data:


```R
plot_PCAs(pca_normalizedX, batchs)
```


    
![png](output_144_0.png)
    


###### PCA of normalized+combat data, without scaling or centering:


```R
plot_PCAs(as.data.table(pca_unscaled_uncenterd_combat4_plate_and_cancerType$x), batchs)
```


    
![png](output_146_0.png)
    


###### PCA of normalized+combat data, per cancer:


```R
class(logCPM_list$`TCGA-BLCA`)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'data.table'</li><li>'data.frame'</li></ol>




```R
class(combat4_plate_and_cancerType)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'matrix'</li><li>'array'</li></ol>




```R
head(combat4_plate_and_cancerType)
```


<table class="dataframe">
<caption>A matrix: 6 × 6981 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>TCGA-BLCA.TCGA-CU-A3KJ-01A-11R-A21D-07</th><th scope=col>TCGA-BLCA.TCGA-K4-A3WU-01B-11R-A23N-07</th><th scope=col>TCGA-BLCA.TCGA-DK-A3IU-01A-11R-A20F-07</th><th scope=col>TCGA-BLCA.TCGA-GV-A40G-01A-11R-A23N-07</th><th scope=col>TCGA-BLCA.TCGA-DK-A3IN-01A-11R-A20F-07</th><th scope=col>TCGA-BLCA.TCGA-SY-A9G0-01A-12R-A38B-07</th><th scope=col>TCGA-BLCA.TCGA-XF-A9SU-01A-31R-A39I-07</th><th scope=col>TCGA-BLCA.TCGA-UY-A78O-01A-12R-A33J-07</th><th scope=col>TCGA-BLCA.TCGA-GD-A2C5-01A-12R-A180-07</th><th scope=col>TCGA-BLCA.TCGA-G2-A2EK-01A-22R-A18C-07</th><th scope=col>⋯</th><th scope=col>TCGA-UCEC.TCGA-E6-A1LX-01A-11R-A14D-07</th><th scope=col>TCGA-UCEC.TCGA-FI-A3PV-01A-11R-A22K-07</th><th scope=col>TCGA-UCEC.TCGA-BG-A0MT-01A-11R-A104-07</th><th scope=col>TCGA-UCEC.TCGA-EY-A549-01A-11R-A27V-07</th><th scope=col>TCGA-UCEC.TCGA-D1-A1O7-01A-11R-A14D-07</th><th scope=col>TCGA-UCEC.TCGA-AJ-A3NF-01A-11R-A22K-07</th><th scope=col>TCGA-UCEC.TCGA-BK-A0C9-01A-11R-A00V-07</th><th scope=col>TCGA-UCEC.TCGA-AJ-A5DW-01A-11R-A27V-07</th><th scope=col>TCGA-UCEC.TCGA-AP-A05J-01A-11R-A00V-07</th><th scope=col>TCGA-UCEC.TCGA-D1-A103-01A-11R-A10J-07</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ENSG00000000003.15</th><td> 5.240923</td><td> 7.2829385</td><td> 4.080372</td><td> 6.673164</td><td> 5.958337</td><td> 5.177666</td><td> 6.179322</td><td> 7.010650</td><td> 6.616479</td><td> 7.096513</td><td>⋯</td><td> 4.841792</td><td> 6.408050</td><td> 6.291513</td><td> 6.016058</td><td> 6.377131</td><td> 6.614543</td><td> 5.244857</td><td> 6.919921</td><td> 6.07314351</td><td> 5.755877</td></tr>
	<tr><th scope=row>ENSG00000000005.6</th><td>-3.824657</td><td>-0.3473206</td><td>-3.279875</td><td>-2.511472</td><td>-3.895275</td><td>-1.358143</td><td>-2.385796</td><td>-3.388396</td><td>-3.229519</td><td>-3.757873</td><td>⋯</td><td>-3.267034</td><td>-4.490457</td><td>-2.576900</td><td>-3.168572</td><td>-3.400200</td><td>-3.290543</td><td>-3.435438</td><td>-1.692732</td><td>-0.08378953</td><td>-3.408239</td></tr>
	<tr><th scope=row>ENSG00000000419.13</th><td> 5.424674</td><td> 4.7260632</td><td> 5.258813</td><td> 5.317270</td><td> 4.289103</td><td> 5.061591</td><td> 5.211746</td><td> 5.706218</td><td> 5.137266</td><td> 5.252700</td><td>⋯</td><td> 5.205149</td><td> 5.504302</td><td> 5.224252</td><td> 5.331553</td><td> 4.754983</td><td> 5.266489</td><td> 4.022501</td><td> 5.646033</td><td> 4.59463864</td><td> 4.676241</td></tr>
	<tr><th scope=row>ENSG00000000457.14</th><td> 3.633643</td><td> 4.4720291</td><td> 3.688951</td><td> 4.032470</td><td> 3.759512</td><td> 3.759937</td><td> 2.685845</td><td> 3.600412</td><td> 4.433219</td><td> 4.158163</td><td>⋯</td><td> 3.675643</td><td> 3.714771</td><td> 4.280032</td><td> 4.656581</td><td> 3.137272</td><td> 3.868885</td><td> 3.923686</td><td> 4.196688</td><td> 3.07401623</td><td> 3.748381</td></tr>
	<tr><th scope=row>ENSG00000000460.17</th><td> 3.306594</td><td> 3.4591858</td><td> 2.836821</td><td> 2.404045</td><td> 3.645799</td><td> 2.402227</td><td> 1.985382</td><td> 2.432520</td><td> 3.332333</td><td> 2.295960</td><td>⋯</td><td> 2.886295</td><td> 2.708516</td><td> 3.543275</td><td> 2.128827</td><td> 2.370223</td><td> 3.050370</td><td> 2.978555</td><td> 3.674334</td><td> 2.53976206</td><td> 2.554597</td></tr>
	<tr><th scope=row>ENSG00000000938.13</th><td> 2.401667</td><td> 3.3609718</td><td> 5.237980</td><td> 1.887517</td><td> 2.881078</td><td> 2.961823</td><td> 3.019278</td><td> 2.208354</td><td> 2.428211</td><td> 2.545943</td><td>⋯</td><td> 1.987255</td><td> 1.250138</td><td> 2.842889</td><td> 3.571543</td><td> 3.372089</td><td> 3.983571</td><td> 3.724744</td><td> 2.745188</td><td> 1.83194884</td><td> 3.279948</td></tr>
</tbody>
</table>



###### PCA of normalized+combat data per cancer type:


```R
#preparation of combat4_plate_and_cancerType_list for the pca plot per cancer
combat4_plate_and_cancerType_list<-combat4_plate_and_cancerType %>% t %>% as.data.table %>% mutate(project_id = batchs$project_id)%>% dplyr::group_split(project_id, .keep = FALSE)
names(combat4_plate_and_cancerType_list)<- names(countdata_list)
combat4_plate_and_cancerType_list<-lapply(combat4_plate_and_cancerType_list,t)
combat4_plate_and_cancerType_list<-lapply(names(countdata_list), function(n){
    colnames(combat4_plate_and_cancerType_list[[n]]) <- colnames(countdata_list[[n]][,-1])
    return(combat4_plate_and_cancerType_list[[n]])
    })
names(combat4_plate_and_cancerType_list)<- names(countdata_list)
```


```R
pca_logCPM_combat <- prcomp_list(combat4_plate_and_cancerType_list,batchs) %>% do.call(what = rbind, args= .)
plot_pca_per_cancer(pca_logCPM_combat, color_by= pca_logCPM_combat$sample_type)+ 
    scale_color_manual(values = setNames(batchs$c_by_type, batchs$sample_type))
```


    
![png](output_153_0.png)
    


### Other exploratory plots using different packages for Principal component analysis (Just or tests: to be removed)


```R
res.pca <- PCA(t(combat4_plate_and_cancerType), graph = FALSE, ncp = 4, scale.unit = F)
```


```R
res.pca_befor <- PCA(t(merged_normalized_data[,-1]), graph = FALSE, ncp = 4, scale.unit = F)
```


```R
#before:
fviz_pca_ind(res.pca_befor, geom = "point",col.ind = batchs$project_id, addEllipses = T)
```


    
![png](output_157_0.png)
    



```R
fviz_pca_ind(res.pca, geom = "point",col.ind = batchs$project_id, addEllipses = T)
```


    
![png](output_158_0.png)
    



```R
#before
fviz_pca_ind(res.pca_befor, geom = "point",col.ind = batchs$groups, addEllipses = T)
```


    
![png](output_159_0.png)
    



```R
fviz_pca_ind(res.pca, geom = "point",col.ind = batchs$groups, addEllipses = T)
```


    
![png](output_160_0.png)
    



```R
mypca = mixOmics::pca(t(combat4_plate_and_cancerType), ncomp = 5, center = FALSE, scale = FALSE)
```


```R
plot(mypca)
```


    
![png](output_162_0.png)
    



```R
p1<-plotIndiv(mypca, comp = 1:2, 
          #ind.names = batchs$participant,
          pch= batchs$groups,
          group = batchs$project_id,
          # graphical parameters
          col = unique(batchs$c_by_cancer), style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Cancer type", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)
```


    
![png](output_163_0.png)
    



```R
p2<-plotIndiv(mypca, comp = 1:2, 
          #ind.names = batchs$participant,
          pch= batchs$groups,
          group = batchs$sample_type,
          # graphical parameters
          col = unique(batchs$c_by_type), style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Center", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)
```


    
![png](output_164_0.png)
    



```R
p3<-plotIndiv(mypca, comp = 1:2, 
          #ind.names = batchs$participant,
          pch= batchs$groups,
          group = batchs$plate,
          # graphical parameters
          col = unique(batchs$c_by_plate), 
              style = "ggplot2", 
          legend = F, legend.position = "right", 
          legend.title.pch = "Plate", ellipse = F, 
          ellipse.level = 0.95, centroid = FALSE)
```


    
![png](output_165_0.png)
    



```R
p4<-plotIndiv(mypca, comp = 1:2, 
          #ind.names = batchs$participant,
          pch= batchs$groups,
          group = batchs$tss,
          # graphical parameters
          col = unique(batchs$c_by_tss), style = "ggplot2", 
          legend = F, legend.position = "right", 
          legend.title = "Center", ellipse = F, 
          ellipse.level = 0.95, centroid = FALSE)
```


    
![png](output_166_0.png)
    


## Data preprocessing for breast cancer study

### Data selection


```R
#check if plat6e used in BRCA are used in other cancer in these data 
```


```R
platesPerCancers <- batchs %>% dplyr::select(c(plate,project_id)) %>% unique %>% group_by(plate) %>% 
    filter(duplicated(plate) | duplicated(plate, fromLast = T))%>% 
    mutate(usedInCancers = n()) %>% arrange(desc(usedInCancers), plate,project_id, ) 
head(platesPerCancers)
```


<table class="dataframe">
<caption>A grouped_df: 6 × 3</caption>
<thead>
	<tr><th scope=col>plate</th><th scope=col>project_id</th><th scope=col>usedInCancers</th></tr>
	<tr><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>A180</td><td>TCGA-BLCA</td><td>6</td></tr>
	<tr><td>A180</td><td>TCGA-BRCA</td><td>6</td></tr>
	<tr><td>A180</td><td>TCGA-COAD</td><td>6</td></tr>
	<tr><td>A180</td><td>TCGA-LIHC</td><td>6</td></tr>
	<tr><td>A180</td><td>TCGA-THCA</td><td>6</td></tr>
	<tr><td>A180</td><td>TCGA-UCEC</td><td>6</td></tr>
</tbody>
</table>




```R
message("The number of plate that are used for more than one type of cancer is : ",
platesPerCancers %>% dplyr::select(c(plate, usedInCancers)) %>% unique %>% dim %>% .[1], " Out of a totle of :",
batchs$plate %>% unique %>% length)
```

    The number of plate that are used for more than one type of cancer is : 78 Out of a totle of :178
    



```R
message("The number of plate that are used for BRCA and more than one type of other cancer is : ",
platesPerCancers %>% filter(project_id=="TCGA-BRCA") %>% unique %>% dim %>% .[1], " Out of a totle of :",
batchs %>% filter(project_id=="TCGA-BRCA") %>% .$plate %>% unique %>% length)
```

    The number of plate that are used for BRCA and more than one type of other cancer is : 33 Out of a totle of :45
    



```R
# Filter the data keeping only  paired data 113 samples 
brca_count <- countdata_list$`TCGA-BRCA`
brca_batchs <- batchs %>% filter(project_id =="TCGA-BRCA")
paired_patients <- brca_batchs %>% dplyr::select(c(participant, groups)) %>% unique() %>% filter(duplicated(participant)) %>% .$participant
```


```R
brca_batchs_paired <- brca_batchs %>% filter(participant %in% paired_patients)
```


```R
brca_batchs_paired %>% filter(participant %in% as.character(brca_batchs_paired %>% filter(sample_type =="Metastatic") %>% .$participant) &
                             groups =="Tumor") %>% dplyr::select( barcode, sample_type, participant) %>% arrange(participant)
```


<table class="dataframe">
<caption>A data.table: 6 × 3</caption>
<thead>
	<tr><th></th><th scope=col>barcode</th><th scope=col>sample_type</th><th scope=col>participant</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>TCGA-E2-A15K-01A-11R-A12P-07</th><td>TCGA-E2-A15K-01A-11R-A12P-07</td><td>Primary Tumor</td><td>A15K</td></tr>
	<tr><th scope=row>TCGA-E2-A15K-06A-11R-A12P-07</th><td>TCGA-E2-A15K-06A-11R-A12P-07</td><td>Metastatic   </td><td>A15K</td></tr>
	<tr><th scope=row>TCGA-BH-A18V-01A-11R-A12D-07</th><td>TCGA-BH-A18V-01A-11R-A12D-07</td><td>Primary Tumor</td><td>A18V</td></tr>
	<tr><th scope=row>TCGA-BH-A18V-06A-11R-A213-07</th><td>TCGA-BH-A18V-06A-11R-A213-07</td><td>Metastatic   </td><td>A18V</td></tr>
	<tr><th scope=row>TCGA-BH-A1FE-01A-11R-A13Q-07</th><td>TCGA-BH-A1FE-01A-11R-A13Q-07</td><td>Primary Tumor</td><td>A1FE</td></tr>
	<tr><th scope=row>TCGA-BH-A1FE-06A-11R-A213-07</th><td>TCGA-BH-A1FE-06A-11R-A213-07</td><td>Metastatic   </td><td>A1FE</td></tr>
</tbody>
</table>



>As we can see in the table above, the metastatic patients also have samples from the primary tumor. **I will remove these samples with a metastatic profile**, as they introduce noise in the data (it was showed by a PCA), and there are only three such samples. We will keep other for these patients with primary tumor type


```R
brca_batchs_paired <- brca_batchs_paired %>% filter(sample_type !="Metastatic")
```


```R
brca_count_paired <- brca_count %>% dplyr::select(all_of(c("rn",brca_batchs_paired$barcode)))
```


```R
brca_batchs_paired <-droplevels(brca_batchs_paired)
```


```R
brca_batchs_paired$sample_type %>%table
```


    .
          Primary Tumor Solid Tissue Normal 
                    119                 113 



```R
mypca_brca = mixOmics::pca(t(brca_count_paired[,-1]), ncomp = 5, center = FALSE, scale = FALSE)
```


```R
plot(mypca_brca)
```


    
![png](output_182_0.png)
    


>PC1 captures the majority of variance (>60%), indicating a strong underlying factor influencing the dataset. PC2 explains around 10%, with subsequent components (PC3 and beyond) contributing progressively less. This suggests that most of the variability in the BRCA cancer data can be captured by the first two principal components, while further PCs add minimal additional variance.


```R
p1<-plotIndiv(mypca_brca, comp = 1:2, col.per.group = levels(brca_batchs_paired$c_by_group), ind.names = 
              FALSE,
          group = brca_batchs_paired$groups,
          # graphical parameters
              style = "ggplot2",
          legend = TRUE, legend.position = "right", 
          legend.title = "Cancer type", 
          legend.title.pch = 'Sample types',              
          legend.pch = FALSE,  # Show point shape in the legend
          ellipse = TRUE, 
          ellipse.level = 0.95)
```


    
![png](output_184_0.png)
    



```R
p1<-plotIndiv(mypca_brca, comp = 1:2, col.per.group = levels(brca_batchs_paired$c_by_plate),
          pch= brca_batchs_paired$groups,
          group = brca_batchs_paired$plate,
          # graphical parameters
              style = "ggplot2",
          legend = TRUE, legend.position = "right", 
          legend.title = "Cancer type", 
          legend.title.pch = 'Sample types',
          legend.pch = FALSE,  # Show point shape in the legend
          ellipse = TRUE, 
          ellipse.level = 0.95)
```


    
![png](output_185_0.png)
    



```R
pca_unscaled_uncenterd_brca_count <- prcomp(t(brca_count_paired[,-1]), scale = F, center = F)
plot_PCAs(as.data.table(pca_unscaled_uncenterd_brca_count$x), brca_batchs_paired)
```


    
![png](output_186_0.png)
    


### Quality control
#### Normalization TMM


```R
# Normalization
dge_brca <- do_normalization(brca_count_paired,brca_batchs_paired,method ="TMM", group_col = "sample_type") 
# Log-transformation
logCPM_brca <- cpm(dge_brca, log=TRUE)
```


```R
rownames(logCPM_brca)<- dge_brca$rn
```


```R
logCPM_brca <- as.data.table(logCPM_brca, keep.rownames = T)
```

#### Outliers detection


```R
# Boxplot of Randomly selected gene expression across sample types (Tumor vs Normal)
# (I have changed the seed and run it several time ...)
set.seed(261090) # For reproducibility
selected_genes <- sample(logCPM_brca$rn, 6)
# Reshape for ggplot (convert to long format)
logCPM_brca_long <- melt(logCPM_brca,id.vars = "rn", variable.name = "Sample", value.name = "Expression") %>%
    mutate(sample_type= brca_batchs_paired$sample_type[Sample])

ggplot(logCPM_brca_long %>% filter(rn %in% selected_genes) , aes(x = rn, y = Expression, fill = sample_type)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 8, outlier.size = 1) +
  labs(title = "Boxplot of Gene Expression Across Sample Types (Tumor vs Normal)",
       x = "Gene",
       y = "Log-Transformed Expression") +
  facet_wrap(~rn, scales = "free_x") +  # Divide the plot into Tumor and Normal panels
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate gene names for readability
```


    
![png](output_192_0.png)
    


#### Trunction at the 99th percentile: a mild way of reducing outliers impact

Truncation at the 99th percentile can partially address outlier removal, but it’s not strictly considered the same as formal outlier removal.
To remove outliers from the data, we can use the Interquartile Range (IQR) method or other approaches to filter out values considered outliers. Typically, outliers are defined as values that are above Q3+1.5×IQRQ3+1.5×IQR or below Q1−1.5×IQRQ1−1.5×IQR, where Q1 and Q3 are the first and third quartiles, respectively, and IQR is the interquartile range (Q3−Q1Q3−Q1).Generaly, outlier removal technics involves completely excluding values that are significantly different from the rest of the data (using statistical methods or visual inspection).
In other hand, truncation at the 99th percentile means capping all values above the 99th percentile to the value at the 99th percentile, which reduces the influence of extreme values. It does not remove the data points entirely, but instead limits their impact on further analysis.
For this breast cancer data, we will **truncate the values at the 99th percentile** as a mild way of reducing the impact of extreme values, without treating them as outliers or eliminating them as the IQR method does.


```R
logCPM_brca_99 <- truncate_at_percentile(logCPM_brca[, -1], truncation_percentil = 0.99)
rownames(logCPM_brca_99) <- logCPM_brca$rn
logCPM_brca_99 <- as.data.table(logCPM_brca_99, keep.rownames = T)
```


```R
# Reshape for ggplot (convert to long format)
logCPM_brca_99_long <- melt(logCPM_brca_99,id.vars = "rn", variable.name = "Sample", value.name = "Expression") %>%
    mutate(sample_type= brca_batchs_paired$sample_type[Sample])
# Create a boxplot for the same selected genes divided by sample type after truncation
ggplot(logCPM_brca_99_long %>% filter(rn %in% selected_genes) , aes(x = sample_type, y = Expression, fill = sample_type)) +
  geom_boxplot(outlier.shape = NA) +  # No outliers displayed in the plot
  labs(title = "Boxplot of Selected Genes (6) - Normal vs Tumor (Truncated at 99th Percentile)",
       x = "Sample Group",
       y = "Log-Transformed CPM") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ rn)  # Separate plots for each gene
```


    
![png](output_196_0.png)
    



```R
pca_unscaled_uncenterd_brca_normalized99 <- prcomp(t(logCPM_brca_99[,-1]), scale = F, center = F)
plot_PCAs(as.data.table(pca_unscaled_uncenterd_brca_normalized99$x), brca_batchs_paired)
```


    
![png](output_197_0.png)
    


The figure shows four PCA plots:
- **(A)** Colored by cancer type, all points are black because we are only analyzing BRCA cancer.
- **(B)** Colored by plate, with no clear batch effect easily detectable by eye.
- **(C)** Colored by TSS, showing no clear batch effect, but it can be observed that the 8 outlier samples all come from the same tissue site.
- **(D)** The PCA plot shows that tumor and control samples are well clustered together, indicating a clear separation between the groups. However, 8 tumor samples were detected as outliers in PC1, positioned far from the main data cluster (above -1000), suggesting potential data anomalies or batch effects that may require further investigation.


```R
wiered_samples <- rownames(pca_unscaled_uncenterd_brca_normalized99$x)[which(pca_unscaled_uncenterd_brca_normalized99$x[,1] > (-1000))]
wiered_samples
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'TCGA-A7-A13E-01B-06R-A277-07'</li><li>'TCGA-A7-A13E-01A-11R-A277-07'</li><li>'TCGA-A7-A0DB-01A-11R-A277-07'</li><li>'TCGA-A7-A0DB-01C-02R-A277-07'</li><li>'TCGA-A7-A13G-01A-11R-A13Q-07'</li><li>'TCGA-A7-A13G-01B-04R-A22O-07'</li><li>'TCGA-A7-A0DC-01B-04R-A22O-07'</li><li>'TCGA-A7-A0DC-01A-11R-A00Z-07'</li></ol>




```R
# lets check their behavior in the normalized data 
logCPM_brca %>% dplyr::select(all_of(wiered_samples)) %>% head
```


<table class="dataframe">
<caption>A data.table: 6 × 8</caption>
<thead>
	<tr><th scope=col>TCGA-A7-A13E-01B-06R-A277-07</th><th scope=col>TCGA-A7-A13E-01A-11R-A277-07</th><th scope=col>TCGA-A7-A0DB-01A-11R-A277-07</th><th scope=col>TCGA-A7-A0DB-01C-02R-A277-07</th><th scope=col>TCGA-A7-A13G-01A-11R-A13Q-07</th><th scope=col>TCGA-A7-A13G-01B-04R-A22O-07</th><th scope=col>TCGA-A7-A0DC-01B-04R-A22O-07</th><th scope=col>TCGA-A7-A0DC-01A-11R-A00Z-07</th></tr>
	<tr><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td> 4.699873</td><td> 5.0921496</td><td>5.7039145</td><td> 4.5904405</td><td> 5.591442</td><td> 3.93900229</td><td> 3.807059</td><td> 6.014709</td></tr>
	<tr><td>-1.101221</td><td>-0.7961921</td><td>0.1991546</td><td>-0.2220211</td><td>-1.376363</td><td>-0.06550389</td><td>-1.853107</td><td>-1.966826</td></tr>
	<tr><td> 4.189540</td><td> 5.4159087</td><td>4.7442972</td><td> 3.2622698</td><td> 4.421612</td><td> 3.05446263</td><td> 3.837242</td><td> 5.584697</td></tr>
	<tr><td> 4.332984</td><td> 5.1623523</td><td>5.3688905</td><td> 5.0701664</td><td> 6.247269</td><td> 5.67310459</td><td> 4.282318</td><td> 4.735253</td></tr>
	<tr><td> 3.862453</td><td> 4.3591536</td><td>4.0085753</td><td> 3.9825512</td><td> 4.201083</td><td> 3.51975195</td><td> 2.733943</td><td> 3.650863</td></tr>
	<tr><td> 1.416184</td><td> 1.7480988</td><td>4.3290317</td><td> 1.2218983</td><td> 1.490362</td><td> 2.16044596</td><td> 2.170087</td><td> 0.718840</td></tr>
</tbody>
</table>




```R
# Create a boxplot for the same selected genes divided by sample type after truncation
ggplot(logCPM_brca_99_long %>% filter(rn %in% selected_genes & sample_type == "Primary Tumor") %>% mutate(weird = Sample %in% wiered_samples) ,
       aes(x = sample_type, y = Expression, fill = weird)) +
  geom_boxplot(outlier.shape = NA) +  # No outliers displayed in the plot
  labs(title = "Boxplot of randomly Selected Genes (6) - weird Tumor vs Tumor (Truncated at 99th Percentile)",
       x = "Sample Group",
       y = "Log-Transformed CPM") +
  theme_minimal() +
  facet_wrap(~ rn)  # Separate plots for each gene
```


    
![png](output_201_0.png)
    


**⚠️ Note:**
We have noticed earlier that some samples are duplicated. Let's revisit and check if there is an intersection between the duplicated tumor samples and these weird tumor samples.


```R
brca_batchs_paired %>% filter(groups =="Tumor" ) %>% 
    filter(duplicated(participant)| duplicated(participant,fromLast = T)) %>% arrange(participant) %>% .$participant %>% unique
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>A0DB</li><li>A0DC</li><li>A13E</li><li>A13G</li></ol>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'A0AU'</li><li>'A0AY'</li><li>'A0AZ'</li><li>'A0B3'</li><li>'A0B5'</li><li>'A0B7'</li><li>'A0B8'</li><li>'A0BA'</li><li>'A0BC'</li><li>'A0BJ'</li><li>'A0BM'</li><li>'A0BQ'</li><li>'A0BS'</li><li>'A0BT'</li><li>'A0BV'</li><li>'A0BW'</li><li>'A0BZ'</li><li>'A0C0'</li><li>'A0C3'</li><li>'A0CE'</li><li>'A0CH'</li><li>'A0D9'</li><li>'A0DB'</li><li>'A0DC'</li><li>'A0DD'</li><li>'A0DG'</li><li>'A0DH'</li><li>'A0DK'</li><li>'A0DL'</li><li>'A0DO'</li><li>'A0DP'</li><li>'A0DQ'</li><li>'A0DT'</li><li>'A0DV'</li><li>'A0DZ'</li><li>'A0E0'</li><li>'A0E1'</li><li>'A0H5'</li><li>'A0H7'</li><li>'A0H9'</li><li>'A0HA'</li><li>'A0HK'</li><li>'A13E'</li><li>'A13F'</li><li>'A13G'</li><li>'A153'</li><li>'A158'</li><li>'A15I'</li><li>'A15K'</li><li>'A15M'</li><li>'A18J'</li><li>'A18K'</li><li>'A18L'</li><li>'A18M'</li><li>'A18N'</li><li>'A18P'</li><li>'A18Q'</li><li>'A18R'</li><li>'A18S'</li><li>'A18U'</li><li>'A18V'</li><li>'A1BC'</li><li>'A1EN'</li><li>'A1EO'</li><li>'A1ET'</li><li>'A1EU'</li><li>'A1EV'</li><li>'A1EW'</li><li>'A1F0'</li><li>'A1F2'</li><li>'A1F6'</li><li>'A1F8'</li><li>'A1FB'</li><li>'A1FC'</li><li>'A1FD'</li><li>'A1FE'</li><li>'A1FG'</li><li>'A1FH'</li><li>'A1FJ'</li><li>'A1FM'</li><li>'A1FN'</li><li>'A1FR'</li><li>'A1FU'</li><li>'A1IG'</li><li>'A1L7'</li><li>'A1LB'</li><li>'A1LH'</li><li>'A1LS'</li><li>'A1N4'</li><li>'A1N5'</li><li>'A1N6'</li><li>'A1N9'</li><li>'A1NA'</li><li>'A1ND'</li><li>'A1NF'</li><li>'A1NG'</li><li>'A1R7'</li><li>'A1RB'</li><li>'A1RC'</li><li>'A1RD'</li><li>'A1RF'</li><li>'A1RH'</li><li>'A1RI'</li><li>'A203'</li><li>'A204'</li><li>'A208'</li><li>'A209'</li><li>'A23H'</li><li>'A2C8'</li><li>'A2C9'</li><li>'A2FB'</li><li>'A2FF'</li><li>'A2FM'</li></ol>
</details>



```R
wiered_participants <-brca_batchs_paired %>% filter(barcode %in% wiered_samples) %>% arrange(participant) %>% .$participant %>% unique
wiered_participants
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>A0DB</li><li>A0DC</li><li>A13E</li><li>A13G</li></ol>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'A0AU'</li><li>'A0AY'</li><li>'A0AZ'</li><li>'A0B3'</li><li>'A0B5'</li><li>'A0B7'</li><li>'A0B8'</li><li>'A0BA'</li><li>'A0BC'</li><li>'A0BJ'</li><li>'A0BM'</li><li>'A0BQ'</li><li>'A0BS'</li><li>'A0BT'</li><li>'A0BV'</li><li>'A0BW'</li><li>'A0BZ'</li><li>'A0C0'</li><li>'A0C3'</li><li>'A0CE'</li><li>'A0CH'</li><li>'A0D9'</li><li>'A0DB'</li><li>'A0DC'</li><li>'A0DD'</li><li>'A0DG'</li><li>'A0DH'</li><li>'A0DK'</li><li>'A0DL'</li><li>'A0DO'</li><li>'A0DP'</li><li>'A0DQ'</li><li>'A0DT'</li><li>'A0DV'</li><li>'A0DZ'</li><li>'A0E0'</li><li>'A0E1'</li><li>'A0H5'</li><li>'A0H7'</li><li>'A0H9'</li><li>'A0HA'</li><li>'A0HK'</li><li>'A13E'</li><li>'A13F'</li><li>'A13G'</li><li>'A153'</li><li>'A158'</li><li>'A15I'</li><li>'A15K'</li><li>'A15M'</li><li>'A18J'</li><li>'A18K'</li><li>'A18L'</li><li>'A18M'</li><li>'A18N'</li><li>'A18P'</li><li>'A18Q'</li><li>'A18R'</li><li>'A18S'</li><li>'A18U'</li><li>'A18V'</li><li>'A1BC'</li><li>'A1EN'</li><li>'A1EO'</li><li>'A1ET'</li><li>'A1EU'</li><li>'A1EV'</li><li>'A1EW'</li><li>'A1F0'</li><li>'A1F2'</li><li>'A1F6'</li><li>'A1F8'</li><li>'A1FB'</li><li>'A1FC'</li><li>'A1FD'</li><li>'A1FE'</li><li>'A1FG'</li><li>'A1FH'</li><li>'A1FJ'</li><li>'A1FM'</li><li>'A1FN'</li><li>'A1FR'</li><li>'A1FU'</li><li>'A1IG'</li><li>'A1L7'</li><li>'A1LB'</li><li>'A1LH'</li><li>'A1LS'</li><li>'A1N4'</li><li>'A1N5'</li><li>'A1N6'</li><li>'A1N9'</li><li>'A1NA'</li><li>'A1ND'</li><li>'A1NF'</li><li>'A1NG'</li><li>'A1R7'</li><li>'A1RB'</li><li>'A1RC'</li><li>'A1RD'</li><li>'A1RF'</li><li>'A1RH'</li><li>'A1RI'</li><li>'A203'</li><li>'A204'</li><li>'A208'</li><li>'A209'</li><li>'A23H'</li><li>'A2C8'</li><li>'A2C9'</li><li>'A2FB'</li><li>'A2FF'</li><li>'A2FM'</li></ol>
</details>


>We can clearly see that these weird samples come from participants who provided two or more samples from the same tumor tissue sites.


```R
brca_batchs_paired$tss %>% table
```


    .
     A7  AC  BH  E2  E9  GI 
     22   8 146  22  30   4 


>The comun thing between these samples that are coming from A7 tss , according to  [TCGA](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tissue-source-site-codes) is Christiana Healthcare source site and from the study named "Breast invasive carcinoma", with a BCR: NCH.

**Hypothesis about these samples:** The issue may be with these specific participants, and we might need to remove their corresponding control samples as well. Could this be due to a potential batch effect? Let’s take a closer look at these participants.


```R
brca_batchs_paired %>%filter(participant %in%  (brca_batchs_paired %>% filter(barcode %in% wiered_samples) %>% .$participant %>% unique)) %>%
    arrange(plate, participant,sampleVial,portionAnalyte )%>% mutate(weird = barcode %in% wiered_samples) %>% dplyr::select(participant, sampleVial,portionAnalyte,plate, groups, weird)
```


<table class="dataframe">
<caption>A data.table: 14 × 6</caption>
<thead>
	<tr><th></th><th scope=col>participant</th><th scope=col>sampleVial</th><th scope=col>portionAnalyte</th><th scope=col>plate</th><th scope=col>groups</th><th scope=col>weird</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;lgl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>TCGA-A7-A0DB-01A-11R-A00Z-07</th><td>A0DB</td><td>01A</td><td>11R</td><td>A00Z</td><td>Tumor </td><td>FALSE</td></tr>
	<tr><th scope=row>TCGA-A7-A0DC-01A-11R-A00Z-07</th><td>A0DC</td><td>01A</td><td>11R</td><td>A00Z</td><td>Tumor </td><td> TRUE</td></tr>
	<tr><th scope=row>TCGA-A7-A0DB-11A-33R-A089-07</th><td>A0DB</td><td>11A</td><td>33R</td><td>A089</td><td>Normal</td><td>FALSE</td></tr>
	<tr><th scope=row>TCGA-A7-A0DC-11A-41R-A089-07</th><td>A0DC</td><td>11A</td><td>41R</td><td>A089</td><td>Normal</td><td>FALSE</td></tr>
	<tr><th scope=row>TCGA-A7-A13E-01A-11R-A12P-07</th><td>A13E</td><td>01A</td><td>11R</td><td>A12P</td><td>Tumor </td><td>FALSE</td></tr>
	<tr><th scope=row>TCGA-A7-A13E-11A-61R-A12P-07</th><td>A13E</td><td>11A</td><td>61R</td><td>A12P</td><td>Normal</td><td>FALSE</td></tr>
	<tr><th scope=row>TCGA-A7-A13G-01A-11R-A13Q-07</th><td>A13G</td><td>01A</td><td>11R</td><td>A13Q</td><td>Tumor </td><td> TRUE</td></tr>
	<tr><th scope=row>TCGA-A7-A13G-11A-51R-A13Q-07</th><td>A13G</td><td>11A</td><td>51R</td><td>A13Q</td><td>Normal</td><td>FALSE</td></tr>
	<tr><th scope=row>TCGA-A7-A0DC-01B-04R-A22O-07</th><td>A0DC</td><td>01B</td><td>04R</td><td>A22O</td><td>Tumor </td><td> TRUE</td></tr>
	<tr><th scope=row>TCGA-A7-A13G-01B-04R-A22O-07</th><td>A13G</td><td>01B</td><td>04R</td><td>A22O</td><td>Tumor </td><td> TRUE</td></tr>
	<tr><th scope=row>TCGA-A7-A0DB-01A-11R-A277-07</th><td>A0DB</td><td>01A</td><td>11R</td><td>A277</td><td>Tumor </td><td> TRUE</td></tr>
	<tr><th scope=row>TCGA-A7-A0DB-01C-02R-A277-07</th><td>A0DB</td><td>01C</td><td>02R</td><td>A277</td><td>Tumor </td><td> TRUE</td></tr>
	<tr><th scope=row>TCGA-A7-A13E-01A-11R-A277-07</th><td>A13E</td><td>01A</td><td>11R</td><td>A277</td><td>Tumor </td><td> TRUE</td></tr>
	<tr><th scope=row>TCGA-A7-A13E-01B-06R-A277-07</th><td>A13E</td><td>01B</td><td>06R</td><td>A277</td><td>Tumor </td><td> TRUE</td></tr>
</tbody>
</table>



#### Possible Problems with outlier samples captured by PCA

The issue likely stems from **batch effects** or **technical variability** associated with specific experimental conditions. All the "weird" samples (marked as `TRUE`) share common features, such as being processed on the same plates (A277, A22O, A13Q). This suggests that the **plate or sample processing conditions** may have introduced **systematic differences** that cause these samples to cluster separately in PCA analysis.

>- The problem is likely **technical artifacts** or **batch effects** tied to specific plates or conditions rather than biological differences.
>- Addressing this may involve **batch correction** methods (e.g., COMBAT) or removing these samples if they cannot be corrected adequately.



```R
brca_batchs_paired %>% arrange(plate) %>% .$plate %>%table
```


    .
    A00Z A056 A084 A089 A115 A12D A12P A137 A13Q A144 A14D A14M A157 A169 A16F A17B 
       7   14    3   22   10   32   30   14   30   12   10    6   18    6    2    4 
    A19E A19W A21T A22O A277 A466 
       1    2    2    2    4    1 


>**Plates 277, A22O  have only these wiered samples!**


```R
platesPerCancers %>% filter(plate %in% c("A277","A22O"))
```


<table class="dataframe">
<caption>A grouped_df: 6 × 3</caption>
<thead>
	<tr><th scope=col>plate</th><th scope=col>project_id</th><th scope=col>usedInCancers</th></tr>
	<tr><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>A277</td><td>TCGA-BLCA</td><td>6</td></tr>
	<tr><td>A277</td><td>TCGA-BRCA</td><td>6</td></tr>
	<tr><td>A277</td><td>TCGA-COAD</td><td>6</td></tr>
	<tr><td>A277</td><td>TCGA-KIRC</td><td>6</td></tr>
	<tr><td>A277</td><td>TCGA-LUAD</td><td>6</td></tr>
	<tr><td>A277</td><td>TCGA-UCEC</td><td>6</td></tr>
</tbody>
</table>



>**"A22O"** plate were used only in BRCA project.

### Data re-selection
Considering that these samples are causing significant noise in our data, we have decided to remove the data from these 4 patients. This decision is based on the need to minimize the impact of these outliers and reanalyze the study. 


```R
wiered_participants
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>A0DB</li><li>A0DC</li><li>A13E</li><li>A13G</li></ol>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'A0AU'</li><li>'A0AY'</li><li>'A0AZ'</li><li>'A0B3'</li><li>'A0B5'</li><li>'A0B7'</li><li>'A0B8'</li><li>'A0BA'</li><li>'A0BC'</li><li>'A0BJ'</li><li>'A0BM'</li><li>'A0BQ'</li><li>'A0BS'</li><li>'A0BT'</li><li>'A0BV'</li><li>'A0BW'</li><li>'A0BZ'</li><li>'A0C0'</li><li>'A0C3'</li><li>'A0CE'</li><li>'A0CH'</li><li>'A0D9'</li><li>'A0DB'</li><li>'A0DC'</li><li>'A0DD'</li><li>'A0DG'</li><li>'A0DH'</li><li>'A0DK'</li><li>'A0DL'</li><li>'A0DO'</li><li>'A0DP'</li><li>'A0DQ'</li><li>'A0DT'</li><li>'A0DV'</li><li>'A0DZ'</li><li>'A0E0'</li><li>'A0E1'</li><li>'A0H5'</li><li>'A0H7'</li><li>'A0H9'</li><li>'A0HA'</li><li>'A0HK'</li><li>'A13E'</li><li>'A13F'</li><li>'A13G'</li><li>'A153'</li><li>'A158'</li><li>'A15I'</li><li>'A15K'</li><li>'A15M'</li><li>'A18J'</li><li>'A18K'</li><li>'A18L'</li><li>'A18M'</li><li>'A18N'</li><li>'A18P'</li><li>'A18Q'</li><li>'A18R'</li><li>'A18S'</li><li>'A18U'</li><li>'A18V'</li><li>'A1BC'</li><li>'A1EN'</li><li>'A1EO'</li><li>'A1ET'</li><li>'A1EU'</li><li>'A1EV'</li><li>'A1EW'</li><li>'A1F0'</li><li>'A1F2'</li><li>'A1F6'</li><li>'A1F8'</li><li>'A1FB'</li><li>'A1FC'</li><li>'A1FD'</li><li>'A1FE'</li><li>'A1FG'</li><li>'A1FH'</li><li>'A1FJ'</li><li>'A1FM'</li><li>'A1FN'</li><li>'A1FR'</li><li>'A1FU'</li><li>'A1IG'</li><li>'A1L7'</li><li>'A1LB'</li><li>'A1LH'</li><li>'A1LS'</li><li>'A1N4'</li><li>'A1N5'</li><li>'A1N6'</li><li>'A1N9'</li><li>'A1NA'</li><li>'A1ND'</li><li>'A1NF'</li><li>'A1NG'</li><li>'A1R7'</li><li>'A1RB'</li><li>'A1RC'</li><li>'A1RD'</li><li>'A1RF'</li><li>'A1RH'</li><li>'A1RI'</li><li>'A203'</li><li>'A204'</li><li>'A208'</li><li>'A209'</li><li>'A23H'</li><li>'A2C8'</li><li>'A2C9'</li><li>'A2FB'</li><li>'A2FF'</li><li>'A2FM'</li></ol>
</details>


After removing these participants, the dataset will consist of:

- Primary Tumor: 109
- Solid Tissue Normal: 109


```R
brca_batchs_paired %>% filter(!participant %in% wiered_participants) %>% .$sample_type %>% table
```


    .
          Primary Tumor Solid Tissue Normal 
                    109                 109 



```R
brca_batchs_paired <- brca_batchs_paired %>% filter(!participant %in% wiered_participants)
brca_batchs_paired <-droplevels(brca_batchs_paired)
brca_count_paired <- brca_count %>% dplyr::select(all_of(c("rn",brca_batchs_paired$barcode)))
```


```R
dim(brca_count_paired[,-1])
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>58329</li><li>218</li></ol>




```R
dim(brca_batchs_paired)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>218</li><li>16</li></ol>




```R
pca_unscaled_uncenterd_brca_count_paired <- prcomp(t(brca_count_paired[,-1]), scale = F, center = F)
plot_PCAs(as.data.table(pca_unscaled_uncenterd_brca_count_paired$x), brca_batchs_paired)
```


    
![png](output_222_0.png)
    


### Normalization TMM and Log-transformation


```R
# Normalization
dge_brca <- do_normalization(brca_count_paired,brca_batchs_paired,method ="TMM", group_col = "sample_type") 
# Log-transformation
logCPM_brca <- cpm(dge_brca, log=TRUE)
rownames(logCPM_brca)<- dge_brca$rn
logCPM_brca <- as.data.table(logCPM_brca, keep.rownames = T)
```

### Trunction at the 99th percentile


```R
logCPM_brca_99 <- truncate_at_percentile(logCPM_brca[, -1], truncation_percentil = 0.99)
rownames(logCPM_brca_99) <- logCPM_brca$rn
logCPM_brca_99 <- as.data.table(logCPM_brca_99, keep.rownames = T)
```


```R
# Reshape for ggplot (convert to long format)
logCPM_brca_99_long <- melt(logCPM_brca_99,id.vars = "rn", variable.name = "Sample", value.name = "Expression") %>%
    mutate(sample_type= brca_batchs_paired$sample_type[Sample])
# Create a boxplot for the same selected genes divided by sample type after truncation
ggplot(logCPM_brca_99_long %>% filter(rn %in% selected_genes) , aes(x = sample_type, y = Expression, fill = sample_type)) +
  geom_boxplot(outlier.shape = NA) +  # No outliers displayed in the plot
  labs(title = "Boxplot of Selected Genes (6) - Normal vs Tumor (Truncated at 99th Percentile)",
       x = "Sample Group",
       y = "Log-Transformed CPM") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ rn)  # Separate plots for each gene
```


    
![png](output_227_0.png)
    


### Principal component analysis


```R
pca_unscaled_uncenterd_brca_normalized99 <- prcomp(t(logCPM_brca_99[,-1]), scale = F, center = F)
plot_PCAs(as.data.table(pca_unscaled_uncenterd_brca_normalized99$x), brca_batchs_paired)
```


    
![png](output_229_0.png)
    



```R
res.pca.brca <- PCA(t(logCPM_brca_99[,-1]), graph = FALSE, ncp = 4, scale.unit = F)
```


```R
fviz_pca_ind(res.pca.brca, geom = "point",  pointsize = 3, addEllipses = T,
             col.ind = brca_batchs_paired$sample_type)
```


    
![png](output_231_0.png)
    



```R
fviz_pca_ind(res.pca.brca, geom = "point",  pointsize = 3, addEllipses =F, ellipse.level=0.95, fill.ind = brca_batchs_paired$sample_type,
             col.ind = brca_batchs_paired$plate)
```


    
![png](output_232_0.png)
    



```R
fviz_pca_ind(res.pca.brca, geom = "point",  pointsize = 3, addEllipses =F, ellipse.level=0.95, fill.ind = brca_batchs_paired$sample_type, 
             col.ind = brca_batchs_paired$tss)
```


    
![png](output_233_0.png)
    



```R
# batch effect if there is : I dont concider that they are batch effect to remove ! but let have a look if data get better 
```

### Batch effect removal using Combat


```R
my_mod_brca <-model.matrix(~as.factor(groups), data=brca_batchs_paired)
# Apply ComBat for plate BE
brca_combat4plate <- sva::ComBat(dat = logCPM_brca_99[,-1],
                 mod = my_mod_brca,
                  batch = brca_batchs_paired$plate, 
                  par.prior = TRUE, 
                  prior.plots = FALSE)
```

    Found 23527 genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.


    Using the 'mean only' version of ComBat
    
    Found20batches
    
    Note: one batch has only one sample, setting mean.only=TRUE
    
    Adjusting for1covariate(s) or covariate level(s)
    
    Standardizing Data across genes
    
    Fitting L/S model and finding priors
    
    Finding parametric adjustments
    
    Adjusting the Data
    
    



```R
rownames(brca_combat4plate)<-logCPM_brca_99$rn
brca_combat4plate<-as.data.table(brca_combat4plate, keep.rownames = T)
```


```R
pca_unscaled_uncenterd_brca_combat4plate <- prcomp(t(brca_combat4plate[,-1]), scale = F, center = F)
plot_PCAs(as.data.table(pca_unscaled_uncenterd_brca_combat4plate$x), brca_batchs_paired)
```


    
![png](output_238_0.png)
    



```R
res.pca.brca.combat <- PCA(t(brca_combat4plate[,-1]), graph = FALSE, ncp = 4, scale.unit = F)
fviz_pca_ind(res.pca.brca.combat, geom = "point",  pointsize = 3, addEllipses = T,
             col.ind = brca_batchs_paired$sample_type)
```


    
![png](output_239_0.png)
    



```R
mypca_brca_combat = mixOmics::pca(t(brca_combat4plate[,-1]), ncomp = 5, center = FALSE, scale = FALSE)
plot(mypca_brca_combat)
plotIndiv(mypca_brca_combat, comp = 1:2, col.per.group = levels(brca_batchs_paired$c_by_group), ind.names = 
              FALSE,
          group = brca_batchs_paired$groups,
          # graphical parameters
              style = "ggplot2",
          legend = TRUE, legend.position = "right", 
          legend.title = "Cancer type", 
          legend.title.pch = 'Sample types',              
          legend.pch = FALSE,  # Show point shape in the legend
          ellipse = TRUE, 
          ellipse.level = 0.95)
```


    
![png](output_240_0.png)
    



    
![png](output_240_1.png)
    


>The PCA plot shows that PC1 explains 97% of the variance in the data, indicating that the first principal component captures most of the variability. This suggests that the major source of variation in our dataset is well represented by PC1, potentially reflecting a key biological or technical factor in the data.

### save data for metabopathia


```R
# Final Check : if the meta data and Gene expression data has the same samples:
all(colnames(brca_combat4plate[,-1]) %in% brca_batchs_paired$barcode) # are all in ?
```


TRUE



```R
fwrite(x = brca_combat4plate, file = file.path(my_dir, paste0("BRCA_119Nx119T_paired_normalized_trun99_combat_data_", version_of_data,".tsv")),quote = F, append = F, sep = "\t", row.names = F, col.names = T, verbose = F)
```


```R
fwrite(x = brca_batchs_paired, file = file.path(my_dir, paste0("BRCA_119Nx119T_paired_metadata_", version_of_data,".tsv")),quote = F, append = F, sep = "\t", row.names = F, col.names = T, verbose = F)
```


```R
# The end !
```
