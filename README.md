# Metabopathia

**Metabopathia** is a computational method designed to model signal transduction through cellular pathways by integrating transcriptomic and metabolic data. It builds on the Canonical Circuit Activity Analysis ([CCAA](https://pubmed.ncbi.nlm.nih.gov/28042959/)), an iterative algorithm that computes signal intensities through networks composed of nodes representing proteins, metabolites, and cellular functions or phenotypes. These networks include both signaling and metabolic pathways, currently utilizing pathway maps from databases like KEGG.
## Table of Contents  
- [Key Concepts](#KeyConcepts)    
- [Metabopathia approach](#MetabopathiaAproach)
- [Quick start](#QuickStart)
    - [Setting up the environment](#SettingUptheEnvironment)
    - [Getting started](#GettingStarted)
    - [Breast cancer case study](#BRCACaseStudy)
        - [Overview](#Overview)
        - [Data preprocessing story: From TCGA repository to Metabopathia input dataset](#preprocessing)
        - [Results and discussion](#ResultsDiscussion) 
- [Future enhancements](#FutureEnhancements)

  
<a name="KeyConcepts"> </a>   
    
### Key Concepts

- **Biological Network**: A representation of biological entities (such as genes, proteins, and metabolites) and their interactions within a living system. These networks capture how molecules interact to carry out cellular functions.

- **Pathways**: A sequence of molecular events or interactions that lead to a specific cellular outcome. Pathways can involve signaling cascades (signaling pathways) or metabolic processes (metabolic pathways), describing how cells respond to signals or perform biochemical reactions.

- **Nodes in Biological Pathways**: In a biological network, nodes represent molecular entities, such as genes, proteins, or metabolites. Each node is a key player in the pathway, contributing to the overall biological process or signal transduction.

- **Interactions in Biological Pathways**: The connections between nodes, known as edges, represent molecular interactions, such as activation, inhibition, or binding events. These interactions define how signals or molecular changes propagate through the network, leading to functional outcomes.
- diagrams maps pathway network : Here I have to be clear about different terminology 


<a name="MetabopathiaAproach"> </a>      

## Metabopathia aproach

The aim of Metabopathia is to offer a novel approach to multi-omics data integration for pathway activity analysis. It reduces the complexity of entire pathways by breaking them down into sub-pathways or circuits (segments that have only one final node, known as the effector). The activity of each individual node is then calculated depending on the type of molecular component involved in signal transduction. These components include proteins and metabolites, with gene expression levels and metabolite activity used as proxies for protein presence and metabolite concentration, respectively. Unlike enrichment-based methods, this algorithm distinguishes between activation and inhibition interactions when inferring signal propagation towards the effector. The final nodes, or effectors, are annotated with their corresponding cellular functions.

To summarize, the algorithm calculates signal intensities as they propagate through the network to effector nodes—the final protein nodes in each sub-pathway or circuit—annotated with specific cellular functions and phenotypic outcomes.

This mechanistic approach enables researchers to better understand the functional impact of gene expression and metabolite dynamics in biological systems, providing a more detailed and accurate representation of cellular signaling and metabolic interactions.

This repository contains the full implementation of Metabopathia, along with an example study using breast cancer data from The Cancer Genome Atlas (TCGA). Metabopathia is currently under development as a web server, accessible [here](http://hipathia.babelomics.org/metabopathia_dev/).

<a name="QuickStart"> </a>       

## Quick Start

<a name="SettingUptheEnvironment"> </a>        

### Setting up the Environment

Before running the code, follow these steps to set up the environment:

1. **Install R version** 4.3.1 on your system (or a newer version).

2.  **[Optional for Linux users]**:
   If you have multiple R versions installed on your machine, open a terminal and set the `RSTUDIO_WHICH_R` environment variable to specify the R version for RStudio. Use the following command:

    ```bash
    export RSTUDIO_WHICH_R=/opt/R/4.3.1/bin/R
    ```
3. **Install required Dependencies**:
The required dependencies are:  ```"BiocManager", "hipathia", "igraph", "SummarizedExperiment", "Matrix", "preprocessCore", "KEGGREST", "optparse"```. You can either install these manually or clone the repository and run the `00_prep_env.R` script to install all the required dependencies and packages.
      
4. **R Session Information**: For reproducibility and troubleshooting, here is the R session information used in this project:

```R
R version 4.3.1 (2023-06-16)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.6 LTS

Matrix products: default
BLAS:   /opt/R/4.3.1/lib/R/lib/libRblas.so 
LAPACK: /opt/R/4.3.1/lib/R/lib/libRlapack.so;  LAPACK version 3.11.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=es_ES.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=es_ES.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=es_ES.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=es_ES.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/Madrid
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
  [1] tidyselect_1.2.1              dplyr_1.1.4                   blob_1.2.4                   
  [4] filelock_1.0.3                Biostrings_2.68.1             bitops_1.0-7                 
  [7] fastmap_1.1.1                 RCurl_1.98-1.14               BiocFileCache_2.8.0          
 [10] promises_1.3.0                digest_0.6.35                 mime_0.12                    
 [13] lifecycle_1.0.4               KEGGREST_1.40.1               interactiveDisplayBase_1.38.0
 [16] RSQLite_2.3.6                 magrittr_2.0.3                compiler_4.3.1               
 [19] rlang_1.1.3                   tools_4.3.1                   igraph_2.0.3                 
 [22] utf8_1.2.4                    yaml_2.3.8                    ggsignif_0.6.4               
 [25] S4Arrays_1.0.6                htmlwidgets_1.6.4             bit_4.0.5                    
 [28] curl_5.2.1                    DelayedArray_0.26.7           plyr_1.8.9                   
 [31] abind_1.4-5                   purrr_1.0.2                   BiocGenerics_0.46.0          
 [34] grid_4.3.1                    stats4_4.3.1                  preprocessCore_1.62.1        
 [37] fansi_1.0.6                   ggpubr_0.6.0                  xtable_1.8-4                 
 [40] colorspace_2.1-0              ggplot2_3.5.1                 scales_1.3.0                 
 [43] MultiAssayExperiment_1.26.0   SummarizedExperiment_1.30.2   cli_3.6.2                    
 [46] crayon_1.5.2                  generics_0.1.3                rstudioapi_0.16.0            
 [49] reshape2_1.4.4                httr_1.4.7                    visNetwork_2.1.2             
 [52] DBI_1.2.2                     cachem_1.0.8                  stringr_1.5.1                
 [55] zlibbioc_1.46.0               AnnotationDbi_1.62.2          BiocManager_1.30.23          
 [58] XVector_0.40.0                matrixStats_1.3.0             vctrs_0.6.5                  
 [61] Matrix_1.6-5                  carData_3.0-5                 jsonlite_1.8.8               
 [64] car_3.1-2                     IRanges_2.34.1                S4Vectors_0.38.2             
 [67] rstatix_0.7.2                 bit64_4.0.5                   hipathia_3.0.2               
 [70] tidyr_1.3.1                   servr_0.30                    glue_1.7.0                   
 [73] stringi_1.8.3                 gtable_0.3.5                  BiocVersion_3.17.1           
 [76] later_1.3.2                   GenomeInfoDb_1.36.4           GenomicRanges_1.52.1         
 [79] munsell_0.5.1                 tibble_3.2.1                  pillar_1.9.0                 
 [82] rappdirs_0.3.3                htmltools_0.5.8.1             GenomeInfoDbData_1.2.10      
 [85] R6_2.5.1                      dbplyr_2.3.4                  shiny_1.8.1.1                
 [88] Biobase_2.60.0                lattice_0.22-6                MetBrewer_0.2.0              
 [91] AnnotationHub_3.8.0           backports_1.4.1               png_0.1-8                    
 [94] broom_1.0.5                   memoise_2.0.1                 httpuv_1.6.15                
 [97] Rcpp_1.0.12                   xfun_0.43                     MatrixGenerics_1.12.3        
[100] pkgconfig_2.0.3
```
### <a name="GettingStarted"> </a> Getting Started

Follow these steps to launch the case study from the command line:

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/kinzaR/metabopathia.git
   cd metabopathia 
   ```
2. **Install dependencies:**
   Execute the `00_prep_env.R` file to install the required dependencies and packages. Use the following command in your R environment:
    ```R
    source("00_prep_env.R")
    ```
   This ensures that RStudio uses the correct R version when launching.
4. **Run the case study**
    ```R
    ./01_main.R --example
    ```
5. **Example Command**
For help: 
   ```bash
   ./01_main.R -h 
   ```
   or 
   ```bash
   Rscript 01_main.R -h 
   ```
4. **Parameters**
    ```bash
    $ ./01_main.R --help
    Usage: ./01_main.R [options]


    Options:
        -s SPE, --spe=SPE
            Species variable. Allowed choices: 'hsa', 'mmu', 'rno'. (default: hsa)

        -v, --verbose
            Enable verbose mode.

        -e EXP_FILE, --exp_file=EXP_FILE
            Path to the expression file.

        -m MET_FILE, --met_file=MET_FILE
            Path to the metabolomics concentration file.

        -t MET_TYPE, --met_type=MET_TYPE
            Allowed values: 
                  inferred: Infer production of metabolites from Metabolizer (doi: 10.1038/s41540-019-0087-2).
                  perturbations: Study perturbance in metabolite concentration.
                  concentration_matrix: Using metabolomics concentration matrix.

        -d DESIGN_FILE, --design_file=DESIGN_FILE
            Path to the design file.

        --group1=GROUP1
            Label of the first group.

        --group2=GROUP2
            Label of the second group to be compared (Reference condition).

        -p PATHWAYS_LIST, --pathways_list=PATHWAYS_LIST
            Vector of the IDs of the pathways to load. By default, all available pathways are loaded. Example: '04014,04015'.

        --paired
            Boolean, whether the samples to be compared are paired. If TRUE, function wilcoxsign_test from package coin is used. If FALSE, function wilcox.test from package stats is used.

        --decompose
            Boolean, whether to compute the values for the decomposed subpathways. By default, effector subpathways are computed.

        --design_type=DESIGN_TYPE
            Type of design. Allowed values: 'categorical' or 'continuous'. Default is 'categorical'.

        --adjust
            Boolean, whether to adjust the p.value with Benjamini-Hochberg FDR method. Default is TRUE.

        --conf.level=CONF.LEVEL
            Level of significance. By default 0.05.

        --difexp
            Boolean, whether to perform differential expression analysis.

        --GO.terms
            Boolean, whether to compute functional analysis with Gene Ontology terms.

        --uni.terms
            Boolean, whether to compute functional analysis with Uniprot keywords.

        --custom.terms=CUSTOM.TERMS
            Path to a file containing a data.frame with the custom annotation of the genes to the functions. First column are gene symbols, second column the functions.

        --analysis=ANALYSIS
            Type of analysis. Allowed values:'overlay','ORA', 'compare', 'predictor_test', 'predictor_train', 'variant_interpreter', 'drug_repurposing'. Default is 'compare'.


                          Differential Signaling 'compare': Provides an estimation of significant cell signaling activity changes across different conditions.

                          Predictor 'predictor': Allows you to train a prediction-test and test it with different data.

                          Variant Functional Interpretation 'variant_interpreter': Provides an estimation of the potential impact of genomic variation on cell signaling and cell functionality.

                          Drug Repurposing 'drug_repurposing': drug repurposing.

        --output_folder=OUTPUT_FOLDER
            Output folder path. Default is 'tmp'.

        --hipathia
            Enable calculation using HiPathia for result comparison between MetaboPathia and HiPathia.

        --example
            Load variables from the example config file (src/example1.R).

        --ready_mgi
            Load ready mgi pre-processed.

        -h, --help
            Show this help message and exit

    ```

<a name="BRCACaseStudy"> </a>       

## Breast cancer case study

<a name="Overview"> </a>         

### Overview

Description od the case study and its objectives....

<a name="Dataset"> </a>         

### Dataset 

<a name="preprocessing"> </a>      

### preprocessing story: From TCGA repository to Metabopathia input data

Mention details about the dataset you recommend for the case study...

### Metabolite Inference

A total of **15 metabolites** were inferred using the `infer_met_from_metabolizer()` function from the Metabolizer package. The identified metabolites are `C00002, C00022, C00026, C00042, C00122, C00130, C00195, C00334, C00410, C00762, C00788, C01598, C01673, C01780, C02465`.
According to [MetaboAnalyst 5.0 platform ](https://www.metaboanalyst.ca/Secure/process/NameMapView.xhtml)  [\[Ref\]](https://doi.org/10.1093/nar/gkp356), the annotation of these inferred metabolites across HMDB, SMILES, and PubChem databases is:

| Query   |                Match                | HMDB         | PubChem | KEGG    | SMILES                                                                                                        |
|---------|-------------------------------------|--------------|--------|---------|----------------------------------------------------------------------------------------------------------------|
| C00002  | Adenosine triphosphate             | HMDB0000538  | 5957   | C00002  | NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O                               |
| C00022  | Pyruvic acid                       | HMDB0000243  | 1060   | C00022  | CC(=O)C(O)=O                                                                                               |
| C00026  | Oxoglutaric acid                   | HMDB0000208  | 51     | C00026  | OC(=O)CCC(=O)C(O)=O                                                                                       |
| C00042  | Succinic acid                      | HMDB0000254  | 1110   | C00042  | OC(=O)CCC(O)=O                                                                                             |
| C00122  | Fumaric acid                       | HMDB0000134  | 444972 | C00122  | OC(=O)\\C=C\\C(O)=O                                                                                        |
| C00130  | Inosinic acid                      | HMDB0000175  | 8582   | C00130  | O[C@@H]1[C@@H](COP(O)(O)=O)O[C@H]([C@@H]1O)N1C=NC2=C1N=CNC2=O                                             |
| C00195  | Cer(d18:1/12:0)                    | HMDB0004947  | 5283562| C00195  | [H][C@@](CO)(NC(=O)CCCCCCCCCCC)[C@H](O)\\C=C\\CCCCCCCCCCCCC                                               |
| C00334  | gamma-Aminobutyric acid            | HMDB0000112  | 119    | C00334  | NCCCC(O)=O                                                                                                  |
| C00410  | Progesterone                       | HMDB0001830  | 5994   | C00410  | [H][C@@]12CC[C@H](C(C)=O)[C@@]1(C)CC[C@@]1([H])[C@@]2([H])CCC2=CC(=O)CC[C@]12C                           |
| C00762  | Cortisone                          | HMDB0002802  | 222786 | C00762  | [H][C@@]12CC[C@](O)(C(=O)CO)[C@@]1(C)CC(=O)[C@@]1([H])[C@@]2([H])CCC2=CC(=O)CC[C@]12C                   |
| C00788  | Epinephrine                        | HMDB0000068  | 5816   | C00788  | CNC[C@H](O)C1=CC(O)=C(O)C=C1                                                                                |
| C01598  | Melatonin                          | HMDB0001389  | 896    | C01598  | COC1=CC2=C(NC=C2CCNC(C)=O)C=C1                                                                             |
| C01673  | Calcitriol                         | HMDB0001903  | 5280453| C01673  | C[C@H](CCCC(C)(C)O)[C@@]1([H])CC[C@@]2([H])\\C(CCC[C@]12C)=C\\C=C1\\C[C@@H](O)C[C@H](O)C1=C              |
| C01780  | Aldosterone                        | HMDB0000037  | 5839   | C01780  | [H][C@@]12CC[C@H](C(=O)CO)[C@]1(C[C@H](O)[C@@]1([H])[C@@]2([H])CCC2=CC(=O)CC[C@]12C)C=O                    |
| C02465  | Liothyronine                       | HMDB0000265  | 5920   | C02465  | N[C@@H](CC1=CC(I)=C(OC2=CC(I)=C(O)C=C2)C(I)=C1)C(O)=O                                                   |

[Here you will find the downloaded CSV from the MetaboAnalyst platform](https://github.com/kinzaR/metabopathia/blob/dev/supplementary_files/brca_caseStudy/inferred_metabolite_brca_metabolizer_v2.csv).


<a name="ResultsDiscussion"> </a>         

### Results and discussion



<a name="FutureEnhancements"> </a>         

## Future enhancements

Several features are planned to expand and refine Metabopathia:
### Amplified Accuracy Through Omics Synergy
Ongoing efforts aim to refine the method’s ability to handle raw metabolic concentrations, improving the interaction modeling between metabolites and proteins. This enhancement will balance the complexity of molecular mechanisms with computational simplicity to reflect biological reality more accurately.
### Extensibility Enhancements
**Integration of Additional Biological Databases:** Future versions will incorporate resources such as SIGNOR to extend the biological knowledge from several databases. Tools like [signor2Hipathia](https://github.com/kinzaR/signor2Hipathia) will enable the integration of these networks to enrich our mechanistic modeling results. This will help avoid bias towards a single database, and the diversity of resources will improve the quality of the results. *Additionally, the stability of results must be ensured across different databases.*

### Interoperability and Versatility
  The widely used graphical standard, [Systems Biology Graphical Notation (SBGN)](https://sbgn.github.io/), divides the graphical languages for representing cellular processes, interactions, and biological networks into three families:
  - Process Descriptions (PD),
  - Entity Relationships (ER),
  - Activity Flows (AF).  
  Metabopathia, the extended version of [Hipathia](https://doi.org/10.18632/oncotarget.14107), currently accepts only Activity Flow diagrams, while others utilize PD and ER maps. These two formats are not straightforward to integrate with our mechanistic activity modeling approaches: [Hipathia](10.1016/j.csbj.2021.05.022), [Cov-Hipathia](10.1186/s13040-021-00234-1), and Metabopathia. Since PD maps are often large, detailed, and complex, there is a need for simplified illustrations. Adding a module to Metabopathia that parses these maps into AF networks will enhance interoperability with approaches that use other network notations and scenarios where biological knowledge is represented in PD and ER notation languages. This module will take advantage of [CaSQ](https://casq.readthedocs.io/en/stable/_modules/casq/celldesigner2qual.html), a tool that converts Process Description networks to SBML-qual with strict semantics. Then, these simplified interaction format will be easy to adapt to the Metabopathia tool.

  A proof of concept using [CaSQ](https://casq.readthedocs.io/en/stable/_modules/casq/celldesigner2qual.html) with [Hipathia](10.1016/j.csbj.2021.05.022) to ensure interoperability was cited in our previous community work. Below is a schema of the use of Hipathia within [an ecosystem developed by the COVID-19 Disease Map community](https://doi.org/10.15252/msb.202110387), ensuring interoperability between all tools in the Disease Map community's work.


<p align="center">
    <img src="https://github.com/user-attachments/assets/f81d0b4e-5d57-45ee-9228-a9d69efd8521" alt="https://www.embopress.org/cms/10.15252/msb.202110387/asset/d583d912-f6bb-4c4f-a38c-a5edc9fdac03/assets/graphic/msb202110387-fig-0001-m.png"/><em>Figure 1. Ecosystem of the COVID‐19 Disease Map Community (https://doi.org/10.15252/msb.202110387) </em>
</p>

