# metabopathia

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

  
<a name="KeyConcepts"/>

### Key Concepts

- **Biological Network**: A representation of biological entities (such as genes, proteins, and metabolites) and their interactions within a living system. These networks capture how molecules interact to carry out cellular functions.

- **Pathways**: A sequence of molecular events or interactions that lead to a specific cellular outcome. Pathways can involve signaling cascades (signaling pathways) or metabolic processes (metabolic pathways), describing how cells respond to signals or perform biochemical reactions.

- **Nodes in Biological Pathways**: In a biological network, nodes represent molecular entities, such as genes, proteins, or metabolites. Each node is a key player in the pathway, contributing to the overall biological process or signal transduction.

- **Interactions in Biological Pathways**: The connections between nodes, known as edges, represent molecular interactions, such as activation, inhibition, or binding events. These interactions define how signals or molecular changes propagate through the network, leading to functional outcomes.
- diagrams maps pathway network : Here I have to be clear about different terminology 


<a name="MetabopathiaAproach"/>    

## Metabopathia aproach

The aim of Metabopathia is to offer a novel approach to multi-omics data integration for pathway activity analysis. It reduces the complexity of entire pathways by breaking them down into sub-pathways or circuits (segments that have only one final node, known as the effector). The activity of each individual node is then calculated depending on the type of molecular component involved in signal transduction. These components include proteins and metabolites, with gene expression levels and metabolite activity used as proxies for protein presence and metabolite concentration, respectively. Unlike enrichment-based methods, this algorithm distinguishes between activation and inhibition interactions when inferring signal propagation towards the effector. The final nodes, or effectors, are annotated with their corresponding cellular functions.

To summarize, the algorithm calculates signal intensities as they propagate through the network to effector nodes—the final protein nodes in each sub-pathway or circuit—annotated with specific cellular functions and phenotypic outcomes.

This mechanistic approach enables researchers to better understand the functional impact of gene expression and metabolite dynamics in biological systems, providing a more detailed and accurate representation of cellular signaling and metabolic interactions.

This repository contains the full implementation of Metabopathia, along with an example study using breast cancer data from The Cancer Genome Atlas (TCGA). Metabopathia is currently under development as a web server, accessible [here](http://hipathia.babelomics.org/metabopathia_dev/).

<a name="QuickStart"/>  

## Quick Start

<a name="SettingUptheEnvironment"/>  

### Setting up the Environment

Before running the code, follow these steps to set up the environment:

1. Install R version 4.3.1 on your system (or a newer version).

2.  **Optional for Linux users: If you have multiple R versions installed on your machine, open a terminal and set the `RSTUDIO_WHICH_R` environment variable to specify the R version for RStudio. Use the following command:**

    ```bash
    export RSTUDIO_WHICH_R=/opt/R/4.3.1/bin/R
    ```

    This ensures that RStudio uses the correct R version when launching.


3. Execute the `00_prep_env.R` file to install the required dependencies and packages. Use the following command in your R environment:

    ```R
    source("00_prep_env.R")
    ```

<a name="GettingStarted"/>  

### Getting Started

Follow these steps to launch the case study from the command line:

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/kinzaR/metabopathia.git
   cd metabopathia 
   ```
2. **Install dependencies:**
    ```R
    source("00_prep_env.R")
    ```
3. **Run the case study**
    ```R
    ./01_main.R --example
    ```
4. **Example Command**
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

<a name="BRCACaseStudy"/>   

## Breast cancer case study

<a name="Overview"/>   

### Overview

Description od the case study and its objectives....

<a name="Dataset"/>   

### Dataset 

<a name="preprocessing"/>   

### preprocessing story: From TCGA repository to Metabopathia input data

Mention details about the dataset you recommend for the case study...


<a name="ResultsDiscussion"/>   

### Results and discussion



<a name="FutureEnhancements"/>   

## Future enhancements

Several features are planned to expand and refine Metabopathia:

### Extensibility Enhancements

**Integration of Additional Biological Databases:** Future versions will incorporate resources such as SIGNOR to extend the biological knowledge from several databases. Tools like [signor2Hipathia](https://github.com/kinzaR/signor2Hipathia) will enable the integration of these networks to enrich our mechanistic modeling results. This will help avoid bias towards a single database, and the diversity of resources will improve the quality of the results. *Additionally, the stability of results must be ensured across different databases.*

### Interoperability and Versatility
  The widely used graphical standard, [Systems Biology Graphical Notation (SBGN)](https://sbgn.github.io/), divides the graphical languages for representing cellular processes, interactions, and biological networks into three families:
  - Process Descriptions (PD),
  - Entity Relationships (ER),
  - Activity Flows (AF).  
  Metabopathia, the extended version of [Hipathia](https://doi.org/10.18632/oncotarget.14107), currently accepts only Activity Flow diagrams, while others utilize PD and ER maps. These two formats are not straightforward to integrate with our mechanistic activity modeling approaches: [Hipathia](10.1016/j.csbj.2021.05.022), [Cov-Hipathia](10.1186/s13040-021-00234-1), and Metabopathia. Since PD maps are often large, detailed, and complex, there is a need for simplified illustrations. Adding a module to Metabopathia that parses these maps into AF networks will enhance interoperability with approaches that use other network notations and scenarios where biological knowledge is represented in PD and ER notation languages. This module will take advantage of [CaSQ](https://casq.readthedocs.io/en/stable/_modules/casq/celldesigner2qual.html), a tool that converts Process Description networks to SBML-qual with strict semantics. Then, these simplified interaction format will be easy to adapt to the Metabopathia tool.  
  A proof of concept using CaSQ to ensure interoperability was cited in our previous community work. Below is a schema of the use of Hipathia within [an ecosystem developed by the COVID-19 Disease Map community](https://doi.org/10.15252/msb.202110387), ensuring interoperability between all tools in the Disease Map community's work.

[http://www.embopress.org/cms/10.15252/msb.202110387/asset/d583d912-f6bb-4c4f-a38c-a5edc9fdac03/assets/graphic/msb202110387-fig-0001-m.png](http://www.embopress.org/cms/10.15252/msb.202110387/asset/d583d912-f6bb-4c4f-a38c-a5edc9fdac03/assets/graphic/msb202110387-fig-0001-m.png)

### Amplified Accuracy Through Omics Synergy
Ongoing efforts aim to refine the method’s ability to handle raw metabolic concentrations, improving the interaction modeling between metabolites and proteins. This enhancement will balance the complexity of molecular mechanisms with computational simplicity to reflect biological reality more accurately.
