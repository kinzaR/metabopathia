# metabopathia
## Background

**Metabopathia** is a computational method designed to model signal transduction through cellular pathways by integrating transcriptomic and metabolomic data. It builds on the Canonical Circuit Activity Analysis, an iterative algorithm that computes signal intensity through networks composed of nodes representing proteins, metabolites, and cellular functions or phenotypes. These networks include both signaling and metabolic pathways, currently utilizing maps from databases like KEGG, with plans to incorporate additional resources such as SIGNOR. Tools like [signor2Hipathia](https://github.com/kinzaR/signor2Hipathia) are used to integrate these networks for enhanced mechanistic modeling.

## Introduction

Metabopathia provides a novel approach to multi-omics data integration, distinguishing it from enrichment-based methods by accurately differentiating between activation and inhibition of pathways. It takes into account expression levels and metabolic activity (for genes and metabolites, respectively), calculating signal intensities as they propagate through the network to effector nodes, which are annotated with specific pathway functions. This approach allows researchers to explore the functional impact of gene expression and metabolite dynamics in a more mechanistically grounded manner.

This repository contains the full implementation of Metabopathia, along with an example study using breast cancer data from The Cancer Genome Atlas (TCGA). The method is currently in development as a web server, accessible [here](http://hipathia.babelomics.org/metabopathia_dev/).


## Setting up the Environment

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

## Case Study 1:

### Overview
Description od the case study and its objectives....

### Proposed Dataset
Mention details about the dataset you recommend for the case study...

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
