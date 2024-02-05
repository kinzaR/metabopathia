# metabopathia
## Abstract (is a quick abstract to re-write )
Metabopathia is a method for computation of signal transduction along signaling pathway and metabolic pathways from transcriptomic and metabolomic data. This method is based on the Canonical Circuit Activity Analysis method which is an iterative algorithm that compute the signal intensity passing through the nodes of a network. This network are composed by proteins and metabolites: signaling pathways and metabolic pathways maps.
The method is taking into account the level of expression, the metabolic concentration (for genes and metabolites respectively) and the intensity of the signal arriving to each node. Unlike the enreachment-based methods, this method is able to deffer between signal that reduce or increase (activation or inhibition).It is providing an other approach of omics data integration for functional analysis allowing to compute the signal arriving to the effector nodes with are annotated by a function for each pathway (:S bad werritten). 
In summary, Metabopathia presents a refined approach to unraveling the intricate interplay of gene expression and metabolite dynamics within cellular pathways.

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
    ./01_main.R -v -e data_examples/Dystrophic_epidermolysis_bullosa/counts_TMM_normalization.tsv \
    -m data_examples/Dystrophic_epidermolysis_bullosa/metabolite_suero.tsv \
    -d  data_examples/Dystrophic_epidermolysis_bullosa/integration_design.tsv \
    --group1 visit1 --group2 control
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

