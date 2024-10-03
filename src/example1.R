opt$spe<-"hsa"
opt$verbose <-TRUE
# opt$exp_file <-"data_examples/Dystrophic_epidermolysis_bullosa/counts_TMM_normalization.tsv"
# opt$met_file <-"data_examples/Dystrophic_epidermolysis_bullosa/metabolite_suero.tsv"
# opt$design_file <-"data_examples/Dystrophic_epidermolysis_bullosa/integration_design.tsv"
# opt$group1 <-"visit1"
# opt$group2 <-"control"
# opt$paired <-TRUE
opt$decompose <-FALSE
opt$design_type <-"categorical"
opt$adjust <-TRUE
opt$difexp <-TRUE
opt$GO.terms <-FALSE
opt$uni.terms <-FALSE
opt$analysis <-"compare"
opt$output_folder <- "tmp"


# ################# example GSE207088
# opt$exp_file <-"data_examples/GSE207088/toBeUsed/RNA_hNESCs_TPM.tsv"
# opt$exp_file <-"data_examples/GSE207088/toBeUsed/RNA_hNESCs_TMM.tsv"
# opt$met_file <-"data_examples/GSE207088/toBeUsed/polar_metab_normalized_annotated.tsv"
# opt$design_file <-"data_examples/GSE207088/toBeUsed/design.tsv"
# opt$group1 <-"IPD"
# opt$group2 <-"CTRL"
# opt$paired <-FALSE# 
################# example for the web 
# opt$exp_file <-"data_examples/BRCA_ER/brca_data_example_ERposneg.tsv"
# opt$design_file <-"data_examples/BRCA_ER/brca_designmatrix_ERposneg.tsv"
# opt$group1 <-"Negative"
# opt$group2 <-"Positive"
################# example GSE207088
opt$exp_file <-"data_examples/metabolizer_as_DS/brca_example1_40_exp.txt"
# opt$met_file <-"data_examples/metabolizer_as_DS/inferedmetabolic_data.tsv"
opt$design_file <-"data_examples/metabolizer_as_DS/brca_example1_40_design.txt"
opt$group1 <-"Tumor"
opt$group2 <-"Normal"
opt$paired <-FALSE
opt$hipathia <- TRUE
# opt$analysis <-"overlay"
opt$analysis <-"compare"
