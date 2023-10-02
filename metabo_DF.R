cat("Welcome to MetaboPathia")
##
# MetaboPathia (drived from Hipathia method) is a method for the computation of signal transduction along signaling pathways
# not only from transcriptomic data but also adding the layer of metabolomic data. 
# As hipathia, it is also based on an iterative algorithm which is able to
# compute the signal intensity passing through the nodes of a network by taking into account
# the level of expression of each gene/metabolite and the intensity of the signal arriving to it. 
# It also provides a new approach to functional analysis allowing to compute the signal arriving to the
# functions annotated to each pathway.
##
# library(igraph, warn.conflicts = F)
# library(ggplot2)
# library("e1071", warn.conflicts = F)
# suppressPackageStartupMessages(library(R.utils, warn.conflicts = F))
library(hipathia, warn.conflicts = F)

################################################################################
# #### INPUT DATA (from web) ## START
# args <- commandArgs(trailingOnly = F, asValues = T,excludeEnvVars = F)
# codebase <- paste0(dirname(normalizePath(args[["file"]])),"/")
# 
# # species
# species <- args[['species']]
# 
# # Genomics: gene variant List file (optional/if there is variants)
# variant_file <- args[['variant_file']]
# 
# # Transcriptomics/protemoics: expression
# ## custom exp file
# custom <- args[["custom"]] # boolean
# exp_file <- args[['exp_file']]
# ## specific tissue from healthi Gtex
# # this will be removed because I will launch it in parallel
# # in this code I will use either custom or gtex tissue
# # TissueList <- args[['TissueList']]
# Tissue <- args[['tissue']]
# 
# # metabolomics
# met_file <- args[['met_file']]
# 
# # param
# unadjusted <- args[["unadjusted"]] # do not adjust pvalue
# paired <- args[["paired"]] # do not adjust pvalue
# 
# #design: here I have to think which scenario will we allow
# design_type <- args[['design_type']]
# design_file <- args[['design_file']]
# cond1 <- args[['cond1']]
# cond2 <- args[['cond2']]
# 
# # method parameters
# pathways_list <- args[['pathways_list']]
# decompose <- args[['decompose']]
# difexp <- args[['difexp']]
# #difexp <- TRUE  overwritten / but removed now because has to be false in prediction test
# 
# # functional analysis
# go <- args[['go']]
# uniprot <- args[['uniprot']]
# 
# # output
# output_folder <- paste(args[['output_folder']], collapse=" ")
# 
# verbose <- args[['verbose']]
# report <- args[['report']]
# #### INPUT DATA (from web) ## END
################################################################################
## EXAMPLE 1: brca_fake_integration
codebase <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(codebase)
species <- "hsa"
# here I will load only easy pathways 
source("filter_path.R")
pathways_list <- get_easy_pathways(hipathia::load_pathways(species))
exp_file <- "data_examples/brca_fake_integration/exp_toy_brca_data_fake.tsv"
design_file <- "data_examples/brca_fake_integration/metabo_brca_design.tsv"
met_file <-"data_examples/brca_fake_integration/metabo_brca_data.tsv"
cond1 <- "TUMOR"
cond2 <- "NORMAL"
paired <- T # not sure !
decompose <- F
difexp <- T
go <- T
uniprot <- T
output_folder <- "test"
verbose <- T

## EXAMPLE 2
# codebase <- "."
# species <- "hsa"
# exp_file <- "pathivar2015/1000gexampledata/GD462RPKM.txt"
# vcf <- "pathivar2015/1000gexampledata/1000g_all_chr_pathigraphgenes.vcf"
# design_file <- "pathivar2015/1000gexampledata/1000gexp_design.txt"
# cond1 <- "GBR"
# cond2 <- "YRI"
# annovar_dir <- "~/annovar_for_pathivar2015/"
# output_folder <- "/mnt/data2/ngs2/web_pretty/prettyvar/test1"

## EXAMPLE 3
# codebase <- "."
# species <- "hsa"
# pathways_list <- "04210,04010"
# exp_file <- "pathivar2015/1000gexampledata/GD462RPKM.txt"
# design_file <- "pathivar2015/1000gexampledata/1000gexp_design.txt"
# cond1 <- "GBR"
# cond2 <- "YRI"
# output_folder <- "/mnt/data2/ngs2/web_pretty/prettyvar/test1"

## EXAMPLE 4
# codebase <- "."
# species <- "hsa"
# pathways_list <- "04210,04010"
# exp_file <- "examples//brca_example1_40__exp.txt"
# design_type <- "categorical"
# design_file <- "examples//brca_example1_40__design.txt"
# cond1 <- "Tumor"
# cond2 <- "Normal"
# decompose <- F
# difexp <- T
# go <- F
# uniprot <- T
# output_folder <- "/mnt/data2/ngs2/web_pretty/last_test2/simple_2"

## CHECK PARAMETERS

#D if(!exists("report") & is.null(report)) report <- F

if(!exists("species") || is.null(species)) species <- "hsa"

if(!exists("decompose") || is.null(decompose)) decompose <- F

if(!exists("go") || is.null(go)) go <- F
if(!exists("uniprot") || is.null(uniprot)) uniprot <- F
if(!exists("difexp") || is.null(difexp)) difexp <- F
if(!exists("analysis") || is.null(analysis)) analysis <- "compare"# or integration .... 
if(!exists("unadjusted") || is.null(unadjusted)) unadjusted <- F
if(!exists("paired") || is.null(paired)) {
  paired <- F
} else {
  paired <- T
}

if(is.null(cond1) | is.null(cond2)) {
  cond1 <- "Tumor"
  cond2 <- "Normal"
}
if(exists("exp_file") && !is.null(exp_file)){
  dataset <- "gex"
} else stop("ERROR: The expression file matrix is not uploaded. Please make sure to upload the expression data before proceeding.")
if(!exists("met_file") || is.null(met_file)){
  warning("Metabolomics profile is missing.")
}
if(!exists("verbose") || is.null(verbose)) verbose <- F

if(!exists("design_type") || is.null(design_type)) design_type <- "categorical"
if(!exists("pathways_list")) pathways_list<-NULL
if(verbose==T){
  cat("\nWELCOME TO METABOPATHIA Alpha-version:\n")
  cat("\nDetected parameters: \n")
  cat("\t\tcodebase: ",codebase,"\n")
  cat("\t\toutput_folder: ",output_folder,"\n")
  cat("\t\tspecies: ",species,"\n")
  cat("\t>>>Expression:\n")
  cat("\t\texp_file: ",exp_file,"\n")
  cat("\t\texpression: ",dataset,"\n")
  cat("\t>>>Metabolomics:\n")
  cat("\t\tmet_file: ",met_file,"\n")
  cat("\t>>>Experimental design:\n")
  cat("\t\tanalysis: ",analysis,"\n")
  cat("\t\tdesign_type: ",design_type,"\n")
  cat("\t\tdesign_file: ",design_file,"\n")
  cat("\t\tcond1: ",cond1,"\n")
  cat("\t\tcond2: ",cond2,"\n")
  cat("\t>>>Method parameters:\n")
  cat("\t\tpathways_list: ",ifelse(is.null(pathways_list), "All", pathways_list),"\n")
  cat("\t\tdecompose: ",decompose,"\n")
  cat("\t\tdifexp: ",difexp,"\n")
  cat("\t>>>Functional analysis:\n")
  cat("\t\tgo: ",go,"\n")
  cat("\t\tuniprot: ",uniprot,"\n")
  # cat("\n\n")
  # cat("\t\treport: ",report,"\n")
  # cat("\n\n")
}

#### PREPARE DATA
if (!dir.exists(output_folder)) dir.create(output_folder)
status <- function(value){
  write(value,file=paste0(output_folder,"/status.txt"))
}
status("0")
# suppress warning messages generated by various functions or operations
# options(warn=-1)
# load pathways
pathways <- hipathia::load_pathways(species = species, pathways_list = pathways_list)
# ## check which pathways are included:
# hipathia::get_pathways_list(pathways)
status("4")
#### LOAD DATA
cat("Loading data...\n")
# load expression
if(dataset=="gex"){
  exp <- read.table(exp_file,header=T,sep="\t",stringsAsFactors=F,row.names = 1,comment.char="",check.names=F)
  exp <- hipathia::translate_data(as.matrix(exp), species=species, verbose=T)
}
# load metabolite concentrations
metabo_data <- read.table(met_file,header=T,sep="\t",stringsAsFactors=F,row.names = 1,comment.char="",check.names=F)
# HERE: I have to translate_data() from metabo names to KEGG IDs

# load design
des <- read.table(design_file,header=F,stringsAsFactors=F)
#paired data has to be ordred # HERE I am assuming that they are ordered forthe fake example
  #if(paired==F){
colnames(des) <- c("sample",c("group","value")[(design_type == "continuous")+1])
  #} else {
  #  colnames(des) <- c("sample",c("group","value")[(design_type == "continuous")+1], "donor")
  #}
rownames(des) <- des$sample
# Filter only the cond1 and cond2
if(design_type == "categorical"){
  des <- des[ des$group==cond1 | des$group==cond2, ]
}
# Make sure that the exp data is the same in the design file
sel_samples <- intersect(colnames(exp),rownames(des))
if(length(sel_samples)==0){
  stop("ERROR: No intersection between samples names in the Expression matrix and the Design matrix; please check your input data ")
}
# at least 3 complete pairs of observation
if(length(sel_samples)<3){
  stop("ERROR: Not enough samples in the Expression matrix (at least 3 complete pairs); please check your input data ")
}
exp <- exp[,sel_samples]
des <- des[sel_samples,]
#paired data suppose t o be ordred
#  if(paired==T & FALSE){
#   sample_order <- order(des$donor)
#    ordered_samples <- des[sample_order,"sample"]
#    exp <- exp[,ordered_samples]
#    des <- des[ordered_samples,]
#  }

#### PREPROCESS DATA
status("20")

#### RUN
cat("Propagating signaling...\n")
### here I have to check pathways metab-nodes
names(pathways)
# # and maybe  merge data 
# hidata <- hipathia::hipathia(exp, pathways, uni.terms = TRUE, GO.terms = TRUE,
#                    decompose = FALSE, verbose=TRUE)
# source from file 
source("utils.R")
metabo_pathways <- add_metabolite_to_mgi(pathways)
# test if the max is dif from min 
# normalize both  data as a block ?! 
# ensure that data from different sources are on a comparable scale 
# they contribute equally to downstream analysis?
genes_vals <- normalize_data(exp, by_quantiles = FALSE, 
                             by_gene = FALSE, percentil = FALSE)
metabo_vals <- normalize_data(as.matrix(metabo_data), by_quantiles = FALSE, 
                             by_gene = FALSE, percentil = FALSE)
metdata <- metabopathia(genes_vals, metabo_vals, metabo_pathways, uni.terms = TRUE, GO.terms = TRUE,
                   decompose = FALSE, verbose=TRUE)

status("50")

## DIFFERENTIAL EXPRESSION

h <- NULL
if(difexp==T &  design_type != "continuous"){
  difcolor <- node.color.per.differential.expression(results, fpathigraphs, cond1, cond2, des)
  h <- summarize.atts(list(difcolor), c("diff.exp"))
}

## LOAD FUNCTIONS
entrez2hgnc <- read.table( paste0(codebase, "/files/annotations/", species, "/entrez_hgnc_", species, ".annot"),header=F,sep="\t",stringsAsFactors=F)
if(go==T | uniprot==T){
  if(go==T){
    go_annot <- load.annot.file( paste0(codebase, "/files/annotations/", species, "/go_bp_", species, ".annot"))
    go.vals <- prettyfuns(results, fpathigraphs, go_annot, entrez2hgnc)
  }
  if(uniprot==T){
    uniprot_annot <- load.annot.file( paste0(codebase, "/files/annotations/", species, "/uniprot_keywords_", species, "__biological_process.annot"))
    fun.vals <- prettyfuns(results, fpathigraphs, uniprot_annot, entrez2hgnc)
  }
}


#### ANALYSIS

if(decompose==T){
  path.vals <- results$all$path.vals
} else {
  path.vals <- results$all$effector.path.vals
}
n <- ncol(path.vals)

# wilcoxon
cat("Comparing groups...\n")

## ANALYSIS

switch(analysis,
       compare = {
         #dir.create(paste0(output_folder,"/report/"))
         if( design_type == "continuous"){
           write_output_cor <- function(vals,comp,pca_mod,prefix){
             write.table(vals,file=paste0(prefix,"_vals.txt"),sep="\t",quote=F,row.names=T,col.names=T)
             out_comp <- comp[,c("path","UP/DOWN","correlation","p.value","FDRp.value")]
             colnames(out_comp) <- gsub("path","path/term",colnames(out_comp))
             out_comp <- out_comp[order(out_comp$p.value,decreasing=F),]
             write.table(out_comp,file=paste0(prefix,"_significance.txt"),sep="\t",quote=F,row.names=F,col.names=T)
           }
           # correlation
           cat("Computing correlation...\n")
           print("path.vals")
           print(path.vals)
           wt <- do.cor(path.vals, des,adjust=(unadjusted==F))
           wt$path <- gsub("\\*","",get.pretty.name.path(pathnames = rownames(wt),pathigraph=fpathigraphs))
           write_output_cor(path.vals,wt,pca_mod,paste0(output_folder,"/paths"))
           if(go==T){
             print("go.vals")
             print(go.vals)
             go_wt <- do.cor(go.vals,des,adjust=(unadjusted==F))
             go_wt$path <- rownames(go_wt)
             write_output_cor(go.vals,go_wt,go_pca_mod,paste0(output_folder,"/go"))
           }
           if(uniprot==T){
             print("fun.vals")
             print(fun.vals)
             fun_wt <- do.cor(fun.vals, des,adjust=(unadjusted==F))
             fun_wt$path <- rownames(fun_wt)
             write_output_cor(fun.vals,fun_wt,fun_pca_mod,paste0(output_folder,"/uniprot"))
           }
         } else {
           write_output <- function(vals,comp,pca_mod,prefix,colors=c("red","blue"),sample_colors=NULL){
             write.table(vals,file=paste0(prefix,"_vals.txt"),sep="\t",quote=F,row.names=T,col.names=T)
             
             out_comp <- comp[,c("path","UP/DOWN","statistic","p.value","FDRp.value", "Fold_Change", "logFC")]
             colnames(out_comp) <- gsub("path","path/term",colnames(out_comp))
             out_comp <- out_comp[order(out_comp$p.value,decreasing=F),]
             write.table(out_comp,file=paste0(prefix,"_significance.txt"),sep="\t",quote=F,row.names=F,col.names=T)
             
             png(paste0(prefix,"_heatmap.png"),width=800,height=800)
             #variable.clust=T
             plot.heatmap(vals,sample_type = group.value,sample_colors=sample_colors,variable.clust=T)
             dev.off()
             cat("\n\npng hecha sin problema \n\n")
             
             png(paste0(prefix,"_pca.png"),width=800,height=600)
             plot.multiple.pca(pca_mod,1:4,colors=colors,legend = names(sample_colors), legend_colors = sample_colors,cex=1.2)
             dev.off()
           }
           # design
           group.value <- des[colnames(exp),"group"]
           colors <- c("#204253", "#ea8c39")[(group.value==cond2)+1]
           sample_colors <- c("#204253", "#ea8c39")
           names(sample_colors) <- c(cond1,cond2)
           
           #random data generation
           #path_vals<-matrix(rexp(40), 10)
           #colnames(path_vals)<-c("A","B","C","D")
           #info<-c("caso","caso","control","control")
           #names(info)<-colnames(path.vals)
           
           #function application
           #print("path.vals")
           #print(path.vals)
           
           info<-group.value
           names(info)<-colnames(path.vals)
           #print("info")
           #print(info)
           #print(cond1)
           #print(cond2)
           fc<-get_FC(path.vals,info,cond1,cond2)
           print("fc")
           print(fc)
           
           wt <- do.wilcox(path.vals, group.value, g1=cond1, g2=cond2,adjust=(unadjusted==F),paired=paired)
           wt <- cbind(wt, fc[rownames(wt),])
           print("wt")
           print(head(wt))
           wt$path <- gsub("\\*","",get.pretty.name.path(pathnames = rownames(wt),pathigraph=fpathigraphs))
           rank <- order(wt$p.value,decreasing=F)
           pca_mod <- do.pca(path.vals[rank[1:(min(20,n))],],cor=(dataset=="gex"))
           write_output(path.vals,wt,pca_mod,paste0(output_folder,"/paths"),colors=colors,sample_colors=sample_colors)
           if(go==T){
             #fold change:
             go_fc<-get_FC(go.vals,info,cond1,cond2)
             print("go_fc")
             print(go_fc)
             
             go_wt <- do.wilcox(go.vals, group.value, g1=cond1, g2=cond2,adjust=(unadjusted==F),paired=paired)
             print("rownames(go_wt)")
             print(rownames(go_wt))
             go_wt <- cbind(go_wt, go_fc[rownames(go_wt),])
             print("go_wt")
             print(head(go_wt))
             go_wt$path <- rownames(go_wt)
             go_rank <- order(go_wt$p.value,decreasing=F)
             ngo <- min(c(20,ncol(go.vals),nrow(go.vals)))
             go_pca_mod <- do.pca(go.vals[go_rank[1:ngo],],cor=(dataset=="gex"))
             write_output(go.vals,go_wt,go_pca_mod,paste0(output_folder,"/go"),colors=colors,sample_colors=sample_colors)
           }
           if(uniprot==T){
             #fold change:
             fun_fc<-get_FC(fun.vals,info,cond1,cond2)
             print("fun_fc")
             print(fun_fc)
             
             fun_wt <- do.wilcox(fun.vals, group.value, g1=cond1, g2=cond2,adjust=(unadjusted==F),paired=paired)
             fun_wt <- cbind(fun_wt, fun_fc[rownames(fun_wt),])
             print("fun_wt")
             print(head(fun_wt))
             fun_wt$path <- rownames(fun_wt)
             fun_rank <- order(fun_wt$p.value,decreasing=F)
             nfun <- min(c(20,ncol(fun.vals),nrow(fun.vals)))
             fun_pca_mod <- do.pca(fun.vals[fun_rank[1:nfun],],cor=(dataset=="gex"))
             write_output(fun.vals,fun_wt,fun_pca_mod,paste0(output_folder,"/uniprot"),colors=colors,sample_colors=sample_colors)
           }
         }
         status("70")
         save.results(results,wt,fpathigraphs,output_folder,effector=(decompose==F),moreatts=h, entrez2hgnc=entrez2hgnc)
         if(report==T){
           create.html.report2(fpathigraphs,wt,codebase,output_folder,effector=(decompose==F),template_name = "index_template.html",output_name = "index.html",clean_out_folder = F)
         }
       },
       train = {
         if(n<=20) {k=3} else if (n<=30) {k=5} else {k=10}
         
         #save.image(file=paste0(output_folder,"/antesGetPredictionModel.RData"))
         cat("Generating prediction model...\n")
         class0 <- cond2
         class1 <- cond1
         mod <- get.prediction.model(path.vals, des, design_type, k, class0, output_folder, filter_paths, paste0(codebase,"files/",species,"_circuit_names.tsv"))
         mod$design_type <- design_type
         mod$pathways_list <- gsub(species, "", pathways_list)
         mod$des <- des
         mod$class0 <- class0
         mod$class1 <- class1
         
         status("75")
         
         cat("Calculating prediction model statistics...\n")
         mod_stats <- get.predmod.stats(mod, design_type)
         #save.image(file=paste0(output_folder,"/despues.RData"))
         
         cat("Optimal hyperparameters: One random split analysis (20% test)...\n")
         # output saved to output_folder/split
         one_split_analysis(mod, design_type, output_folder,class0,class1)
         # save
         save.pred.res(
           mod,
           mod_stats,
           codebase,
           path.vals,
           fpathigraphs,
           output_folder,
           filter_paths,
           effector=(decompose==F)
         )
         
       },
       test = {
         cat("Predicting...\n")
         class0<-model$class0
         class1 <- model$class1
         test_preds <- predict.newdataset(model$model, as.data.frame(t(path.vals)), type=design_type,
                                          decision = T, proba = T, class0=class0)
         status("75")
         cat("Calculating prediction model statistics...\n")
         #test_stats <- get.stats(test_preds$pred_res$prediction, des$group, NULL, design_type, class0)
         
         #save.prednewdataset.res(test_preds, test_stats, path.vals, fpathigraphs, output_folder, class0, class1)
         save.prednewdataset.res(test_preds, path.vals, fpathigraphs, output_folder, class0, class1)
       }
)


status("95")

################## REPORT

cat("Creating report...\n")
switch(species,
       hsa = {speciesTitle <- "Human (Homo sapiens)"},
       mmu = {speciesTitle <- "Mouse (Mus musculus)"},
       rno = {speciesTitle <- "Rat (Rattus norvegicus)"})

results <- init.report(analysis)

results <- add.section(results,"Input parameters",0)

if(dataset=="gex"){
  results <- add.param(results,"Expression file",basename(exp_file),1)
}
if(analysis != "test"){
  results <- add.param(results,"Design file",basename(design_file),1)
  if( design_type != "continuous"){
    results <- add.param(results,"Comparison",paste0(cond1," vs ",cond2),1)
    
    results <- add.param(results,"Paired analysis",ifelse(paired, "Yes", "No"),1)
  } else {
    results <- add.param(results,"Comparison","Correlation with a continuous variable",1)
  }
  #results <- add.param(results,"Design file",paste0(basename(design_file)))
}
results <- add.param(results,"Species",speciesTitle,1)
#results <- add.param(results,"Decomposed paths",decompose)



switch(analysis,
       compare = {
         results <- add.section(results,"Pathways",0)
         results <- add.pathwayViewer(results, "session", "pathways", 1)
         results <- add.section(results,"Circuit values",0)
         results <- add.download(results,"Circuit values","paths_vals.txt",1)
         if( design_type != "continuous"){
           results <- add.image(results,"Heatmap","paths_heatmap.png",1)
           results <- add.image(results,"PCA","paths_pca.png",1)
         }
         results <- add.table(results,"Circuit significance","paths_significance.txt",1)
         
         #results <- add.html(results,"network.html")
         
         if(go==T | uniprot==T){
           results <- add.section(results,"Function based analysis",0)
           if(go==T){
             results <- add.download(results,"GO terms values","go_vals.txt",1)
             if( design_type != "continuous"){
               results <- add.image(results,"Heatmap","go_heatmap.png",1)
               results <- add.image(results,"PCA","go_pca.png",1)
             }
             results <- add.table(results,"GO term significance","go_significance.txt",1)
           }
           if(uniprot==T){
             results <- add.download(results,"Uniprot keywords values","uniprot_vals.txt",1)
             if( design_type != "continuous"){
               results <- add.image(results,"Heatmap","uniprot_heatmap.png",1)
               results <- add.image(results,"PCA","uniprot_pca.png",1)
             }
             results <- add.table(results,"Uniprot keyword significance","uniprot_significance.txt",1)
           }
         }
       },
       train = {
         results <- add.section(results,"Circuit values",0)
         results <- add.download(results,"Circuit values","paths_vals.txt",1)
         results <- add.section(results, "Model training and evaluation",0)
         #k-fold cross validation: k
         results <- add.param(results,paste0(k,"-fold cross-validation"),k,1)
         #La gráfica de svm.performance.heatmap.png
         #results <- add.image(results,"SVM performance heatmap","svm.performance.heatmap.png",1)
         #Optimal hyperparameters: One random split analysis (20% test)
         #results <- add.param(results,"Optimal hyperparameters", "One random split analysis (20% test)",1)
         #Las gráficas y el output en .TXT del directorio split en output_folder.
         #el archivo model.stats
         #results <- add.table(results,"Test model statistics","split/test_model_stats.txt",1)
         results <- add.param(results,"Validation of typical split","",1)
         results <- add.image(results,
                              "Split train precision and recall",
                              "split/split_Train_pr.png",1)
         results <- add.image(results,
                              "Split train receiver operating characteristic",
                              "split/split_Train_roc.png",1)
         #el archivo model.stats
         #results <- add.table(results,"Test model statistics","split/test_model_stats.txt",1)
         results <- add.image(results,
                              "Split test precision and recall",
                              "split/split_Test_pr.png",1)
         results <- add.image(results,
                              "Split test receiver operating characteristic",
                              "split/split_Test_roc.png",1)
         results <- add.section(results,"Probability distribution",0)
         results <- add.image(results,
                              "Test probability boxplot",
                              "split/test_probability_boxplot.png",1)
         #results <- add.section(results,"Prediction model",0)
         #results <- add.param(results,"K-fold cross-validation",k,1)
         #results <- add.table(results,"Model statistics","model_stats.txt",1)
         if(filter_paths) {
           #results <- add.table(results,"Selected features","filtered_features_cfs.txt",page_size = 14,1)
           results <- add.table(results,"The most relevant circuits along with their interaction sign:","relevances.tsv",page_size = 20,1)
           #results <- add.html(results,"network.html",1)  sifs4CellMaps
           #results <- add.pathwayViewer(results, "session", "pathways", 1)
         }
       },
       test = {
         results <- add.section(results,"Circuit values",0)
         results <- add.download(results,"Circuit values","paths_vals.txt",1)
         results <- add.section(results,"Prediction model",0)
         results <- add.table(results,"Prediction results","prediction_results.txt",1)
         #results <- add.download(results,"Model statistics","prediction_stats.txt",1) 
         results <- add.table(results,"Model statistics","prediction_stats.txt",1)   
       }
)
write(render.xml(results),file=paste0(output_folder,"/report.xml"))
#unlink(paste0(output_folder,"/sifs4CellMaps"),recursive = T)

if(exists("wt")){
  wt$status <- wt$"UP/DOWN"
  wt$has_changed <- T
  path_json <- create.path.info(wt,metaginfo$pathigraphs)
  write(path_json,file=paste0(output_folder,"/pathways/path_info.json"))
}

#cp_command <- paste0("cp -r ",codebase,"/report-files/pretty_legend.png ",output_folder)
#system(cp_command)


cat("[Finished]\n")

status("100")




# # SAVE subgraph info
# subgraphs <- unlist(lapply(metaginfo$pathigraph,function(x) x$effector.subgraphs),recursive = F)
# subgraph_nodes <- lapply(subgraphs,function(subgraph){(cbind(node=V(subgraph)$name,label=V(subgraph)$label,entrezs=sapply(V(subgraph)$genesList,function(x){x<-setdiff(x,"/");paste(x,collapse=",")})))})
# subgraph_nodes2 <- as.data.frame(do.call("rbind",mapply(USE.NAMES = T,function(x,y) (cbind(path=rep(x,nrow(y)),y)),names(subgraph_nodes),subgraph_nodes)))
# subgraph_nodes2$path <- gsub("\\.","__",subgraph_nodes2$path)
# rownames(subgraph_nodes2) <- NULL
# write.table(subgraph_nodes2,file="path_nodes.txt",row.names=F,col.names=T,quote=F,sep="\t")





