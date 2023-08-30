cat("Welcome to MetaboPathia")
# library(igraph, warn.conflicts = F)
# library("e1071", warn.conflicts = F)
# suppressPackageStartupMessages(library(R.utils, warn.conflicts = F))
library(hipathia, warn.conflicts = F)
#### INPUT DATA
# args <- commandArgs(trailingOnly = F, asValues = T,excludeEnvVars = F)
codebase <- paste0(dirname(normalizePath(args[["file"]])),"/")

# species
species <- args[['species']]

# expression
exp_file <- args[['exp_file']]
dataset <- args[["exp_db"]] # gex or ea or hpa
tissue <- args[["tissue"]]  #if ea or hpa selected as dataset

# variants
vcf <- args[["variant_file"]]
ref_genome <- args[["ref_genome"]] # hg19 or hg38
inheritance_pattern <- args[["inheritance_pattern"]] # recessive or dominant
compound_hetero <- args[["compound_hetero"]]
sift <- args[["sift"]]
polyphen <- args[["polyphen"]]
phastcons <- args[["phastcons"]]
annovar_dir <- args[["annovar_dir"]]
annovar_db <- args[["annovar_db"]]
vcftools_dir <- args[["vcftools_dir"]]

# analysis
analysis <- args[["analysis"]] # compare, train, test
# do not adjust pvalue
if(is.null(args[["adjusted"]])){
  unadjusted <- TRUE
}else unadjusted <- FALSE
paired <- args[["paired"]] # do not adjust pvalue

# design
design_type <- args[['design_type']] #categorical, continuous or mutation_effect
design_file <- args[['design_file']]
cond1 <- args[['cond1']]
cond2 <- args[['cond2']]

# method parameters
pathways_list <- args[['pathways_list']]
decompose <- args[['decompose']]
difexp <- args[['difexp']]



# train
Rfilter_paths <- args[['filter_paths']] # TRUE or FALSE
# test
model_file <- args[['model_file']]

# functional analysis
go <- args[['go']]
uniprot <- args[['uniprot']]

# output
output_folder <- paste(args[['output_folder']], collapse=" ")

verbose <- args[['verbose']]
report <- args[['report']]

# 
# #/mnt/httpd/bioinfo/opencga-0.7/tools/pathways/differential-signaling.sh --difexp  --ref_genome hg19 --pathways_list 04014,04015,04010,04012,04310,04330,04340,04350,04390,04370,04630,04064,04668,04066,04068,04020,04071,04024,04022,04151,04152,04150,04110,04114,04210,04115,04510,04520,04530,04540,04611,04620,04621,04622,04650,04660,04662,04664,04666,04670,04062,04910,04922,04920,03320,04912,04915,04914,04921,04919,04916,04261,04270,04722,05200,05231,05202,05205 --cond1 YRI --design_type categorical --sift  --output_folder /mnt/httpd/bioinfo/opencga-0.7/sessions/jobs/J_qA6tOTxLmC/ --inheritance_pattern recessive --species hsa --exp_file /mnt/httpd/bioinfo/opencga-0.7/sessions/users/test1/projects/23837/23838/sample_1000g_GBR_vs_YRI__exp.txt --design_file /mnt/httpd/bioinfo/opencga-0.7/sessions/users/test1/projects/23837/23838/sample_1000g_GBR_vs_YRI__design.txt --cond2 GBR --variant_file /mnt/httpd/bioinfo/opencga-0.7/sessions/users/test1/projects/23837/23838/sample_1000g_GBR_vs_YRI.vcf
# hipathia_home <- "."
# go <- NULL
# uniprot <- NULL
# analysis <- "compare"
# tissue <- NULL
# difexp <- T
# ref_genome <- "hg19"
# pathways_list <- "04014,04015,04010,04012,04310,04330,04340,04350,04390,04370,04630,04064,04668,04066,04068,04020,04071,04024,04022,04151,04152,04150,04110,04114,04210,04115,04510,04520,04530,04540,04611,04620,04621,04622,04650,04660,04662,04664,04666,04670,04062,04910,04922,04920,03320,04912,04915,04914,04921,04919,04916,04261,04270,04722,05200,05231,05202,05205"
# cond1 <- "YRI"
# design_type <- "categorical" 
# sift <- T
# output_folder <- "/mnt/data2/ngs2/web_pretty/last_test2/1000g_mutation_effect"
# inheritance_pattern <- "recessive"
# species <- "hsa"
# compound_hetero <- NULL
# exp_file <- "~/code/prettyways/examples/sample_1000g_GBR_vs_YRI__exp.txt"
# design_file <- "~/code/prettyways/examples/sample_1000g_GBR_vs_YRI__design.txt"
# cond2 <- "GBR"
# vcf <- "~/code/prettyways/examples/sample_1000g_GBR_vs_YRI.vcf"
# annovar_dir <- "~/annovar_for_pathivar2015/"

## EXAMPLE 1
# hipathia_home <- "."
# species <- "hsa"
# exp_file <- "examples//brca_example1_40__exp.txt"
# design_file <- "examples//brca_example1_40__design.txt"
# cond1 <- "Tumor"
# cond2 <- "Normal"
# decompose <- F
# difexp <- T
# go <- T
# uniprot <- T
# output_folder <- "/mnt/data2/ngs2/web_pretty/test5000"

## EXAMPLE 2
# hipathia_home <- "."
# species <- "hsa"
# exp_file <- "pathivar2015/1000gexampledata/GD462RPKM.txt"
# vcf <- "pathivar2015/1000gexampledata/1000g_all_chr_pathigraphgenes.vcf"
# design_file <- "pathivar2015/1000gexampledata/1000gexp_design.txt"
# cond1 <- "GBR"
# cond2 <- "YRI"
# annovar_dir <- "~/annovar_for_pathivar2015/"
# output_folder <- "/mnt/data2/ngs2/web_pretty/prettyvar/test1"

## EXAMPLE 3
# hipathia_home <- "."
# species <- "hsa"
# pathways_list <- "04210,04010"
# exp_file <- "pathivar2015/1000gexampledata/GD462RPKM.txt"
# design_file <- "pathivar2015/1000gexampledata/1000gexp_design.txt"
# cond1 <- "GBR"
# cond2 <- "YRI"
# output_folder <- "/mnt/data2/ngs2/web_pretty/prettyvar/test1"

## EXAMPLE 4
# hipathia_home <- "."
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

if(is.null(report)) report <- F

if(is.null(species)) species <- "hsa"

if(is.null(decompose)) decompose <- FALSE

if(is.null(go)) go <- F
if(is.null(uniprot)) uniprot <- F
if(is.null(difexp)) difexp <- "difexp"
#if(is.null(difexp)) difexp <- T
if(is.null(analysis)) analysis <- "compare"
if(is.null(unadjusted)) unadjusted <- F
if(is.null(paired)) {
  paired <- F
} else {
  paired <- T
}

if(is.null(cond1) | is.null(cond2)) {
  cond1 <- "Tumor"
  cond2 <- "Normal"
}
if(!is.null(exp_file)){
  dataset <- "gex"
} else if(is.null(dataset)) stop("ERROR: no expression is defined")

# variant defaults
if(is.null(ref_genome)) ref_genome <- "hg19"
if(is.null(inheritance_pattern)) inheritance_pattern <- "dominant"

if(is.null(sift) & is.null(polyphen) & is.null(phastcons)) {
  consequence <- c("SIFT","PolyPhen")
} else {
  consequence <- c()
  if(!is.null(sift)) consequence <- c(consequence,"SIFT")
  if(!is.null(polyphen)) consequence <- c(consequence,"PolyPhen")
  if(!is.null(phastcons)) consequence <- c(consequence,"PhastCons")  
}

if(is.null(compound_hetero)) compound_hetero <- T
if(dataset!="gex" & is.null(vcf)) stop("ERROR: vcf or expression file must be defined")

#Modified By kinza
#if(is.null(Rfilter_paths)) filter_paths <- F
#if(!is.null(Rfilter_paths)) filter_paths <- T
# Hardcoded 
filter_paths <- F

if(is.null(verbose)) verbose <- F

if(is.null(design_type)) design_type <- "categorical"

if(verbose==T){
  cat("\nWELCOME TO hiPathia 2!!!\n")
  cat("\nDetected parameters: \n")
  cat("\t\thipathia_home: ",hipathia_home,"\n")
  cat("\t\toutput_folder: ",output_folder,"\n")
  cat("\t\tspecies: ",species,"\n")
  cat("\t>>>Expression:\n")
  cat("\t\texp_file: ",exp_file,"\n")
  cat("\t\texpression: ",dataset,"\n")
  cat("\t\ttissue: ",tissue,"\n")
  cat("\t>>>Variants:\n")
  cat("\t\tvcf: ",vcf,"\n")
  cat("\t\tref_genome: ",ref_genome,"\n")
  cat("\t\tinheritance_pattern: ",inheritance_pattern,"\n")
  cat("\t\tcompound_hetero: ",compound_hetero,"\n")
  cat("\t\tconsequence: ",consequence,"\n")
  cat("\t\tannovar_dir: ",annovar_dir,"\n")
  cat("\t\tannovar_DB: ",annovar_db,"\n")
  cat("\t\tvcftools_dir: ",vcftools_dir,"\n")
  cat("\t>>>Experimental design:\n")
  cat("\t\tanalysis: ",analysis,"\n")
  cat("\t\tdesign_type: ",design_type,"\n")
  cat("\t\tdesign_file: ",design_file,"\n")
  cat("\t\tcond1: ",cond1,"\n")
  cat("\t\tcond2: ",cond2,"\n")
  cat("\t>>>Method parameters:\n")
  cat("\t\tpathways_list: ",pathways_list,"\n")
  cat("\t\tdecompose: ",decompose,"\n")
  cat("\t\tdifexp: ",difexp,"\n")
  cat("\t>>>Train and test:\n")
  cat("\t\tRECIEVED filter_paths: ",Rfilter_paths,"\n")
  cat("\t\tfilter_paths: ",filter_paths,"\n")
  cat("\t\tmodel_file: ",model_file,"\n")
  cat("\t>>>Functional analysis:\n")
  cat("\t\tgo: ",go,"\n")
  cat("\t\tuniprot: ",uniprot,"\n")
  cat("\n\n")
}

#### PREPARE DATA

#source(paste0(hipathia_home,"/prettyways.R"))
#source(paste0(hipathia_home,"/stats.R"))
#source(paste0(hipathia_home,"/functions.r"))
source(paste0(hipathia_home,"/report.r"))

cat("sourcing VCF and Prediction requirement ... \n")
source(paste0(hipathia_home,"/vcf2genes.R"))
source(paste0(hipathia_home,"/vcf2othertools.R"))
source(paste0(hipathia_home,"/predict.r"))
cat("Ready to call the hipathia2 package")
cat("\n\n")

dir.create(output_folder, showWarnings = FALSE)
status <- function(value){
  write(value,file=paste0(output_folder,"/status.txt"))
}
status("0")

# load and filter graphs
cat("all pathways loading begin...\n\n")
options(warn=-1)
pathways <- load_pathways(species = species)
cat("all pathways loading Done\n\n")
#pathways_list = c("hsa03320", "hsa04012")

#load(paste0(hipathia_home,"/files/meta_graph_info_", species, ".RData"))

if(analysis == "test"){
  load(model_file)
  pathways_list <- model$pathways_list
  design_type <- model$design_type
  if(!is.null(vcf) & !is.null(model$ref_genome)){
    ref_genome <- model$ref_genome
    consequence <- model$consequence
    inheritance_pattern <- model$inheritance_pattern
    compound_hetero <- model$compound_hetero
  }
}
if(!is.null(pathways_list)){
  pathways_list <- gsub(" ","",unlist(strsplit(pathways_list,",")))
  pathways_list <- paste0(species,pathways_list)
  if(length(setdiff(pathways_list,names(pathways$pathigraphs)))>0){
    #stop("ERROR: detected strange pathway id")
    print("ERROR: detected strange pathway id")
    pathways_list <- intersect(pathways_list,names(pathways$pathigraphs))
  }
  cat("just selected pathways loading begin...\n\n")
  pathways <- load_pathways(species = species,
                            pathways_list = pathways_list)
  cat("just selected pathways loading Done\n\n")
  
  
}
options(warn=0)
status("4")


#### LOAD DATA

cat("Loading data...\n")

# load expression
if(dataset=="gex"){
  exp <- read.table(exp_file,header=T,sep="\t",stringsAsFactors=F,row.names = 1,comment.char="",check.names=F)
  trans_data <- hipathia::translate_data(data.matrix(exp),species)
  exp_data <- normalize_data(trans_data, by_quantiles =(analysis=="train" | analysis=="test"), by_gene = F,percentil = F) 
  
}
if(!is.null(vcf)){
  vcf2othertools(vcf,vcf_outdir=output_folder,ref_genome=ref_genome,annovar_dir=annovar_dir,annovar_db=annovar_db,vcftools_dir=vcftools_dir)
  exp_vcf <- vcf2genes(output_folder,dataset,tissue,inheritance_pattern,compound_hetero,exp_data,consequence,pathigraph.genes=pathigraph.genes)
  if( design_type != "mutation_effect"){
    exp_data <- exp_vcf
  } else {
    colnames(exp_vcf) <- paste0( colnames(exp_vcf), "_vcf")
    exp_data <- cbind(exp_data, exp_vcf)
  }
}

# load design
if( analysis != "test"){
  if(design_type == "mutation_effect"){
    des <- as.data.frame(cbind(colnames(exp_data), rep(c("exp_data", "exp_vcf"), each=(ncol(exp_data)/2))), stringsAsFactors=F)
    #if(paired==F){
    colnames(des) <- c("sample", "group")
    #} else {
    #colnames(des) <- c("sample", "group", "donor")
    #}
    rownames(des) <- des$sample
    cond1 <- "exp_vcf"
    cond2 <- "exp_data"
  }else{
    des <- read.table(design_file,header=F,stringsAsFactors=F)
    #if(paired==F){
    colnames(des) <- c("sample",c("group","value")[(design_type == "continuous")+1])
    #} else {
    #  colnames(des) <- c("sample",c("group","value")[(design_type == "continuous")+1], "donor")
    #}
    rownames(des) <- des$sample
    
    if(design_type == "categorical"){
      des <- des[ des$group==cond1 | des$group==cond2, ]  
    }
    
    sel_samples <- intersect(colnames(exp_data),rownames(des))
    exp_data <- exp_data[,sel_samples]
    des <- des[sel_samples,]
    
    #if(paired==T){
    #  sample_order <- order(des$donor)
    #  ordered_samples <- des[sample_order,"sample"]
    #  exp_data <- exp_data[,ordered_samples]
    #  des <- des[ordered_samples,]
    #}
    print(des)
    print(head(exp_data))
    
  }
} else {
  des <- model$des
}

#### PREPROCESS DATA

cat("Scaling expression values...\n")

status("20")

#### RUN

##############################
# Hipathia
##############################
cat("Propagating signaling...\n")
results <- hipathia(exp_data, pathways, decompose = decompose, verbose = TRUE)


status("50")


## LOAD FUNCTIONS

#if(go==T | uniprot==T){
#  entrez2hgnc <- read.table( paste0(hipathia_home, "/files/annotations/", species, "/entrez_hgnc_", species, ".annot"),header=F,sep="\t",stringsAsFactors=F)
#  if(go==T){
#    go_annot <- load.annot.file( paste0(hipathia_home, "/files/annotations/", species, "/go_bp_", species, ".annot"))
#    go.vals <- prettyfuns(results, fpathigraphs, go_annot, entrez2hgnc)
#  }
#  if(uniprot==T){
#  	uniprot_annots <- get_pathways_annotations(rownames(comp), pathways, "uniprot")
#    uniprot_annot <- load.annot.file( paste0(hipathia_home, "/files/annotations/", species, "/uniprot_keywords_", species, "__biological_process.annot"))
#    fun.vals <- prettyfuns(results, fpathigraphs, uniprot_annot, entrez2hgnc)
#  }
#}


#### ANALYSIS

# Descriptive plots
#---------------------------------
#path_vals <- get_paths_matrix(results)
path_vals <- get_paths_data(results,matrix=TRUE )
# Define groups to compare
sample_group <- des[colnames(path_vals),"group"]



n <- ncol(path_vals)
r <- nrow(path_vals)
# wilcoxon
cat("Comparing groups...\n")

## ANALYSIS

switch(analysis,       
       compare = {
         #dir.create(paste0(output_folder,"/report/"))
         if( design_type == "continuous"){
           # correlation
           cat("Computing correlation...\n")
           
           wt <- do_cor(path_vals, des,adjust=(unadjusted==F))
           path_names <- get_path_names(pathways, rownames(wt))
           #wt <- cbind(path_names, wt)
           
           
           write.table(path_vals,file=paste0(output_folder,"/paths_vals.txt"),sep="\t",quote=F,row.names=T,col.names=T)
           write.table(wt,file=paste0(output_folder,"/paths_significance.txt"),sep="\t",quote=F,row.names=F,col.names=T)
           
           
           if(go==T){
             # GO terms
             #-----------
             go_vals <- quantify_terms(results, pathways, "GO")
             go_wt <- do_cor(go_vals,des,adjust=(unadjusted==F))
             go_wt$path <- rownames(go_wt)
             write_output_cor(go_vals,go_wt,go_pca_mod,paste0(output_folder,"/go"))
           }  
           if(uniprot==T){
             uniprot_annots <- get_pathways_annotations(rownames(comp), pathways, "uniprot")
             # Uniprot Keywords
             #------------------
             uniprot_vals <- quantify_terms(results, pathways, "uniprot")
             
             fun_wt <- do_cor(uniprot_vals, des,adjust=(unadjusted==F))
             fun_wt$path <- rownames(fun_wt)
             write_output_cor(uniprot_vals,fun_wt,fun_pca_mod,paste0(output_folder,"/uniprot"))
           }
           
         } else {
           write_output <- function(vals,comp,pca_mod,prefix,colors=c("red","blue"),sample_colors=NULL){
             if(class(vals)=="SummarizedExperiment"){
               vals<-assays(vals)$path
             }
             print("class(vals)")
             print(class(vals))
             write.table(vals,file=paste0(prefix,"_vals.txt"),sep="\t",quote=F,row.names=T,col.names=T)
             #out_comp <- comp[,c("path","UP/DOWN","statistic","p.value","FDRp.value")]
             #colnames(out_comp) <- gsub("path","path/term",colnames(out_comp))
             #out_comp <- out_comp[order(out_comp$p.value,decreasing=F),]
             #write.table(out_comp,file=paste0(prefix,"_significance.txt"),sep="\t",quote=F,row.names=F,col.names=T)
             write.table(comp,file=paste0(prefix,"_significance.txt"),sep="\t",quote=F,row.names=F,col.names=T)
             ranked_path_vals <- vals[order(comp$p.value, decreasing = FALSE),]
             #cat("ranked_path_vals[comp$FDRp.value < 0.05,]")
             
             #Volcano of Fold change vs Pvalue
             
             #options(bitmapType='cairo')
             #cat("Volcano of Fold change vs Pvalue")
             #EnhancedVolcano(comp,'','Fold_Change','p.value')
             print(head(comp))
             #png(paste0(prefix,"_Volcano.png"),width=800,height=600)
             
             #EnhancedVolcano(comp, lab = comp$path_names, x = 'logFC',y = 'p.value',xlim = c(-4, 4), FCcutoff = 1)
             #EnhancedVolcano(comp,'','Fold_Change','p.value')
             
             #dev.off()
             print("Volcano Done :P")
             
             if(nrow(vals)>=2){
               options(bitmapType='cairo')
               #paste0(output_folder,"/go","_heatmap.png")
               png(paste0(prefix,"_heatmap.png"),width=800,height=800)
               if(length(ranked_path_vals[comp$FDRp.value < 0.05,1])>=2){
                 heatmap_plot(ranked_path_vals[comp$FDRp.value < 0.05,],sample_group,colors = "hipathia", variable_clust = TRUE,labRow=rownames(ranked_path_vals[comp$FDRp.value < 0.05,]), labCol=colnames(ranked_path_vals[comp$FDRp.value < 0.05,]))# Hipathia web + filter by selected groups
               }else{
                 heatmap_plot(ranked_path_vals,sample_group,colors = "hipathia", variable_clust = TRUE,labRow=rownames(ranked_path_vals), labCol=colnames(ranked_path_vals))
               }
               dev.off()
               png(paste0(prefix,"_pca.png"),width=800,height=600)
               if(nrow(vals)==2){
                 multiple_pca_plot(pca_model,sample_group, cex = 3, plot_variance = TRUE, comps = seq_len(2))
               }else{
                 multiple_pca_plot(pca_mod,sample_group, cex = 3, plot_variance = TRUE)
               }
               dev.off()
               
             }else{
               print("nothing to do here :P")
             }  
           }
           
           # design
           group.value <- des[colnames(exp_data),"group"] #50B7AE #B6EBE7
           colors <- c("#204253", "#ea8c39")[(group.value==cond2)+1]
           sample_colors <- c("#204253", "#ea8c39")
           names(sample_colors) <- c(cond1,cond2)
           
           #Fold change
           #Function to be moved to other statFile.r
           #FC and logFC function
           get_FC<-function(path_vals,info,case,control){
             print("info")
             print(info)
             print("path_vals")
             print(path_vals[1,])
             logFC<-mat.or.vec(nrow(path_vals),2)
             rownames(logFC)<-rownames(path_vals)
             colnames(logFC)<-c("Fold_Change","logFC")
             cat("path_vals")
             cat(head(path_vals))
             cat("names(info[info==case])")
             cat(head(names(info[info==case])))
             cat("names(info[info==control])")
             cat(head(names(info[info==control])))
             logFC[,"Fold_Change"]<-rowMeans(path_vals[,names(info[info==case])])/rowMeans(path_vals[,names(info[info==control])])
             logFC[,"logFC"]<-log2(logFC[,"Fold_Change"])
             return(logFC) 
           }
           ####
           info<-group.value
           names(info)<-colnames(path_vals)
           fc<-get_FC(path_vals,info,case=cond1, control=cond2)
           print("fc")
           print(fc)
           
           wt <- do_wilcoxon(path_vals, sample_group, g1=cond1, g2=cond2,adjust=(unadjusted==F),paired=paired)
           #add FC col to wt
           wt <- cbind(wt, fc[rownames(wt),])
           print("wt")
           print(head(wt))
           
           path_names <- get_path_names(pathways, rownames(wt))
           wt <- cbind(path_names, wt)
           
           pathways_summary <- get_pathways_summary(wt, pathways)
           pathways_summary <- data.frame(names = row.names(pathways_summary), pathways_summary)
           
           colnames(pathways_summary)<- c("Pathway_name", "id_pathways", "num_total_paths", "num_significant_paths", "percent_significant_paths", "num_up_paths", "percent_up_paths", "num_down_paths", "percent_down_paths")
           #print(colnames(pathways_summary))
           cat("pathways_summary : \n")
           print(pathways_summary)
           
           head(pathways_summary, n = 15)
           write.table(pathways_summary,file=paste0(output_folder,"/pathways_summary_vals.txt"),sep="\t",quote=F,row.names=F,col.names=T)
           
           #summary pathways to visualize it
           
           pathways_summary_table <- data.frame(
             pathways_summary[,1],
             pathways_summary[,2],
             trunc(pathways_summary[,3]),
             paste0("<paper-progress value='",pathways_summary[,5],"' title='",pathways_summary[,5],"%' style='width:100%' ></paper-progress>"),
             paste0("<paper-progress class='red paper-progress-0'  title='Up: ",pathways_summary[,7],"% Vs Down : ",pathways_summary[,9],"%' value='",pathways_summary[,7],"' secondary-progress='",pathways_summary[,9]+pathways_summary[,7],"'  style='width:100%' ></paper-progress>")
           )
           colnames(pathways_summary_table)<- c("Pathway_name", "Pathway_id", "Total_paths", "Significant_paths",  "Up Vs down")
           write.table(pathways_summary_table,file=paste0(output_folder,"/pathways_summary_table.txt"),sep="\t",quote=F,row.names=F,col.names=T)
           print(head(pathways_summary_table[,3]))
           #end of summary pathways 
           
           
           ranked_path_vals <- path_vals[order(wt$p.value, decreasing = FALSE),]
           rank <- wt[order(wt$p.value, decreasing = FALSE), ]
           table(wt$FDRp.value < 0.05)
           
           # PCA
           #--------------
           # Perform PCA model
           pca_model <- do_pca(ranked_path_vals[1:min(nrow(ranked_path_vals), ncol(ranked_path_vals)),]) # add a  filter by selected groups
           
           #pca_model <- do_pca(ranked_path_vals[1:ncol(ranked_path_vals),])
           #pca_mod <- do_pca(path_vals[rank[1:(min(20,n))],],cor=(dataset=="gex"))      
           
           write_output(path_vals,wt,pca_model,paste0(output_folder,"/paths"),colors=colors,sample_colors=sample_colors)
           if(go==T){
             go_vals <- quantify_terms(results, pathways, "GO", out_matrix=TRUE)        
             go_wt <- do_wilcoxon(go_vals, group.value, g1=cond1, g2=cond2,adjust=(unadjusted==F),paired=paired)
             go_wt$path <- rownames(go_wt)
             go_rank <- order(go_wt$p.value,decreasing=F)
             ngo <- min(c(20,ncol(go_vals),nrow(go_vals)))
             go_pca_mod <- do_pca(go_vals[go_rank[1:ngo],],cor=(dataset=="gex"))
             write_output(go_vals,go_wt,go_pca_mod,paste0(output_folder,"/go"),colors=colors,sample_colors=sample_colors)
             HSGOA<-paths_to_go_ancestor(pathways, wt, go_wt, pval = 0.05)
             write.table(HSGOA,file=paste0(output_folder,"/HSGOA.txt"),sep="\t",quote=F,row.names=F,col.names=T)
             
           }
           if(uniprot==T){
             uniprot_annots <- get_pathways_annotations(rownames(comp), pathways, "uniprot")
             # Uniprot Keywords
             #------------------
             uniprot_vals <- quantify_terms(results, pathways, "uniprot", out_matrix=TRUE)
             
             fun_wt <- do_wilcoxon(uniprot_vals, group.value, g1=cond1, g2=cond2,adjust=(unadjusted==F),paired=paired)
             fun_wt$path <- rownames(fun_wt)
             fun_rank <- order(fun_wt$p.value,decreasing=F)
             nfun <- min(c(20,ncol(uniprot_vals),nrow(uniprot_vals)))
             print("uniprot_vals")
             print(uniprot_vals)
             print("length(uniprot_vals)")
             print(length(uniprot_vals[,1]))
             print("uniprot_vals[fun_rank[1:nfun],]")
             print(uniprot_vals[fun_rank[1:nfun],])
             if(length(uniprot_vals[,1])>1){
               fun_pca_mod <- do_pca(uniprot_vals[fun_rank[1:nfun],],cor=(dataset=="gex"))
             }else{
               fun_pca_mod <-NULL
             }
             write_output(uniprot_vals,fun_wt,fun_pca_mod,paste0(output_folder,"/uniprot"),colors=colors,sample_colors=sample_colors)
           }
           
         }
         
         status("70")
         ##############################
         # Pathway visualization
         ##############################
         
         if(difexp == "difexp"){
           # Define node colors
           colors_de <- node_color_per_de(results, pathways, sample_group, "Tumor","Normal")
           #colors_de_hipathia <- node_color_per_de(results, pathways, sample_group,"Tumor", "Normal", colors = "hipathia")
           create_report(wt, pathways, path=output_folder,output_folder="", node_colors = colors_de)
         }else{
           if(difexp == "difexpuni"){
             # Node colors with Uniprot grouping
             #colors_uni <- node_color_per_de(results, pathways, sample_group, "Tumor", "Normal", group_by = "uniprot")
             #create_report(wt, pathways, path=output_folder,output_folder="", node_colors = colors_uni)
             
             colors_de_uni <- node_color_per_de(results, pathways, sample_group, "Tumor","Normal", group_by = "uniprot")
             create_report(wt, pathways, path=output_folder,output_folder="",node_colors = colors_de_uni, group_by = "uniprot")
             
           }else{
             if(difexp == "difexpgp"){
               # Node colors with genes grouping
               colors_gen  <- node_color_per_de(results, pathways, sample_group, "Tumor", "Normal", group_by = "genes")
               create_report(wt, pathways, path=output_folder,output_folder="", node_colors = colors_gen, group_by = "genes")
             }else{
               if(difexp == "difexpgo"){
                 # Node colors with genes grouping
                 colors_de_go  <- node_color_per_de(results, pathways, sample_group, "Tumor", "Normal", group_by = "GO")
                 create_report(wt, pathways, path=output_folder,output_folder="", node_colors = colors_de_go, group_by = "GO")
               }
             }
           }
         }
         #create_report(wt, pathways, output_folder, node_colors = colors_de)
         #save.results(results,wt,fpathigraphs,output_folder,effector=(decompose==F),moreatts=h)
         #if(report==T){
         #  create.html.report2(fpathigraphs,wt,hipathia_home,output_folder,effector=(decompose==F),template_name = "index_template.html",output_name = "index.html",clean_out_folder = F)
         #}
         
       },
       train = {
         
         if(n<=20) {k=3} else if (n<=30) {k=5} else {k=10}
         
         cat("Generating prediction model...\n")
         mod <- get.prediction.model(as.matrix(path_vals), des, design_type, k, filter_paths)
         mod$design_type <- design_type
         mod$pathways_list <- gsub(species, "", pathways_list)
         mod$des <- des
         #if(!is.null(vcf)){      
         mod$ref_genome <- ref_genome
         mod$consequence <- consequence
         mod$inheritance_pattern <- inheritance_pattern
         mod$compound_hetero <- compound_hetero
         #}
         
         status("75")
         
         cat("Calculating prediction model statistics...\n")
         mod_stats <- get.predmod.stats(mod, design_type)
         
         # save
         save.pred.res(mod,mod_stats,results,path_vals,fpathigraphs,output_folder,filter_paths, effector=(decompose==F))
         
       },
       test = {
         
         cat("Predicting...\n")
         
         predres <- predict.newdataset(model[["model"]], t(path_vals), des, design_type) 
         status("75")
         
         cat("Calculating prediction model statistics...\n")
         #pred_stats <- get.stats(predres[["svm_pred"]], model, design_type)
         
         # save
         save.prednewdataset.res(predres[["pred_res"]],NULL,results,path_vals,fpathigraphs,output_folder)
         
       }
)


status("95")

################## REPORT

cat("Creating report...\n")

results <- init.report(analysis)

results <- add.section(results,"Input parameters",0)

results <- add.param(results,"species",species)
if(dataset=="gex"){
  results <- add.download(results,"Expression file",basename(exp_file),1)
} else {
  results <- add.param(results,"Expression DB",paste0(dataset, " (",tissue,")"))  
}
if(!is.null(vcf)){
  results <- add.download(results,"VCF file",basename(vcf),1)
  results <- add.param(results,"Reference genome",ref_genome)
  results <- add.param(results,"Inheritance pattern",inheritance_pattern)
  results <- add.param(results,"Include compound heterozygotes",compound_hetero)
  results <- add.param(results,"Mutation impact criteria",paste(consequence,collapse=","))  
}
if(design_type!="mutation_effect" & analysis != "test"){
  results <- add.download(results,"Design file",basename(design_file),1)
  if( design_type != "continuous"){
    results <- add.param(results,"Comparison",paste0(cond1," vs ",cond2))
    results <- add.param(results,"Paired analysis",paired)
  } else {
    results <- add.param(results,"Comparison","Correlation with a continuous variable")    
  }    
  results <- add.param(results,"Design file",paste0(basename(design_file)))
}
results <- add.param(results,"Decomposed paths",decompose)

results <- add.section(results,"Path values",0)
results <- add.download(results,"Path values","paths_vals.txt",1)

switch(analysis,
       compare = {
         results <- add.section(results,"Pathways",0)
         results <- add.pathwayViewer(results, "session", "pathways", 1)
         results <- add.table(results,"Pathways summary","pathways_summary_table.txt",1, page_size = 18, originFile="pathways_summary_vals.txt")
         results <- add.table(results,"Path significance","paths_significance.txt",1)
         if( design_type != "continuous"){
           if(file.exists(paste0(output_folder,"/paths_heatmap.png"))){
             results <- add.image(results,"Heatmap","paths_heatmap.png",1)
           }
           if(file.exists(paste0(output_folder,"/paths_pca.png"))){
             results <- add.image(results,"PCA","paths_pca.png",1)
           }
           if(file.exists(paste0(output_folder,"/paths_Volcano.png"))){
             results <- add.image(results,"Volcano","paths_Volcano.png",1)
           }
         }
         #results <- add.table(results,"Path significance","paths_significance.txt",1)
         #results <- add.section(results,"Pathways",1)
         #results <- add.pathwayViewer(results, "session", "pathways", 1)
         
         #results <- add.html(results,"network.html")
         
         if(go==T | uniprot==T){
           results <- add.section(results,"Function based analysis",0) 
           if(go==T){
             results <- add.download(results,"GO terms values","go_vals.txt",1)
             if( design_type != "continuous"){
               if(file.exists(paste0(output_folder,"/go_heatmap.png"))){
                 results <- add.image(results,"Heatmap","go_heatmap.png",1)
               }
               if(file.exists(paste0(output_folder,"/go_pca.png"))){
                 results <- add.image(results,"PCA","go_pca.png",1)
               }
               if(file.exists(paste0(output_folder,"/go_Volcano.png"))){
                 results <- add.image(results,"Volcano","go_Volcano.png",1)
               }
             }
             results <- add.table(results,"GO term significance","go_significance.txt",1)
             results <- add.table(results,"Highest significant GO ancestor (HSGOA)","HSGOA.txt",1)
           }
           if(uniprot==T){
             results <- add.download(results,"Uniprot keywords values","uniprot_vals.txt",1)
             if( design_type != "continuous"){
               #paste0(output_folder,"/go","_heatmap.png")
               if(file.exists(paste0(output_folder,"/uniprot_heatmap.png"))){
                 results <- add.image(results,"Heatmap","uniprot_heatmap.png",1)            
               }
               if(file.exists(paste0(output_folder,"/uniprot_pca.png"))){
                 results <- add.image(results,"PCA","uniprot_pca.png",1)            
               }
               if(file.exists(paste0(output_folder,"/uniprot_Volcano.png"))){
                 results <- add.image(results,"Volcano","uniprot_Volcano.png",1)            
               }
             }
             results <- add.table(results,"Uniprot keyword significance","uniprot_significance.txt",1)
           }
         }
       },
       train = {
         results <- add.section(results,"Prediction model",0) 
         results <- add.param(results,"K-fold cross-validation",k)
         results <- add.table(results,"Model statistics","model_stats.txt",1)
         if(filter_paths) {      
           results <- add.table(results,"Selected features","filtered_features_cfs.txt",page_size = 14,1)
           results <- add.html(results,"network.html",1)
         }
       },
       test = {
         results <- add.section(results,"Prediction model",0) 
         results <- add.table(results,"Prediction results","prediction_results.txt",1)
         #results <- add.download(results,"Model statistics","prediction_stats.txt")    
       }
)

write(render.xml(results),file=paste0(output_folder,"/report.xml"))
#unlink(paste0(output_folder,"/sifs4CellMaps"),recursive = T)

if(exists("wt")){
  wt$status <- wt$"UP/DOWN"
  wt$has_changed <- T
  #path_json <- create_path_info(wt,pathways$pathigraphs)
  #write(path_json,file=paste0(output_folder,"pathway-viewer/pathways/path_info.json"))
  
  cp_command2 <- paste0("cp  -r ",output_folder,"/pathway-viewer/pathways ",output_folder)
  system(cp_command2)
}

#cp_command <- paste0("cp -r ",output_folder,"/pathway-viewer/report_legend.png ",output_folder)
#system(cp_command)


cat("[Finished]\n")

status("100")




# # SAVE subgraph info
# subgraphs <- unlist(lapply(pathways$pathigraph,function(x) x$effector.subgraphs),recursive = F)
# subgraph_nodes <- lapply(subgraphs,function(subgraph){(cbind(node=V(subgraph)$name,label=V(subgraph)$label,entrezs=sapply(V(subgraph)$genesList,function(x){x<-setdiff(x,"/");paste(x,collapse=",")})))})
# subgraph_nodes2 <- as.data.frame(do.call("rbind",mapply(USE.NAMES = T,function(x,y) (cbind(path=rep(x,nrow(y)),y)),names(subgraph_nodes),subgraph_nodes)))
# subgraph_nodes2$path <- gsub("\\.","__",subgraph_nodes2$path)
# rownames(subgraph_nodes2) <- NULL
# write.table(subgraph_nodes2,file="path_nodes.txt",row.names=F,col.names=T,quote=F,sep="\t")





