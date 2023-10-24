rm(list=ls())
# 
# library(preprocessCore)
# library(igraph)
# library(limma)
# library(R.utils,quietly=T,verbose=F,warn.conflicts=F)
# library("gplots")

midas_home <-  "."
species <- "hsa"
exp_file <- "../../data_examples/Dystrophic_epidermolysis_bullosa/counts_TMM_normalization.tsv"
met_file <- "../../data_examples/Dystrophic_epidermolysis_bullosa/metabolite_suero.tsv"
design_file  <- "../../data_examples/Dystrophic_epidermolysis_bullosa/integration_design.tsv"
cond1 <- "control"
cond2 <- "visit1"
output_folder <- "../../test/"
analysis <- "compare" #compare(ModuleActivity), pathact(InsilicoKO), train ot test(Prediction)
design_type <- "categorical" #"categorical", "continuous", "mutation_effect"
paired <- T  # T or F


# To delet
corvar <-  "BRCA1.expression"
methodcor <- "pearson"



if(is.null(design_type)) {design_type <- "categorical"}
if(is.null(methodcor)) {methodcor <- "spearman"}


source("pretty_functions_Sep2016.R")
# load("/KEGG_module/Analysis_May2016/Project_Secondment_ETH (fromcrom16)/KEGG_module_April2016/hsa_module_data_April2016.RData")
load(paste0(species,"_module_data_Dec2016.RData"))


# if check.names=T then TCGA-6D-AA2E-01A changes into TCGA.6D.AA2E.01A.
# aaa <- read.table("/brca_example1_40__exp_HGNC.txt",sep="\t",row.names = 1,header=T)
combat.vals <- as.matrix(read.table(exp_file,sep="\t",header=T,stringsAsFactors=F,check.names=F))
rownames(combat.vals) <- combat.vals[,1]
samplenames <- colnames(combat.vals)[-1]
combat.vals <- as.matrix(t(apply(combat.vals[,-1],1,as.numeric)))
colnames(combat.vals) <- samplenames
print("des:")
print(head(combat.vals))
# xref
if(species %in% c("hsa","mmu","rno")){
  
  load(paste0(midas_home,"/files/xref/",species,"/","xref.rdata"))
  newgenenames <- xref[rownames(combat.vals)]
  missinggenenames <- which(sapply(newgenenames, is.null))
  onetomanygenes <- which(sapply(newgenenames, function(x) length(x)>1))
  rownames(combat.vals) <- newgenenames
  
  for(omg in onetomanygenes){
    addgene <- newgenenames[[omg]]
    addgene_mat <- matrix(data =  rep(combat.vals[omg,],length(addgene)),nrow = length(addgene),byrow = T)
    rownames(addgene_mat) <-   addgene
    combat.vals <- rbind(combat.vals,addgene_mat)
  }
  if(length(missinggenenames)>0){
    combat.vals <- combat.vals[-missinggenenames,]
  }
}


if(design_type!="continuous"){
  cancer_types <- as.matrix(read.table(design_file,sep="\t",header=F,stringsAsFactors=F,check.names=F))
  print("cancer_types")
  print(cancer_types)
  rownames(cancer_types) <- cancer_types[,1]
  print("c(cond1,cond2)")
  print(c(cond1,cond2))
  idx_mycondsamples <- which(cancer_types[,2] %in% c(cond1,cond2))
  cancer_types<- cancer_types[idx_mycondsamples,]
  idx_mycondsamples <- match(cancer_types[,1],colnames(combat.vals ))
  idx_mycondsamplesNA <- which(is.na(idx_mycondsamples))
  if(length(idx_mycondsamplesNA)>0){
    idx_mycondsamples <- idx_mycondsamples[-idx_mycondsamplesNA]
    cancer_types <- cancer_types[-idx_mycondsamplesNA,]
  }
  combat.vals <- combat.vals[,idx_mycondsamples]
}else {
  cancer_types <- as.matrix(read.table(design_file,sep="\t",header=T,stringsAsFactors=F,check.names=F))
  rownames(cancer_types) <- cancer_types[,1]
  cancer_types <- data.frame(cancer_types[,1], as.numeric(cancer_types[,match(corvar,colnames(cancer_types))]),stringsAsFactors = F)
  colnames(cancer_types) <- c("V1",corvar)
  sharedsamples <- intersect(colnames(combat.vals),cancer_types[,1])
  combat.vals <- combat.vals[,match(sharedsamples,colnames(combat.vals))]
  cancer_types <- cancer_types[match(sharedsamples,cancer_types[,1]),]
}
print("cancer_types")
print(cancer_types)

################# All module genes ################# 

all_metabolites <- sort(unique(unlist(sapply(hsa_module_data, function(x) {nodes <- c(x$graphobject$SIF[,1],x$graphobject$SIF[,3])
metabolites <- nodes[c(grep("C",nodes),grep("G",nodes))]
return(metabolites)}))))

all_module_genes_vec <- get.All.module.genes(hsa_module_data)

######### get all module rxns ###########

all_module_rxn_vec <- get.All.module.rxns(hsa_module_data)

########### calculate rxn vals from gene exp ########### 

rxn_gene_mat <- get.rxn.gene.matrix(hsa_module_data)

gene_exp_mat <- get.gene.exp.of.module.genes(combat.vals, all_module_genes_vec, min.exp=0.0001)

rxn_vals_mat <- get.RXNvals.from.genes(all_module_rxn_vec, gene_exp_mat, rxn_gene_mat)
metabolite_matrix <- mat.or.vec(nr = length(all_metabolites), nc = ncol(rxn_vals_mat))
rownames(metabolite_matrix) <- all_metabolites
metabolite_matrix[,] <- 1

# probably "discard_compounds" is not any more important because I will use metabolites with default value 1.
#hsa_module_data <- discard_compounds(hsa_module_data)
results_module_activities <- methyways(hsa_module_data,rxn_vals_mat=rxn_vals_mat,expbased=T,fluxbased=F,verbose = F, default_value=0.5,metabolitematrix = metabolite_matrix)

#### do it as function
results_module_activities_list <- sapply(results_module_activities, function(x){ x[[1]][1] })
results_module_activities_matrix <- do.call("rbind",results_module_activities_list)

write.table(x=results_module_activities_matrix,file=paste0(output_folder,"/module_vals.txt"),quote=F,sep="\t")

idx_NA_discard <- which(apply(results_module_activities_matrix,1, function(x){ all(is.na(x)) }))
if(length(idx_NA_discard)>0){
  results_module_activities_matrix <- results_module_activities_matrix[-idx_NA_discard,]
}

if(1==0){
  results_module_node_vals_list <- sapply(results_module_activities, function(x){ x[[1]][3] })
  results_module_node_vals_matrix <- do.call("rbind",results_module_node_vals_list)
  allnodes <- unique(rownames(results_module_node_vals_matrix))
  allnodes <- allnodes[grep("R",allnodes)]
  results_module_node_vals_matrix <- results_module_node_vals_matrix[match(allnodes, rownames(results_module_node_vals_matrix)),]
}

######################### diff nodes and modules ######################### 
if(design_type!="continuous"){
  # difexp.nod <- compute.difexp(vals = results_module_activities_matrix, control.string=cond1, disease.string=cond2, experimental.design=cancer_types)
  difexp.nod <- do.wilcox(sel.vals = results_module_activities_matrix,group.value = cancer_types,g1 = cond1,g2 = cond2,paired = F,verbose = T)
  updown <- sapply(difexp.nod$statistic, function(x){if(x <0){"down"}else if(x > 0){"up"}else{"none"}})
  updown <- sapply(difexp.nod$statistic, function(x){if(x <0){"down"}else if(x > 0){"up"}else{"none"}})
  # col <- get.colors.from.pval(updown, difexp.nod$p.value)
  # difexp.nod$Status <- updown
  names(updown) <- rownames(results_module_activities_matrix)
}else{
  cor_res <- do.cor(sel.vals = results_module_activities_matrix,design = cancer_types,adjust = T,methodcor = "spearman")
  difexp.nod <- cor_res[,c(3,1,4)]
  updown <- sapply(difexp.nod$correlation, function(x){if(x <0){"down"}else if(x > 0){"up"}else{"none"}})
  updown <- sapply(difexp.nod$correlation, function(x){if(x <0){"down"}else if(x > 0){"up"}else{"none"}})
  names(updown) <- rownames(results_module_activities_matrix)
}

# module_info <- read.delim(paste0(midas_home,"/files/","95_module_list_w_4_sub_category.csv"),sep = "\t",stringsAsFactors=F)
# module_description_vec<- c()
# for(de in module_info$Description.Name){  
#   desc_split <- strsplit(de,", ")[[1]]
#   idx_no <- grep("=>",desc_split)  
#   if(length(idx_no)>0){
#     module_description_vec <- c(module_description_vec,paste0(desc_split[-idx_no],collapse=", "))
#   }else{module_description_vec <- c(module_description_vec,paste0(desc_split,collapse=", "))}
# }
# 
# rownames(difexp.nod)[rownames(difexp.nod)=="R00385"] <- "M00085_2"
# rownames(difexp.nod)[rownames(difexp.nod)=="R07162"] <- "M00085_1"
# rownames(difexp.nod)[rownames(difexp.nod)=="R07809"] <- "M00079_1"
# rownames(difexp.nod)[rownames(difexp.nod)=="R07810"] <- "M00079_2"
# rownames(difexp.nod)[rownames(difexp.nod)=="R04355"] <- "M00082_1"
# rownames(difexp.nod)[rownames(difexp.nod)=="R10707"] <- "M00082_2"

load(paste0(midas_home,"/files/moduleinfo.RData"))
load(paste0(midas_home,"/files/pathinfo.RData"))

difexp.nod$ModuleID <- rownames(difexp.nod)
difexp.nod$ModuleDescription <-  moduleinfo$ModuleDescription[match(difexp.nod$ModuleID ,moduleinfo$ModuleID)]
difexp.nod$StartMetaboliteID <- moduleinfo$StartMetaboliteID[match(difexp.nod$ModuleID ,moduleinfo$ModuleID)]
difexp.nod$StartMetaboliteName <- moduleinfo$StartMetaboliteName[match(difexp.nod$ModuleID ,moduleinfo$ModuleID)]
difexp.nod$EndMetaboliteID <- moduleinfo$EndMetaboliteID[match(difexp.nod$ModuleID ,moduleinfo$ModuleID)]
difexp.nod$EndMetaboliteName <- moduleinfo$EndMetaboliteName[match(difexp.nod$ModuleID ,moduleinfo$ModuleID)]
module_name_link<- sapply(difexp.nod$ModuleID, function(x)strsplit(x,"_")[[1]][1])
difexp.nod$KEGGlink <- paste0("http://www.genome.jp/kegg-bin/show_module?map=",module_name_link,"&org=",species,"&select_org=",species)
difexp.nod$Pathway <- path_mat[,3][match(difexp.nod$ModuleID ,path_mat[,1])]
if(design_type!="continuous"){
  #difexp.nod <- difexp.nod[,c(6,7,13,1,2,3,5,8,9,10,11,12)]
  difexp.nod <- difexp.nod[,c(5,6,12,1,2,3,4,7,8,9,10,11)]
  difexp.nod <- difexp.nod[order(difexp.nod$adj.p.value,decreasing = F),]
}else{
  difexp.nod <- difexp.nod[,c(4,5,11,1,2,3,6,7,8,9,10)]
  difexp.nod <- difexp.nod[order(abs(difexp.nod$correlation),decreasing = T),]
}

write.table(x=difexp.nod,file=paste0(output_folder,"/module_stats.txt"),quote=F,sep="\t",row.names = F)

if(design_type!="continuous"){
  #nodeCols <- node.color.per.differential.expression(results_module_node_vals_matrix,control.string=cond1, disease.string=cond2, experimental.design=cancer_types)
  nodeCols <- node.color.per.differential.expression(rxn_vals_mat,control.string=cond1, disease.string=cond2, experimental.design=cancer_types)
}else{
  nodeCols <- c()
}

source(paste0(midas_home,"/createJson.R"))
dir.create(paste0(output_folder,"/pathways/"),showWarnings=F)

# pathFrameModule <- get.pathFrame(hsa_module_data="/KEGG_module/Analysis_May2016/Project_Secondment_ETH (fromcrom16)/KEGG_module_April2016/hsa_module_data_April2016.RData")
pathFrameModule <- get.pathFrame(hsa_module_data=paste0(midas_home,"/files/midasspecies/",species,"/",species,"_module_data_Dec2016.RData"),midas_home,species)

createModuleJson(pathFrame=pathFrameModule$pathFrame, PathModuleMatrix=pathFrameModule$PathModuleMatrix, updown = updown,nodeCols = nodeCols,difexp.nod = difexp.nod, rxn_gene_mat=NULL,SIFdir=paste0(midas_home,"/files/midasspecies/",species,"/SIF/pathways/"),saveDir=paste0(output_folder,"/pathways/"))

######################## plots ######################## 
if(design_type!="continuous"){
  zero_var <- which(apply(results_module_activities_matrix,1,var)==0)
  if(length(zero_var)>0){
    results_module_activities_matrix <- results_module_activities_matrix[-zero_var,]
  }
  
  print(dim(results_module_activities_matrix))
  print(length(cancer_types[,2]))
  
  #png(paste0(output_folder,"/module_pca.png"),width=800,height=600)
  png(paste0(output_folder,"/module_pca.png"),width = 10, height = 10, units = 'in', res = 600)
  pca_mod <- do.pca.2(results_module_activities_matrix)
  my_colors <-  cancer_types[,2]
  my_colors[my_colors==cond1] <- "chocolate"
  my_colors[my_colors==cond2] <- "aquamarine4"
  plot.multiple.pca(pca_mod,1:4,colors=my_colors,legend = c(cond1,cond2), legend_colors = c("chocolate","aquamarine4"),cex = 1.2)
  dev.off()
  
  png(paste0(output_folder,"/module_heatmap.png"),width = 10, height = 10, units = 'in', res = 600)
  path.vals <- t(apply(results_module_activities_matrix,1,function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))))
  
  cexColX = 0.25 * sqrt(100/ncol(results_module_activities_matrix)+10)
  cexRowX = 0.22 * sqrt(100/nrow(results_module_activities_matrix)+10)
  idx_labRow <- match(rownames(results_module_activities_matrix),difexp.nod$ModuleID)
  labRowX = paste0(difexp.nod$ModuleID[idx_labRow],": ",difexp.nod$ModuleDescription[idx_labRow])
  sample_chr_length <- max(nchar(colnames(results_module_activities_matrix)))
  if(sample_chr_length>0){sample_chr_length<-30}
  #heatmap(path.vals,margins = c(8+(sample_chr_length*0.2),10),labRow = labRowX, labCol = colnames(results_module_activities_matrix), scale="none",Rowv=NA,Colv=T,ColSideColors = my_colors,col=colorRampPalette(c("blue","gray","red"))(256)) 
  heatmap.2(path.vals,margins = c(8+(sample_chr_length*0.2),10),cexCol = cexColX, cexRow = cexRowX ,labRow = labRowX, 
            labCol = colnames(results_module_activities_matrix), scale="none",Rowv=NA,Colv=T,ColSideColors = my_colors,
            col=colorRampPalette(c("blue","gray","red"))(256), trace="none", density="none", keysize = 1) 
  legend(x = -0.01, y=0.88,legend= c(cond1, cond2) ,col=c("chocolate", "aquamarine4"),pch=15,lwd=2,xpd=T,cex=1,border=NA,pt.cex=1,bty = "n")
  
  dev.off()
}

############# Report ###########

source(paste0(midas_home,"/report.r"))

report <- init.report("module-activity")
report <- add.section(report,"Input parameters",0)
report <- add.param(report,"species",species)
report <- add.param(report,"Expression file",basename(exp_file))
report <- add.param(report,"Design file",basename(design_file))
if(design_type!="continuous"){
  report <- add.param(report,"Comparison",paste0(cond1," vs ", cond2))
}else{
  report <- add.param(report,key = "Correlation ",value = "with a continuous variable")
}

report <- add.section(report,"Module Values",0)
report <- add.download(report,"Module Values","module_vals.txt")

if(design_type!="continuous"){
  report <- add.image(report,"Heatmap","module_heatmap.png",width = 800,level = 1)
  report <- add.image(report,"PCA","module_pca.png",width = 800,level = 1)
}

report <- add.table(report,"Path significance","module_stats.txt",1)
report <- add.section(report,"Pathways",1)
report <- add.pathwayViewer(report, "session", "pathways",1)

write(render.xml(report),file=paste0(output_folder,"/report.xml"))

