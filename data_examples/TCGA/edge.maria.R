if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
source("https://bioconductor.org/biocLite.R")

BiocManager::install("edgeR")
BiocManager::install("sva")
BiocManager::install("hipathia", version = "3.8")
BiocManager::install("ComplexHeatmap", version = "3.8")

library(edgeR)
library(hipathia)
library(preprocessCore)
library(biomaRt)
library(sva)

DEBUG <- TRUE
NORMALIZATION <- "limma"    #puede ser "combat" o "norm"


###############  C O N S T A N T S  ###############
WD <- paste0("/mnt/sshfs/gattaca1/mnt/lustre/scratch/CBRA/collaborations/AMartinMontalvo/results/old/hipathia2/analysis_proteinsInAllExperiments_", NORMALIZATION, "/")
DATA_DIR <- "/mnt/sshfs/gattaca1/mnt/lustre/scratch/CBRA/collaborations/AMartinMontalvo/results/old/hipathia2/"
DATA_FILE <- paste0(DATA_DIR, "merge.csv")
DESIGN_FILE <- paste0(DATA_DIR, "design_limma.txt")
HIPATHIA_DESIGN <- paste0(DATA_DIR, "design.txt")
CASES <- c( 
  c("hfd", "std"), 
  c("hfdsb","stdsb"), 
  c("stdsb", "std"),
  c("hfdsb", "hfd"))


###############  I N I T  ###############
#establece directorio de trabajo
dir.create(WD, showWarnings = FALSE, recursive = TRUE)
setwd(WD)

#hacer listas para guardar resultados
de<-list()#lista para resultados de expresion diferencial
ds<-list()#lista para resultados de señalizacion diferencial

#read expression files
exp<-read.table(DATA_FILE, header=T, sep="\t", stringsAsFactors=F, row.names=1)

#read and parse sample info files
clin<-read.table(DESIGN_FILE, header=T, sep="\t", stringsAsFactors = T)
#clin$groups<-as.factor(paste(clin$diet,clin$drug,sep="."))

data  <- as.matrix(exp)
block <- factor(clin$experiment)
group <- factor(clin$groups,levels = levels(clin$groups))


#para la conversión de identificadores
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))



###############  N O R M A L I Z E  ###############
##remove low-count for differential expression and normalize data
y<-DGEList(counts=data,group=group)
y<-calcNormFactors(y)
keep<-rowSums(cpm(y))>3
y <- y[keep, , keep.lib.sizes=T]
data_norm<-cpm(y,log=T,prior.count = 3) #Para sacar los cpms si quieres usarlos a parte
r<-as.data.frame(data_norm)
my_names<-getBM(filters = "entrezgene", attributes= c("mgi_symbol","entrezgene"),values=rownames(r),mart= mart)
r$id <- rownames(r)
r<-merge(my_names,r, by.y="id",by.x="entrezgene")
write.table(r, "data_norm", sep = '\t', quote = FALSE, row.names = F)


if (NORMALIZATION == "combat"){
  #normalize2
  data_combat<-data
  plotMDS(data_combat, col=ifelse(clin$experiment=="iTRAQ1", "blue","red"))
  var.data <- apply(data_combat, 1, var)
  data <- data_combat[which(var.data != 0 ),]  # eliminamos entras con 0 varianza, porque si no peta ComBat
  plotMDS(data, col=ifelse(clin$experiment=="iTRAQ1", "blue","red"))
  my_DGEList <- DGEList(counts=data_combat, group=group)
  #my_mod = model.matrix(~1, data=clin)
  my_mod = model.matrix(~group, data=clin)
  All_cancer_norm <- calcNormFactors(my_DGEList) #, method="upperquartile")
  my_data = cpm(All_cancer_norm, log=TRUE, prior.count=2)
  combat <- ComBat(dat=my_data, batch=block, mod=my_mod, par.prior=TRUE, prior.plots=TRUE)
  plotMDS(combat, col=ifelse(clin$experiment=="iTRAQ1", "blue","red"))
  r<-as.data.frame(combat)
  my_names<-getBM(filters = "entrezgene", attributes= c("mgi_symbol","entrezgene"),values=rownames(r),mart= mart)
  r$id <- rownames(r)
  r<-merge(my_names,r, by.y="id",by.x="entrezgene")
  write.table(r, "data_combat", sep = '\t', quote = FALSE, row.names = F)
}



######## nueva normalización
#normalize3
data2<-data
plotMDS(data2, col=ifelse(clin$experiment=="iTRAQ1", "blue","red"))
data_norm <- removeBatchEffect(cpm(data2, log=T, prior.count=3),batch=factor(clin$experiment))
plotMDS(data_norm, col=ifelse(clin$experiment=="iTRAQ1", "blue","red"))
r<-as.data.frame(data_norm)
my_names<-getBM(filters = "entrezgene", attributes= c("mgi_symbol","entrezgene"),values=rownames(r),mart= mart)
r$id <- rownames(r)
r<-merge(my_names,r, by.y="id",by.x="entrezgene")
write.table(r, "data_norm2", sep = '\t', quote = FALSE, row.names = F)




#design<-model.matrix(~group)
design<-model.matrix(~0+group+block)
#vwts <- voomWithQualityWeights(data_norm, design=design, plot=TRUE)  # RECOMENDABLE!!!!!!!!!!!!!!!!
colnames(design)<-c(levels(clin$groups), "iTRAQ2")
cont.matrix <- makeContrasts(
  SBinSTD=stdsb-std,              # genes DE en STD cuando se aplica SB
  SBinHFD=hfdsb-hfd,              # genes DE en HFD cuando se aplica SB
  HFDvsSTD=hfd-std,
  HFDSBvsSTDSB=hfdsb-stdsb,
  diff=(hfdsb-hfd)-(stdsb-std),   # genes que responden diferencialmente frente a la droga cuando se compara HFD y STD (diferencia de diferencias)
  levels=design)





y<-data_norm
y<-estimateDisp(y,design)
fit<-glmQLFit(y,design)
qlf <- glmQLFTest(fit, contrast = cont.matrix[,"SBinSTD"])
r<-as.data.frame(topTags(qlf,adjust.method = "fdr",sort.by = "p.value" ,p.value = 0.05,n=Inf))
my_names<-getBM(filters = "entrezgene", attributes= c("mgi_symbol","entrezgene"),values=rownames(r),mart= mart)
r$id <- rownames(r)
r<-merge(my_names,r, by.y="id",by.x="entrezgene")
write.csv(r, "DE.SBinSTD", row.names = F)


qlf <- glmQLFTest(fit, contrast = cont.matrix[,"SBinHFD"])
r<-as.data.frame(topTags(qlf,adjust.method = "fdr",sort.by = "p.value" ,p.value = 0.05,n=Inf))
my_names<-getBM(filters = "entrezgene", attributes= c("mgi_symbol","entrezgene"),values=rownames(r),mart= mart)
r$id <- rownames(r)
r<-merge(my_names,r, by.y="id",by.x="entrezgene")
write.csv(r, "DE.SBinHFD", row.names = F)

qlf <- glmQLFTest(fit, contrast = cont.matrix[,"diff"])
r<-as.data.frame(topTags(qlf,adjust.method = "fdr",sort.by = "p.value" ,p.value = 0.05,n=Inf))
my_names<-getBM(filters = "entrezgene", attributes= c("mgi_symbol","entrezgene"),values=rownames(r),mart= mart)
r$id <- rownames(r)
r<-merge(my_names,r, by.y="id",by.x="entrezgene")
write.csv(r, "DE.diff", row.names = F)


y<-estimateDisp(y,design)
fit<-glmQLFit(y,design)
qlf <- glmQLFTest(fit, contrast = cont.matrix[,"HFDvsSTD"])
r<-as.data.frame(topTags(qlf,adjust.method = "fdr",sort.by = "p.value" ,p.value = 0.05,n=Inf))
my_names<-getBM(filters = "entrezgene", attributes= c("mgi_symbol","entrezgene"),values=rownames(r),mart= mart)
r$id <- rownames(r)
r<-merge(my_names,r, by.y="id",by.x="entrezgene")
write.csv(r, "DE.HFDvsSTD", row.names = F)


y<-estimateDisp(y,design)
fit<-glmQLFit(y,design)
qlf <- glmQLFTest(fit, contrast = cont.matrix[,"HFDSBvsSTDSB"])
r<-as.data.frame(topTags(qlf,adjust.method = "fdr",sort.by = "p.value" ,p.value = 0.05,n=Inf))
my_names<-getBM(filters = "entrezgene", attributes= c("mgi_symbol","entrezgene"),values=rownames(r),mart= mart)
r$id <- rownames(r)
r<-merge(my_names,r, by.y="id",by.x="entrezgene")
write.csv(r, "DE.HFDSBvsSTDSB", row.names = F)




#p_values<-topTags(qlf,adjust.method = "fdr",sort.by = "p.value" ,p.value = 1,n=Inf)
#sig<-as.data.frame(sig)
#p_values<-as.data.frame(p_values)









#########################
###  H I P H A T I A  ###
#########################
dd<-read.table(HIPATHIA_DESIGN, header=T, check.names=T,  stringsAsFactors=F, row.names = 1) 

pathways <- load_pathways(species = "mmu")
#ds.exp <- translate_matrix(y, "mmu")



if (NORMALIZATION == "combat"){
  hipdata <- normalize_data(combat, by_quantiles = TRUE)

} else {
  hipdata <- normalize_data(data_norm, by_quantiles = TRUE)
}
r<-as.data.frame(hipdata)
my_names<-getBM(filters = "entrezgene", attributes= c("mgi_symbol","entrezgene"),values=rownames(r),mart= mart)
r$id <- rownames(r)
r<-merge(my_names,r, by.y="id",by.x="entrezgene")
write.table(r, paste0("hipdata_", NORMALIZATION), sep = '\t', quote = FALSE, row.names = F)

results <- hipathia(hipdata, pathways, verbose=F)
comp <- get_paths_matrix(results)
path.vals <- comp

dd <- dd[match(colnames(comp),rownames(dd)),,drop=F]
all(colnames(comp)==rownames(dd))
sample.group <- dd[colnames(comp),"group"]


#fit<-lmFit(comp,design)
#fit <- eBayes(fit)
#plotSA(fit)
#ds_res<-topTable(fit, number = Inf, coef = c("hfd","std"))
#ds_res$path_name <- get_path_names(pathways, rownames(ds_res))

#heatmap_plot(path.vals, 
#             sample_type=sample.group, 
#             colors="hipathia", 
#             variable_clust = TRUE,
#             sample_clust = TRUE,
#             main="Genes")

# Uniprot Keywords
#uniprot.vals <- quantify_terms(results, pathways, "uniprot")
# plot
#heatmap_plot(uniprot.vals, 
#             sample_type=sample.group, 
#             variable_clust = TRUE, 
#             main="Uniprot")

# GO terms
go.vals <- quantify_terms(results, pathways, "GO")
#plot heatmap
#heatmap_plot(go.vals, 
#             sample_type=sample.group, 
#             colors = "hipathia", 
#             variable_clust = TRUE,
#             main="GO terms")


for (caseIndex in seq(1,length(CASES),2)) {
  g1 <- CASES[caseIndex]
  g2 <- CASES[caseIndex+1]
  study <- paste0(g1, "_", g2)
  pdf(file=paste0("plots_",study,".pdf"))
  
  # select most relevant
  comp <- do_wilcoxon(path.vals, sample.group, g1 = g1, g2 = g2, adjust = TRUE)
  path.names <- get_path_names(pathways, rownames(comp))
  comp <- cbind(path.names, comp)
  conf <- max(comp[comp$p.value<=0.05,]$FDRp.value) + 0.0001  
  if (conf>1){conf<-1}
  
  pathways.summary <- get_pathways_summary(comp, pathways, conf=conf)
  pathways.selected.list <- pathways.summary[pathways.summary$num_significant_paths>=1,]$id_pathways
  pathways.selected <-  load_pathways(species="mmu", pathways_list = pathways.selected.list)
  write.csv(pathways.summary[pathways.summary$num_significant_paths>=1,], paste0("pathways_summary_",study, ".csv"))
  head(pathways.summary, n = 15)
  
  #ranked.path.vals <- path.vals[order(comp$p.value, decreasing = FALSE),]
  signif.path.vals <- path.vals[comp$p.value<0.05,]
  write.csv(signif.path.vals, paste0("modules_signif_vals_", study, ".csv"))
  signif.comp <- comp[comp$p.value<0.05, ]
  write.csv(signif.comp, paste0("modules_signif_pval_", study, ".csv"))
  table(comp$p.value<0.05)
  
  
  
  
  # Plot heatmap for most relevant paths
  #heatmap_plot(path.vals[comp$p.value<0.05,],
  #             sample.group, 
  #             colors="hipathia", 
  #             variable_clust = TRUE,
  #             sample_clust = TRUE,
  #             main=paste("D.E. pathways", study, "(p<=0.05)"),
  #             legend = FALSE,
  #             labRow = rownames(path.vals) )
  heatmap_plot(signif.path.vals,
               sample.group, 
               colors="hipathia", 
               variable_clust = TRUE,
               sample_clust = TRUE,
               main=paste("Cellular signaling", study, "(p<=0.05)"),
               legend = FALSE)#,
               #labRow = rownames(path.vals))
  
  # Plot heatmap for most relevant paths
  #heatmap_plot(ranked.path.vals[1:15,], sample_type=sample.group, variable_clust = TRUE)
  #heatmap_plot(ranked.path.vals[1:99,], 
  #             sample.group, 
  #             colors="hipathia", 
  #             variable_clust = TRUE,
  #             sample_clust = TRUE,
  #             main=paste0("Top pathways ", study),
  #             legend = FALSE)
  
  
  
  
  # Visualization
  #------------------------
  # Visualize comparison
  #colors.de <-node_color_per_de(results, pathways.selected, sample.group, g1, g2, conf=conf)
  
  #pathway_comparison_plot(comp, metaginfo = pathways, pathway = "mmu05014", colors = "hipathia")
  #colors.de.hipathia <- node_color_per_de(results, pathways, sample.group, g1, g2, colors = "hipathia", conf = conf)
  #pathway_comparison_plot(comp, metaginfo = pathways, pathway = "mmu03320", node_colors = colors.de.hipathia, colors = "hipathia")
  
  # Visualize comparison in server
  #create_report(comp, pathways, "save_noColors")
  create_report(comp, pathways.selected, paste0("report_",study), conf=conf)
  #visualize_report("save_colors")
  

  # PCA
  pca.model <- do_pca(signif.path.vals[1:min(ncol(signif.path.vals), nrow(signif.path.vals)),])
  # plot PCA
  #pca_plot(pca.model, sample.group)
  multiple_pca_plot(pca.model, 
                    sample.group, 
                    cex=3, 
                    plot_variance = TRUE, 
                    main=paste("Cellular signaling",study))
  
  
  
  # Decomposed pathways
  #pathways.decomposed <- load_pathways(species = "mmu", pathways_list = pathways.selected.list)
  #results.decomposed <- hipathia(exp.data, pathways.decomposed, verbose=TRUE, decompose = TRUE)
  #dec.vals <- get_paths_matrix(results.decomposed)
  
  #comp.dec <- do_wilcoxon(dec.vals, sample.group, g1=g1, g2=g2)
  #dec.names <- get_path_names(pathways.decomposed, rownames(comp.dec))
  #comp.dec <- cbind(dec.names, comp.dec)
  
  #create_report(comp.dec, pathways.decomposed, paste0("report_decomposed_",study), node_colors=colors.de)
  #visualize_report("save_decomposed")
  
  
  # Function analysis
  #---------------------------------
  #comp.uni <- do_wilcoxon(uniprot.vals, sample.group, g1=g1, g2=g2)
  #ranked.uni.comp <- comp.uni[order(comp.uni$p.value, decreasing = FALSE),]
  #ranked.uni.vals <- uniprot.vals[order(comp.uni$p.value, decreasing = FALSE),]
  #write.csv(ranked.uni.comp, paste0("ranked_uniprot_", study, ".csv"))
  #write.csv(ranked.uni.vals, paste0("ranked_uniprot_vals_", study, ".csv"))
  
  #heatmap_plot(ranked.go.vals, sample_type=sample.group, variable_clust = TRUE)
#  heatmap_plot(ranked.uni.vals[ranked.uni.comp$p.value<=0.05,], 
#               sample_type=sample.group, 
#               variable_clust = TRUE,
#               main=paste("Uniprot functions", study, "(p<=0.05)"),
#               legend = FALSE,
#               labRow = rownames(ranked.uni.vals) )
  
  
  
  
  # Select most relevant GO terms
  comp.go <- do_wilcoxon(go.vals, sample.group, g1=g1, g2=g2)
  signif.go.comp <- comp.go[comp.go$p.value<0.05,]
  signif.go.vals <- go.vals[comp.go$p.value<0.05,]
  write.csv(signif.go.comp, paste0("signif_GO_pvals_", study, ".csv"))
  write.csv(signif.go.vals, paste0("signif_GO_vals_", study, ".csv"))
  
  table(comp.go$p.value<0.05)
  
  
  #heatmap_plot(ranked.go.vals, sample_type=sample.group, variable_clust = TRUE)
  heatmap_plot(signif.go.vals, 
               sample_type=sample.group, 
               variable_clust = TRUE,
               main=paste("Cellular functions",study,"(p<=0.05)"),
               legend = FALSE)#, 
               #labRow = rownames(signif.go.vals))
  
  
  
  # PCA
  pca.model <- do_pca(signif.go.vals[1:min(ncol(signif.go.vals), nrow(signif.go.vals)),])
  # plot PCA
  #pca_plot(pca.model, sample.group)
  multiple_pca_plot(pca.model, 
                    sample.group, 
                    cex=3, 
                    plot_variance = TRUE, 
                    main=paste("Cellular functions", study))
  
  
  #heatmap_plot(ranked.go.vals[comp.go$p.value<0.05  & comp.go$FDRp.value<1,], sample.group, colors="hipathia", variable_clust = TRUE,main="Significant GO terms")
  
  dev.off()  
}




###################
###  P L O T S  ###
###################

#PANEL A - Valores de expresión absolutos (normalizados) para proteínas DE en al menos una comparación (STDSB/STD, HFDSB/HFD, HFD/STD, HFDSB/STDSB)
WD <- paste0("/mnt/sshfs/gattaca1/mnt/lustre/scratch/CBRA/collaborations/AMartinMontalvo/results/analysis_proteinsInAllExperiments_combat/")
setwd(WD)

DE.HFDvsSTD<-read.delim("diffExpr/DE.HFDvsSTD", sep=",")
DE.HFDSBvsSTDSb<-read.delim("diffExpr/DE.HFDSBvsSTDSB", sep=",")
DE.HFDSBvsHFD<-read.delim("diffExpr/DE.SBinHFD", sep=",")
DE.STDSBvsSTD<-read.delim("diffExpr/DE.SBinSTD", sep=",")

# data_combat <- read.delim(paste0(WD,"data/data_combat"))
# #data_combat.dedup <- unique(data_combat[,c(1,3:ncol(data_combat))])
# data_combat.dedup <- unique(data_combat[,2:ncol(data_combat)])
# data_combat.dedup <- data_combat.dedup[!data_combat.dedup$mgi_symbol == "",]
# rownames(data_combat.dedup) <- data_combat.dedup$mgi_symbol
# data_combat.dedup <- data_combat.dedup[,2:ncol(data_combat.dedup)]
# data_combat.dedup <- data_combat.dedup[, c(1,4,9,12,
#                                            3,7,11,15,
#                                            2,6,10,14,                      5
#                                            ,8,13,16)]
# 
# heatmap(as.matrix(data_combat.dedup), Colv=NA, scale = "row", col=palette)
# 
# 
# 
# 
# #data        = data.matrix(data_combat.dedup)
# data <- data.matrix(data_combat.dedup[rownames(data_combat.dedup) %in% unique(DE.HFDvsSTD$mgi_symbol),])
# distance    = dist(data, method = "canberra")
# cluster     = hclust(distance, method="ward.D2")
# dendrogram  = as.dendrogram(cluster)
# Rowv        = rowMeans(data, na.rm = T)
# dendrogram  = reorder(dendrogram, Rowv)
# reorderfun = function(d,w) { d }
# palette <- colorRampPalette(c('#ff0000','#ffff00','#0000ff'))(256)
# heatmap(data, scale="row", col=palette, Rowv=dendrogram, reorderfun=reorderfun)
# #
# 
# 
# 
# 
# DE <- merge(DE.HFDvsSTD, DE.HFDSBvsSTDSb, by="mgi_symbol", all=T)
# rownames(DE) <- DE$mgi_symbol
# DE$logCPM.x[isNA(DE$logCPM.x)] <- 0
# DE$logCPM.y[isNA(DE$logCPM.y)] <- 0
# DE <- unique(DE[,c("logCPM.x", "logCPM.y")])
# palette <- colorRampPalette(c('#ff0000','#ffff00','#0000ff'))(256)
# heatmap(as.matrix(DE), Colv=NA, scale = "none")
# 
# 
# 
# data <- data.matrix(combat[rownames(combat) %in% unique(DE.HFDvsSTD$entrezgene),])
# distance    = dist(data, method = "canberra")
# cluster     = hclust(distance, method="ward.D2")
# dendrogram  = as.dendrogram(cluster)
# Rowv        = rowMeans(data, na.rm = T)
# dendrogram  = reorder(dendrogram, Rowv)
# reorderfun = function(d,w) { d }
# palette <- colorRampPalette(c('#ff0000','#ffff00','#0000ff'))(256)
# order <- c(1,4,9,12,
#            3,7,11,15,
#            2,6,10,14,
#            5,8,13,16)
# heatmap(data[, order], col=palette, Colv=NA, scale = "row", Rowv=dendrogram, reorderfun=reorderfun)



data <- data.matrix(data_norm[rownames(data_norm) %in% unique(DE.HFDvsSTD$entrezgene),])
distance    = dist(data, method = "canberra")
cluster     = hclust(distance, method="ward.D2")
dendrogram  = as.dendrogram(cluster)
Rowv        = rowMeans(data, na.rm = T)
dendrogram  = reorder(dendrogram, Rowv)
reorderfun = function(d,w) { d }
palette <- colorRampPalette(c('#ff0000','#ffff00','#0000ff'))(256)
order <- c(1,4,9,12,
           3,7,11,15,
           2,6,10,14,
           5,8,13,16)
heatmap(data[, order], col=palette, Colv=NA, scale = "row", Rowv=dendrogram, reorderfun=reorderfun, labRow=NA)




data <- data.matrix(data_norm[rownames(data_norm) %in% unique(DE.HFDSBvsSTDSb$entrezgene),])
distance    = dist(data, method = "canberra")
cluster     = hclust(distance, method="ward.D2")
dendrogram  = as.dendrogram(cluster)
Rowv        = rowMeans(data, na.rm = T)
dendrogram  = reorder(dendrogram, Rowv)
reorderfun = function(d,w) { d }
palette <- colorRampPalette(c('#ff0000','#ffff00','#0000ff'))(256)
order <- c(1,4,9,12,
           3,7,11,15,
           2,6,10,14,
           5,8,13,16)
heatmap(data[, order], col=palette, Colv=NA, scale = "row", Rowv=dendrogram, reorderfun=reorderfun)





data <- data.matrix(data_norm[rownames(data_norm) %in% unique(DE.HFDSBvsHFD$entrezgene),])
distance    = dist(data, method = "canberra")
cluster     = hclust(distance, method="ward.D2")
dendrogram  = as.dendrogram(cluster)
Rowv        = rowMeans(data, na.rm = T)
dendrogram  = reorder(dendrogram, Rowv)
reorderfun = function(d,w) { d }
palette <- colorRampPalette(c('#ff0000','#ffff00','#0000ff'))(256)
order <- c(1,4,9,12,
           3,7,11,15,
           2,6,10,14,
           5,8,13,16)
heatmap(data[, order], col=palette, Colv=NA, scale = "row", Rowv=dendrogram, reorderfun=reorderfun)
