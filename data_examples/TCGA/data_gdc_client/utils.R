library(edgeR)
library(sva)
do_combat <- function(exp_data){
  #Combat 
  data4combat<-assay(exp_data) #dim fom BLCA: 60660 X 431
  # Colors for sample types:
  data_frame % > % mutate(Price_status=case_when(Price >= 500000 & Price <= 900000 ~ "Average", Price > 900000 ~ "High", TRUE ~ "Low")) 
  colors <- data.frame("id" = exp_data$barcode,
    "type" = exp_data$sample_type) %>% 
    mutate(color= case_when(type == "Solid Tissue Normal" ~ "#009E73", 
                            type == "Primary Solid Tumor" | type == "Primary Tumor" ~ "#E69F00",
                            type == "Recurrent Solid Tumor" | type == "Recurrent Tumor" ~ "#D55E00",
                            type == "Primary Blood Derived Cancer" ~ "#56B4E9",
                            type == "Recurrent Blood Derived Cancer" ~ "#0072B2",
                            type == "Metastatic" ~ "#CC79A7"))
  plotMDS(data4combat, col=colors$color)
  var.data <- apply(data4combat, 1, var)
  data <- data4combat[which(var.data != 0 ),]  # Maria and Mrin said: eliminamos entras con 0 varianza, porque si no peta ComBat
  plotMDS(data, col=colors$color)
  group <- factor(exp_data$sample_type)
  my_DGEList <- DGEList(counts=data4combat, group=group)
  my_mod = model.matrix(~group, data=colData(exp_data))
  All_cancer_norm <- calcNormFactors(my_DGEList) #, method="upperquartile")
  my_data = cpm(All_cancer_norm, log=TRUE, prior.count=2) # prior.count =1 to have log 1 = 0 then the meaning will not change maybe ?
  # unique(exp_data$preservation_method)
  # unique(exp_data$icd_10_code)
  blocks<- do.call(rbind,lapply(exp_data$barcode, function(bc){
    center <- strsplit(bc,"-")[[1]][7]
    plate <-  strsplit(bc,"-")[[1]][6]
    tss <-  strsplit(bc,"-")[[1]][2]
    return(c(center, plate, tss))
    }
  ))%>% as.data.frame(stringsAsFactors = T)
  colnames(blocks) <- c("center", "plate", "tss") 
  rownames(blocks)<- exp_data$barcode
  blocks$group <-factor (exp_data$sample_type)
  # other colors 
  rownames(colors)<- colnames(exp_data)
  # by center 
  cols_center<-setNames(rainbow(length(levels(blocks$center))),levels(blocks$center))
  colors <- colors %>%  mutate(by_center = cols_center[strsplit(id,"-")[[1]][7]])
  plotMDS(my_data, col= colors$by_center)
  # by plate
  cols_plate<-setNames(rainbow(length(levels(blocks$plate))),levels(blocks$plate))
  colors <- colors %>% rowwise()%>% mutate(by_plate = cols_plate[strsplit(id,"-")[[1]][6]])
  plotMDS(my_data, col= colors$by_plate)
  # by tss
  cols_tss<-setNames(rainbow(length(levels(blocks$tss))),levels(blocks$tss))
  colors <- colors %>%  mutate(by_tss = cols_tss[strsplit(id,"-")[[1]][2]])
  plotMDS(my_data, col= colors$by_tss)
  
  # combat <- ComBat(dat=my_data, batch=blocks$plate, mod=blocks$group, par.prior=TRUE, prior.plots=TRUE)
  combat <- ComBat(dat=my_data, batch=blocks$plate, mod=blocks$group)
  # recheck here !
  plotMDS(combat, col=colors$by_plate)
  plotMDS(combat, col=colors$color)
  
  # combat <- ComBat(dat=my_data, batch=blocks$plate, mod=blocks$group, par.prior=TRUE, prior.plots=TRUE)
  combat <- ComBat(dat=my_data, batch=blocks$tss, mod=blocks$group)
  # recheck here !
  plotMDS(combat, col=colors$by_tss)
  plotMDS(combat, col=colors$color)
  
  # combat <- ComBat(dat=my_data, batch=blocks$plate, mod=blocks$group, par.prior=TRUE, prior.plots=TRUE)
  combat <- ComBat(dat=my_data, batch=blocks$center, mod=blocks$group) # in the curent example BLCA cancer :  only one center no sene here !
  # recheck here !
  plotMDS(combat, col=colors$by_center)
  plotMDS(combat, col=colors$color)
  
  r<-as.data.frame(combat)
  my_names<-getBM(filters = "entrezgene", attributes= c("mgi_symbol","entrezgene"),values=rownames(r),mart= mart)
  r$id <- rownames(r)
  r<-merge(my_names,r, by.y="id",by.x="entrezgene")
  write.table(r, "data_combat", sep = '\t', quote = FALSE, row.names = F) # write RData!
}
## Function to load Rdata from main folder:
load_data <- function(cancer_id, path2dir) {
  file_path <- file.path(path2dir, cancer_id, paste0("counts_", cancer_id, ".RData"))
  if (file.exists(file_path)) {
   load(file_path)
    return(data)
  } else {
    stop(paste("File not found for cancer:", cancer))
  }
}

