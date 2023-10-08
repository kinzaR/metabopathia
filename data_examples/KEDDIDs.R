library(httr)
library(xml2)
ampollas_raw <- read.csv("data_examples/Dystrophic_epidermolysis_bullosa/Carlos_leon_results_Ampollas_LC_V3.csv")

getKEGGID <- function(exact_mass,formula,mol_weight_min,mol_weight_max) {
  # Define the KEGG API base URL
  kegg_api_url <- "http://rest.kegg.jp/find/compound"
  
  # Make an HTTP GET request to the KEGG API
  # response <- GET(paste0(kegg_api_url, "/", trunc(exact_mass * 100) / 100,"/exact_mass/",
  response <- GET(paste0(kegg_api_url, "/",
                         gsub("\\[|\\]|\\+", "", formula) ,"/formula/",
                         mol_weight_min,"-",mol_weight_max,"/mol_weight"))
  
  # Check if the response is successful
  if (http_type(response) == "text/html") {
    return("NA")
  }
  keggId <- content(response, as = "text") %>% gsub("cpd:(\\w+).*", "\\1", .)
  
  if(nchar(content(response, as = "text"))>5) {
    print(response)
    cat(exact_mass," - ",content(response, as = "text"),"\n")
  }
  return(keggId)
}

# n<- 230.04554
# exact_mass<-trunc(n * 100) / 100
# formula<-"[C10H10ClN3Na]+"
# mol_weight_min<-325
# mol_weight_max<-337.5

keggIDs <- c()
for (i in seq(length(ampollas_raw$formula))) {
  keggIDs[i] <- getKEGGID(ampollas_raw$mass[i],ampollas_raw$formula[i],ampollas_raw$start[i],ampollas_raw$end[i])
}

# Add a new column for KEGG IDs
metabolites_df$kegg_id <- sapply(c(ampollas_raw$mass,ampollas_raw$formula,ampollas_raw$start,ampollas_raw$end), getKEGGID)
metabolites_df$kegg_id <- sapply(metabolites_df$formula, getKEGGID2)

# Print the updated data frame with KEGG IDs
print(metabolites_df)


kegg_id<-"C00681" 
kegg_names<-function(kegg_id){
  url <- paste("http://rest.kegg.jp/get/", kegg_id, sep="")
  response <- GET(url)
  # Check if the request was successful
  if (http_status(response)$category == "Success") {
    # Parse the response to extract metabolite name
    kegg_info <- content(response, "text")
    lines <- strsplit(kegg_info, "\n")[[1]]
    
    # Extract the metabolite name
    metabolite_name <- NULL
    if(any(startsWith(kegg_id,c("C","D")))){
      for (line in lines) {
        if (startsWith(line, "NAME")) {
          metabolite_name <- gsub("NAME +", "", line)
          break
        }
      }
    } else if(startsWith(kegg_id,"G")){
      for (line in lines) {
        if (startsWith(line, "COMPOSITION")) {
          metabolite_name <- gsub("COMPOSITION +", "", line)
          break
        }
      }
    }
    for (line in lines) {
      if (startsWith(line, "FORMULA")) {
        metabolite_formula <- gsub("FORMULA +", "", line)
        break
      }
    }
    
    # Print the metabolite name
    if (!is.null(metabolite_name)) {
      metabolite_name<-sub(";$", "", metabolite_name)
      # print(metabolite_name)
    } else {
      cat("Metabolite name not found for", kegg_id, "\n")
      metabolite_name <- NA
    }
  } else {
    cat("Failed to fetch KEGG information for", kegg_id, "\n")
    metabolite_name <- "NI"
  }
  return(c(kegg_id,metabolite_name,metabolite_formula))
}
metabo<-sapply(metabo_pathways$all.metabolite, kegg_names)
metab_df <- as.data.frame(metabo) %>% t
colnames(metab_df)<- c("KEGGID","NAME","FORMULA")
metab_df<-as.data.frame(metab_df)
write.table(x = metab_df, file = "metabolite_120_mgi.tsv", quote = F,sep = "\t",row.names = F)
class(metab_df)

# compare between ampollas raw formula
library(dplyr)
get_formula_from_data<-function(ampollas_raw){
  ampollas_clean_formul<- ampollas_raw
  ampollas_clean_formul[!grepl(";",ampollas_clean_formul)] <- sub("^(.*?)_.*$","\\1",ampollas_raw[!grepl(";",ampollas_raw)])
  ampollas_clean_formul[grepl("; ",ampollas_clean_formul)] <- sapply(strsplit(ampollas_clean_formul[grepl("; ",ampollas_clean_formul)]," "),
                                                                              function(xx) {
                                                                                x <- sub(" $", "", xx)
                                                                                formula <-grep("_M\\+H(;)?$", x, value = TRUE)
                                                                                if(isEmpty(formula)) {
                                                                                  formula <-grep("_M\\+NH4(;)?$", x, value = TRUE)
                                                                                  if(isEmpty(formula)){
                                                                                    formula <-grep("_M\\+Na(;)?$", x, value = TRUE)
                                                                                  }
                                                                                }
                                                                                if(!isEmpty(formula)) {
                                                                                  formula <- sub("^(.*?)_.*$","\\1",formula)
                                                                                  if(length(formula)==1) return(formula)
                                                                                  if(length(formula)!=1 && n_distinct(formula)){
                                                                                    return(formula[1])
                                                                                  }else{
                                                                                    cat("\n",formula," has differents formulas options, need a manual feedback!! \n")
                                                                                  }
                                                                                }else{
                                                                                  cat("\n",x," has neither _M+H, _M+Na nor M-NH4+\nThe original formula was: \n")
                                                                                  print(ampollas_raw[ampollas_raw==paste(xx, collapse = " "),c(1:8)])
                                                                                  return(NA)
                                                                                }
                                                                              }) %>%unlist()
  return(ampollas_clean_formul)
}


ampolla_formulas<-get_formula_from_data(ampollas_raw$anot)
#suero
suero_data <- read.csv("data_examples/Dystrophic_epidermolysis_bullosa/Carlos_leon_hilic_positivev2 (1) (2).csv")
suero_data_neg <- read.csv("data_examples/Dystrophic_epidermolysis_bullosa/Carlos_leon_hilic_negativev2 (2).csv")
suero_formulas<-get_formula_from_data(suero_data$anot)
# suero_formulas<-get_formula_from_data(suero_data_neg$anot)

ampollas_raw$anot<-ampolla_formulas
suero_data$anot<-suero_formulas

intersect(metab_df$FORMULA, ampolla_formulas) %>%length()
intersect(metab_df$FORMULA, suero_formulas) %>%length()

data_for_metabo_amp<-ampollas_raw[ampollas_raw$anot %in%metab_df$FORMULA, ]
data_for_metabosuero<-suero_data[suero_data$anot %in%metab_df$FORMULA, ]

length(unique(data_for_metabo_amp$anot))
length(unique(data_for_metabosuero$anot))

amp_met<-metab_df[metab_df$FORMULA%in%data_for_metabo_amp$anot,1] %>%unique()
sue_met<-metab_df[metab_df$FORMULA%in%data_for_metabosuero$anot,1] %>%unique()
# p <- metabo_pathways$pathigraphs$hsa04923
library(purrr)
met_per_path_amp<-lapply(metabo_pathways$pathigraphs, function(p){
  commun <- intersect(V(p$graph)$metaboID , amp_met)
  if(length(commun)>0)  return(c(all(V(p$graph)$metaboID[!is.na(V(p$graph)$metaboID)] %in% amp_met),
                                 100*length(commun)/length(V(p$graph)$metaboID[!is.na(V(p$graph)$metaboID)]),
                                 paste(commun,collapse = ","),
                                 paste(unique(V(p$graph)$metaboID[!is.na(V(p$graph)$metaboID)]),collapse = ",")))
}) %>% purrr::compact(.)
met_per_path_suer<-lapply(metabo_pathways$pathigraphs, function(p){
  commun <- intersect(V(p$graph)$metaboID , sue_met)
  if(length(commun)>0)  return(c(all(V(p$graph)$metaboID[!is.na(V(p$graph)$metaboID)]%in% sue_met),
                                 100*length(commun)/length(V(p$graph)$metaboID[!is.na(V(p$graph)$metaboID)]),
                                 paste(commun,collapse = ","),
                                 paste(unique(V(p$graph)$metaboID[!is.na(V(p$graph)$metaboID)]),collapse = ",")))
}) %>% purrr::compact(.)

inter_ampolla_data <- do.call(rbind,met_per_path_amp)%>%as.data.frame()
inter_suero_data <- do.call(rbind,met_per_path_suer)%>%as.data.frame()
colnames(inter_ampolla_data)<- c("fully covered", "percentage of coverage %", "detected metabolites","metabolite in pathway")
colnames(inter_suero_data)<- c("fully covered", "percentage of coverage %", "detected metabolites","metabolite in pathway")

write.table(x = inter_ampolla_data, file = "data_examples/Dystrophic_epidermolysis_bullosa/intersection_Kegg_ampolla.tsv", quote = F, sep = "\t", row.names = T, col.names = T)
write.table(x = inter_suero_data, file = "data_examples/Dystrophic_epidermolysis_bullosa/intersection_Kegg_suero.tsv", quote = F, sep = "\t", row.names = T, col.names = T)

inter_suero_data[as.logical(inter_suero_data$`fully covered`),]
inter_ampolla_data[as.logical(inter_ampolla_data$`fully covered`),]
