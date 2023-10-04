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
metab_df$FORMULA %in% ampollas_raw$anot
found<-c()
for (met in metab_df$FORMULA) {
  if(length(metab_df[grepl(met, ampollas_raw$anot),1])>0) {
      found<- c(found,met)
      cat("**",met," was found in ",ampollas_raw$anot,"\n")
    }
}
found_matr<-matrix(found, ncol = 2)
