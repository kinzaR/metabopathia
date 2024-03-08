## Selecting and Identifying MGI pathways without any complication. 
# Simple pathways are :
#   - Have only entrezID or NA in the genesList
# - If there is metabolites: only one per metaboliteID; no + & no hsa: 

# library(hipathia)
# species <- "hsa"
# p <- hipathia::load_pathways(species = species)
# pgs <- p$pathigraphs
library("stringr")
library("igraph")
get_easy_pathways<- function(p){
  pgs <- p$pathigraphs
  easy_paths<-c()
  for(i in pgs){
    graph_att <-do.call(what = cbind, args = vertex_attr(i$graph)) %>% as.data.frame(.)
    # condition 1: check the genesList is it is only entrezId or / or NA
    #this just to test manualy
    # for (j in seq_along(graph_att$genesList)) {
    #   l<-graph_att$genesList[[j]]
    #   if(!(is.null(l) || all(is.na(l)) || all(!is.na(suppressWarnings(as.numeric(l))) | l=="/")))
    #   cat(j," : ",l," ==> ",is.null(l) || all(is.na(l)) || all(!is.na(suppressWarnings(as.numeric(l))) | l=="/"),"\n")
    # }
    cond1 <-sapply(graph_att$genesList, function(l){
      is.null(l) || all(is.na(l)) || all(!is.na(suppressWarnings(as.numeric(l))) | l=="/")
    }) %>% all()
    # if(!cond1){
    #   cat("Cond1 not satisfied in : ", i$path.id)
    # }
    # condition 2: check metabolite
    detected_metabolites <- graph_att$shape=="circle" & !is.na(str_extract(graph_att$tooltip, "(?<=href=http://www.kegg.jp/dbget-bin/www_bget\\?)[^>]+"))
    detected_metabolites <- graph_att[which(detected_metabolites),"tooltip"] %>% str_extract(., "(?<=href=http://www.kegg.jp/dbget-bin/www_bget\\?)[^>]+")
    # check if there is a CID + CID
    cond2 <- all(!grepl("\\+", detected_metabolites))
    # if(!cond2){
    #   cat("Cond2 not satisfied in : ", i$path.id)
    #   print(detected_metabolites)
    # }
    if(cond1 && cond2) easy_paths <- c(easy_paths,i$path.id)
  }
  cat("We will work with ",length(easy_paths), " out of ", length(p$pathigraphs))
  return(easy_paths)
}
