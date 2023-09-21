#
species <- "hsa"
library(KEGGREST)
# Convert metabolite names to KEGG IDs
kegg_ids_pubchem <- KEGGREST::keggConv(target = "compound", source = "pubchem")
kegg_ids_chebi <- KEGGREST::keggConv(target = "compound", source = "chebi")
library(hipathia)
p <- hipathia::load_pathways(species = species)
pgs <- p$pathigraphs
#unique(unlist
kegg_ids_pathNodes_df <- lapply(pgs, function(x) {
  annot$id <- x$path.id
  tooltip<-V(x$graph)$tooltip[which(V(x$graph)$shape == 
                                      "circle")]
  annot <- vertex_attr(x$graph) %>% as_data_frame(.) %>%  .[.$shape == "circle",c("label","tooltip")]
  # Extract the kegg id
  annot$nodeID <- str_extract(annot$tooltip, "(?<=href=http://www.kegg.jp/dbget-bin/www_bget\\?)[^>]+")
  annot$path.id <- x$path.id
  annot$path.name<- x$path.name
  y <- annot$tooltip
  names(y) <- annot$label
  return(annot)
  # return(y)
})
kegg_ids_pathNodes_df2<-do.call(rbind, kegg_ids_pathNodes_df)
write.table(kegg_ids_pathNodes_df2, file = "metabolite_146KEGGpathways.tsv", quote = F, append = F, sep = "\t", row.names = F, col.names = T)

filtered_df <- kegg_ids_pathNodes_df2 %>%
  filter(grepl("\\+", nodeID))
write.table(filtered_df, file = "metabolite_146KEGGpathways_withComplexity.tsv", quote = F, append = F, sep = "\t", row.names = F, col.names = T)

#only relevant pathways c("hsa04064", "hsa04064","hsa05205","hsa05231")
rlv <-c("hsa04064", "hsa04520","hsa05205","hsa05231")
kegg_ids_pathNodes_df2_rlv <- kegg_ids_pathNodes_df2[which(kegg_ids_pathNodes_df2$path.id %in% rlv),]
View(kegg_ids_pathNodes_df2_rlv)
View(as_data_frame(vertex_attr(p$pathigraphs$hsa05205$graph)))
# les relev
rlv2 <-c("hsa04390","hsa04724","hsa04919","hsa04971","hsa05132","hsa05164","hsa05168")
kegg_ids_pathNodes_df2_rlv2 <- kegg_ids_pathNodes_df2[which(kegg_ids_pathNodes_df2$path.id %in% rlv2),]
View(kegg_ids_pathNodes_df2_rlv2)
View(as_data_frame(vertex_attr(p$pathigraphs$hsa04390$graph)))

kegg_ids_pathNodes2<-unname(kegg_ids_pathNodes) %>% unlist(.)

hipathia::DApathway("hsa04724", p)
