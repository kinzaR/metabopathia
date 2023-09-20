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
kegg_ids_pathNodes <- sapply(pgs, function(x) {
  tooltip<-V(x$graph)$tooltip[which(V(x$graph)$shape == 
                                      "circle")]
  annot <- vertex_attr(x$graph) %>% as_data_frame(.) %>%  .[.$shape == "circle",c("label","tooltip")]
  # Extract the kegg id
  annot$tooltip <- str_extract(annot$tooltip, "(?<=href=http://www.kegg.jp/dbget-bin/www_bget\\?)[^>]+")
  y <- annot$tooltip
  names(y) <- annot$label
  return(y)
})
kegg_ids_pathNodes2<-unname(kegg_ids_pathNodes) %>% unlist(.)
