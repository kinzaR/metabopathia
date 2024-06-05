# nodes.vals=nodes.vals; subgraph=graphobject[[1]]; ininodes=graphobject[[2]]; endnode=graphobject[[3]]
# method="pond"; maxnum = 100; tol = 0.000001; divide=F; response.tol = 0
##############
path.value<-function( nodes.vals, subgraph, ininodes, endnode, method="pond", maxnum = 100, tol = 0.000001, divide=F, response.tol = 0 ){
  # Initialize lists
  ready <- ininodes
  # cat(ready)
  processed <- list()
  
  # Initialize node values
  node.signal <- matrix(NA, ncol=ncol(nodes.vals), nrow = length(V(subgraph)), dimnames = list(V(subgraph)$name, colnames(nodes.vals)))
  # head(node.signal)
  endnode.signal.dif <- 10
  # cat(endnode.signal.dif)
  
  num <- 0 
  while( length(ready) > 0 && num <= maxnum){
    num <- num + 1
    # cat(num)
    actnode <- ready[[1]]
    # cat(actnode)
    old.signal <- node.signal[actnode,]
    # cat(old.signal)
    
    # Compute node signal 
    if(divide && actnode != endnode){
      nfol <- length(incident(subgraph, actnode, mode="out"))      
      # cat(nfol)
    }else{
      nfol <- 1
      # cat(nfol)
    }
    
    ###### added by Cankut
    # node.signal[actnode,] <- compute.node.signal2(actnode, nodes.vals[actnode,], node.signal, subgraph, method, response.tol) / nfol
    # head(node.signal)
    aaa<-compute.node.signal2(actnode, nodes.vals[actnode,], node.signal, subgraph, method, response.tol) / nfol
    node.signal[actnode,] <- as.numeric(aaa)
    ######
    
    # Transmit signal
    nextnodes <- get.edgelist(subgraph)[incident(subgraph, actnode, mode="out"),2]
    dif <- old.signal - node.signal[actnode,]
    cat(actnode,"==",node.signal[actnode,])
    
    if(actnode==endnode){
      if(!all(is.na(dif)))
        endnode.signal.dif <- c(endnode.signal.dif, sqrt(sum(dif^2)))
      #num <- num+1
    }
    
    if(all(is.na(old.signal)) || endnode.signal.dif > tol )
      ready <- unique(c(ready, nextnodes))
    ready <- ready[-1]
  }
  return(list(node.signal[endnode,], endnode.signal.dif, nodes.vals))
}


###
# node.val<-nodes.vals[actnode,]
# subgraph

compute.node.signal2<-function(actnode, node.val, node.signal, subgraph, method="pond", response.tol = 0){
  
  incis <- incident(subgraph, actnode, mode="in")  
  # cat(incis)
  
  if(length(incis)==0){    
    signal <- rep(1, length(node.val))
    # cat(signal)
  } else {    
    
    # get activators and inhibitors signal
    prevs <- get.edgelist(subgraph)[incis,1]  
    input_signals <- node.signal[prevs,,drop=F]
    nas <- is.na(input_signals[,1])
    prevs <- prevs[!nas]
    incis <- incis[!nas]
    input_signals <- input_signals[!nas,,drop=F]
    typeincis <- E(subgraph)$relation[incis]
    activators <- typeincis==1
    nactivators <- sum(activators)
    inhibitors <- typeincis==-1
    ninhibitors <- sum(inhibitors)
    activator_signals <- input_signals[activators,,drop=F]
    inhibitor_signals <- input_signals[inhibitors,,drop=F]
    
    if( method == "sum"){
      s1 <- prettyifelse(nactivators>0, colSums(activator_signals), rep(1,length(node.val)))
      s2 <- prettyifelse(ninhibitors>0, colSums(inhibitor_signals), rep(0,length(node.val)))
      signal <- s1-s2
    }
    else if( method == "pond"){
      s1 <- prettyifelse(nactivators>0, apply(1- activator_signals, 2, prod), rep(0,length(node.val)))
      s2 <- prettyifelse(ninhibitors>0, apply(1- inhibitor_signals, 2, prod), rep(1,length(node.val)))      
      signal <- (1-s1)*s2
    }
    else if( method == "min"){
      s1 <- prettyifelse(nactivators>0, apply(activator_signals,2,min), rep(1,length(node.val)))
      s2 <- prettyifelse(ninhibitors>0, 1-apply(inhibitor_signals,2,max), rep(1,length(node.val)))
      signal <- s1*s2
    }
    else {
      stop("Unknown propagation rule")
    } 
    
    # If signal too low, signal do not propagate
    if(sum(nas) == 0 && signal < response.tol)
      signal <- rep(0,length(node.val))
    
  }
  
  signal[signal>1] <- 1
  signal[signal<0] <- 0
  signal <- signal*node.val    
  
  return(signal)
}

###
prettyifelse<-function(test,una,olaotra){
  if(test){
    return(una)
  } else {
    return(olaotra)
  }
}


############


############ for igraph v0.07
moduleSIFandGraph_07<-function(module_file_path,modulename,saveSIFpath,save_modulename=NULL, mystructure=NULL){
  
  if(is.null(save_modulename)){save_modulename<-modulename}
  if(!is.null(mystructure)){ modulename_w_in_out  <- modulename; modulename  <- strsplit(modulename,"_")[[1]][1]}
  
  library(KEGGgraph)
  library(graph)
  library(igraph)
  
  module<-read.delim(paste0(module_file_path,"/",modulename,".tsv"),sep="\t",stringsAsFactors=F,header=F)
  if(nrow(module)>1){
    
    module_pathways<-read.delim(paste0(module_file_path,"/",modulename,"pathways.txt"),sep="\t",stringsAsFactors=F,header=F)
    
    SIF_tmp1<-mat.or.vec(nr=nrow(module),nc=3)
    SIF_tmp2<-mat.or.vec(nr=nrow(module),nc=3)
    
    for(i in 1:nrow(module)){
      
      # met1 + met2 | reax1,reax2
      first_node<-strsplit(module$V1[i],split=" ")[[1]]
      reaction<-first_node[(grep("R",first_node))]
      compoundinputC<-first_node[(grep("C",first_node))]
      compoundinputG<-first_node[(grep("G",first_node))]
      second_node<-strsplit(module$V2[i],split=" ")[[1]]
      compoundoutputC<-second_node[(grep("C",second_node))]
      compoundoutputG<-second_node[(grep("G",second_node))]
      
      
      if(length(compoundinputC)>0){
        SIF_tmp1[i,1]<-paste0(compoundinputC,collapse="+")
      }else{SIF_tmp1[i,1]<-paste0(compoundinputG,collapse="+")}
      SIF_tmp1[i,2]<-"activation"
      SIF_tmp1[i,3]<-paste0(reaction,collapse="+")
      SIF_tmp2[i,1]<-paste0(reaction,collapse="+")
      SIF_tmp2[i,2]<-"activation"
      if(length(compoundoutputC)>0){
        SIF_tmp2[i,3]<-paste0(compoundoutputC,collapse="+")
      }else{SIF_tmp2[i,3]<-paste0(compoundoutputG,collapse="+")}
    }
    
    SIF<-rbind(SIF_tmp1,SIF_tmp2)
    # here can be fixed example M00009 "C00024+C00036" initial node now only C00024
    #inidnode<-strsplit(SIF[1,1],"[+]")[[1]][1]
    if(is.null(mystructure)){
      inidnode<-SIF[1,1]
      endnode<-SIF[nrow(SIF),3]
      cycle <- NULL
    }else{
      inidnode <- mystructure[[grep(modulename_w_in_out,names(mystructure))]]$init_nodes
      endnode <- mystructure[[grep(modulename_w_in_out,names(mystructure))]]$out_node
      cycle <- mystructure[[grep(modulename_w_in_out,names(mystructure))]]$cycle
      allnodes <- mystructure[[grep(modulename_w_in_out,names(mystructure))]]$nodes
    }
    #empty_parts<-unique(c(which(SIF[,1]==""),which(SIF[,3]=="")))
    #if(length(empty_parts)>0){SIF<-SIF[-empty_parts,] }
    
    compoundwithplus1<-unique(c(grep("[+]C",SIF[,1]),grep("[+]G",SIF[,1])))
    if(length(compoundwithplus1)>0){
      addings_right1<-sapply(compoundwithplus1,function(x)  {y<-strsplit(SIF[x,1],"[+]")[[1]]; z<-rep(SIF[x,3],length(y))})
      addings_left1<-sapply(compoundwithplus1,function(x)  {y<-strsplit(SIF[x,1],"[+]")[[1]]; return(y)})
      SIF_discard<-sapply(compoundwithplus1,function(x)  {x})
      SIF<-SIF[-SIF_discard,]
      addings_mat1<-mat.or.vec(nr=length(c(addings_left1)),nc=3)
      addings_mat1[,1]<-addings_left1
      addings_mat1[,2]<-rep("activation",length(addings_left1))
      addings_mat1[,3]<-addings_right1
      SIF<-rbind(SIF[-nrow(SIF),],addings_mat1,SIF[nrow(SIF),])  # takes in the middle of SIF to keep end node as orginal
    }
    
    compoundwithplus2<-unique(c(grep("[+]C",SIF[,3]),grep("[+]G",SIF[,3])))
    if(length(compoundwithplus2)>0){
      addings_right2<-sapply(compoundwithplus2,function(x)  {y<-strsplit(SIF[x,3],"[+]")[[1]]; z<-rep(SIF[x,1],length(y))})
      addings_left2<-sapply(compoundwithplus2,function(x)  {y<-strsplit(SIF[x,3],"[+]")[[1]]; return(y)})
      SIF_discard<-sapply(compoundwithplus2,function(x)  {x})
      SIF<-SIF[-SIF_discard,]
      addings_mat2<-mat.or.vec(nr=length(c(addings_left2)),nc=3)
      addings_mat2[,1]<-addings_right2
      addings_mat2[,2]<-rep("activation",length(addings_left2))
      addings_mat2[,3]<-addings_left2
      SIF<-rbind(SIF[-nrow(SIF),],addings_mat2,SIF[nrow(SIF),])
    }
    
    if(exists("allnodes")){
      idx_allnodes  <- unique(c(which(SIF[,1] %in% allnodes),  which(SIF[,2] %in% allnodes)))
      SIF <- SIF[idx_allnodes,]  
    }
    
    df <- SIF  
    # SIF TO igraph
    g<-graph.data.frame(df[,c(1,3)],directed=T,vertices=unique(c(df[,1],df[,3])))
    
    # deal with problems such M00029
    outdegree <- igraph::degree(graph=g,v=V(g),mode="out")
    out_nodes <-  V(g)$name[which(outdegree==0)]
    if(length(out_nodes)>1){
      discard_out_nodes <- out_nodes[out_nodes!=endnode]
      idx_don_left <- which(df[,1] %in%  discard_out_nodes)
      idx_don_right <- which(df[,3] %in%  discard_out_nodes)
      
      idx_don_all <- c(idx_don_left,idx_don_right)
      df <- df[-idx_don_all,]
    }
    g<-graph.data.frame(df[,c(1,3)],directed=T,vertices=unique(c(df[,1],df[,3])))
    g<-set.edge.attribute(g, "relation", value=rep(1,length(E(g))))
    
    if(!is.null(mystructure)){  
      write.table(df,paste0(saveSIFpath,"/",save_modulename,"SIF.txt"),sep="\t",row.names=F,quote=F,col.names=F)
    }
    
    return(list("igraph"=g,"initial_node"=inidnode,"end_node"=endnode, "cycle"=cycle,"SIF"=SIF,"modulePATH"=module_pathways,"module_name"=paste0(modulename,"_",endnode)))
  }else{return("short module is discarded")}
}



############ for igraph v1.0
moduleSIFandGraph<-function(module_file_path,modulename,saveSIFpath,save_modulename=NULL, mystructure=NULL){
  
  if(is.null(save_modulename)){save_modulename<-modulename}
  if(!is.null(mystructure)){ modulename_w_in_out  <- modulename; modulename  <- strsplit(modulename,"_")[[1]][1]}
  
  library(KEGGgraph)
  library(graph)
  library(igraph)
  
  module<-read.delim(paste0(module_file_path,"/",modulename,".tsv"),sep="\t",stringsAsFactors=F,header=F)
  if(nrow(module)>1){
    
    module_pathways<-read.delim(paste0(module_file_path,"/",modulename,"pathways.txt"),sep="\t",stringsAsFactors=F,header=F)
    
    SIF_tmp1<-mat.or.vec(nr=nrow(module),nc=3)
    SIF_tmp2<-mat.or.vec(nr=nrow(module),nc=3)
    
    for(i in 1:nrow(module)){
      
      # met1 + met2 | reax1,reax2
      first_node<-strsplit(module$V1[i],split=" ")[[1]]
      reaction<-first_node[(grep("R",first_node))]
      compoundinputC<-first_node[(grep("C",first_node))]
      compoundinputG<-first_node[(grep("G",first_node))]
      second_node<-strsplit(module$V2[i],split=" ")[[1]]
      compoundoutputC<-second_node[(grep("C",second_node))]
      compoundoutputG<-second_node[(grep("G",second_node))]
      
      
      if(length(compoundinputC)>0){
        SIF_tmp1[i,1]<-paste0(compoundinputC,collapse="+")
      }else{SIF_tmp1[i,1]<-paste0(compoundinputG,collapse="+")}
      SIF_tmp1[i,2]<-"activation"
      SIF_tmp1[i,3]<-paste0(reaction,collapse="+")
      SIF_tmp2[i,1]<-paste0(reaction,collapse="+")
      SIF_tmp2[i,2]<-"activation"
      if(length(compoundoutputC)>0){
        SIF_tmp2[i,3]<-paste0(compoundoutputC,collapse="+")
      }else{SIF_tmp2[i,3]<-paste0(compoundoutputG,collapse="+")}
    }
    
    SIF<-rbind(SIF_tmp1,SIF_tmp2)
    # here can be fixed example M00009 "C00024+C00036" initial node now only C00024
    #inidnode<-strsplit(SIF[1,1],"[+]")[[1]][1]
    #empty_parts<-unique(c(which(SIF[,1]==""),which(SIF[,3]=="")))
    #if(length(empty_parts)>0){SIF<-SIF[-empty_parts,] }
    
    compoundwithplus1<-unique(c(grep("[+]C",SIF[,1]),grep("[+]G",SIF[,1])))
    if(length(compoundwithplus1)>0){
      addings_right1<-sapply(compoundwithplus1,function(x)  {y<-strsplit(SIF[x,1],"[+]")[[1]]; z<-rep(SIF[x,3],length(y))})
      addings_left1<-sapply(compoundwithplus1,function(x)  {y<-strsplit(SIF[x,1],"[+]")[[1]]; return(y)})
      SIF_discard<-sapply(compoundwithplus1,function(x)  {x})
      SIF<-SIF[-SIF_discard,]
      addings_mat1<-mat.or.vec(nr=length(c(addings_left1)),nc=3)
      addings_mat1[,1]<-addings_left1
      addings_mat1[,2]<-rep("activation",length(addings_left1))
      addings_mat1[,3]<-addings_right1
      SIF<-rbind(SIF[-nrow(SIF),],addings_mat1,SIF[nrow(SIF),])  # takes in the middle of SIF to keep end node as orginal
    }
    
    compoundwithplus2<-unique(c(grep("[+]C",SIF[,3]),grep("[+]G",SIF[,3])))
    if(length(compoundwithplus2)>0){
      addings_right2<-sapply(compoundwithplus2,function(x)  {y<-strsplit(SIF[x,3],"[+]")[[1]]; z<-rep(SIF[x,1],length(y))})
      addings_left2<-sapply(compoundwithplus2,function(x)  {y<-strsplit(SIF[x,3],"[+]")[[1]]; return(y)})
      SIF_discard<-sapply(compoundwithplus2,function(x)  {x})
      SIF<-SIF[-SIF_discard,]
      addings_mat2<-mat.or.vec(nr=length(c(addings_left2)),nc=3)
      addings_mat2[,1]<-addings_right2
      addings_mat2[,2]<-rep("activation",length(addings_left2))
      addings_mat2[,3]<-addings_left2
      SIF<-rbind(SIF[-nrow(SIF),],addings_mat2,SIF[nrow(SIF),])
    }
    
    if(is.null(mystructure)){
      inidnode<-SIF[1,1]
      endnode<-SIF[nrow(SIF),3]
      cycle <- NULL
    }else{
      inidnode <- mystructure[[grep(modulename_w_in_out,names(mystructure))]]$init_nodes
      endnode <- mystructure[[grep(modulename_w_in_out,names(mystructure))]]$out_node
      cycle <- mystructure[[grep(modulename_w_in_out,names(mystructure))]]$cycle
      allnodes <- mystructure[[grep(modulename_w_in_out,names(mystructure))]]$nodes
      idx_allnodes  <- unique(c(which(SIF[,1] %in% allnodes),  which(SIF[,2] %in% allnodes)))
      SIF <- SIF[idx_allnodes,]  
    }
    
    df <- as.data.frame(SIF[,c(1,3)],stringsAsFactors = F)
    # SIF TO igraph
    g<-graph.data.frame(d = df,directed=T,vertices=unique(c(df[,1],df[,2])))
    
    # deal with problems such M00029
    outdegree <- igraph::degree(graph=g,v=V(g),mode="out")
    out_nodes <-  V(g)$name[which(outdegree==0)]
    if(length(out_nodes)>1 & !is.null(mystructure)){
      discard_out_nodes <- out_nodes[out_nodes!=endnode]
      idx_don_left <- which(df[,1] %in%  discard_out_nodes)
      idx_don_right <- which(df[,2] %in%  discard_out_nodes)
      
      idx_don_all <- c(idx_don_left,idx_don_right)
      df <- df[-idx_don_all,]
    }
    g<-graph.data.frame(d = df,directed=T,vertices=unique(c(df[,1],df[,2])))
    g<-set.edge.attribute(g, "relation", value=rep(1,length(E(g))))
    
    if(!is.null(mystructure)){  
      write.table(df,paste0(saveSIFpath,"/",save_modulename,"SIF.txt"),sep="\t",row.names=F,quote=F,col.names=F)
    }
    
    return(list("igraph"=g,"initial_node"=inidnode,"end_node"=endnode, "cycle"=cycle,"SIF"=SIF,"modulePATH"=module_pathways,"module_name"=paste0(modulename,"_",endnode)))
  }else{return("short module is discarded")}
}






###
# pathID="hsa00330"
# KGMLsave_dir="/home/cankut/Desktop/Metabolic pathyway/KEGG_module/KEGG_module_data/KGML/"
# SIFsave_dir="/home/cankut/Desktop/"
# KGML_local_dir="/home/cankut/Desktop/Metabolic pathyway/KEGG_module/moduleWebApp/hsa/KGML/"


KEGG_Met_path_parse<-function(KGML_local_dir=NULL, pathID, KGMLsave_dir, SIFsave_dir,organism,biomaRt_data){
  print(pathID)
  dir.create(SIFsave_dir,showWarnings=F,recursive=T)
  if(!is.null(KGML_local_dir)){
    toyKGML<-paste0(KGML_local_dir,"/",pathID,".xml")
    if(!file.exists(toyKGML)){
      toyKGML<- paste0("http://rest.kegg.jp/get/",pathID,"/kgml")
      #toyKGML<-paste0("http://www.kegg.jp/kegg-bin/download?entry=",pathID, "&format=kgml")
      download.file(toyKGML,destfile=paste0(KGMLsave_dir,"/",pathID,".xml")) 
    }
  }else{
    toyKGML<- paste0("http://rest.kegg.jp/get/",pathID,"/kgml")
    #toyKGML<-paste0("http://www.kegg.jp/kegg-bin/download?entry=",pathID, "&format=kgml")
    download.file(toyKGML,destfile=paste0(KGMLsave_dir,"/",pathID,".xml"))
  }
  
  pathway <- parseKGML(toyKGML)
  
  entryID<-sapply(pathway@nodes,function(x) x@entryID)
  name<-unlist(sapply(pathway@nodes,function(x) paste0(x@name,collapse=" ")))
  type<-sapply(pathway@nodes,function(x) x@type)
  link<-sapply(pathway@nodes,function(x) x@link)
  reaction<-sapply(pathway@nodes,function(x) x@reaction)
  # "rn:R01230 rn:R01231" bunla ilgili sorun va
  # df icindede bundan var R04560,R06975 
  # reaction<-gsub("rn:","",reaction)
  graphics_name<-unlist(sapply(pathway@nodes,function(x) x@graphics@name,simplify=F))
  graphics_X<-sapply(pathway@nodes,function(x) x@graphics@x)
  graphics_Y<-sapply(pathway@nodes,function(x) x@graphics@y)
  
  # if it is (the shape) a "line" it is smthg like mmu00510 (part below)
  graphics_shape<-sapply(pathway@nodes,function(x) x@graphics@type)
  graphics_width<-sapply(pathway@nodes,function(x) x@graphics@width)
  graphics_height<-sapply(pathway@nodes,function(x) x@graphics@height)
  graphics_color<-sapply(pathway@nodes,function(x) x@graphics@bgcolor)
  
  direction <- sapply(pathway@reactions,function(x) x@type)
  names(direction) <- sapply(pathway@reactions,function(x) x@name)
  idx_direction <- match(names(direction),reaction)
  direction_attr <- rep(NA,length(reaction))
  direction_attr[idx_direction[which(!is.na(idx_direction))]] <- direction[which(!is.na(idx_direction))]
  
  KEGG_met_path_node<-data.frame(entryID,name,type,reaction,graphics_name,graphics_color,graphics_shape,graphics_width,"x"=graphics_X,"y"=graphics_Y*-1,direction_attr,stringsAsFactors=F)
  
  #     # Reaction1->subtype(compound,metabolite)->Reaction2
  #     edge_type<-sapply(pathway@edges,function(x) x@type)
  #     edge_entry1<-sapply(pathway@edges,function(x) x@entry1ID)
  #     edge_entry2<-sapply(pathway@edges,function(x) x@entry2ID)
  #     #edge_subtype<-sapply(pathway@edges,function(x) x@subtype$subtype@value)
  #     edge_subtype<-sapply(pathway@edges, function(x){ ifelse(is.null(x@subtype$subtype),NA,x@subtype$subtype@value)})
  
  #     KEGG_met_path_SIF_tmp1<-mat.or.vec(nr=length(edge_entry2),nc=3)
  #     KEGG_met_path_SIF_tmp2<-mat.or.vec(nr=length(edge_entry2),nc=3)
  #     
  #     for(i in 1:length(edge_type)){
  #       
  #       KEGG_met_path_SIF_tmp1[i,1]<-edge_entry1[i]
  #       KEGG_met_path_SIF_tmp1[i,2]<-edge_type[i]
  #       KEGG_met_path_SIF_tmp1[i,3]<-edge_subtype[i]
  #       
  #       KEGG_met_path_SIF_tmp2[i,1]<-edge_subtype[i]
  #       KEGG_met_path_SIF_tmp2[i,2]<-edge_type[i]
  #       KEGG_met_path_SIF_tmp2[i,3]<-edge_entry2[i]
  #       
  #     }
  #     KEGG_met_path_SIF<-rbind(KEGG_met_path_SIF_tmp1,KEGG_met_path_SIF_tmp2)
  #     rm(KEGG_met_path_SIF_tmp1,KEGG_met_path_SIF_tmp2)
  #     
  #     KEGG_met_path_SIF[,1]<- unlist(sapply(KEGG_met_path_SIF[,1],function(x)  graphics_name[which(entryID %in% x)]))
  #     KEGG_met_path_SIF[,3]<- unlist(sapply(KEGG_met_path_SIF[,3],function(x)  graphics_name[which(entryID %in% x)]))
  
  ## Reaction genes April2016
  entrezIDs <- sapply(name, function(x) {
    isgene<-grep(paste0(organism,":"),x)
    if(length(isgene)>0){
      #my_ID<-as.character(sapply(strsplit(x," ")[[1]],function(y) strsplit(y,":")[[1]][2]))
      my_ID<-gsub(" ",";",gsub(paste0(organism,":"),"",x))
    }else{
      my_ID<-"isNOTGene"    
    } 
    return(my_ID)
  })
  
  KEGG_met_path_node$hgnc_symbols<-NA
  KEGG_met_path_node$entrezIDs<-entrezIDs
  
  if(1==0){
    # I add this part because one reaction can be catalyzed by various enzyme. check hsa00010 R01786 1) 3098;3099;3101;80201 2) 2645 (April2016)
    # R01600 and R01786 have same genes. xml has info about R01600 but if I search it on KEGG image ther is no result(http://www.genome.jp/kegg-bin/show_pathway?hsa00010+R01600), WHY? 
    rxns_for_entrez<-unique(KEGG_met_path_node$reaction[!is.na(KEGG_met_path_node$reaction)])  
    for(rr in rxns_for_entrez){ 
      idx_rr <- which(KEGG_met_path_node$reaction %in% rr)
      entrezID_graph <- unique(entrezIDs[idx_rr[which(entrezIDs[idx_rr]!="isNOTGene")]])
      if(length(entrezID_graph)>0){
        KEGG_met_path_node$entrezIDs[idx_rr] <- paste0(entrezID_graph,collapse = ";")
      }
    }
  }
  
  ## !!!! there are multiplied rows, they come from 1 gene different enzymes(E.C. numbers)
  # e.g. http://www.genome.jp/kegg-bin/show_pathway?scale=1.0&query=BPGM&map=hsa00010&scale=&auto_image=&show_description=hide&multi_query=
  # maybe I have to generate SIF with EC numbers
  # KEGG_met_path_SIF<-unique(KEGG_met_path_SIF)
  ################################################
  
  ########### prepare SIF File ########### 
  
  reaction_name <- sapply(pathway@reactions,function(x) x@name)
  reaction_name_substrate <- list()
  for(i in 1:length(pathway@reactions)){ reaction_name_substrate[[i]] <- pathway@reactions[[i]]@substrateName }
  reaction_name_product <- list()
  for(i in 1:length(pathway@reactions)){ reaction_name_product[[i]] <- pathway@reactions[[i]]@productName }
  
  SIF_left <- c()
  SIF_rigft <- c() 
  
  idx_push <- 1
  for(s in 1:length(reaction_name)){  
    
    for(l in 1:length(reaction_name_substrate[[s]])){
      SIF_left[idx_push] <- reaction_name_substrate[[s]][l]
      SIF_rigft[idx_push] <- reaction_name[s]
      idx_push <- idx_push + 1
      if(direction[s]=="reversible"){
        SIF_left[idx_push] <- reaction_name[s]
        SIF_rigft[idx_push] <- reaction_name_substrate[[s]][l]
        idx_push <- idx_push + 1
      }
    }
  }
  
  
  for(s in 1:length(reaction_name)){  
    
    for(l in 1:length(reaction_name_product[[s]])){
      SIF_left[idx_push] <-  reaction_name[s]
      SIF_rigft[idx_push] <- reaction_name_product[[s]][l]
      idx_push <- idx_push + 1
      if(direction[s]=="reversible"){
        SIF_left[idx_push] <-  reaction_name_product[[s]][l]
        SIF_rigft[idx_push] <- reaction_name[s]
        idx_push <- idx_push + 1
      }
    }
  }
  
  isSIF   <- data.frame(SIF_left,"to",SIF_rigft,stringsAsFactors=F)
  
  ########### prepare Node Attribute File ########### 
  
  idx_attr_gene <-which(type=="gene")
  rxnID <- reaction[idx_attr_gene]
  rxnX <- graphics_X[idx_attr_gene]
  rxnY <- graphics_Y[idx_attr_gene]
  rxnEntry <- entryID[idx_attr_gene]
  rxnGN <- graphics_name[idx_attr_gene]
  rxnGS <- graphics_shape[idx_attr_gene]
  rxnEntrez <- KEGG_met_path_node$entrezIDs[idx_attr_gene]
  rxnDirection <- direction_attr[idx_attr_gene]
  
  idx_attr_cpd <- which(type=="compound")
  cpdID <- name[idx_attr_cpd]
  cpdX <- graphics_X[idx_attr_cpd]
  cpdY <- graphics_Y[idx_attr_cpd]
  cpdEntry <- entryID[idx_attr_cpd]
  cpdGN <- graphics_name[idx_attr_cpd]
  cpdGS <- graphics_shape[idx_attr_cpd]
  cpdEntrez <- KEGG_met_path_node$entrezIDs[idx_attr_cpd]
  cpdDirection <- direction_attr[idx_attr_cpd]
  
  attrName <- c(rxnID,cpdID)
  attrX <- c(rxnX,cpdX)
  attrY <- c(rxnY,cpdY)
  attrEntry <- c(rxnEntry,cpdEntry)
  attrGN <- c(rxnGN,cpdGN)
  attrShape <- c(rxnGS,cpdGS)
  attrEntrez <- c(rxnEntrez,cpdEntrez)
  attrType <- c(rep("rxn",length(idx_attr_gene)),rep("cpd",length(idx_attr_cpd)))
  attrDirection <- c(rxnDirection,cpdDirection)
  
  
  ## if there is a pathways like M00510, can have some nodes without coordinates. See below M00510
  idx_no_coord <- unique(c(which(is.na(attrX)), which(is.na(attrY))))
  if(length(idx_no_coord)!=0){
    attrName <-  attrName[-idx_no_coord]
    attrX <- attrX[-idx_no_coord]
    attrY <- attrY[-idx_no_coord]
    attrShape <- attrShape[-idx_no_coord]
    attrEntry <- attrEntry[-idx_no_coord]
    attrGN <- attrGN[-idx_no_coord]
    attrEntrez <- attrEntrez[-idx_no_coord]
    attrType <- attrType[-idx_no_coord]
    attrDirection <- attrDirection[-idx_no_coord]
  }
  
  isAttr <- data.frame(attrName,"x"=attrX,"y"=attrY*-1,attrShape, attrEntry,attrGN, attrEntrez, attrType, attrDirection, stringsAsFactors=F)
  
  ## add missing node information #
  ## e.g cpd:C00085 in hsa00520 has no information it has id="100" and C05345 also has id="100". Strange node.
  
  nodesList <- unique(c(isSIF$SIF_left,isSIF$SIF_rigft))
  idx_missing_node_info <- match(nodesList,isAttr$attrName)
  missing_node_names <- nodesList[which(is.na(idx_missing_node_info))]
  
  if(length(missing_node_names)>0){
    
    for(m in 1:length(missing_node_names)){
      
      node_left <- isSIF[which(isSIF==missing_node_names[m],arr.ind=T)[,1],"SIF_left"]
      node_right <- isSIF[which(isSIF==missing_node_names[m],arr.ind=T)[,1],"SIF_rigft"]
      
      node_LR <- unique(c(node_left,node_right))
      node_LR <- node_LR[node_LR!=missing_node_names[m]]
      missing_node_x<- mean(isAttr$x[match(node_LR,isAttr$attrName)])
      missing_node_y<- mean(isAttr$y[match(node_LR,isAttr$attrName)])
      
      index_add_missing <- nrow(isAttr)
      if(length(grep("cpd",missing_node_names[m]))==1){
        missing_attrShape <- "circle"
        missing_attrEntry <- paste0("missing_node",m)
        missing_attrGN <- missing_node_names[m]
        missing_attrEntrez <- "isNOTGene"
        missing_attrType <- "cpd"
        missing_attrDirection <- NA
      }else if(length(grep("gl",missing_node_names[m]))==1){   
        missing_attrShape <- "circle"
        missing_attrEntry <- paste0("missing_node",m)
        missing_attrGN <- missing_node_names[m]
        missing_attrEntrez <- "isNOTGene"
        missing_attrType <- "cpd"
        missing_attrDirection <- NA
      }else if(length(grep("rn",missing_node_names[m]))==1){
        missing_attrShape <- "rectangle" 
        missing_attrEntry <- paste0("missing_node",m)
        missing_attrGN <- missing_node_names[m]
        missing_attrEntrez <- NA
        missing_attrType <- "rxn"
        missing_attrDirection <- "undirected"
      }else(cat("unknown node type"))
      isAttr[index_add_missing,] <- c(missing_node_names[m],missing_node_x,missing_node_y,missing_attrShape,missing_attrEntry,missing_attrGN,missing_attrEntrez,missing_attrType,missing_attrDirection)
      
    }
  }
  
  ##
  
  
  if(1==0){
    
    ## Convert geneIDs
    # this was part is not any more fuctional. I added code above. The geneIds are taken from name instead of graphics name April2016
    load(biomaRt_data)
    hgnc_symbols<-c()
    entrezIDs<-c()
    
    for(gene in 1:nrow(KEGG_met_path_node)){
      if(KEGG_met_path_node$type[gene]=="gene"){
        geneX<-KEGG_met_path_node$graphics_name[gene]
        geneY<-sapply(strsplit(gsub("[..]","",geneX),split="[,] "),function(x) as.vector(x))
        x<-1
        while(x<=length(geneY)){
          entrezID<-as.character()
          trans.table <-getBM(attributes=c("rgd_symbol","entrezgene"), filters = "rgd_symbol", values = as.vector(geneY)[x], mart = ensembl_rat)
          entrezID<-trans.table[1,2]
          hgnc_symbol<-trans.table[1,1]
          if(trans.table$entrezgene>0 && !is.na(trans.table$entrezgene) && nrow(trans.table)>0){
            break
          }
          x<-x+1
        }
        if(is.na(entrezID)){
          entrezID<-"NotConverted"
          hgnc_symbol<-as.vector(geneY)[1]
        }else{entrezID<-entrezID
              hgnc_symbol<-hgnc_symbol
        }
      }else{entrezID<-"isNOTGene"; hgnc_symbol<-"isNOTGene"}
      hgnc_symbols<-c(hgnc_symbols,hgnc_symbol)
      entrezIDs<-c(entrezIDs,entrezID)
    }
    KEGG_met_path_node$hgnc_symbols<-hgnc_symbols
    KEGG_met_path_node$entrezIDs<-entrezIDs
  }
  
  
  ## corrdinates are not conpatible with cellmaps because same gene, compound can be found several times inside the KGML
  write.table(KEGG_met_path_node,paste0(SIFsave_dir,"/",pathID,"_node.txt"), sep="\t", row.names=F, quote=F, col.names=T)
  write.table(isSIF,paste0(SIFsave_dir,"/",pathID,"_SIF.txt"), sep="\t", row.names=F, quote=F, col.names=F)
  write.table(isAttr,paste0(SIFsave_dir,"/",pathID,"_attr.txt"), sep="\t", row.names=F, quote=F, col.names=T)
  
  # This changed as above because 2 different reaction can be controled by same gene. In this case we miss gene expression data for some Rxns  
  #   conversionEntrez<-unique(KEGG_met_path_node$entrezIDs[KEGG_met_path_node$entrezIDs!="isNOTGene"])
  #   idx_conversionReaction<-sapply(conversionEntrez,function(x) which(KEGG_met_path_node$entrezIDs %in% x)[1])
  #   conversionReaction<-KEGG_met_path_node$reaction[idx_conversionReaction]
  #   conversionFrame<-data.frame(EntrezID=conversionEntrez,ReaxID=conversionReaction,stringsAsFactors=F)
  #   
  
  discard<-unique(c(which(is.na(KEGG_met_path_node$reaction)),which(KEGG_met_path_node$entrezIDs=="isNOTGene")))
  conversionReaction<-unique(KEGG_met_path_node$reaction[-discard])
  KEGG_met_path_node_sub<- KEGG_met_path_node[-discard,]
  #idx_conversionReaction<-sapply(conversionReaction,function(x) which(KEGG_met_path_node$reaction %in% x)[1])
  #conversionEntrez<-KEGG_met_path_node_sub$entrezIDs[idx_conversionReaction]
  #conversionFrame<-data.frame(EntrezIDs=conversionEntrez,ReaxID=gsub("rn:","",conversionReaction),stringsAsFactors=F)
  idx_conversionReaction<-sapply(conversionReaction,function(x) which(KEGG_met_path_node_sub$reaction %in% x))
  
  if(is.matrix(idx_conversionReaction)){ idx_conversionReaction<-idx_conversionReaction[1,] }
  
  entrezsCR <- c()
  for(icr in names(idx_conversionReaction)){
    entrezsCR <- c(entrezsCR,paste(KEGG_met_path_node_sub$entrezIDs[idx_conversionReaction[[icr]]],collapse = ";"))
  }
  #conversionEntrez<-KEGG_met_path_node_sub$entrezIDs[idx_conversionReaction]
  conversionFrame<-data.frame(EntrezIDs=entrezsCR, ReaxID=gsub("rn:","",conversionReaction),stringsAsFactors=F)
  
  # there are some rxn ids which are duplicated. for that reason I discard "isNOTGene" again. e.g. "R00014" in hsa00010" 
  discard2<-which(conversionFrame$EntrezIDs=="isNOTGene")
  if(length(discard2)>0){ conversionFrame<- conversionFrame[-discard2,]}
  
  return(list("KEGG_met_path_node_info"=KEGG_met_path_node,"Reax2Entrez"=conversionFrame))
}




#graphobject<-hsa_module_data$M00036$graphobject
#KEGG_met_path_node<-hsa_module_data$M00036$KEGG_met_path_node
# rxn_vals_mat<-GeneExp
# default_value=0.5

## simplify node_vals function
# this function checks and ignores missing data(gene exp or node val)
# does not require prior function such as add missing data
node_vals<-function(graphobject,KEGG_met_path_node,rxn_vals_mat,default_value){
  
  # graphobject is output of moduleSIFandGraph
  # KEGG_met_path_node is output of KEGG_Met_path_parse
  # rxn_vals_mat is nrow=n rxns ncol= m samples nxm matrix
  
  if(!is.null(graphobject$reduced_SIF)){SIF<-graphobject$reduced_SIF}else{SIF<-graphobject$SIF}
  
  entries<-unique(c(SIF[,1],SIF[,3]))
  nodes.vals.total<-data.frame(name=seq(1:length(entries)))
  
  idxleft<-grep("R",SIF[,1])
  SIFleft<-SIF[,1][idxleft]
  idxright<-grep("R",SIF[,3])
  SIFright<-SIF[,3][idxright]
  reactions<-unique(c(SIFleft,SIFright))
  
  y<-sapply(reactions,function(x) {strsplit(x,split=",")[[1]]})
  if(all(reactions == y)){y<-as.list(y)}
  # fake data
  # y[1]<-list(c("R1+R2","R3","R4")) #c)
  # y[2]<-list(c("R5+R6")) # a
  # y[3]<-list(c("R7","R8")) #h
  # y[4]<-list(c("R9")) #e
  # y[5]<-list(c("R19+R20+R21","R22+R23")) #b
  # y[6]<-list(c("R24+R25+R26","R27+R28","R29","R30")) #d
  
  y1 <- sapply(y,function(node) {
    
    #complex_enzyme_rxn
    complex_enzyme_rxn <-node[grep("[+]",node)]
    
    if(length(complex_enzyme_rxn)>0){
      node_comp<- sapply(complex_enzyme_rxn,function(comp_node) {
        
        r_comp <- strsplit(comp_node,"[+]")[[1]]
        r_comp <- r_comp[r_comp %in% rownames(rxn_vals_mat)]
        if(length(r_comp)>1){
          calc_comp <- apply(rxn_vals_mat[which(rownames(rxn_vals_mat) %in% r_comp),],2,min)
        }else if(length(r_comp)==1){
          calc_comp <- rxn_vals_mat[which(rownames(rxn_vals_mat) %in% r_comp),]
        }else{ calc_comp <- NULL }
        return(calc_comp)    
      })
      
      # HERE: Kinza change this class assesment !
      # if(class(node_comp)=="list"){
      if(is.list(node_comp)){
        node_comp <- node_comp[sapply(node_comp,function(x) !is.null(x))]
        node_comp  <- do.call("rbind",node_comp)
      }else{ node_comp <- t(node_comp) }
      
      
      #isoenzymes_rxn
      isoenzymes_rxn <- node[grep(invert = T,"[+]",node)]
      
      if(length(isoenzymes_rxn)>0){
        isoenzymes_rxn <- isoenzymes_rxn[isoenzymes_rxn %in% rownames(rxn_vals_mat)]
        node_iso <- rxn_vals_mat[which(rownames(rxn_vals_mat) %in% isoenzymes_rxn),]
        # HERE: Kinza change class assesment :
        # if(class(node_iso)=="numeric") { node_iso <- t(as.matrix(node_iso))}
        if(any(class(node_iso)=="numeric")) { node_iso <- t(as.matrix(node_iso)) }
        # if(isNumeric(node_iso)) { node_iso <- t(as.matrix(node_iso)) }
        if(nrow(node_iso)==0){ node_iso<-NULL }
      }else{node_iso<-NULL}
      
      comp_iso_rxn<-rbind(node_comp,node_iso)
      if(!is.null(comp_iso_rxn)){
        calculated_node_val <- apply(comp_iso_rxn,2,max)
      }else{calculated_node_val<-default_value}
      
    }else{
      
      node <- node[which(node %in% rownames(rxn_vals_mat))]
      idx_iso_node<-which(rownames(rxn_vals_mat) %in% node)
      
      if(length(idx_iso_node)>1){
        node_iso <- rxn_vals_mat[idx_iso_node,]
        calculated_node_val <- apply(node_iso,2,max) 
      }else if(length(idx_iso_node)==1){
        calculated_node_val <- rxn_vals_mat[idx_iso_node,]
      }else{ calculated_node_val<-default_value }
      
    }
    return(calculated_node_val) 
  } )
  # HERE : Kinza change this text to be compatible with the new versio of R
  # if(class(y1)=="matrix"){
  if(is.matrix(y1)){
    node.vals<-t(y1)
  }else if(class(y1)=="list"){  # HERE . is list instead?
    node.vals<-do.call(what="rbind",y1)
  }else{ node.vals<-mat.or.vec(nr = length(y1) ,nc = ncol(rxn_vals_mat))
         colnames(node.vals) <- colnames(rxn_vals_mat)
         rownames(node.vals) <- names(y1)
         #node.vals[,1] <- default_value
         #node.vals[,2] <- default_value
         node.vals[,] <- default_value
         message(sprintf("Warning for %s: All node vals are %s",graphobject$module_name,default_value)) 
  }
  print("node.vals:")
  print(node.vals)
  if(length(is.na(node.vals[,1]))>0){
    node.vals[is.na(node.vals[,1]),] <- default_value
  }
  
  # nodes.vals.total <- t(apply(node.vals, 1, function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))))
  # idx_NAs_node_vals <- unique(c(which(is.na(nodes.vals.total),arr.ind = T)[,1],which(nodes.vals.total=="NaN",arr.ind = T)[,1]))
  # nodes.vals.total[idx_NAs_node_vals,]<-0.5
  
  nodes.vals.total <- node.vals
  
  return(nodes.vals.total)
  
}


# graphobject<-hsa_module_data$M00141$graphobject
# KEGG_met_path_node<-hsa_module_data$M00141$KEGG_met_path_node
# rxn_vals_mat<-GeneExp
# default_value=0.5

## simplify node_vals function
# this function checks and ignores missing data(gene exp or node val)
# does not require prior function such as add missing data
node_vals_from_gex<-function(graphobject,KEGG_met_path_node,gene_exp_mat,default_value=0.5){
  
  # graphobject is output of moduleSIFandGraph
  # KEGG_met_path_node is output of KEGG_Met_path_parse
  # gene_exp_mat is nrow=genes ncol=1m samples sammple nxm matrix
  
  if(!is.null(graphobject$reduced_SIF)){SIF<-graphobject$reduced_SIF}else{SIF<-graphobject$SIF}
  
  entries<-unique(c(SIF[,1],SIF[,3]))
  nodes.vals.total<-data.frame(name=seq(1:length(entries)))
  
  for(sample in 1:ncol(gene_exp_mat)){
    
    idxleft<-grep("R",SIF[,1])
    SIFleft<-SIF[,1][idxleft]
    idxright<-grep("R",SIF[,3])
    SIFright<-SIF[,3][idxright]
    reactions<-unique(c(SIFleft,SIFright))
    
    ## default data for non Reaction node entries. like compound. metabolite, glycan ...etc   # !!!!!!!!!! it seems I do not need these part
    SIFleft_nonreax<-SIF[,1][-idxleft]
    SIFright_nonreax<-SIF[,3][-idxright]
    non_reaction_enrties<-unique(c(SIFleft_nonreax,SIFright_nonreax))
    non_reaction_enrties_frame<-as.data.frame(mat.or.vec(nr=length(non_reaction_enrties),nc=ncol(gene_exp_mat)),stringsAsFactors=F)
    if(nrow(non_reaction_enrties_frame)>=1){non_reaction_enrties_frame[,]<-default_value} # sebep: hic compound yoksa bos frame donuyo ve burada hata verio. 
    #non_reaction_enrties_frame[,]<-runif(nrow(non_reaction_enrties_frame)*ncol(non_reaction_enrties_frame), 0,1) # as an example random generated data
    ###
    
    y<-sapply(reactions,function(x) {strsplit(x,split=",")[[1]]})
    if(all(reactions == y)){y<-as.list(y)}
    # fake data
    # y[1]<-list(c("R1+R2","R3","R4")) #c)
    # y[2]<-list(c("R5+R6")) # a
    # y[3]<-list(c("R7","R8")) #h
    # y[4]<-list(c("R9")) #e
    # y[5]<-list(c("R19+R20+R21","R22+R23")) #b
    # y[6]<-list(c("R24+R25+R26","R27+R28","R29","R30")) #d
    
    
    
    for(m in 1:length(y)){
      
      #x<-as.vector(y[[3]])
      x<-as.vector(y[[m]])
      if_complex<-grep("[+]R",x)
      
      # a)
      if(length(x)==1 && length(if_complex)>0 ){
        
        complex<-unlist(strsplit(x,split="[+]"))
        
        my_expression_data<-list()
        for(comp in complex){ my_entezIDs<-KEGG_met_path_node[[2]]$EntrezID[which(KEGG_met_path_node[[2]]$ReaxID %in% comp)]; #print(my_entezIDs);
                              if(length(my_entezIDs)>0){
                                my_gene<-which(rownames(gene_exp_mat) %in% my_entezIDs)
                                if(length(my_gene)>0){  
                                  my_expression_data[paste0(comp,collapse="+")]<-as.numeric(gene_exp_mat[my_gene,sample]) }}}
        
        if(length(my_expression_data)>0){
          ccc<-min(unlist(my_expression_data))
        }else{ccc<-default_value 
              #print(sprintf("NODE:%s has default value as %s",paste0(x,collapse=","),default_value))
        }
        #print("# a")
        #print(ccc)
        y[[m]] <- as.numeric(ccc)      
        
        # b)  
      }else if(length(x)>1 && length(if_complex)>1 && length(if_complex)==length(x)){
        
        x<-x[if_complex]
        aaa<-sapply(x,function(xs)  strsplit(xs,split="[+]"))
        
        my_expression_data<-list()
        for(comp in aaa){ my_entezIDs<-KEGG_met_path_node[[2]]$EntrezID[which(KEGG_met_path_node[[2]]$ReaxID %in% comp)]; #print(my_entezIDs);
                          if(length(my_entezIDs)>0){
                            my_gene<-which(rownames(gene_exp_mat) %in% my_entezIDs)
                            if(length(my_gene)>0){  
                              my_expression_data[paste0(comp,collapse="+")]<-as.numeric(gene_exp_mat[my_gene,sample]) }}}
        
        if(length(my_expression_data)>0){
          bbb<-sapply(my_expression_data,function(xs) min(xs))    
          ccc<- max(bbb)
        }else{ccc<-default_value
              #print(sprintf("NODE:%s has default value as %s",paste0(x,collapse=","),default_value))
        }
        #print("# b")
        #print(ccc)
        y[[m]] <- as.numeric(ccc)      
        
        # c)
      }else if(length(x)>1 && length(if_complex)==1 && length(if_complex)<length(x)){ 
        
        x1<-x[if_complex]
        aaa1<-sapply(x1,function(xs)  strsplit(xs,split="[+]"))
        
        my_expression_data<-list()
        for(comp in aaa1){ my_entezIDs<-KEGG_met_path_node[[2]]$EntrezID[which(KEGG_met_path_node[[2]]$ReaxID %in% comp)]; #print(my_entezIDs);
                           if(length(my_entezIDs)>0){
                             my_gene<-which(rownames(gene_exp_mat) %in% my_entezIDs)
                             if(length(my_gene)>0){  
                               my_expression_data[paste0(comp,collapse="+")]<-as.numeric(gene_exp_mat[my_gene,sample]) }}}
        
        if(length(my_expression_data)>0){
          bbb1<-sapply(my_expression_data,function(xs) min(xs))
        }else{bbb1<-list()}
        
        x2<-x[-if_complex]
        
        my_expression_data2<-list()
        for(comp in x2){ my_entezIDs<-KEGG_met_path_node[[2]]$EntrezID[which(KEGG_met_path_node[[2]]$ReaxID %in% comp)]; #print(my_entezIDs);
                         if(length(my_entezIDs)>0){
                           my_gene<-which(rownames(gene_exp_mat) %in% my_entezIDs)
                           if(length(my_gene)>0){  
                             my_expression_data2[comp]<-as.numeric(gene_exp_mat[my_gene,sample]) }}}
        
        #for(comp in x2){ my_entezIDs<-KEGG_met_path_node$entrezIDs[which(KEGG_met_path_node$reaction  %in% comp)[1]];  my_expression_data2[comp]<-as.numeric(gene_exp_mat[which(rownames(gene_exp_mat) %in% my_entezIDs),1]) }
        
        
        
        if(length(my_expression_data2)>0 && length(bbb1)>0){
          my_expression_data2<-max(unlist(my_expression_data2))
          ccc<-max(my_expression_data2,bbb1)
        }else if(length(my_expression_data2)==0 && length(bbb1)>0){
          ccc<-max(bbb1)
        }else if(length(my_expression_data2)>0 && length(bbb1)==0){
          my_expression_data2<-max(unlist(my_expression_data2))
          ccc<-max(my_expression_data2)
        }else{ccc<-default_value
              #print(sprintf("NODE:%s has default value as %s",paste0(x,collapse=","),default_value))
        }
        
        #print("# c")
        #print(ccc)
        y[[m]] <- as.numeric(ccc)      
        
        # d)  
      }else if(length(x)>1 && length(if_complex)>1 && length(if_complex)<length(x)){
        
        x1<-x[if_complex]
        aaa1<-sapply(x1,function(xs)  strsplit(xs,split="[+]"))
        
        my_expression_data<-list()
        for(comp in aaa1){ my_entezIDs<-KEGG_met_path_node[[2]]$EntrezID[which(KEGG_met_path_node[[2]]$ReaxID %in% comp)]; #print(comp);print(my_entezIDs);
                           if(length(my_entezIDs)>0){
                             my_gene<-which(rownames(gene_exp_mat) %in% my_entezIDs)
                             if(length(my_gene)>0){  
                               my_expression_data[paste0(comp,collapse="+")]<-as.numeric(gene_exp_mat[my_gene,sample]) }}}
        
        if(length(my_expression_data)>0){
          bbb1<-sapply(my_expression_data,function(xs) min(xs))
        }else{bbb1<-list()}
        
        x2<-x[-if_complex]
        
        my_expression_data2<-list()
        for(comp in x2){ my_entezIDs<-KEGG_met_path_node[[2]]$EntrezID[which(KEGG_met_path_node[[2]]$ReaxID %in% comp)]; #print(my_entezIDs);
                         if(length(my_entezIDs)>0){
                           my_gene<-which(rownames(gene_exp_mat) %in% my_entezIDs)
                           if(length(my_gene)>0){  
                             my_expression_data2[comp]<-as.numeric(gene_exp_mat[my_gene,sample]) }}}
        
        
        
        if(length(my_expression_data2)>0 && length(bbb1)>0){
          ccc<-max(my_expression_data2,bbb1)
        }else if(length(my_expression_data2)==0 && length(bbb1)>0){
          ccc<-max(bbb1)
        }else if(length(my_expression_data2)>0 && length(bbb1)==0){
          ccc<-max(my_expression_data2)
        }else{ccc<-default_value
              #print(sprintf("NODE:%s has default value as %s",paste0(x,collapse=","),default_value))
        }
        
        #print("# d")
        #print(ccc)
        y[[m]] <- as.numeric(ccc)      
        
        # e)
      }else if(length(x)==1 && length(if_complex)==0){
        
        my_expression_data<-list()
        for(comp in x){ my_entezIDs<-KEGG_met_path_node[[2]]$EntrezID[which(KEGG_met_path_node[[2]]$ReaxID %in% comp)]; #print(comp);print(my_entezIDs);
                        if(length(my_entezIDs)>0){
                          my_gene<-which(rownames(gene_exp_mat) %in% my_entezIDs)
                          if(length(my_gene)>0){  
                            my_expression_data[comp]<-as.numeric(gene_exp_mat[my_gene,sample]) }}}
        
        if(length(my_expression_data)>0){
          ccc<-unlist(my_expression_data)
        }else{ccc<-default_value 
              #print(sprintf("NODE:%s has default value as %s",paste0(x,collapse=","),default_value))
        }
        #print("# e")
        #print(ccc)
        y[[m]] <- as.numeric(ccc)      
        
        # h)  
      }else if(length(x)>1 && length(if_complex)==0){
        
        my_expression_data<-list()
        for(comp in x){ my_entezIDs<-KEGG_met_path_node[[2]]$EntrezID[which(KEGG_met_path_node[[2]]$ReaxID %in% comp)]; #print(comp);#print(my_entezIDs);
                        if(length(my_entezIDs)>0){
                          my_gene<-which(rownames(gene_exp_mat) %in% my_entezIDs)
                          if(length(my_gene)>0){
                            my_expression_data[comp]<-as.numeric(gene_exp_mat[my_gene,sample]) }}}
        
        
        if(length(my_expression_data)>0){
          ccc<-max(unlist(my_expression_data))
        }else{ccc<-default_value
              #print(sprintf("NODE:%s has default value as %s",paste0(x,collapse=","),default_value))
        }
        
        #print("# h")
        #print(ccc)
        y[[m]] <- as.numeric(ccc)   
        
      }else{print("unkown Reaction status")}
      
    }
    
    node.vals<-do.call(what="rbind",y)
    colnames(node.vals)<-colnames(gene_exp_mat)[sample]
    if(nrow(non_reaction_enrties_frame)>=1){ # if there are compounds in the sif do else do not
      non_reaction_enrties_frame_sub_df<-as.data.frame(non_reaction_enrties_frame[,sample])
      colnames(non_reaction_enrties_frame_sub_df)<-colnames(gene_exp_mat)[sample]
      rownames(non_reaction_enrties_frame_sub_df)<-non_reaction_enrties
      node.vals<-rbind(node.vals,non_reaction_enrties_frame_sub_df)
    }
    rownames(nodes.vals.total)<-c(rownames(node.vals))
    nodes.vals.total<- cbind(nodes.vals.total,node.vals)
    
    #     print(head(nodes.vals.total))
    #     print(head(node.vals))
  }
  nodes.vals.total<-nodes.vals.total[,-1]
  
  nodes.vals.total <- t(apply(nodes.vals.total, 1, function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))))
  # I think default values are returned as NaN
  
  if(length(which(nodes.vals.total=="NaN"))>0){
    for(na in unique(which(nodes.vals.total=="NaN",arr.ind=T)[,1])){ nodes.vals.total[na,]<-default_value}
  }
  
  return(nodes.vals.total)
  
}


node_vals_for_flux<-function(graphobject,KEGG_met_path_node,flux_dist,default_value=0.5){
  
  # graphobject is output of moduleSIFandGraph
  # KEGG_met_path_node is output of KEGG_Met_path_parse
  # flux_dist is nrow=genes ncol=1m samples sammple nxm matrix
  
  if(!is.null(graphobject$reduced_SIF)){SIF<-graphobject$reduced_SIF}else{SIF<-graphobject$SIF}
  
  entries<-unique(c(SIF[,1],SIF[,3]))
  nodes.vals.total<-data.frame(name=seq(1:length(entries)))
  
  for(sample in 1:ncol(flux_dist)){
    
    idxleft<-grep("R",SIF[,1])
    SIFleft<-SIF[,1][idxleft]
    idxright<-grep("R",SIF[,3])
    SIFright<-SIF[,3][idxright]
    reactions<-unique(c(SIFleft,SIFright))
    
    ## default data for non Reaction node entries. like compound. metabolite, glycan ...etc   # !!!!!!!!!! it seems I do not need these part
    SIFleft_nonreax<-SIF[,1][-idxleft]
    SIFright_nonreax<-SIF[,3][-idxright]
    non_reaction_enrties<-unique(c(SIFleft_nonreax,SIFright_nonreax))
    non_reaction_enrties_frame<-as.data.frame(mat.or.vec(nr=length(non_reaction_enrties),nc=ncol(flux_dist)),stringsAsFactors=F)
    if(nrow(non_reaction_enrties_frame)>=1){non_reaction_enrties_frame[,]<-default_value} # sebep: hic compound yoksa bos frame donuyo ve burada hata verio. 
    #non_reaction_enrties_frame[,]<-runif(nrow(non_reaction_enrties_frame)*ncol(non_reaction_enrties_frame), 0,1) # as an example random generated data
    ###
    
    y<-sapply(reactions,function(x) {strsplit(x,split=",")[[1]]})
    if(all(reactions == y)){y<-as.list(y)}
    # fake data
    # y[1]<-list(c("R1+R2","R3","R4")) #c)
    # y[2]<-list(c("R5+R6")) # a
    # y[3]<-list(c("R7","R8")) #h
    # y[4]<-list(c("R9")) #e
    # y[5]<-list(c("R19+R20+R21","R22+R23")) #b
    # y[6]<-list(c("R24+R25+R26","R27+R28","R29","R30")) #d
    
    
    
    for(m in 1:length(y)){
      
      #x<-as.vector(y[[3]])
      x<-as.vector(y[[m]])
      if_complex<-grep("[+]R",x)
      
      # a)
      if(length(x)==1 && length(if_complex)>0 ){
        
        complex<-unlist(strsplit(x,split="[+]"))
        
        my_expression_data<-list()
        for(comp in complex){ my_entezIDs<-KeggBIGGIids$Vid[unique(which(KeggBIGGIids$KEGG.ID %in% comp)[1])]; print(my_entezIDs)
                              if(length(my_entezIDs)>0 && !is.na(my_entezIDs) && my_entezIDs!=""){
                                my_expression_data[paste0(comp,collapse="+")]<-as.numeric(flux_dist[which(rownames(flux_dist) %in% my_entezIDs),sample]) }}
        
        if(length(my_expression_data)>0){
          ccc<-min(unlist(my_expression_data))
        }else{ccc<-default_value; print(sprintf("NODE:%s has default value as %s",paste0(x,collapse=","),default_value))}
        print("# a")
        print(ccc)
        y[[m]] <- as.numeric(ccc)      
        
        # b)  
      }else if(length(x)>1 && length(if_complex)>1 && length(if_complex)==length(x)){
        
        x<-x[if_complex]
        aaa<-sapply(x,function(xs)  strsplit(xs,split="[+]"))
        
        my_expression_data<-list()
        for(comp in aaa){ my_entezIDs<-KeggBIGGIids$Vid[which(KeggBIGGIids$KEGG.ID %in% comp)[1]]; print(my_entezIDs);
                          if(length(my_entezIDs)>0 && !is.na(my_entezIDs) && my_entezIDs!=""){
                            my_expression_data[paste0(comp,collapse="+")]<-as.numeric(flux_dist[which(rownames(flux_dist) %in% my_entezIDs),sample]) }}
        
        if(length(my_expression_data)>0){
          bbb<-sapply(my_expression_data,function(xs) min(xs))    
          ccc<- max(bbb)
        }else{ccc<-default_value; print(sprintf("NODE:%s has default value as %s",paste0(x,collapse=","),default_value))}
        print("# b")
        print(ccc)
        y[[m]] <- as.numeric(ccc)      
        
        # c)
      }else if(length(x)>1 && length(if_complex)==1 && length(if_complex)<length(x)){ 
        
        x1<-x[if_complex]
        aaa1<-sapply(x1,function(xs)  strsplit(xs,split="[+]"))
        
        my_expression_data<-list()
        for(comp in aaa1){ my_entezIDs<-KeggBIGGIids$Vid[which(KeggBIGGIids$KEGG.ID %in% comp)[1]]; print(my_entezIDs);
                           if(length(my_entezIDs)>0 && !is.na(my_entezIDs) && my_entezIDs!=""){
                             my_expression_data[paste0(comp,collapse="+")]<-as.numeric(flux_dist[which(rownames(flux_dist) %in% my_entezIDs),sample]) }}
        
        if(length(my_expression_data)>0){
          bbb1<-sapply(my_expression_data,function(xs) min(xs))
        }else{bbb1<-list()}
        
        x2<-x[-if_complex]
        
        my_expression_data2<-list()
        for(comp in x2){ my_entezIDs<-KeggBIGGIids$Vid[which(KeggBIGGIids$KEGG.ID %in% comp)[1]]; print(my_entezIDs);
                         if(length(my_entezIDs)>0 && !is.na(my_entezIDs) && my_entezIDs!=""){
                           my_expression_data2[comp]<-as.numeric(flux_dist[which(rownames(flux_dist) %in% my_entezIDs),sample]) }}
        
        #for(comp in x2){ my_entezIDs<-KEGG_met_path_node$entrezIDs[which(KEGG_met_path_node$reaction  %in% comp)[1]];  my_expression_data2[comp]<-as.numeric(flux_dist[which(rownames(flux_dist) %in% my_entezIDs),1]) }
        
        
        
        if(length(my_expression_data2)>0 && length(bbb1)>0){
          my_expression_data2<-max(unlist(my_expression_data2))
          ccc<-max(my_expression_data2,bbb1)
        }else if(length(my_expression_data2)==0 && length(bbb1)>0){
          ccc<-max(bbb1)
        }else if(length(my_expression_data2)>0 && length(bbb1)==0){
          my_expression_data2<-max(unlist(my_expression_data2))
          ccc<-max(my_expression_data2)
        }else{ccc<-default_value; print(sprintf("NODE:%s has default value as %s",paste0(x,collapse=","),default_value))}
        
        print("# c")
        print(ccc)
        y[[m]] <- as.numeric(ccc)      
        
        # d)  
      }else if(length(x)>1 && length(if_complex)>1 && length(if_complex)<length(x)){
        
        x1<-x[if_complex]
        aaa1<-sapply(x1,function(xs)  strsplit(xs,split="[+]"))
        
        my_expression_data<-list()
        for(comp in aaa1){ my_entezIDs<-KeggBIGGIids$Vid[which(KeggBIGGIids$KEGG.ID %in% comp)[1]]; print(comp);print(my_entezIDs);
                           if(length(my_entezIDs)>0 && !is.na(my_entezIDs) && my_entezIDs!=""){
                             my_expression_data[paste0(comp,collapse="+")]<-as.numeric(flux_dist[which(rownames(flux_dist) %in% my_entezIDs),sample]) }}
        
        if(length(my_expression_data)>0){
          bbb1<-sapply(my_expression_data,function(xs) min(xs))
        }else{bbb1<-list()}
        
        x2<-x[-if_complex]
        
        my_expression_data2<-list()
        for(comp in x2){ my_entezIDs<-KeggBIGGIids$Vid[which(KeggBIGGIids$KEGG.ID %in% comp)[1]]; print(my_entezIDs);
                         if(length(my_entezIDs)>0 && !is.na(my_entezIDs) && my_entezIDs!=""){
                           my_expression_data2[comp]<-as.numeric(flux_dist[which(rownames(flux_dist) %in% my_entezIDs),sample]) }}
        
        
        
        if(length(my_expression_data2)>0 && length(bbb1)>0){
          ccc<-max(my_expression_data2,bbb1)
        }else if(length(my_expression_data2)==0 && length(bbb1)>0){
          ccc<-max(bbb1)
        }else if(length(my_expression_data2)>0 && length(bbb1)==0){
          ccc<-max(my_expression_data2)
        }else{ccc<-default_value; print(sprintf("NODE:%s has default value as %s",paste0(x,collapse=","),default_value))}
        
        print("# d")
        print(ccc)
        y[[m]] <- as.numeric(ccc)      
        
        # e)
      }else if(length(x)==1 && length(if_complex)==0){
        
        my_expression_data<-list()
        for(comp in x){ my_entezIDs<-KeggBIGGIids$Vid[which(KeggBIGGIids$KEGG.ID %in% comp)[1]]; print(comp);print(my_entezIDs);
                        if(length(my_entezIDs)>0 && !is.na(my_entezIDs) && my_entezIDs!=""){
                          my_expression_data[comp]<-as.numeric(flux_dist[which(rownames(flux_dist) %in% my_entezIDs),sample]) }}
        
        if(length(my_expression_data)>0){
          ccc<-unlist(my_expression_data)
        }else{ccc<-default_value; print(sprintf("NODE:%s has default value as %s",paste0(x,collapse=","),default_value))}
        print("# e")
        print(ccc)
        y[[m]] <- as.numeric(ccc)      
        
        # h)  
      }else if(length(x)>1 && length(if_complex)==0){
        
        my_expression_data<-list()
        for(comp in x){ my_entezIDs<-KeggBIGGIids$Vid[which(KeggBIGGIids$KEGG.ID %in% comp)[1]]; print(comp);print(my_entezIDs);
                        if(length(my_entezIDs)>0 && !is.na(my_entezIDs) && my_entezIDs!=""){
                          my_expression_data[comp]<-as.numeric(flux_dist[which(rownames(flux_dist) %in% my_entezIDs),sample]) }}
        
        
        if(length(my_expression_data)>0){
          ccc<-max(unlist(my_expression_data))
        }else{ccc<-default_value; print(sprintf("NODE:%s has default value as %s",paste0(x,collapse=","),default_value))}
        
        print("# h")
        print(ccc)
        y[[m]] <- as.numeric(ccc)   
        
      }else{print("unkown Reaction status")}
      
    }
    
    node.vals<-do.call(what="rbind",y)
    colnames(node.vals)<-colnames(flux_dist)[sample]
    if(nrow(non_reaction_enrties_frame)>=1){ # if there are compounds in the sif do else do not
      non_reaction_enrties_frame_sub_df<-as.data.frame(non_reaction_enrties_frame[,sample])
      colnames(non_reaction_enrties_frame_sub_df)<-colnames(flux_dist)[sample]
      rownames(non_reaction_enrties_frame_sub_df)<-non_reaction_enrties
      node.vals<-rbind(node.vals,non_reaction_enrties_frame_sub_df)
    }
    rownames(nodes.vals.total)<-c(rownames(node.vals))
    nodes.vals.total<- cbind(nodes.vals.total,node.vals)
    
    #     print(head(nodes.vals.total))
    #     print(head(node.vals))
  }
  nodes.vals.total<-nodes.vals.total[,-1]
  
  nodes.vals.total <- t(apply(nodes.vals.total, 1, function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))))
  # I think default values are returned as NaN
  
  if(length(which(nodes.vals.total=="NaN"))>0){
    for(na in unique(which(nodes.vals.total=="NaN",arr.ind=T)[,1])){ nodes.vals.total[na,]<-default_value}
  }
  
  return(nodes.vals.total)
  print(sample)
}


# SIF<-graphobject[[4]]
# initial_node<-graphobject[[2]]
# end_node<-graphobject[[3]]
# measured_metabolites=c("C00118","C00231")

reduce_SIF<-function(SIF,initial_node,end_node,measured_metabolites=NULL){
  
  keep_metabolites<-c()
  if(length(which(SIF[,1]==initial_node))>1){
    keep_metabolites<-c(keep_metabolites,initial_node)
  }
  if(length(which(SIF[,3]==end_node))>1){
    keep_metabolites<-c(keep_metabolites,end_node)
  }
  if(length(which(SIF[,1]==initial_node))==1 && length(which(SIF[,3]==end_node))==1){
    keep_metabolites<-c()}   
  
  #!!! # both is.null(keep_metabolites) only X null only y null seklinde  yap
  if(is.null(keep_metabolites)){
    reaction_left<-grep("R",SIF[,1])
    reaction_right<-grep("R",SIF[,3])
    SIFnew<-SIF[c(reaction_left,reaction_right),]
  }else{
    reaction_left<-grep("R",SIF[,1])
    reaction_right<-grep("R",SIF[,3])
    measured_left<-which(SIF[,1] %in% measured_metabolites)
    measured_right<-which(SIF[,3] %in% measured_metabolites)
    compound_left<-which(SIF[,1] %in% keep_metabolites)
    compound_right<-which(SIF[,3] %in% keep_metabolites)
    
    for_newSIF<-sort(unique(c(reaction_left,reaction_right)))
    for_newSIF<-for_newSIF[!for_newSIF %in% unique(c(compound_left,compound_right,measured_left,measured_right))]
    
    SIFnew<-SIF[unique(c(for_newSIF)),]
    SIFnew<-SIFnew[order(SIFnew[,1],decreasing=T),]
    SIFnew_later_add<-SIF[unique(c(compound_left,compound_right,measured_left,measured_right)),]
  }
  
  SIF_only_reaction<-mat.or.vec(nr=1,nc=3)
  
  for(rx in 1:length(grep("R",SIFnew[,1]))){
    
    left<-SIFnew[rx,1]
    if(length(which(measured_metabolites %in% SIFnew[rx,3]))>=1){
      right<-SIFnew[rx,3]
    }else{   
      right<-SIFnew[which(SIFnew[,1]==SIFnew[rx,3]),3]
    }
    print(left)
    print(right)
    
    #print(sprintf("%s to %s",left,right))
    if(length(left)==1 && length(right)==1){
      
      to_add<-mat.or.vec(nr=1,nc=3)
      to_add[1,1]<-left
      to_add[1,2]<-"activation"
      to_add[1,3]<-right
      print("x")
      print(to_add)
      
    }else if(length(left)>1 && length(right)==1){
      
      to_add<-mat.or.vec(nr=length(left),nc=3)
      to_add[1:length(left),1]<-left
      to_add[1:length(left),2]<-"activation"
      to_add[1:length(left),3]<-rep(right,times=length(left))
      print("y")
      print(to_add)
      
    }else if(length(left)==1 && length(right)>1){
      
      to_add<-mat.or.vec(nr=length(right),nc=3)
      to_add[1:length(right),1]<-rep(left,times=length(right))
      to_add[1:length(right),2]<-"activation"
      to_add[1:length(right),3]<-right
      print("z")
      print(to_add)
      
    }else if(length(left)==1 && length(right)==0){
      
      to_add<-mat.or.vec(nr=1,nc=3)
      to_add[1,1]<-left
      to_add[1,2]<-"activation"
      to_add[1,3]<-end_node
      print("w")
      print(to_add)
      
    }else{print("there is an error in this SIF check it")}
    
    SIF_only_reaction<-rbind(SIF_only_reaction,to_add)
    
  }
  SIF_only_reaction<-rbind(SIF_only_reaction,SIFnew_later_add)
  SIF_only_reaction<-unique(SIF_only_reaction[-1,])
  return(SIF_only_reaction)
  
}


# updated version of the reduce_SIF 
# measured_metabolites<-c(initial_node,end_node)

reduce_SIF2<-function(SIF,end_node=NULL,measured_metabolites=NULL, verbose=F){
  
  # I could not remember for what  
  #   keep_metabolites<-c()
  #   if(length(which(SIF[,1]==initial_node))>1){
  #     keep_metabolites<-c(keep_metabolites,initial_node)
  #   }
  #   if(length(which(SIF[,3]==end_node))>1){
  #     keep_metabolites<-c(keep_metabolites,end_node)
  #   }
  #   if(length(which(SIF[,1]==initial_node))==1 && length(which(SIF[,3]==end_node))==1){
  #     keep_metabolites<-c()}   
  #   
  #   
  #   #!!! # both is.null(keep_metabolites) only X null only y null seklinde  yap
  #   if(is.null(keep_metabolites)){
  #     reaction_left<-grep("R",SIF[,1])
  #     reaction_right<-grep("R",SIF[,3])
  #     SIFnew<-SIF[c(reaction_left,reaction_right),]
  #   }else{
  #      reaction_left<-grep("R",SIF[,1])
  #      reaction_right<-grep("R",SIF[,3])
  #      measured_left<-which(SIF[,1] %in% measured_metabolites)
  #      measured_right<-which(SIF[,3] %in% measured_metabolites)
  #     compound_left<-which(SIF[,1] %in% keep_metabolites)
  #     compound_right<-which(SIF[,3] %in% keep_metabolites)
  #     
  #     for_newSIF<-sort(unique(c(reaction_left,reaction_right)))
  #     for_newSIF<-for_newSIF[!for_newSIF %in% unique(c(compound_left,compound_right,measured_left,measured_right))]
  #     
  #     SIFnew<-SIF[unique(c(for_newSIF)),]
  #     SIFnew<-SIFnew[order(SIFnew[,1],decreasing=T),]
  #     SIFnew_later_add<-SIF[unique(c(compound_left,compound_right,measured_left,measured_right)),]
  #   }
  
  # this is same as below inside the if(is.null(keep_metabolites)){
  reaction_left<-grep("R",SIF[,1])
  reaction_right<-grep("R",SIF[,3])
  SIFnew<-SIF[c(reaction_left,reaction_right),]
  
  SIF_only_reaction<-mat.or.vec(nr=1,nc=3)
  
  for(rx in 1:length(grep("R",SIFnew[,1]))){
    
    left<-SIFnew[rx,1]
    if(length(which(measured_metabolites %in% SIFnew[rx,3]))>=1){
      right<-SIFnew[rx,3]
    }else{   
      right<-SIFnew[which(SIFnew[,1]==SIFnew[rx,3]),3]
    }
    #print(left)
    #print(right)
    
    #print(sprintf("%s to %s",left,right))
    if(length(left)==1 && length(right)==1){
      
      to_add<-mat.or.vec(nr=1,nc=3)
      to_add[1,1]<-left
      to_add[1,2]<-"activation"
      to_add[1,3]<-right
      if(verbose){ print("x"); print(to_add) }
      do_bind<-T
      
    }else if(length(left)>1 && length(right)==1){
      
      to_add<-mat.or.vec(nr=length(left),nc=3)
      to_add[1:length(left),1]<-left
      to_add[1:length(left),2]<-"activation"
      to_add[1:length(left),3]<-rep(right,times=length(left))
      if(verbose){ print("y"); print(to_add)}
      do_bind<-T
      
    }else if(length(left)==1 && length(right)>1){
      
      to_add<-mat.or.vec(nr=length(right),nc=3)
      to_add[1:length(right),1]<-rep(left,times=length(right))
      to_add[1:length(right),2]<-"activation"
      to_add[1:length(right),3]<-right
      if(verbose){print("z"); print(to_add)}
      do_bind<-T
      
    }else if(length(left)==1 && length(right)==0){
      
      
      #       to_add<-mat.or.vec(nr=1,nc=3)
      #       to_add[1,1]<-left
      #       to_add[1,2]<-"activation"
      #       to_add[1,3]<-end_node
      #       print("w")
      #       print(to_add)
      do_bind<-F
      next
      
    }else{  if(verbose){ print("there is an error in this SIF check it") }}
    
    if(do_bind){ SIF_only_reaction<-rbind(SIF_only_reaction,to_add)}
    
  }
  #  SIF_only_reaction<-rbind(SIF_only_reaction,SIFnew_later_add)
  SIF_only_reaction<-unique(SIF_only_reaction[-1,])
  
  return(SIF_only_reaction)
  
}



reaction_and_its_substrates<-function(SIF){
  
  reaction_right<-grep("R",SIF[,3])
  right_rxns_uniq<-unique(SIF[reaction_right,3])
  rxn_and_its_substrate_metabolite<-mat.or.vec(nr=length(right_rxns_uniq),nc=2)
  colnames(rxn_and_its_substrate_metabolite)<-c("reaction","substrate")  
  rxn_and_its_substrate_metabolite<-as.data.frame(rxn_and_its_substrate_metabolite)
  
  for(r_rxn in 1:length(right_rxns_uniq)){
    
    rxn_for_subs<-right_rxns_uniq[r_rxn]
    subs_of_rxn<-SIF[,1][which(SIF[,3] %in% rxn_for_subs)]
    rxn_and_its_substrate_metabolite$substrate[r_rxn]<- paste(subs_of_rxn,collapse=";")
    rxn_and_its_substrate_metabolite$reaction[r_rxn]<- rxn_for_subs
  }
  
  return(rxn_and_its_substrate_metabolite)
}



discard_compounds<-function(hsa_module_data){
  
  hsa_module_data_with_reduced <- lapply(hsa_module_data,function(x){ 
    # print(x$graphobject$module_name)
    gobj <- x$graphobject
    SIF_reduced<-reduce_SIF2(gobj$SIF,measured_metabolites=NULL) 
    rxn_and_substrate<-reaction_and_its_substrates(gobj$SIF)
    
    if(is.null(dim(SIF_reduced))){ 
      SIF_reduced<-as.data.frame(matrix(SIF_reduced,ncol=3),stringsAsFactors=F) 
    }else{ SIF_reduced<-as.data.frame(SIF_reduced,stringsAsFactors=F)  }                       
    
    g_reduced<-graph.data.frame(SIF_reduced[,c(1,3)],directed=T,vertices=unique(c(SIF_reduced[,1],SIF_reduced[,3])))
    g_reduced<-set.edge.attribute(g_reduced, "relation", value=rep(1,length(E(g_reduced)))) # indicates activation for compute.node.signal2
    
    x$graphobject[["reduced_SIF"]]<-SIF_reduced
    x$graphobject[["reduced_igraph"]]<-g_reduced
    x$graphobject[["reduced_initial_node"]]<- SIF_reduced[1,1] 
    #!!! here the end node selected arbitrary, to skip actnode==endnode from path.value function
    if(SIF_reduced[nrow(SIF_reduced),3]!=SIF_reduced[1,1]){ # checks if it is a cycle or not
      x$graphobject[["reduced_end_node"]]<- gobj$SIF[,1][which(gobj$SIF[,3] %in% gobj$end_node)]
    }else{
      x$graphobject[["reduced_end_node"]]<- SIF_reduced[nrow(SIF_reduced)-1,3]
    }
    
    x$graphobject[["rxn_and_substrate"]]<-rxn_and_substrate
    
    return(x)
  }
  )
  # message("model is reduced (metabolites are discarded from model)")
  return(hsa_module_data_with_reduced)  
  
}



sif2graph<-function(sif){
  
  #if the input is a character it shoud be the name ot the sif file
  #otherwise a matrix in the sif format
  if (is.vector(sif) && (typeof(sif) == "character")){
    sif = read.table(sif) 
  }
  
  # build the unique vertices from the column 1 and 3 of the SIF file
  vertices = unique(c(as.character(sif[,1]), as.character(sif[,3])))
  # some aliases
  v1 = sif[,1]
  v2 = sif[,3]
  edges = as.numeric(sif[,2])
  
  l = length(vertices) - 1
  g <- new("graphNEL", nodes=vertices, edgemode="directed")
  #weights = rep(1, l)
  weights = edges
  for (i in 1:length(v1)){
    g <- addEdge(as.character(v1[i]), as.character(v2[i]), g, weights[i])
  }
  return(g)
  
}

getmodulePahtID<-function(graphobject){
  k<-1
  pathID<-graphobject$modulePATH[k,1]
  toyKGML<-paste0("http://www.kegg.jp/kegg-bin/download?entry=",pathID,"&format=kgml")
  while(getURL(toyKGML)=="\n"){
    pathID<-graphobject$modulePATH[k,1]
    if(is.na(pathID))
      break
    toyKGML<-paste0("http://www.kegg.jp/kegg-bin/download?entry=",pathID,"&format=kgml")
    k<-k+1
  }
  return(pathID)
}


curate_module<-function(modulename,inits,ends,module_file_path,saveSIFpath,curation_list,save_modulename){
  
  # curation_list<-"/home/cankut/Desktop/Metabolic pathyway/KEGG_module/Restful_based/manual_curation_list.csv"
  # module_file_path<-"/home/cankut/Desktop/Metabolic pathyway/KEGG_module/Restful_based/parsed_txt/"
  # saveSIFpath<-"/home/cankut/Desktop/Metabolic pathyway/KEGG_module/Restful_based/SIF/"
  # modulename<-"M00044"
  # inits<-c("C00082")
  # ends<-c("C00164","C00122")
  # curated_name<-"M00044_1"
  
  diff_vec<-setdiff(curation_list$new_name[curation_list$new_name!=""],unique(curation_list$new_name[curation_list$new_name!=""]))
  if(length(diff_vec)==0){
    
    # trim leading and flanking blank spaces
    inits<-strsplit(inits,"[+]")[[1]]
    inits<-gsub("^\\s+|\\s+$", "", inits)
    ends<-strsplit(ends,"[+]")[[1]]
    ends<-gsub("^\\s+|\\s+$", "", ends)
    
    updated_modules<-list()
    
    graphobject<-moduleSIFandGraph(module_file_path=module_file_path,
                                   modulename=modulename,  #"M00014"
                                   saveSIFpath=saveSIFpath,save_modulename=save_modulename)
    
    #plot(graphobject[[1]],vertex.size=1,edge.arrow.size=1)
    
    if(is.null(ends)){ends<-graphobject$end_node}
    
    SIF<-graphobject$SIF
    
    idx_init_maxtrix<-which(SIF[,1] %in% inits) 
    multi_initial_nodels<-c()
    
    if(length(idx_init_maxtrix)>1){
      
      initial_matrix<-mat.or.vec(nr=length(inits),nc=3)
      for(i in 1:length(inits)){
        initial_matrix[i,1]<-inits[i]
        initial_matrix[i,2]<-"activation"
        initial_matrix[i,3]<- SIF[idx_init_maxtrix[i],3]
      }
      SIF<-SIF[-idx_init_maxtrix,]
      multi_initial_nodels<-inits
      
      #       for(i in 1:length(inits)){
      #         initial_matrix[i,1]<-pseudo_node
      #         initial_matrix[i,2]<-"activation"
      #         initial_matrix[i,3]<- inits[i]
      #       }
    }else{
      initial_matrix<-mat.or.vec(nr=length(inits),nc=3)
      for(i in 1:length(inits)){
        initial_matrix[i,1]<-paste(inits,collapse="+")
        initial_matrix[i,2]<-"activation"
        initial_matrix[i,3]<- SIF[idx_init_maxtrix[i],3]
      }
      SIF<-SIF[-idx_init_maxtrix,]
    }
    
    idx_end_maxtrix<-which(SIF[,3] %in% ends) 
    new_SIF<-rbind(initial_matrix,SIF[-idx_end_maxtrix,],SIF[idx_end_maxtrix,])
    new_SIF <- as.data.frame(new_SIF,stringsAsFactors=F)
    g<-graph.data.frame(new_SIF[,c(1,3)],directed=T,vertices=unique(c(new_SIF[,1],new_SIF[,3])))
    g<-set.edge.attribute(g, "relation", value=rep(1,length(E(g))))
    plot(g,vertex.size=1,edge.arrow.size=1)
    
    write.table(new_SIF,paste0(saveSIFpath,"/",save_modulename,"SIF.txt"),sep="\t",row.names=F,quote=F,col.names=F)
    
    if(!is.null(multi_initial_nodels)){
      updated_modules<-list("igraph"=g,"initial_node"=multi_initial_nodels,"end_node"=new_SIF[nrow(new_SIF),3],"SIF"=as.matrix(new_SIF),"modulePATH"=graphobject$modulePATH,"module_name"=save_modulename)
      
    }else{
      updated_modules<-list("igraph"=g,"initial_node"=new_SIF[1,1],"end_node"=new_SIF[nrow(new_SIF),3],"SIF"=as.matrix(new_SIF),"modulePATH"=graphobject$modulePATH,"module_name"=save_modulename)
    }
    return(updated_modules)
  }else{message("There are duplicated new_module names. Check your curation_list file! ")}
}



KEGG_module_and_path_func<-function(hsa_module_data){
  KEGG_module_and_path=data.frame(Module="",KGML="",Paths="",module_link="",start_molecule="",end_molecule="",stringsAsFactors=F)
  for(m in 1:length(hsa_module_data)){
    graphobject<-hsa_module_data[[m]]$graphobject
    pathID<-getmodulePahtID(graphobject)
    KEGG_module_and_path[m,"Module"] <- graphobject$module_name
    if(is.na(pathID)){pathID<-"KGML is not exist"}else{pathID<-pathID}
    KEGG_module_and_path[m,"KGML"] <- pathID
    KEGG_module_and_path[m,"Paths"] <- paste(graphobject$modulePATH$V1,collapse=";")
    KEGG_module_and_path[m,"module_link"] <- paste0("http://www.genome.jp/kegg-bin/show_module?map=",graphobject$module_name,"&org=hsa&select_org=hsa")
    KEGG_module_and_path[m,"start_molecule"] <- paste(graphobject$initial_node,collapse="+")
    KEGG_module_and_path[m,"end_molecule"] <- graphobject$end_node
  }
  return(KEGG_module_and_path)
}


# g<-graph.data.frame(SIF_only_reaction[,c(1,3)],directed=T,vertices=unique(c(SIF_only_reaction[,1],SIF_only_reaction[,3])))
# plot(g)

###### biomart data

### hsa
# library("biomaRt")
# #listMarts(host="useast.ensembl.org")
# ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="useast.ensembl.org",dataset="hsapiens_gene_ensembl")
# #head(ensembl@attributes,100)
# #grep("symbol",ensembl@filters$name)
# save(ensembl,file="/home/cankut/Desktop/Metabolic pathyway/KEGG_module/biomaRt.Rdata")

### rat
# library("biomaRt")
# listMarts(host="jul2015.archive.ensembl.org")
# ensembl_rat = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="rnorvegicus_gene_ensembl", host = "jul2015.archive.ensembl.org")
# save(ensembl_rat,file="/home/cankut/Desktop/Metabolic pathyway/KEGG_module/retreat_MLPM/biomaRt_rat.Rdata")


## buna cok gerek duymayabilirim cunku default value nodes.valdede ekleniyo
add.missing.genes<-function(hsa_module_data,rxn_vals_mat,default_value=0.5){
  
  genes_and_compound_per_module<-sapply(hsa_module_data, function(x) x$KEGG_met_path_node$Reax2Entrez$EntrezID)
  genes_and_compound_per_module<-sort(unique(unlist(genes_and_compound_per_module)))
  genes_and_compound_per_module<-genes_and_compound_per_module[!genes_and_compound_per_module=="NotConverted"]
  
  missing_info<- genes_and_compound_per_module[!genes_and_compound_per_module %in% rownames(rxn_vals_mat)]
  missing_info_mat<-mat.or.vec(nr=length(missing_info),nc=ncol(rxn_vals_mat))
  colnames(missing_info_mat)<-colnames(rxn_vals_mat)
  rownames(missing_info_mat)<-missing_info
  missing_info_mat[,]<-0.5  # default vaule for missing genes
  rxn_vals_mat<-rbind(rxn_vals_mat,missing_info_mat)
  
  return(rxn_vals_mat)
  
}


module_stats_table<-function(hsa_module_data="/home/cankut/Desktop/Metabolic pathyway/KEGG_module/Restful_based/hsa_module_data.RData",
                             KEGG_compound_info="/home/cankut/Desktop/Metabolic pathyway/KEGG_module/KEGG_compound/KEGG_compound_info.RData",
                             output="/home/cankut/Desktop/Metabolic pathyway/KEGG_module/Restful_based/KEGG_module_and_path_table.txt"){
  
  load(hsa_module_data)
  load(KEGG_compound_info)
  KEGG_module_and_path<-KEGG_module_and_path_func(hsa_module_data)
  
  for(c in 1:length(KEGG_module_and_path$start_molecule)){ 
    compound<-KEGG_module_and_path$start_molecule[c]
    comps<-strsplit(compound,"[+]")[[1]]  
    if(length(comps)>1){
      my_complex_comps<-c()
      for(c2 in 1:length(comps)){
        my_complex_comps<-c(my_complex_comps,compound_info[[comps[c2]]][1])
      }
      KEGG_module_and_path$start_molecule[c]<- paste0(compound,": ", paste(my_complex_comps,collapse=" + "))
    }else{ 
      KEGG_module_and_path$start_molecule[c] <- paste0(compound,": ",compound_info[[comps]][1])
    }
  }
  
  for(c in 1:length(KEGG_module_and_path$end_molecule)){ 
    compound<-KEGG_module_and_path$end_molecule[c]
    comps<-strsplit(compound,"[+]")[[1]]  
    if(length(comps)>1){
      my_complex_comps<-c()
      for(c2 in 1:length(comps)){
        my_complex_comps<-c(my_complex_comps,compound_info[[comps[c2]]][1])
      }
      KEGG_module_and_path$end_molecule[c]<- paste0(compound,": ", paste(my_complex_comps,collapse=" + "))
    }else{ 
      KEGG_module_and_path$end_molecule[c] <- paste0(compound,": ",compound_info[[comps]][1])
    }
  }
  
  write.table(KEGG_module_and_path,output,sep="\t",row.names=F,col.names=T,quote=F)
  cat("table saved to",output)
}



methyways<-function(hsa_module_data,rxn_vals_mat,expbased=T,fluxbased=F,flux_dist, verbose=T, default_value=0.5, metabolitematrix=NULL){
  
  res<-lapply(hsa_module_data,function(x){ 
    
    
    graphobject<-x$graphobject
    KEGG_met_path_node<-x$KEGG_met_path_node
    
    if(expbased || !fluxbased ){
      nodes.vals<-node_vals(graphobject,KEGG_met_path_node=KEGG_met_path_node,rxn_vals_mat,default_value)
    }else{
      nodes.vals<-node_vals_for_flux(graphobject,KEGG_met_path_node,flux_dist)
    }
    
    if(!is.null(graphobject$reduced_igraph))
    {
      subgraph<-graphobject$reduced_igraph
      ininodes<-graphobject$reduced_initial_node
      endnode<-graphobject$reduced_end_node
    }else{
      subgraph<-graphobject$igraph
      ininodes<-graphobject$initial_node
      endnode<-graphobject$end_node 
    }
    
    if(!is.null(metabolitematrix)){
      nodes.vals <- rbind(nodes.vals, metabolitematrix)
    }
    
    path.val<-path.value(nodes.vals=nodes.vals, subgraph ,ininodes,endnode)
    if(verbose){ cat(graphobject$module_name,"...DONE \n") }
    return(list(path.val))
    
  })
  
  return(res)
}

get.All.module.genes <- function(hsa_module_data){
  
  # hsa_module_data: list contains module and their features. Such as hsa_module_data_April2016.RData
  
  all_module_genes <-sapply(hsa_module_data, function(x) {
    module_elements<-unique(as.vector(x$graphobject$SIF[,c(1,3)]))
    rxns<-as.vector(unlist(sapply(module_elements[grep("R",module_elements)],function(x) {strsplit(x,split = ",")})))
    rxns<-as.vector(unlist(sapply(rxns,function(x) {strsplit(x,split = "[+]")})))
    #module_genes<-x$KEGG_met_path_node$Reax2Entrez$EntrezID[which(x$KEGG_met_path_node$Reax2Entrez$ReaxID %in% rxns)]
    
    Reax2Entrez <- x$KEGG_met_path_node$Reax2Entrez
    to_split <- grep(" R",Reax2Entrez$ReaxID)
    
    if(length(to_split)>0){
      for(ts in to_split){
        idx_push <- nrow(Reax2Entrez) + 1
        rxn_to_split <- Reax2Entrez$ReaxID[ts]
        item_to_push <- strsplit(rxn_to_split," ")[[1]]
        Reax2Entrez[c(idx_push:(idx_push+length(item_to_push)-1)),1] <- rep(Reax2Entrez$EntrezIDs[ts],length(item_to_push)) 
        Reax2Entrez[c(idx_push:(idx_push+length(item_to_push)-1)),2] <- item_to_push
      }
      Reax2Entrez <- Reax2Entrez[-to_split,]
    }
    module_genes<-Reax2Entrez$EntrezID[which(Reax2Entrez$ReaxID %in% rxns)]
    
    return(module_genes)
  }
  )
  
  all_module_genes<-unique(unlist(all_module_genes))
  
  all_module_genes_vec<-c()
  for(amg in all_module_genes){  all_module_genes_vec<-c(all_module_genes_vec, strsplit(amg,";")[[1]])}
  all_module_genes_vec<-unique(all_module_genes_vec)
  
  return(all_module_genes_vec)
  
}


get.All.module.rxns <- function(hsa_module_data){
  
  # hsa_module_data: list contains module and their features. Such as hsa_module_data_April2016.RData
  
  all_module_rxn_vec <- sapply(hsa_module_data, function(x) {
    module_elements<-unique(as.vector(x$graphobject$SIF[,c(1,3)]))
    rxns<-as.vector(unlist(sapply(module_elements[grep("R",module_elements)],function(x) {strsplit(x,split = ",")})))
    rxns<-as.vector(unlist(sapply(rxns,function(x) {strsplit(x,split = "[+]")})))
    return(rxns)
  }
  )
  
  all_module_rxn_vec<-unique(unlist(all_module_rxn_vec))
  return(all_module_rxn_vec)
  
}


get.rxn.gene.matrix <- function(hsa_module_data){
  
  # hsa_module_data: list contains module and their features. Such as hsa_module_data_April2016.RData
  
  rxn_gene_mat<-mat.or.vec(nc = 2, nr = 1)
  colnames(rxn_gene_mat) <- c("EntrezIDs","ReaxID")
  
  for(mat in 1:length(hsa_module_data)){
    rxn_gene_mat <- rbind(rxn_gene_mat,hsa_module_data[[mat]]$KEGG_met_path_node$Reax2Entrez)
  }
  
  rxn_gene_mat<-rxn_gene_mat[-1,]
  rxn_gene_mat<-unique(rxn_gene_mat)
  
  Reax2Entrez <- rxn_gene_mat
  to_split <- grep(" R",Reax2Entrez$ReaxID)
  
  if(length(to_split)>0){
    for(ts in to_split){
      idx_push <- nrow(Reax2Entrez) + 1
      rxn_to_split <- Reax2Entrez$ReaxID[ts]
      item_to_push <- strsplit(rxn_to_split," ")[[1]]
      Reax2Entrez[c(idx_push:(idx_push+length(item_to_push)-1)),1] <- rep(Reax2Entrez$EntrezIDs[ts],length(item_to_push)) 
      Reax2Entrez[c(idx_push:(idx_push+length(item_to_push)-1)),2] <- item_to_push
    }
    Reax2Entrez <- Reax2Entrez[-to_split,]
  }
  
  return(Reax2Entrez)  
}



get.gene.exp.of.module.genes <- function(combat.vals, all_module_genes_vec, min.exp=0.001, onesample=F){
  
  # combat.vals is normalized and batch corrected RNASeq (gene expression) data 
  # all_module_genes_vec: output of "get.All.module.genes"
  # min.exp lets us propagate flux/signal. This is in case RXN = 1 gene = 0.001
  
  idx_sellect_genes <- match(all_module_genes_vec,rownames(combat.vals))
  idx_NA_genes <- which(is.na(idx_sellect_genes))
  if(length(idx_NA_genes)>0){
    idx_sellect_genes <- idx_sellect_genes[-idx_NA_genes]
  }
  
  gene_exp_mat<-as.matrix(combat.vals[idx_sellect_genes,])
  
  # combat returns negative values which can be replaced by 0
  gene_exp_mat[gene_exp_mat<0]<-0
  if(onesample){
    gene_exp_mat <- (gene_exp_mat-min(gene_exp_mat,na.rm=T))/(max(gene_exp_mat,na.rm=T)-min(gene_exp_mat,na.rm=T)) 
  }
  else{
    gene_exp_mat <- t(apply(gene_exp_mat, 1, function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))))   
  }
  gene_exp_mat[is.na(gene_exp_mat)]<-0
  gene_exp_mat[gene_exp_mat==0] <- min.exp
  return(gene_exp_mat)
}




get.RXNvals.from.genes <- function(all_module_rxn_vec,gene_exp_mat,rxn_gene_mat){
  
  # all_module_rxn_vec: output of "get.All.module.rxns"
  # gene_exp_mat: output of "get.gene.exp.of.module.genes"
  
  rxn_vals_mat <- as.matrix(mat.or.vec(nr = length(all_module_rxn_vec) ,nc = ncol(gene_exp_mat)))
  colnames(rxn_vals_mat) <- colnames(gene_exp_mat)
  rownames(rxn_vals_mat) <- all_module_rxn_vec
  
  for(r in 1:length(all_module_rxn_vec)){
    
    # use grep to deal with M00015: R01251, hsa00330: R01251 R01248
    #r_genes_1 <- rxn_gene_mat$EntrezIDs[which(rxn_gene_mat$ReaxID %in% all_module_rxn_vec[r])]
    r_genes_1 <- rxn_gene_mat$EntrezIDs[grep(all_module_rxn_vec[r],rxn_gene_mat$ReaxID)]
    
    # I add this because there are some module reactions which are without gene info e.g. R02189 in M00001 catalyzed by EC 2.7.1.63 is not exist  for homosapiens
    if(length(r_genes_1)>0){
      
      r_genes_2 <- unique(unlist(sapply(r_genes_1,strsplit,";")))
      
      r_genes_3 <- r_genes_2[r_genes_2 %in% rownames(gene_exp_mat)]
      
      r_genes_3_vals <- gene_exp_mat[which(rownames(gene_exp_mat) %in% r_genes_3),]
      
      
      if(is.vector(r_genes_3_vals)){
        r_genes_3_vals <- r_genes_3_vals
      }else if(is.matrix(r_genes_3_vals)){
        #r_genes_3_vals <- quantile(r_genes_3_vals,0.9)
        r_genes_3_vals <- t(apply(r_genes_3_vals,2, function(x){quantile(x,0.9)}))
      }else{ r_genes_3_vals<-NA }
      
    }else{
      r_genes_3_vals<-NA
    }
    
    if(1==0){ 
      if(any(is.na(r_genes_3_vals))){cat("Rxn: ",all_module_rxn_vec[r]," has NA \n")}
    }
    rxn_vals_mat[r,] <- r_genes_3_vals
    
    if((r %% 50)==0){
      cat(r,"...DONE\n")
    }
  }
  
  rxns_missing_val <- which(apply(rxn_vals_mat,1,function(x) all(is.na(x))))
  # there is no gene expression info to calculate these reactions
  cat(paste0("number reactions with missing value: ", length(rxns_missing_val)))
  
  rxn_vals_mat <- rxn_vals_mat[-rxns_missing_val,]
  return(rxn_vals_mat)
}





# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

###################### taken from prettyways repos. on 29 January 2016

# PCA 
do.pca.2 <- function(data){
  fit <- prcomp(t(data))
  fit$scores <- fit$x
  return(fit)
}


do.pca <- function(data,cor=F){
  fit <- princomp(t(data), cor=cor)
  fit$var <- fit$sdev^2
  fit$explain_var <- fit$var/sum(fit$var)
  fit$acum_explain_var <- cumsum(fit$explain_var)
  return(fit)
}


plot.pca.variance <- function(fit,thresh=0,acum=F){
  if(acum==F){
    barplot(fit$explain_var[ fit$explain_var>thresh],ylab="explain variance",xlab="",las=2,cex.names=0.5,ylim=c(0,1))
  } else {
    barplot(fit$acum_explain_var[ fit$acum_explain_var<(1-thresh)],ylab="acum explain variance",xlab="",las=2,cex.names=0.5,ylim=c(0,1))
  }
}

plot.pca <- function(fit,cp1=1,cp2=2,colors=NULL,legend=NULL,legend_colors=NULL,cex=0.5,pch=20,text=F){
  if(is.null(colors)) colors <- "black"
  cpv1 <- fit$scores[,cp1]
  cpv2 <- fit$scores[,cp2]
  if(text==F){
    plot(cpv1,cpv2,xlab=paste("pc",cp1),ylab=paste("pc",cp2),col=colors,pch=pch,cex=cex)
  } else {
    plot(cpv1,cpv2,xlab=paste("pc",cp1),ylab=paste("pc",cp2),type="n")
    text(cpv1,cpv2,labels = rownames(fit$scores),col=colors,pch=pch,cex=cex)
  }
  if(!is.null(legend)){
    legend("left",legend=legend,col=legend_colors,pch=pch,xpd=T,cex=cex,border=0)
  }
}

plot.multiple.pca <- function(fit,comps=1:3,colors,plot.variance=F,legend=NULL,legend_colors=NULL,cex=0.9,pch=20,text=F){
  
  if(ncol(fit$scores)<=2){
    
    warning(paste0("There are not enough Principal Components. Argument [comp] is set as 1:",ncol(fit$scores)))
    comps <- 1:ncol(fit$scores)
  }
  
  combs <- combn(comps,2)
  ncombs <- ncol(combs)
  nn <- ncombs
  if(!is.null(legend)) nn <- nn + 1
  if(plot.variance==T) nn <- nn + 1  
  
  nr <- floor(sqrt(nn))
  nc <- ceiling((nn)/nr)
  par(mfrow=c(nr,nc))
  for(i in 1:(ncombs)){
    plot.pca(fit,cp1=combs[1,i],cp2=combs[2,i],colors=colors,cex=cex,pch=pch,text=text)
  }
  if(!is.null(legend)){
    plot(1,type="n",axes=F,xlab="",ylab="")
    legend("center",legend=legend,col=legend_colors,pch=15,lwd=2,xpd=T,cex=1.5,border=NA,pt.cex=1.2)
  }
  if(plot.variance==T){
    plot.pca.variance(fit,acum=T,thresh=0.1)
  }
  par(mfrow=c(1,1))
}


plot.heatmap <- function(path.vals,sample_type,colors=c("blue","gray","red"),sample.clust=T,variable.clust=F,labRow=NULL, labCol=NULL, sample_colors=NULL){  
  if(sample.clust==F){
    colv <- NA
  } else {
    colv <- T
  }
  if(variable.clust==F){
    rowv <- NA
  } else {
    rowv <- T
  }
  if(is.null(labRow)){
    labRow <- rownames(path.vals)
  }
  if(is.null(labCol)){
    labCol <- colnames(path.vals)
  }
  if(is.null(sample_colors)){
    sample_colors <- rainbow(length(unique(sample_type)))
    names(sample_colors) <- unique(sample_type)
  }
  path.vals <- t(apply(path.vals,1,function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))))
  #heatmap(path.vals,margins = par("mar")[1:2],labRow = labRow, labCol = labCol, scale="none",Rowv=rowv,Colv=colv,ColSideColors = sample_colors[sample_type],col=colorRampPalette(colors)(256))  
  heatmap(path.vals,margins = c(8,8),labRow = labRow, labCol = labCol, scale="none",Rowv=rowv,Colv=colv,ColSideColors = sample_colors[sample_type],col=colorRampPalette(colors)(256))  
}



do.wilcox <- function(sel.vals, group.value, g1, g2, paired=F, verbose=T ){
  
  if(all(colnames(sel.vals)!=group.value[,1])){
    group.value <- group.value[match(colnames(sel.vals),group.value[,1]),]
  }
  
  g1_indexes <- which(group.value[,2]==g1)
  g2_indexes <- which(group.value[,2]==g2)
  
  stat.vals <- calculate.wilcox.test(sel.vals,g2_indexes,g1_indexes,paired=paired)
  
  return(stat.vals)
  
}


calculate.wilcox.test <- function( data, control, disease, paired){
  if(paired == TRUE){
    dat <- apply(data, 1, wilcoxsign.test.fun, control, disease)
    testData <- do.call("rbind",dat)
    fdrData <- p.adjust(testData[,1], method = "fdr")
    data2 <- data.frame(testData, fdrData, stringsAsFactors=F)
  }else{
    testData <- do.call("rbind",apply(data, 1, wilcox.test.fun, control, disease, paired))
    fdrData <- p.adjust(testData[,1], method = "fdr")
    # Standardize statistic
    m <- length(control)*length(disease)/2
    sigma <- sqrt(length(control)*length(disease)*(length(control)+length(disease)+1)/12)
    z <- (testData[,3]-m)/sigma
    data2 <- data.frame(testData[,1:2], z, fdrData,stringsAsFactors=F)
    data2[which(data2$pvalue==1 & data2$class == "0" & data2$fdrData == 1), "z"] <- 0
  }
  colnames(data2) <- c("p.value", "UP/DOWN", "statistic", "adj.p.value")
  data2[data2$statistic>0,"UP/DOWN"] <- "UP"
  data2[data2$statistic<0,"UP/DOWN"] <- "DOWN"
  data2[data2$statistic==0,"UP/DOWN"] <- "No Change"
  data2 <- data2[,c(3,1,4,2)]
  return(data2)
}


wilcoxsign.test.fun <- function(x, control, disease){
  r <- try(wilcoxsign_test( as.numeric(x[disease]) ~ as.numeric(x[control])))
  if (class(r) == "try-error"){
    pvalue <- 1
    class <- "0"
    stat <- 0
  } else{
    pvalue <- pvalue(r)
    stat <- statistic(r, "standardized")[1]
    if(stat > 0){
      class <- "UP"
    }else if (stat < 0){
      class <- "DOWN"
    }else{
      class <- "0"
    }
  }
  result <- data.frame(pvalue=pvalue, class=class, stat=stat, stringsAsFactors=F)
  return(result)
}


wilcox.test.fun <- function(x, control, disease, paired){
  r <- try(wilcox.test(x=as.numeric(x[disease]), y=as.numeric(x[control]) , conf.int=TRUE, alternative="two.sided", paired=paired))
  result <- wilcox.data.frame(r)
  return(result)
}

wilcox.data.frame <- function(wilcox){
  if (class(wilcox) == "try-error"){
    pvalue <- 1
    class <- "0"
    stat <- 0
  } else{
    pvalue <- wilcox$p.value
    esti <- wilcox$estimate
    stat <- wilcox$statistic[[1]]
    if (esti < 0){
      class <- "DOWN" ## regarding DISEASE
    } else if (esti > 0){
      class <- "UP" ## regarding DISEASE
    } else if (esti == 0){
      if (wilcox$conf.int[1] == 0){
        class <- "UP"
      } else if (wilcox$conf.int[2] == 0){
        class <- "DOWN"
      } else{
        class <- 0
      }
    }
  }
  result <- data.frame(pvalue, class, stat,stringsAsFactors=F)
  return(result)
}


do.cor <- function(sel.vals, design,adjust=T,methodcor="spearman"){
  
  data <- sel.vals[,design[,1]]
  
  testData <- do.call("rbind",apply(data, 1, function(x) {cor.test.fun(x,values=as.numeric(design[,2]),method=methodcor)}))
  if(adjust==T){
    fdrData <- p.adjust(testData[,1], method = "fdr")
  } else {
    fdrData <- testData[,1]
  }
  data2 <- data.frame(testData[,1:3], fdrData,stringsAsFactors=F)
  
  colnames(data2) <- c("p.value", "UP/DOWN", "correlation", "adj.p.value")
  data2[data2$statistic>0,"UP/DOWN"] <- "UP"
  data2[data2$statistic<0,"UP/DOWN"] <- "DOWN"
  return(data2)
}

cor.test.fun <- function(x, values,method){
  r <- try(cor.test(x=as.numeric(x), y=as.numeric(values), method = methodcor))
  result <- cor.data.frame(r)
  return(result)
}

cor.data.frame <- function(wilcox){
  if (class(wilcox) == "try-error" | is.na( wilcox$p.value )){
    pvalue <- 1
    class <- "0"
    esti <- 0
  } else{
    pvalue <- wilcox$p.value
    esti <- wilcox$estimate[1]
    if (esti < 0){
      class <- "DOWN" ## regarding DISEASE
    } else if (esti > 0){
      class <- "UP" ## regarding DISEASE
    } else if (esti == 0){
      if (wilcox$conf.int[1] == 0){
        class <- "UP"
      } else if (wilcox$conf.int[2] == 0){
        class <- "DOWN"
      } else{
        class <- 0
      }
    }
  }
  result <- data.frame(pvalue, class, esti, stringsAsFactors=F)
  return(result)
}



compute.difexp <- function(vals, control.string, disease.string, experimental.design){
  control <- experimental.design[which(experimental.design[,2]==control.string),1]
  disease <- experimental.design[which(experimental.design[,2]==disease.string),1]
  sgroup <- experimental.design[experimental.design[,2] == control.string | experimental.design[,2]== disease.string,2]
  names(sgroup) <- experimental.design[experimental.design[,2] == control.string | experimental.design[,2]== disease.string,1]
  design <- cbind(sgroup==control.string, sgroup==disease.string)+0
  colnames(design) <- c("grupo1", "grupo2")
  
  fit <- lmFit(vals[,c(control, disease)], design[c(control, disease),])
  cont.matrix <- makeContrasts(grupo1-grupo2, levels=design[c(control, disease),])
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  padj <- p.adjust(fit2$"p.value")
  result <- data.frame(statistic=as.numeric(fit2$t), p.value=as.numeric(fit2$p.value), adj.p.value=padj, laterality=as.factor(fit2$t>0))
  rownames(result) <- rownames(fit2)
  
  return(result)
}



get.colors.from.pval <- function(updown, pvals, up.col="#da1f1f", down.col="#1f9cda", no.col="white", both.col="#959595",conf=0.05){
  colors <- sapply(1:length(updown), function(i){
    if(!is.na(pvals[i]) && pvals[i] <= conf){
      trans <- (1-18*pvals[i])
      if(is.na(updown[i])){
        return(no.col)
      }else if(updown[i] == "up"){
        cc <- colorRamp(c(no.col, up.col))(trans)/255
        return(rgb(cc[1], cc[2], cc[3]))
      }else if(updown[i] == "down" ){
        cc <- colorRamp(c(no.col, down.col))(trans)/255
        return(rgb(cc[1], cc[2], cc[3]))
      }else if(updown[i] == "none" ){
        cc <- colorRamp(c(no.col, both.col))(trans)/255
        return(rgb(cc[1], cc[2], cc[3]))
      }
    }else{
      return(no.col)
    }
  })
  return(colors)
}


node.color.per.differential.expression <- function(results_module_node_vals_matrix, control.string, disease.string, experimental.design){ 
  
  difexp.nod <- compute.difexp(vals = results_module_node_vals_matrix, control.string, disease.string, experimental.design=cancer_types)
  updown <- sapply(difexp.nod$statistic, function(x){if(x <0){"down"}else if(x > 0){"up"}else{"none"}})
  col <- get.colors.from.pval(updown, difexp.nod$p.value)
  names(col) <- rownames(results_module_node_vals_matrix)
  # Add function colors
  # toadd <- V(fpgs[[path]]$graph)$name[!V(fpgs[[path]]$graph)$name %in% rownames(difexp.nod)]
  # coltoadd <- rep("white", length(toadd))
  # names(coltoadd) <- toadd
  # col <- c(col, coltoadd)
  
  name.vec <- names(col)
  node.names <-strsplit(gsub("[+]"," ",gsub("[,]"," ",names(col))),split=" ") 
  length.node.names <- sapply(node.names,length) 
  idx.length <- which(length.node.names>1)
  
  if(length(idx.length)>0){
    
    for(i in 1:length(idx.length)){
      
      for(j in 1:length.node.names[idx.length[i]]){
        idx_push <- length(col)+1 
        col[idx_push] <- col[i]
        name.vec[idx_push] <- node.names[[idx.length[i]]][j]
      }
    }
    
    col <- col[-idx.length]
    names(col) <- name.vec[-idx.length]
  }
  return(col)
  
}



path.value<-function( nodes.vals, subgraph, ininodes, endnode, method="pond", maxnum = 100, tol = 0.000001, divide=F, response.tol = 0 ){
  
  # Initialize lists
  ready <- ininodes
  processed <- list()
  
  # Initialize node values
  node.signal <- matrix(NA, ncol=ncol(nodes.vals), nrow = length(V(subgraph)), dimnames = list(V(subgraph)$name, colnames(nodes.vals)))
  endnode.signal.dif <- 10
  
  num <- 0
  reached_last <- F
  while( length(ready) > 0 && num <= maxnum){
    num <- num + 1
    actnode <- ready[[1]]
    old.signal <- node.signal[actnode,]
    
    # Compute node signal
    if(divide && actnode != endnode){
      nfol <- length(incident(subgraph, actnode, mode="out"))
    }else{
      nfol <- 1
    }
    # print(nfol)
    node.signal[actnode,] <- compute.node.signal2(actnode, nodes.vals[actnode,], node.signal, subgraph, method, response.tol) / nfol
    
    # Transmit signal
    nextnodes <- get.edgelist(subgraph)[incident(subgraph, actnode, mode="out"),2]
    dif <- old.signal - node.signal[actnode,]
    
    if(actnode==endnode){
      reached_last <- T
      if(!all(is.na(dif)))
        endnode.signal.dif <- c(endnode.signal.dif, sqrt(sum(dif^2)))
      #num <- num+1
    }
    if(all(is.na(old.signal)) || endnode.signal.dif[length(endnode.signal.dif)] > tol )
      ready <- unique(c(ready, nextnodes))
    ready <- ready[-1]
  }
  if(reached_last==F){
    endnode.signal.dif <- NA
  }
  return(list(node.signal[endnode,], endnode.signal.dif,nodes.vals))
}

compute.node.signal2<-function(actnode, node.val, node.signal, subgraph, method="pond", response.tol = 0){
  incis <- incident(subgraph, actnode, mode="in")
  
  if(length(incis)==0){
    signal <- rep(1, length(node.val))
    
  } else {
    
    # get activators and inhibitors signal
    prevs <- get.edgelist(subgraph)[incis,1]
    input_signals <- node.signal[prevs,,drop=F]
    nas <- is.na(input_signals[,1])
    prevs <- prevs[!nas]
    incis <- incis[!nas]
    input_signals <- input_signals[!nas,,drop=F]
    typeincis <- E(subgraph)$relation[incis]
    activators <- typeincis==1
    nactivators <- sum(activators)
    inhibitors <- typeincis==-1
    ninhibitors <- sum(inhibitors)
    activator_signals <- input_signals[activators,,drop=F]
    inhibitor_signals <- input_signals[inhibitors,,drop=F]
    
    if( method == "sum"){
      s1 <- prettyifelse(nactivators>0, colSums(activator_signals), rep(1,length(node.val)))
      s2 <- prettyifelse(ninhibitors>0, colSums(inhibitor_signals), rep(0,length(node.val)))
      signal <- s1-s2
    }
    else if( method == "pond"){
      s1 <- prettyifelse(nactivators>0, apply(1- activator_signals, 2, prod), rep(0,length(node.val)))
      s2 <- prettyifelse(ninhibitors>0, apply(1- inhibitor_signals, 2, prod), rep(1,length(node.val)))
      signal <- (1-s1)*s2
    }
    else if( method == "min"){
      s1 <- prettyifelse(nactivators>0, apply(activator_signals,2,min), rep(1,length(node.val)))
      s2 <- prettyifelse(ninhibitors>0, 1-apply(inhibitor_signals,2,max), rep(1,length(node.val)))
      signal <- s1*s2
    }
    else {
      stop("Unknown propagation rule")
    }
    # If signal too low, signal do not propagate
    # HERE: kinza changed this to force same dimentiality for the doble &&
    # if(sum(nas) == 0 && signal < response.tol)
    if(sum(nas) == 0 && any(signal < response.tol))
      signal <- rep(0,length(node.val))
    
  }
  
  signal[signal>1] <- 1
  signal[signal<0] <- 0
  signal <- signal*node.val
  
  return(signal)
}
