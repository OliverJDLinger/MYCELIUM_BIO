#function to get "all"one" shortest path edges from receptor_lig to intTFs based on weights
shortest_path_edges <- function(s,t,g){
  if (length(s) == 0) stop ('No intermediates found. You may decrease the percentile in order to find intermediates.') # decrease percentile cutoff
  if (length(t) == 0) stop ('No non-terminal differentially expressed TFs found. You may decrease the cutoff.') # decrease normal cutoff
  paths=(get.all.shortest.paths(g, s, t, mode = c("out"))$res)
  if (length(paths) == 0) stop ('No shortest path found. You may decrease the cutoff in order to find shortest path.') # decrease normal cutoff
  edges=lapply(paths,edge_shortest,g)
  return(edges)
}

#function to get all shortest path edges from receptor_lig to intTFs
shortest_path_edges_all <- function(s,t,g){
  if (length(s) == 0) stop ('No intermediates found. You may decrease the percentile in order to find intermediates.') # decrease percentile cutoff
  if (length(t) == 0) stop ('No non-terminal differentially expressed TFs found. You may decrease the cutoff.') # decrease normal cutoff
  paths=(get.all.shortest.paths(g, s, t, mode = c("out"),weight=NA)$res)
  if (length(paths) == 0) stop ('No shortest path found. You may decrease the cutoff in order to find shortest path.') # decrease normal cutoff
  edges=lapply(paths,edge_shortest,g)
  return(edges)
}

#function for returning all shortest paths
shortest.paths.all <- function(s,t,g){
    if (length(s) == 0) stop ('No intermediates found. You may decrease the percentile in order to find intermediates.') # decrease percentile cutoff
    if (length(t) == 0) stop ('No non-terminal differentially expressed TFs found. You may decrease the cutoff.') # decrease normal cutoff
    paths=(get.all.shortest.paths(g, s, t, mode = c("out"),weight=NA)$res)
    if (length(paths) == 0) stop ('No shortest path found. You may decrease the cutoff in order to find shortest path.') # decrease normal cutoff
    return(paths)
}

#funtion to get paths for each tf
paths.tf <- function(t,s,g) #x is a list of comp score
{
  paths=lapply(s,shortest.paths.all,t,g)
  #edges=lapply(paths,edge_shortest,g)
}

#function for retaining edges
edge_shortest <- function(path, graph)
{
  edges=E(graph, path=path)
}

# Function to trim results for display
#
# The function returns shortlist of best results
#
# @param results all results comuputed by SigHotSpotter_pipeline
# @param active Active (TRUE) or inactive (FALSE)
# @return Shortlisted results
# @export
.trimResults <- function(results, active = TRUE) {

  res_trimmed = results[,1:2]
  if (active){
    res_trimmed <- head(res_trimmed, 10L)
    res_trimmed <- res_trimmed[res_trimmed[,2]>0.5,]
    colnames(res_trimmed) <- c('Active signaling hotspots', 'Compatibility score')
  } else
  {
    res_trimmed <- tail(res_trimmed, 10L)
    res_trimmed <- res_trimmed[res_trimmed[,2]<0.5,]
    colnames(res_trimmed) <- c('Inactive signaling hotspots', 'Compatibility score')
    res_trimmed <- res_trimmed[order(res_trimmed$'Compatibility score'),]
  }
  res_trimmed[,2] = round(res_trimmed[,2],4)
  rownames(res_trimmed) <- NULL
  return(res_trimmed)
}

#making the edge weight non-negative and taking the absolute
prun.int.g <- function(g){
  #E(g)$weight=1-(abs(E(g)$weight))
  #function to delet edges from tf to dummy un the network
  del=incident(g, "Dummy", mode = c("in"))
  g <- delete.edges(g,del)
  g <- set.vertex.attribute(g,"name","Dummy","NICHE")
  return(g)
}

#making the edge weight abs
abs_weight <- function(g){
  E(g)$weight=(abs(E(g)$weight))
  return(g)
}



#function to generate the shortest path network given source, target and a network
shortest_path_network <- function (s,t,g, mst){
  #changing the edge attributes of the integrated network
  #g=non_neg_weight(g)
  #removing the edges from TFs to dummy as it affectes the shortest paths
  #del=as_adj_edge_list(g, mode = c("in"))$Dummy
  #g <- delete.edges(g,del)
  if (mst==T){
    edges_g=lapply(s,shortest_path_edges,t,g)
    net <- subgraph.edges(g, unlist(edges_g), delete.vertices = T)
    net=mst(net)

  } else
  {edges_g=lapply(s,shortest_path_edges,t,g)
  net <- subgraph.edges(g, unlist(edges_g), delete.vertices = T)
  }
  return(net)
}

# function for converting the input graph and the ints to a shortestpat network
# @param a,i active and inactive source nodes or int
# @param g input graph on which shortest path network must be inferred i.e. gintg
# @param t terminal nodes
to_sp_net <- function(a,i,g,t){
  #changing the edge attributes of the integrated network
  #g=non_neg_weight(g)
  #removing the edges from dummy to TFs as it affectes the shortest paths
  #del=as_adj_edge_list(g, mode = c("in"))$Dummy
  #g <- delete.edges(g,del)
  #Shortest path edges
  edges_a=lapply(a,shortest_path_edges,t,g)
  edges_i=lapply(i,shortest_path_edges,t,g)
  edges_d=lapply("NICHE",shortest_path_edges_all,c(i),g)
  #subnetwork from SP edges
  sp_sub_net=subgraph.edges(g, c(unlist(edges_a),unlist(edges_i),unlist(edges_d)), delete.vertices = T)
  sp_sub_net=set_vertex_attr(sp_sub_net, "group", index = c(a), value="a_int")
  sp_sub_net=set_vertex_attr(sp_sub_net, "group", index = c(i), value="i_int")
  sp_sub_net=set_vertex_attr(sp_sub_net, "group", index = c(t), value="tf")
  return(sp_sub_net)
}

# function for converting the input graph and the ints to a shortestpat network
# @param a,i active and inactive source nodes or int
# @param g input graph on which shortest path network must be inferred i.e. gintg
# @param t terminal nodes
to_sp_net_int <- function(s,g,t,deg,non_interface_TFs){
  #changing the edge attributes of the integrated network
  #g=non_neg_weight(g)
  #removing the edges from dummy to TFs as it affectes the shortest paths
  #del=as_adj_edge_list(g, mode = c("in"))$Dummy
  #g <- delete.edges(g,del)
  #Shortest path edges
  t=intersect(t,V(g)$name)
  edges_a=lapply(s,shortest_path_edges_all,t,g)
  edges_d=lapply("NICHE",shortest_path_edges_all,c(s),g)
  #classifying up and downregulated TFs
  up_t=up_down_tfs(g,deg[deg[2]==1,],non_interface_TFs)
  down_t=up_down_tfs(g,deg[deg[2]==-1,],non_interface_TFs)
  #subnetwork from SP edges
  sp_sub_net=subgraph.edges(g, c(unlist(edges_a),unlist(edges_d)), delete.vertices = T)
  sp_sub_net=set_vertex_attr(sp_sub_net, "group", index = c(s), value="int")
  sp_sub_net=set_vertex_attr(sp_sub_net, "group", index = c(up_t[up_t %in% V(sp_sub_net)$name]), value="upregulated")
  sp_sub_net=set_vertex_attr(sp_sub_net, "group", index = c(down_t[down_t %in% V(sp_sub_net)$name]), value="downregulated")
  return(sp_sub_net)
}

# function for changing the visnetwork edge color
# @param visg the visnetwork object
vis.edge.color <- function(visg){
  visg$edges$color[visg$edges$Effect==1]="green"
  visg$edges$color[visg$edges$Effect==-1]="red"
  return(visg)
}

# obsoleted, this function is used from app.R
  # Function to plot the Visnetwork object
  # @param visg the visnetwork object
 vis.net.plot <- function(visg){
    #hierarchy
  visNetwork(visg$nodes,visg$edges) %>% visNodes(visg, shape="box") %>%
  visIgraphLayout(layout = "layout_as_tree",root="NICHE",flip.y = F) %>%
  visEdges(arrows = "to") %>%  visOptions(highlightNearest = list(enabled =TRUE, degree = 1, hover = T), nodesIdSelection = TRUE)  %>%
  visEdges(smooth = T) %>% visGroups(visg, groupname="int", shape="circle", color ="orange") %>%
  visGroups(visg, groupname="upregulated", color = "green",shape="triangle") %>%
  visGroups(visg, groupname="downregulated", color = "red", shape="triangle") %>%
  visPhysics(stabilization = FALSE) %>% visEdges(smooth = FALSE) %>%
  visExport()
}

#unused funtion for plotting the entire network
#' function for plotting a union network
#' @param visg is the visnetwork graph object
#' @param x: is the table of trimmed intermeduates
vis_plot_union <- function(visg,x){
  int=as.vector.factor((x[1])[1:nrow(x),])
  visge=lapply(visg,function(x) x$edges)
  visgn=lapply(visg,function(x) x$nodes)
  visgnj=join_all(visgn,by=c("id"),type="full")
  visgej=join_all(visge,by=c("from","to"),type="full")
  vis_union_A <- list(nodes=visgnj,edges=visgej)
  vis_union_A$nodes$group[vis_union_A$nodes$id %in% int]="int"
  return(vis_union_A)
}

#' Function to identify feed back arc sets
#' @param g:input graph for which feedback edges are to be found
#' @param s sources: set of nodes defined as the source nodes
#' @param t sinks: set of nodes defined as the target nodes 
findFAS <- function (g){
  gReal=g #saving the original graph for later use
  # calculate in and out degree
  V(g)$indegree<-igraph::degree(g, v=V(g), mode = c("in"))
  V(g)$outdegree<-igraph::degree(g, v=V(g), mode = c("out"))
  #define 3 bins s1,s2,d denoting the sources, sinks and intermediates
  s1<-V(g)[indegree==0]$name
  s2<-V(g)[outdegree==0]$name
  i=V(g)[indegree>0&outdegree>0]$name
  
  del=c(s1,s2)
  #delete vertices from graph
  while(length(del)!=0)
  {
    g <- delete.vertices(g,del)
    V(g)$indegree<-igraph::degree(g, v=V(g), mode = c("in"))
    V(g)$outdegree<-igraph::degree(g, v=V(g), mode = c("out"))
    s1<-c(s1,(V(g)[indegree==0])$name)
    s2<-c(s2,(V(g)[outdegree==0])$name)
    i<-c(i,setdiff(V(g)$name,c(s1,s2)))
    del<-c(V(g)[indegree==0]$name,V(g)[outdegree==0]$name)
  }
  i<-setdiff(V(g)$name,c(s1,s2))
  # Once these source and sinks are defined, work with the difference
  while(length(V(g))!=0)
  {
    g <- delete.vertices(g,del)
    #del=character(0)
    V(g)$indegree<-igraph::degree(g, v=V(g), mode = c("in"))
    V(g)$outdegree<-igraph::degree(g, v=V(g), mode = c("out"))
    V(g)$diffdegree<-V(g)$outdegree-V(g)$indegree
    if (length(del)==0) {
      del<-((V(g)[diffdegree==(max(V(g)$diffdegree))])$name)
      #print(del)
      s1<-c(s1,del)
    } else {
      del<-c(V(g)[indegree==0]$name,V(g)[outdegree==0]$name)
      #print(del)
    }
    s1<-c(s1,(V(g)[indegree==0])$name)
    s2<-c((V(g)[outdegree==0])$name,s2)
    nodeOrder<-c(s1,s2)
  }
  # Retrieve the feedback edges
  # Creating an attribute rank for the edge
  g=gReal
  V(g)[nodeOrder]$rank <- seq_along(nodeOrder)
  # Creating an edge attribute
  el=get.edgelist(g)
  # Creating an edge attribute FAS by taking the difference between the rank of two nodes ineach edge. Edges with positive ranks are FAS
  E(g)$fas <- V(g)[el[,1]]$rank - V(g)[el[,2]]$rank
  # Flag FAS with positive scores to be deleted
  fas<-E(g)[E(g)$fas>0]
  gDag <- delete.edges(g,fas)
  results <- list(gDag=gDag,originalgraph=gReal,FAS=fas)
  return(results)
}

#' Yen's K directed shortest paths
#' @return the shortest path as a list of vertices or NULL if there is no path between src and dest
shortest_path <- function(graph, src, dest){
  path <- suppressWarnings(get.shortest.paths(graph, src, dest,mode=c("out"), weights=NA))
  path <- names(path$vpath[[1]])
  if (length(path)==1) NULL else path
} 

#'@return the sum of the weights of all the edges in the given path
path_weight <- function(path, graph) sum(E(graph, path=path)$weight)

#' @returnthe effect of each edge in a path
path_effect <- function(path, graph) (E(graph, path=path)$Effect)

#'@description sorts a list of paths based on the weight of the path
sort_paths <- function(graph, paths) paths[paths %>% sapply(path_weight, graph) %>% order]

#'@description creates a list of edges that should be deleted
find_edges_to_delete <- function(A,i,rootPath){
  edgesToDelete <- NULL
  for (p in A){
    rootPath_p <- p[1:i]
    if (all(rootPath_p == rootPath)){
      edge <- paste(p[i], ifelse(is.na(p[i+1]),p[i],p[i+1]), sep = '|')
      edgesToDelete[length(edgesToDelete)+1] <- edge
    }
  }
  unique(edgesToDelete)
}

#returns the k shortest path from src to dest
#sometimes it will return less than k shortest paths. This occurs when the max possible number of paths are less than k
k_shortest_yen <- function(graph, src, dest, k){
  if (src == dest) stop('src and dest can not be the same (currently)')
  
  #accepted paths
  A <- list(shortest_path(graph, src, dest))
  if (k == 1) return (A)
  #potential paths
  B <- list()
  
  for (k_i in 2:k){
    prev_path <- A[[k_i-1]]
    num_nodes_to_loop <- length(prev_path)-1
    for(i in 1:num_nodes_to_loop){
      spurNode <- prev_path[i]
      rootPath <- prev_path[1:i]
      
      edgesToDelete <- find_edges_to_delete(A, i,rootPath)
      t_g <- delete.edges(graph, edgesToDelete)
      #for (edge in edgesToDelete) t_g <- delete.edges(t_g, edge)
      
      spurPath <- shortest_path(t_g,spurNode, dest)
      
      if (!is.null(spurPath)){
        total_path <- list(c(rootPath[-i], spurPath))
        if (!total_path %in% B) B[length(B)+1] <- total_path
      }
    }
    if (length(B) == 0) break
    B <- sort_paths(graph, B)
    A[k_i] <- B[1]
    B <- B[-1]
  }
  A
}

#' Function for perturbing each receptor/lig and calculating the compatability score for only one phenotype
#' @param r receptor/intermediate to be perturbed, can bea list or a single receptor or intermediate
#' @param g subg before perturbation
#' @param real real simulation data
perturb_receptor <- function(r,g,real,percentile,DE_Genes){
  #delete all receptors except 1
  g=delete.vertices(g,r)
  #To ensure reachability for the Markov chain
  V(g)$degree=igraph::degree(g, v=V(g), mode = c("in"))
  #Select Nodes to be deleted
  del=V(g)[degree==0]
  #delete vertices from graph
  while(length(del)!=0)
  {
    g <- delete.vertices(g,del)
    V(g)$degree=igraph::degree(g, v=V(g), mode = c("in"))
    del=V(g)[degree==0]
  }
  #Same as above but remove nodes with with zero out degree
  V(g)$degree=igraph::degree(g, v=V(g), mode = c("out"))
  #Select Nodes to be deleted
  del=V(g)[degree==0]
  while(length(del)!=0)
  {
    g <- delete.vertices(g,del)
    V(g)$degree=igraph::degree(g, v=V(g), mode = c("out"))
    del=V(g)[degree==0]
  }
  #####TO EXTRACT THE LARGEST STRONGLY CONNECTED COMPONENT
  members <- membership(clusters(g, mode="strong"))
  #l=lapply(unique(members), function (x) induced.subgraph(g, which(members == x)))
  #g=l[[1]]
  SCC <- clusters(g, mode="strong")
  subg <- induced.subgraph(g, which(membership(SCC) == which.max(sizes(SCC))))
  #subg=simplify(subg,edge.attr.comb=list("first"))
  
  Steady_state_true=Markov_chain_stationary_distribution(subg)
  
  ## Retrieves high probability intermediates
  int=high_probability_intermediates(Steady_state_true, intermediates, percentile)
  ## Retrieves high probability intTFs
  tflist<-as.data.frame(tflist$Symbol)
  names(tflist)<-"Gene"
  intTF<-high_probability_intermediates(Steady_state_true, tflist, (percentile-15))
  intTF<-as.data.frame(intTF)
  names(intTF)<-"Gene"
  gintg=integrate_sig_TF(subg,Steady_state_true,DE_Genes, non_interface_TFs, TF_TF_interactions )
  nTF=nonterminal_DE_TFs(gintg,DE_Genes,non_interface_TFs)
  
  #comp_score for each TF individually
  score=lapply(nTF,comp_score_tf, int, gintg)
  #converting the nested list into a matrix whose row sum will give probability of each intermediate
  score_m=(matrix(unlist(score), ncol=length(score), byrow=F))
  score_m_means=as.list(rowMeans(score_m))
  final_score=compatability_score(score_m_means,Steady_state_true,int)
  perturbSSscore <- join_all(list(as.data.frame(intTF),real$SS,Steady_state_true,real),by="Gene",type="inner")
  perturbSSscore<-cbind(perturbSSscore,perturbSSscore[3]/perturbSSscore[2])
  perturbSSscore <- perturbSSscore[order(perturbSSscore$'foldChangePerturbed'),]
  names(perturbSSscore)<-c("Gene","Steady_state_True","Steady_state_perturbed", "foldChangePerturbed")
  ## Computing networks for visualization
  trimmed_score_I=.trimResults(final_score,F)
  trimmed_score_A=.trimResults(final_score,T)
  toiintI=c(as.matrix(trimmed_score_I$`Inactive signaling hotspots`))
  toiintA=c(as.matrix(trimmed_score_A$`Active signaling hotspots`))
  twoi=c(toiintI,toiintA)
  
  #pruning the integrated networks
  gintg.p=prun.int.g(gintg)
  
  #building networks for all intermediates for active signaling hotspots
  sp_int_A <- lapply(toiintA,to_sp_net_int,gintg.p,nTF,DE_Genes,non_interface_TFs)
  
  #building networks for inactive signaling hotspots
  sp_int_I <- lapply(toiintI,to_sp_net_int,gintg.p,nTF,DE_Genes,non_interface_TFs)
  
  #converting to visNetwork
  vis_net_A <- lapply(sp_int_A,toVisNetworkData)
  vis_net_I <- lapply(sp_int_I,toVisNetworkData)
  
  #for edge color
  vis_net_A <- lapply(vis_net_A,vis.edge.color)
  vis_net_I <- lapply(vis_net_I,vis.edge.color)
  ret_value = list(Steady_state_true=Steady_state_true,final_score=final_score,
                   trimmed_score_A=trimmed_score_A,
                   trimmed_score_I=trimmed_score_I,
                   vis_net_A=vis_net_A,
                   vis_net_I=vis_net_I,perturbedScore=perturbSSscore )
  return(ret_value)
}

#' @function to save igraph graph in edge list format with two edge attribute effect and weight
#' @param g: Igraph object
saveIgraph <- function(g){
  el<-as.data.frame(as_edgelist(g,names = T))
  el<-cbind(el,edge_attr(g)$Effect,edge_attr(g)$weight)
  names(el)<-c("Source","Target","Effect","weight")
  return(el)
}

#' @Function for integrating signaling and TF networks, integrates with all immediates downstream TFs of IntTFs.
#' @param g=signaling graph (subg), x steady state vector, deg=differentially expressed genes
integrate_sig_TFGeneral <- function(g,x,deg, non_interface_TFs, TF_TF_interactions )
{
  #delete sccDummy inthe network #only for FAS removal
  #g=delete.vertices(g,"sccDummy")
  el=as_edgelist(g)
  graph_markov=as.data.frame(cbind(el,E(g)$Effect))
  colnames(graph_markov)=c("Source","Target","Effect")
  colnames(deg)=c("Gene","Bool_DEG")
  #' Not joining by overll interface TFs but data sepcific interface TFs
  # non_interface_DE_TFs=join(deg,non_interface_TFs,by=c("Gene"),type="inner")
  #' data specific non-interfaceDETFs
  #' defining trur nonIntDETFs
  #' whatever DEG intersects with subg is true interface TF
  intTF=intersect((as.vector(t(deg[1]))),V(g)$name)
  nonintTF=intersect((as.vector(t(deg[1]))),V(graph.data.frame(TF_TF_interactions))$name)
  #' New non-terminal DETFS
  non_interface_DE_TFs=as.data.frame(setdiff(nonintTF,intTF))
  names(non_interface_DE_TFs)<-"Gene"
  non_interface_DE_TFs<-join(deg,non_interface_DE_TFs,by=c("Gene"),type="inner")
  #Get the phenotype from the non_interface_DE_TFs and map it to the TF-TF interaction network
  colnames(non_interface_DE_TFs)[1]="Target"
  DE_TF_TF_interactions_target=join(TF_TF_interactions,non_interface_DE_TFs,by=c("Target"),type="left")
  DE_TF_TF_interactions_target=na.omit(DE_TF_TF_interactions_target)
  ab=DE_TF_TF_interactions_target
  #ab$Effect=ab$Effect*ab$Phenotype1 #second Phenotype will be opposite of this #must inpust second data for second phenotype
  names(ab)<-NULL
  ab=as.matrix(ab)
  ab[,3]=as.numeric(ab[,3])*as.numeric(ab[,4])
  ab=as.data.frame(ab)
  names(ab)=c("Source","Target","Effect","DEG")
  graph_markov$Effect=as.numeric(as.character(graph_markov$Effect))
  graph_markov=rbind(graph_markov,ab[1:3]) #merging the nTF interaction with appropriate sign Effect with the original graph
  colnames(x)[1] <- "Source"
  ab=join(graph_markov,x,by=c("Source"),type="left",match="first")
  colnames(ab)[3:ncol(ab)]="source"
  colnames(x)[1] <- "Target"
  ab1=join(graph_markov,x,by=c("Target"),type="left",match="first")
  names(ab1) <- NULL
  names(ab) <- NULL
  ab=ab[,4:ncol(ab)]
  ab1=ab1[,4:ncol(ab1)]
  #creating node SS as the edge property
  weight=as.numeric(as.matrix(ab))
  #edge_P=as.data.frame(weight)
  graph_markov=(cbind(graph_markov,weight))
  graph_markov[is.na(graph_markov)] <- 1  #Making TF-TF interactions dependent only on the expression status
  graph_markov$Effect=as.numeric(as.matrix((graph_markov$Effect)))
  g3 <- graph.data.frame(as.data.frame(graph_markov))
  #updating the graph attribute for the adjacency matrix i.e. product SS (weight) and effect
  E(g3)$weight=E(g3)$weight*E(g3)$Effect
  #deleting TF nodes with no indegree
  V(g3)$degree=igraph::degree(g3, v=V(g3), mode = c("in"))
  #Select Nodes to be deleted
  del=V(g3)[degree==0]
  #delete vertices from graph
  while(length(del)!=0)
  {
    g3 <- delete.vertices(g3,del)
    V(g3)$degree=igraph::degree(g3, v=V(g3), mode = c("in"))
    del=V(g3)[degree==0]
  }
  g3
}

#' @Function for integrating signaling and transcriptional networks considereing compatible TF-TF interactions
#' Also here the pathlength from the interfaceTFs will to nonintDETFs will be more than 1
#' non_interface_TFs are those that do not have a indegree from a signaling molecule. This is initially calculated based on the background interactome and TF-TF interactions from metacore
#' All other TFs that either have indegree only from signaling molecules or from both signaling molecules and TFs are interface_TFs. This set of interface_TFs will vary from data to data
#' 
integrateSignalingcompatibleTRN <- function(g,deg,x, non_interface_TFs, TF_TF_interactions){
  unionBackSigTftf <- graph.data.frame(rbind(Background_signaling_interactome,TF_TF_interactions))
  unionBackSigTftf <- delete.vertices(unionBackSigTftf,as.vector(t(non_interface_TFs$Gene)))
  a=as.data.frame(cbind(as_edgelist(unionBackSigTftf),E(unionBackSigTftf)$Effect))
  names(a)<-c("Source","Target","Effect")
  deginNet=intersect((as.vector(t(deg[1]))),V(g)$name)
  allTF=intersect((as.vector(t(deg[1]))),V(graph.data.frame(TF_TF_interactions))$name)
  interfaceDETFs <- intersect(V(g)$name,allTF)
  non_interface_DE_TFs=as.data.frame(setdiff(allTF,deginNet))
  names(non_interface_DE_TFs)<-"Gene"
  non_interface_DE_TFs<-join(deg,non_interface_DE_TFs,by=c("Gene"),type="inner")
  
}

#' @Function to add TF-TF interactions to the backgroundsignaling interaction
#' @param g: background_signalin interactome
#' @param t: TF-TF interactions
#' @return a: Background network with Tf-Tf interactions
addTfIntBack <- function(Background_signaling_interactome,TF_TF_interactions){
  unionBackSigTftf <- graph.data.frame(rbind(Background_signaling_interactome,TF_TF_interactions))
  unionBackSigTftf <- delete.vertices(unionBackSigTftf,as.vector(t(non_interface_TFs$Gene)))
  a=as.data.frame(cbind(as_edgelist(unionBackSigTftf),E(unionBackSigTftf)$Effect))
  names(a)<-c("Source","Target","Effect")
  return(a)
}

#' @Function to find non_interface_DE_TFs
#' @param g: subg
#' @param deg: differentially expressed genes
#' 
  nonIntDeTfs <- function(g,deg){
  deginNet=intersect((as.vector(t(deg[1]))),V(g)$name)
  allTF=intersect((as.vector(t(deg[1]))),V(graph.data.frame(TF_TF_interactions))$name)
  interfaceDETFs <- intersect(V(g)$name,allTF)
  non_interface_DE_TFs=as.data.frame(setdiff(allTF,deginNet))
  names(non_interface_DE_TFs)<-"Gene"
  non_interface_DE_TFs<-join(deg,non_interface_DE_TFs,by=c("Gene"),type="inner")
  }
  
#' @Function add attributes from a dataframe
#' @param graph input graph
#' @param D dataframe with attributes
  SetNodeAttributes<-function (graph, D) 
  {
    if (!is.igraph(graph)) 
      stop("Not a graph object")
    if (!is.data.frame(D)) 
      stop("Not a dataframe")
    for (i in colnames(D)) graph <- set.vertex.attribute(graph, i, value = D[, i])
    return(graph)
  }

#' @Function to build TF-TF network based on expressed genes
#' @param Deg
#' @param input_data
#' @param cutoff
#' @param Steady_state
#' @param cutoff
  buildTrn <- function(deg,input_data,cutoff){
    b=as.data.frame(input_data)
    dd=join(deg[1],b,by=c("Gene"),type="left")
    #filter genes not expressed in more than these many percent pf cells
    b=filter_exp(b,cutoff)
    b=rbind(b,dd)
    b=unique(b)
    #This is to convert chr into numericversion
    b[2:ncol(b)]<-as.data.frame(lapply(b[2:ncol(b)],as.numeric,b[2:ncol(b)]))
    ##Renaming the first column for finding the union
    colnames(b)[1] <- "Source"
    # loading TF-TF interactions
    a<-TF_TF_interactions
    a=graph.data.frame(TF_TF_interactions)
    g=induced.subgraph(a,intersect(as.vector(t(b$Source)),V(a)$name))
    SCC <- clusters(g, mode="weak")
    trng <- induced.subgraph(g, which(membership(SCC) == which.max(sizes(SCC))))
    return(trng)
  }
  
#' @Function to create union of two graphs while preserving the attribute names
#' @param g1, g2 the two input graphs
  union2<-function(g1, g2){
    
    #Internal function that cleans the names of a given attribute
    CleanNames <- function(g, target){
      #get target names
      gNames <- parse(text = (paste0(target,"_attr_names(g)"))) %>% eval 
      #find names that have a "_1" or "_2" at the end
      AttrNeedsCleaning <- grepl("(_\\d)$", gNames )
      #remove the _x ending
      StemName <- gsub("(_\\d)$", "", gNames)
      
      NewnNames <- unique(StemName[AttrNeedsCleaning])
      #replace attribute name for all attributes
      for( i in NewnNames){
        
        attr1 <- parse(text = (paste0(target,"_attr(g,'", paste0(i, "_1"),"')"))) %>% eval
        attr2 <- parse(text = (paste0(target,"_attr(g,'", paste0(i, "_2"),"')"))) %>% eval
        
        g <- parse(text = (paste0("set_",target,"_attr(g, i, value = ifelse(is.na(attr1), attr2, attr1))"))) %>%
          eval
        
        g <- parse(text = (paste0("delete_",target,"_attr(g,'", paste0(i, "_1"),"')"))) %>% eval
        g <- parse(text = (paste0("delete_",target,"_attr(g,'", paste0(i, "_2"),"')"))) %>% eval
        
      }
      
      return(g)
    }
    
    g <- igraph::union(g1, g2) 
    #loop through each attribute type in the graph and clean
    for(i in c("graph", "edge", "vertex")){
      g <- CleanNames(g, i)
    }
    
    return(g)
    
  }

  #' @Function to extract edge weight and node weight, given a path
  #' @param 
  #' @param 
  # Extracting the edge weights or effects
  product_path_weight1<-function(path, graph){
  ew<-as.numeric(E(graph, path=path)$Effect)
  # extracting the node weights
  nw=as.numeric(path$weight)
  ew=prod(ew[ew!=0])
  nw=prod(nw[nw!=0])
  sum_weight=sum(ew,nw)
  if (sum_weight==0) {p_edge_weights=0} else
  {p_edge_weights=prod(ew,nw)}
  x=p_edge_weights
  x=unlist(x)
  x_pos=(x[x>0])
  x_neg=(x[x<0])
  #SPcompatability<-list(compatible=length(x_pos),incompatible=length(x_neg))
  return(list(p_edge_weights=p_edge_weights,compatiblePaths=x_pos,incompatiblePaths=x_neg))
  }
  
#' @Function to create the integrated network
#' @param subg,trng
#' @param steady state, deg
  sigTrnIntegratedNet <- function(subg,trng,ss,deg){
    # union of subg and trng
    intg<-union2(subg,trng)
    freenodes<-as.data.frame(setdiff(V(intg)$name,V(subg)$name))
    names(freenodes)<-"Gene"
    deg1=join(freenodes,deg,by="Gene",type="inner")
    names(deg1)<-c("Gene","weight")
    # combining steady state and deg
      names(deg)<-c("Gene","weight")
      names(ss)<-c("Gene","weight")
    # find those genes that are differentially expressed and also have a SS
      #ssdegOverlaP<-intersect(as.vector(t(ss[1])),as.vector(t(deg[1])))
      #deg1<-deg[ ! deg$Gene %in% ssdegOverlaP, ]
    # remove deg not inthe subg
    #rmDeg<-setdiff(as.vector(t(deg1[1])),as.vector(t(ss[1])))
    #deg1<-deg1[ ! deg1$Gene %in% rmDeg, ]
    
    weig<-rbind(ss,deg1)
    tfexp<-as.data.frame(setdiff(V(intg)$name,as.vector(t(weig$Gene))))
    tfexp=cbind(tfexp,mean(ss$weight)/100)
    names(tfexp)<-c("Gene","weight")
    weig<-rbind(weig,tfexp)
    #making sure the node attributes are mapped correctly
    ll=as.data.frame(V(intg)$name)
    names(ll)="Gene"
    ll=join(ll,weig,by=c("Gene"),type="left")
    intg<-SetNodeAttributes(intg,ll)  
    return(intg)
  }

  #' Filter TFs based on a cutoff (percentage of expression) 
  #'
  #' @param x tfs expression matrix
  #' @param y expression cutoff (in percentage)
  #' @param normData if norm data = true, do not replace value < 1 to 0 because do not make any sense
  #'
  #' @return filtered matrix based on the expression cutoff
  #' 
  filter_exp <- function(x,y,normData=FALSE){ 
    if(!normData){
      x[x < 1] <- 0 #make any TF expression less than 1 as 0
    }
    percent_expressed=((rowSums(x != 0))*100)/ncol(x)
    ip_tfs_percent=cbind(x,percent_expressed)
    ip_expressed=ip_tfs_percent[ip_tfs_percent$percent_expressed>=y,]
    #dropping the percent column
    ip_expressed$percent_expressed<-NULL
    return(ip_expressed)
  }
  
  #' @Function to remove nTFs that cannot be reached
  #' @param g: graph
  #' @param s: source
  #' @param t: target
    RmNonreachNtFs <- function(g,s,t){
    distMat <- distances(g,s,t,mode = "out")
    cs=(colSums(distMat))
    nTF<-names(cs[cs!=Inf])
    return(nTF)
  }
 
  #' @Function to remove incompatible edges in GRN
  #' @param g: graph
  #' @param deg: DEG
    RmInCompEdges <- function(g,deg){
      
    }
    
  #' @Function to assess the significance of comptible shortest paths based on hypergeometric test
  #' @param g: gintg or the integrated network
  #' @param deg: differentially expressed TFs
  #' 
  #' 

  #' @Function to convert the gintg without weight into one with edge weights
  #' @param g: gintg
   convertGintg <- function(s,g,t){
     v=as.data.frame(cbind(V(g)$name,V(g)$weight))
     names(v)=c("target","wei")
     gel=as.data.frame(as_edgelist(g,names = T))
     names(gel)=c("source","target")
     glu=join(gel,v,by="target",type="left",match="all")
     g=graph.data.frame(glu)
     g=set.edge.attribute(g,"weight",index = E(g),value = as.numeric(as.numeric(c((E(g)$wei)))))
     paths=(get.all.shortest.paths(g, s, t, mode = c("out"), weights=NA)$res)
     x=lapply(paths,product_path_weight,g)
     x=unlist(x)
     x_pos=x[x>0]
     x_neg=x[x<0]
     SPcompatability<-list(compatible=length(x_pos),incompatible=length(x_neg))
     return(SPcompatability)
   }