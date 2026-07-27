# Function for equi distant herearchial layout
l=layout_as_tree(g,root="NICHE",flip.y = T)
len=max(l)-min(l)
lev=c(unique(l[,2]))
# finding length of each levels
levLen=lapply(lev,function(lev,l) length((l[l[,2]==lev,]))/2,l)
ele=max(unlist(levLen))
eleLen=lapply(lev,function(lev,l) ((l[l[,2]==lev,])),l)



newl=cbind(l[,1]*1000,l[,2])


if(species == "MOUSE"){
  load(system.file("extdata", "MOUSE_Background_Network_omnipath_reactome_metacore_01042019.RData", package = "SigHotSpotter"))
} else {
  if(species == "HUMAN"){
    load(system.file("extdata", "HUMAN_Background_Network_omnipath_reactome_metacore_01042019.RData", package = "SigHotSpotter"))
  } else {
    stop("Only the following species are supported: 'MOUSE', 'HUMAN'")
  }
}
idata = idata
DE_Genes=DE_Genes
if(invert_DE)
{
  DE_Genes[,2] = -DE_Genes[,2]
}

subg=Data_preprocessing(idata,cutoff,species)

## Calculate stationary distribution of the MC
Steady_state_true=Markov_chain_stationary_distribution(subg)

## Retrieves high probability intermediates
int=high_probability_intermediates(Steady_state_true, intermediates, percentile)
gintg=integrate_sig_TF(subg,Steady_state_true,DE_Genes, non_interface_TFs, TF_TF_interactions )
nTF=nonterminal_DE_TFs(gintg,DE_Genes,non_interface_TFs)

#comp_score for each TF individually
score=lapply(nTF,comp_score_tf, int, gintg)
#converting the nested list into a matrix whose row sum will give probability of each intermediate
score_m=(matrix(unlist(score), ncol=length(score), byrow=F))
score_m_means=as.list(rowMeans(score_m))
final_score=compatability_score(score_m_means,Steady_state_true,int)

## Computing networks for visualization
trimmed_score_I=.trimResults(final_score,F)
trimmed_score_A=.trimResults(final_score,T)
toiintI=c(as.matrix(trimmed_score_I$`Inactive signaling hotspots`))
toiintA=c(as.matrix(trimmed_score_A$`Active signaling hotspots`))
twoi=c(toiintI,toiintA)

#pruning the integrated networks
gintg.p=prun.int.g(gintg)

#building networks for all intermediates for active signaling hotspots
sp_int_A <- lapply(toiintA,to_sp_net_int_1path,gintg.p,nTF,DE_Genes,non_interface_TFs)

#building networks for inactive signaling hotspots
sp_int_I <- lapply(toiintI,to_sp_net_int_1path,gintg.p,nTF,DE_Genes,non_interface_TFs)

#converting to visNetwork
vis_net_A <- lapply(sp_int_A,toVisNetworkData)
vis_net_I <- lapply(sp_int_I,toVisNetworkData)

#for edge color
vis_net_A <- lapply(vis_net_A,vis.edge.color)
vis_net_I <- lapply(vis_net_I,vis.edge.color)
ret_value = list(final_score=final_score,
                 trimmed_score_A=trimmed_score_A,
                 trimmed_score_I=trimmed_score_I,
                 vis_net_A=vis_net_A,
                 vis_net_I=vis_net_I )

to_sp_net_int_1path <- function(s,g,t,deg,non_interface_TFs){
  #changing the edge attributes of the integrated network
  #g=non_neg_weight(g)
  #removing the edges from dummy to TFs as it affectes the shortest paths
  #del=as_adj_edge_list(g, mode = c("in"))$Dummy
  #g <- delete.edges(g,del)
  #Shortest path edges
  t=intersect(t,V(g)$name)
  edges_a=lapply(s,shortest_path_edges,t,g)
  edges_d=lapply("NICHE",shortest_path_edges,c(s),g)
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

#Function for integrating signaling and TF networks, g=signaling graph, x steady state vector, deg=differentially expressed genes
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

#' Define interface_TF
#' 
interface_TFs<-as.data.frame(setdiff(V(graph.data.frame(TF_TF_interactions))$name,as.character(as.vector(non_interface_TFs)$Gene)))
names(interface_TFs)<-"Gene"
#select high probable intTFs
intInterTFS<-high_probability_intermediates(Steady_state_true,interface_TFs,70)
#select those intTFs that are not high probability and delete them
delIntTFs<-setdiff(high_probability_intermediates(Steady_state_true,interface_TFs,0),high_probability_intermediates(Steady_state_true,interface_TFs,75))
gintg.pDel<-delete.vertices(gintg.p,delIntTFs)

#' defining trur nonIntDETFs
#' whatever DEG intersects with subg is true interface TF
intTF=intersect((as.vector(t(deg[1]))),V(subg)$name)
nonintTF=intersect((as.vector(t(deg[1]))),V(graph.data.frame(TF_TF_interactions))$name)
#' New non-terminal DETFS
non_interface_DE_TFs=as.data.frame(setdiff(nonintTF,intTF))
names(non_interface_DE_TFs)<-"Gene"
non_interface_DE_TFs<-join(deg,non_interface_DE_TFs,by=c("Gene"),type="inner")

#' Adding akt1 foxo interactions to the database
aktints=read.table("aktintrctions.txt",sep="\t",header = T)
Background_signaling_interactome=rbind(Background_signaling_interactome,aktints)

#' comparing int perturbation with the original response
map2k1P=perturb_receptor("Map2k1",subg)
itf_ss=merge(interface_TFs,Steady_state_true)
itf_ssP=merge(interface_TFs,map2k1P$Steady_state_true)
merge_its=merge(itf_ss,itf_ssP,by="Gene")
merge_its=cbind(merge_its,merge_its[2]/merge_its[3])


gintg=integrate_sig_TFGeneral(subg,Steady_state_true,DE_Genes, non_interface_TFs, TF_TF_interactions)
nTF=intersect(V(gintg)$name,as.vector(t(non_interface_DE_TFs$Gene)))


 trng=buildTrn(DE_Genes,idata,99)
 gintg=sigTrnIntegratedNet(subg,trng,Steady_state_true,DE_Genes)
 non_interface_DE_TFs <- nonIntDeTfs(subg,DE_Genes)
 nTF=intersect(V(gintg)$name,as.vector(t(non_interface_DE_TFs$Gene)))
 length(nTF)
