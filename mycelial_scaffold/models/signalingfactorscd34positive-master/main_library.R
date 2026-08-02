Data_preprocessing <- function(input_data,cutoff,species,FAS=NULL)
{
  if(species == "MOUSE"){
    #load("C:/Users/olive/OneDrive/Desktop/MYCELIUM_BIO/mycelial_scaffold/models/signalingfactorscd34positive-master/MOUSE_Background_Network_omnipath_reactome_metacore_01042019.RData")
  } else {
    if(species == "HUMAN"){
      #load(system.file("extdata", "HUMAN_Background_Network_omnipath_reactome_metacore_01042019.RData", package = "SigHotSpotter"))
    } else {
      stop("Only the following species are supported: 'MOUSE', 'HUMAN'")
    }
  }
  b=as.data.frame(input_data)
  ##COMVERT CELLS WITH LESS THAT 1 FPKM TO 0
  b[b < 1] <- 0
  ##Add a new gene Dummy with expression value 1, Only works if Dummy is present in the initial network
  b[nrow(b)+1, ] <- c("Dummy", rep(1,(ncol(b)-1)))
  #This is to convert chr into numericversion
  b[2:ncol(b)]<-as.data.frame(lapply(b[2:ncol(b)],as.numeric,b[2:ncol(b)]))
  ##Renaming the first column for finding the union
  colnames(b)[1] <- "Source"
  
  #' Adding akt1 foxo interactions to the database since they were having conflicting evedeice in Omnipath and removed originally
  Background_signaling_interactome=rbind(Background_signaling_interactome,c("Akt1","Foxo1",-1))
  Background_signaling_interactome=rbind(Background_signaling_interactome,c("Akt1","Foxo3",-1))
  Background_signaling_interactome=rbind(Background_signaling_interactome,c("Akt1","Foxo4",-1))
  #adding Igf1r and Akt1 interaction
  Background_signaling_interactome=rbind(Background_signaling_interactome,c("Igf1r","Akt1",1))
  
  
  #loading background signaling interactme
  #' Removing FAS
  if (is.null(FAS)){
    #a=addTfIntBack(Background_signaling_interactome,TF_TF_interactions)
    a=Background_signaling_interactome #in case you want to us the original background interactome
  } else {
    a=Background_signaling_interactomeFAS
  }
  ab=join(a,b,by=c("Source"),type="left",match="first")
  colnames(ab)[3:ncol(ab)]="source"
  colnames(b)[1] <- "Target"
  ab1=join(a,b,by=c("Target"),type="left",match="first")
  names(ab1) <- NULL
  names(ab) <- NULL
  ab=ab[,4:ncol(ab)]
  ab1=ab1[,4:ncol(ab1)]
  ab=as.matrix(ab)
  ab1=as.matrix(ab1)
  ########Elementwise product
  g=ab * ab1
  ########Sum of elementwise product
  sum_product_expression=rowSums(g, na.rm = FALSE, dims = 1)
  g3=cbind(a,sum_product_expression)
  ########Calculation of precentage of cells expressed
  h=rowSums(g != 0)
  percent_expressed=(h*100)/ncol(ab)
  g3=cbind(g3,percent_expressed)
  g3[is.na(g3)]<-0
  #names(g3)=c("V1","V2","Effect","V3","V4")
  #g3[g3$percent_expressed <as.numeric(cutoff) ,3]=0
  #write.table(g3,out, sep="\t",row.names=FALSE)
  ######NETWORK preprocessing
  g <- graph.data.frame(as.data.frame(g3))
  #g=simplify(g,edge.attr.comb=list("first"))
  #DELETE those edges whoese sum is zero
  #del=E(g)[V3==0]
  #DELETE those edges whoese sum is zero (V3) OR whose average expression is lessthan args(3) (V4)
  #convert cutoff to numeric ##NOT DELETING NODES ANYMORE
  del=E(g)[sum_product_expression==0|percent_expressed<as.numeric(cutoff)]
  g <- delete.edges(g,del)
  g<-delete.vertices(g,"Hsp90aa1")
  #SINCE THE TFs AND RECEPTORS ARE ALREADY CONNECTED TO DUMMY, REMOVE ANY NODE THAT HAS ZERO in degree or zero out degree
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
  #SCC <- clusters(g, mode="strong")
  #subg <- induced.subgraph(g, which(membership(SCC) == which.max(sizes(SCC))))
  #subg=simplify(subg,edge.attr.comb=list("first"))
  #subg
  
  if (is.null(FAS)){
    SCC <- clusters(g, mode="strong")
    subg <- induced.subgraph(g, which(membership(SCC) == which.max(sizes(SCC))))
    #subg=simplify(subg,edge.attr.comb=list("first"))
    subg
  } else {
    subg=add_sccDummy(g)
  }
  return(subg)
}

##STATIONARY DISTRIBUTION IN R ITSELF
##make a sparce-matrix from the edgelist
Markov_chain_stationary_distribution <- function(subg) #The function takes the edgelist with probabilitys and computes the SS probability
{
  ####Write this subgraph as edgelist to make normal graph with ids ie.e names=F
  transition_probability=as.data.frame(as_edgelist(subg,names=F))
  transition_probability$probability=paste(E(subg)$sum_product_expression)
  transition_probability[3]=as.numeric(transition_probability[[3]])
  myMatrix = sparseMatrix(i = transition_probability[1:nrow(transition_probability),1], j = transition_probability[1:nrow(transition_probability),2],x = transition_probability[1:nrow(transition_probability),3])
  #Making a stochastic matrix
  myMatrix = (myMatrix)/Matrix::rowSums((myMatrix))
  #write.table(myMatrix,"matrix.txt",sep="\t")
  #print(myMatrix)
  ##Eigen value of the sparse matrix
  #ev=eigen(Matrix::t(myMatrix))
  ##eigen value by Rspectra package
  el=eigs(Matrix::t(myMatrix),1,which="LR")
  ##eigen value that matches 1
  #match(1.0000,Re(round(ev$values, digits = 5)))
  #col=which.min(abs(ev$values - 1.00000000))
  ##STATIONARY DISTRIBUTION
  #SD=(abs(ev$vectors[,col]))/sum(abs(ev$vectors[,col]))
  #Stationary distribution from Rspectra eigen vectors
  SD=(abs(el$vectors))/sum(abs(el$vectors))
  SD=as.data.frame(SD)
  SD
  SD=cbind((as.data.frame(V(subg)$name)),SD)
  SD=as.data.frame(SD)
  colnames(SD)[1] <- "Gene"
  out_SS=paste("Steady_state",sep="")
  colnames(SD)[2] <- out_SS
  SD
}




#----------------------------------------------------------------------------------
# TO CALCULATE THE COMPATABILITY SCORE FOR THE HIGH PROBABILITY INTERMEDIATES
#----------------------------------------------------------------------------------

high_probability_intermediates <- function(x, intermediates, percentile)
{
  intermediates=unique(join(intermediates,x,by=c("Gene"),type="inner"))
  #Selecting top 90 percentile of intermediates with steady state probabilities
  percentile=as.numeric(percentile)/100
  SS_90percentile=as.vector(quantile(intermediates[,2], as.numeric(percentile)))
  ##Shortlisting high-probability intermediates > 90 percentile of SS probability
  int=as.vector((subset(intermediates, intermediates[2] > SS_90percentile , select=c("Gene")))$Gene)
  int
}


#Function for integrating signaling and TF networks, g=signaling graph, x steady state vector, deg=differentially expressed genes
integrate_sig_TF <- function(g,x,deg, non_interface_TFs, TF_TF_interactions )
{
  #delete sccDummy inthe network #only for FAS removal
  #g=delete.vertices(g,"sccDummy")
  el=as_edgelist(g)
  graph_markov=as.data.frame(cbind(el,E(g)$Effect))
  colnames(graph_markov)=c("Source","Target","Effect")
  colnames(deg)=c("Gene","Bool_DEG")
  non_interface_DE_TFs=join(deg,non_interface_TFs,by=c("Gene"),type="inner")
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



#Function of shortlisting non-terminal differentially expressed genes
nonterminal_DE_TFs <- function(g,deg,non_interface_TFs)
{
  colnames(deg)=c("Gene","Bool_DEG")
  non_interface_DE_TFs=join(deg,non_interface_TFs,by=c("Gene"),type="inner") #IF THIS IS ZERO NEED TO ABORT
  #load non-terminal TFs
  nTF=non_interface_DE_TFs[1]
  names(nTF)=NULL
  nTF=as.vector(t(nTF))
  nTF<-intersect(nTF,V(g)$name) #Some TFs must still be missing in the final g3
  if (length(nTF) == 0) stop ('No downstream TF found for the cutoff employed. You may decrease the cutoff in order to find TFs.') # decrease percentile cutoff
  nTF
}

#' Define interface_TF This is done after the data is loaded
load("C:/Users/olive/OneDrive/Desktop/MYCELIUM_BIO/mycelial_scaffold/models/signalingfactorscd34positive-master/MOUSE_Background_Network_omnipath_reactome_metacore_01042019.RData")
interface_TFs<-as.data.frame(setdiff(V(graph.data.frame(TF_TF_interactions))$name,as.character(as.vector(non_interface_TFs)$Gene)))
names(interface_TFs)<-"Gene"

#Function of shortlisting non-terminal differentially expressed genes But without stop command for building networks
#Function for classifying nTF as up or down regulated
up_down_tfs <- function(g,deg,non_interface_TFs)
{
  colnames(deg)=c("Gene","Bool_DEG")
  non_interface_DE_TFs=join(deg,non_interface_TFs,by=c("Gene"),type="inner") #IF THIS IS ZERO NEED TO ABORT
  #load non-terminal TFs
  nTF=non_interface_DE_TFs[1]
  names(nTF)=NULL
  nTF=as.vector(t(nTF))
  nTF<-intersect(nTF,V(g)$name) #Some TFs must still be missing in the final g3
  nTF
}

##Path_weight
#product_path_weight <- function(path, graph) sum(E(graph, path=path)$weight)/length(path)
product_path_weight <- function(path, graph)
{
  edge_weights=E(graph, path=path)$weight
  v_Edge_weight=sum(edge_weights)
  if (v_Edge_weight==0) {p_edge_weights=0} else
  {p_edge_weights=prod(edge_weights[edge_weights!=0])}
}

path_weights <- function(path, graph) (E(graph, path=path)$weight)

##--Function for spliting and taking product of shortest path res file---
##This function takes ONE shortest path (x) and the adjacency matrix (l) as as input, and returns the product of SS probability of the intermediates in that shortest path.


##Function for Compatability score
##This function takes s="source" (one source gene), t="target" (a vector of target gene),g="graph", l="adjacency matrix" as input and finds the shortest paths and passes the argument to spsplit.
##Then it gets the product of intermediates of each path for a source and all its targets and returs its product as final output.

spcal_path_weight <- function(s,t,g){
  if (length(s) == 0) stop ('No intermediates found. You may decrease the percentile in order to find intermediates.') # decrease percentile cutoff
  if (length(t) == 0) stop ('No non-terminal differentially expressed TFs found. You may decrease the cutoff.') # decrease normal cutoff
  paths=(get.all.shortest.paths(g, s, t, mode = c("out"), weights=NA)$res)
  #paths=k_shortest_yen(g,s,t,10)
  if (length(paths) == 0) stop ('No shortest path found. You may decrease the cutoff in order to find shortest path.') # decrease normal cutoff
 # if (length(paths) == 0){
 #   return (0)
 # }
 # To remove NULL paths
  #paths = paths[-which(sapply(paths, is.null))]
  # old path weight estimate
  #weight=lapply(paths,product_path_weight,g)
  # using the updated sig_trn integration
  weightcomp=lapply(paths,product_path_weight1,g)
  weight=lapply(weightcomp, function(xl) xl$p_edge_weights)
  #compatability=lapply(weightcomp, function(xl) xl$SPcompatability)
  cp=c(unlist(lapply(weightcomp, function(xl) xl$compatiblePaths)))
  icp=c(unlist(lapply(weightcomp, function(xl) xl$incompatiblePaths)))
  #s=skewness(unlist(weight))
  s=weight_probability(unlist(weight))
  #s=sum(unlist(weight))
  out=list(weight=s,AllcompatiblePaths=cp,AllincompatiblePaths=icp)
  return(out)
}

#' @Function to calculate activation probability
#' @param x: list of weights and compatability paths obtained from spcal_path_weight
  calActProb <- function(x){
    return(mean(unlist(c(x[1,]))))
    }

#' @Function to calculate significance based on hypergeometric test
#' @param x: list of weights and compatability paths obtained from spcal_path_weight
#' @param m: total number of compatible shortest paths
#' @param n: total number of incompatible shortest paths
  calSig<-function(x,xlist){
    cpSum<-sum(unlist(c(x[2,]))) # q i.e. number of white balls
    icpSum<-sum(unlist(c(x[3,]))) # number of black balls drawn
    totalcp<-sum(unlist(lapply(xlist, function(y) sum(unlist(y[2,]))))) # m total number of white balls in the urn
    totalIcp<-sum(unlist(lapply(xlist, function(y) sum(unlist(y[3,]))))) # n total number of black balls in the urn
    k<-cpSum+icpSum # number of balls drawn from the urn
    if (mean(unlist(c(x[1,])))>0.5) {
      p.value <-  phyper(q=cpSum-1, m=totalcp, n=totalIcp, k=cpSum+icpSum, lower.tail=F)
                         } else {
      p.value <-  phyper(q=cpSum-1, m=totalcp, n=totalIcp, k=cpSum+icpSum, lower.tail=T)
    }
      
    #p.value <-  phyper(q=cpSum-1, m=totalcp, n=totalIcp, k=cpSum+icpSum, lower.tail=T)
    return(p.value)
  }
  
#Function to parallelize the Compatability Score calculations
parallel.function.comp <- function(i){
  mclapply(i,spcal,nTF,g3,l)
}

#FUNCTION FOR SKEWNESS # FROM MOMENTS PACKAGE
skewness <- function (x, na.rm = FALSE)
{
  if (is.matrix(x))
    apply(x, 2, skewness, na.rm = na.rm)
  else if (is.vector(x)) {
    if (na.rm)
      x <- x[!is.na(x)]
    n <- length(x)
    (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
  }
  else if (is.data.frame(x))
    sapply(x, skewness, na.rm = na.rm)
  else skewness(as.vector(x), na.rm = na.rm)
}

weight_probability <- function(x)
{
  x=unlist(x)
  x_pos=x[x>0]
  x_neg=x[x<0]
  x_tot=sum(abs(x_pos),abs(x_neg))
  #probability of the intermediate to be compatible: closer to 1 more compatible it is and closer t zero more incompatible it is
  #p_neg=sum(x_neg)/x_tot
  p_pos=sum(x_pos)/x_tot
}

#function to compare actual probabilities of compatible and incompatible shortest paths
comp_score_tf <- function(t,s,g) #x is a list of comp score
{
  comp_score=lapply(s,spcal_path_weight,t,g)
}

compatability_score <- function(x,y,z,int) #where x is compatability score as list and y is steady state, z is the p-value
{
  x=unlist(x)
  z=unlist(z)
  x=cbind(as.data.frame(int),as.data.frame(unlist(x)),as.data.frame(z))
  colnames(x)=c("Gene", "Activation_probability","pval")
  x=join(x,y,by=c("Gene"),type="inner")
  x=x[order(x$Activation_probability,decreasing = TRUE),]
}
