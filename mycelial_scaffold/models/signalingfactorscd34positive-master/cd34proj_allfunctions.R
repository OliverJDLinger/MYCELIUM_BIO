#' Extract most relevant columns i.e ctrl or EAE and cluster information from "networkAnnotableFile"
#' Funtion to extract relevant cells from the "networkExpressionFile" based on the required cell type
#' @param x: networkAnnotationFile, y: networkExpressionFile, z: is a string representing the clusterName, g: is a string  representing "Ctrl" or "EAE" or "l2" if it is a level 2 clustering
extractClusterExp <- function(x,y,z,g){
  if (g=="l2") {
    #' Extract rownames
    rn=row.names(x[x$ReorderedViolin==z,])
  }
  else {
    rn=row.names(x[x$ReorderedViolin==z&x$Group==g,])}
  clu=y[,c(rn)]
  return(clu)
}

#' Function to combine user defined cell subpopulations
#' @param x:name for the combined cluster, ... list of elements to be clustered
combineClusters <- function(...){
  x=cbind(...)
  return(x)
}

#' Function for differential expression
#' Perform a differential expression analysis
#'
#' @param df data frame with the cells and target or background as a type
#' @param m single cell count matrix
#'
#' @return list of DE genes with padj < 0.05
#'
runDE <- function(df,m){
  #Construction of the summarizedExperiment object
  #Object with the matrix containing the genes expressed in the 3 subpopulations
  se <- SummarizedExperiment(assays = list(counts = m), colData=df) 
  zinb <- zinbFit(se, K=2, epsilon=1e12) 
  seZinb <- zinbwave(se, fitted_model = zinb, K = 2, epsilon=1e12)
  weights <- assay(seZinb, "weights")
  
  #DESeq2
  dds <- DESeqDataSet(seZinb, design = ~ type)
  dds$type <- relevel(dds$type, ref = "background")
  
  #LRT adapt? aux UMI counts - in fact WALD test one-tailed without LFshrinkage (off thanks to LRT as model to run)
  #poscounts adapte? au single-cell
  dds <- DESeq(dds, parallel =TRUE, sfType="poscounts", minmu=1e-6, test = "LRT", full = ~ type, reduced = ~ 1)
  res <- results(dds, alpha = 0.05, lfcThreshold = 0, altHypothesis = "greater", test = "Wald")
  dfRes <- data.frame("Genes"=res@rownames, "pValue"=res@listData$pvalue, "logFC"=res@listData$log2FoldChange, "padj"= res@listData$padj)
  
  return(dfRes$Genes[which(dfRes$padj<0.05)])
}

#' Perform a differential expression analysis DEFAULT
#'
#' @param x count matrix of treatment
#' @param y count matrix of target
#'
#' @return list of DE genes with padj < 0.05
#'
runDE_S <- function(x,y){
  #Construction of the summarizedExperiment object
  #Object with the matrix containing the genes expressed in the 3 subpopulations
  m=cbind(x,y)
  m=round(m)
  #filter cells not expressed in a given number of cells
  m=filter_exp(m,5,F)
  df<-data.frame('Cell'=colnames(m),'type'=rep(NA,length(colnames(m))))
  df$type[df$Cell %in% colnames(x)] <- "target"
  df$type[is.na(df$type)] <- "background"
  
  se <- SummarizedExperiment(assays = list(counts = m), colData=df) 
  zinb <- zinbFit(se, K=2, epsilon=1e12) 
  seZinb <- zinbwave(se, fitted_model = zinb, K = 2, epsilon=1e12)
  weights <- assay(seZinb, "weights")
  
  #DESeq2
  dds <- DESeqDataSet(seZinb, design = ~ type)
  dds$type <- relevel(dds$type, ref = "background")
  
  #LRT adapt? aux UMI counts - in fact WALD test one-tailed without LFshrinkage (off thanks to LRT as model to run)
  #poscounts adapte? au single-cell
  dds <- DESeq(dds, parallel =F, sfType="poscounts", minmu=1e-6, test = "LRT", full = ~ type, reduced = ~ 1)
  res <- results(dds, alpha = 0.05, lfcThreshold = 0, altHypothesis = "greater", test = "Wald")
  dfRes <- data.frame("Genes"=res@rownames, "pValue"=res@listData$pvalue, "logFC"=res@listData$log2FoldChange, "padj"= res@listData$padj)
  
  return(dfRes$Genes[which(dfRes$padj<0.05)])
}

pvalue_cal<- function(idx,a,b)
{
  l=t.test(a[idx,2:ncol(a)],b[idx,2:ncol(b)])
  mean=as.numeric(l$estimate)
  fc=log2(mean[1]/mean[2])
  list(as.data.frame(cbind(idx, mean[1], mean[2], fc, l$p.value)))
}

output_pval <- function(pval_list,idx)
{
  ll=join_all((pval_list),by=c("idx"),type="full")
  lol=cbind(t(t(names(idx))),ll[2:5])
  names(lol)=c("Gene","mean_population1","mean_population2","Log_Fold_change","p_value")
  adj_p=p.adjust(lol$`p_value`, method = "fdr")
  lol=cbind(lol,adj_p)
}

DEG_pvalue_singlecell <- function(a,b)
{
  idx <- seq_len(nrow(a))
  names(idx) <- (a[,1])
  ttest_p <- sapply(idx,pvalue_cal,a,b)
  pp=output_pval(ttest_p,idx)
  names(pp)=c("Gene","mean_population1","mean_population2","Log_Fold_change","p_value","adj_p")
  return(pp)
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
  x[((rowSums(x != 0))*100)/ncol(x)>y,]
}

#' Perform a differential expression analysis
#'
#' @param df data frame with the cells and target or background as a type
#' @param m single cell count matrix
#'
#' @return list of DE genes with padj < 0.05
#'
runDE <- function(df,m){
  #Construction of the summarizedExperiment object
  #Object with the matrix containing the genes expressed in the 3 subpopulations
  se <- SummarizedExperiment(assays = list(counts = m), colData=df) 
  zinb <- zinbFit(se, K=2, epsilon=1e12) 
  seZinb <- zinbwave(se, fitted_model = zinb, K = 2, epsilon=1e12)
  weights <- assay(seZinb, "weights")
  
  #DESeq2
  dds <- DESeqDataSet(seZinb, design = ~ type)
  dds$type <- relevel(dds$type, ref = "background")
  
  #LRT adapt? aux UMI counts - in fact WALD test one-tailed without LFshrinkage (off thanks to LRT as model to run)
  #poscounts adapte? au single-cell
  dds <- DESeq(dds, parallel =TRUE, sfType="poscounts", minmu=1e-6, test = "LRT", full = ~ type, reduced = ~ 1)
  res <- results(dds, alpha = 0.05, lfcThreshold = 0, altHypothesis = "greater", test = "Wald")
  dfRes <- data.frame("Genes"=res@rownames, "pValue"=res@listData$pvalue, "logFC"=res@listData$log2FoldChange, "padj"= res@listData$padj)
  
  return(dfRes$Genes[which(dfRes$padj<0.05)])
}

#' Perform a differential expression analysis
#'
#' @param x count matrix of treatment
#' @param y count matrix of target
#'
#' @return list of DE genes with padj < 0.05
#'
runDE_S <- function(x,y){
  #Construction of the summarizedExperiment object
  #Object with the matrix containing the genes expressed in the 3 subpopulations
  m=cbind(x,y)
  m=round(m)
  #filter genes not expressed in a given number of cells
  m=filter_exp(m,5,F)
  #filter cells with zero coults
  m = m[,colSums(m) > 0]
  #filter only Tfs
  m=filter_tfs_r(m)
  df<-data.frame('Cell'=colnames(m),'type'=rep(NA,length(colnames(m))))
  df$type[df$Cell %in% colnames(x)] <- "target"
  df$type[is.na(df$type)] <- "background"
  
  se <- SummarizedExperiment(assays = list(counts = m), colData=df) 
  zinb <- zinbFit(se, K=2, epsilon=1e12) 
  seZinb <- zinbwave(se, fitted_model = zinb, K = 2, epsilon=1e12)
  weights <- assay(seZinb, "weights")
  
  #DESeq2
  dds <- DESeqDataSet(seZinb, design = ~ type)
  dds$type <- relevel(dds$type, ref = "background")
  
  #LRT adapt? aux UMI counts - in fact WALD test one-tailed without LFshrinkage (off thanks to LRT as model to run)
  #poscounts adapte? au single-cell
  dds <- DESeq(dds, parallel =F, sfType="poscounts", minmu=1e-6, test = "LRT", full = ~ type, reduced = ~ 1)
  res <- results(dds, alpha = 0.05, lfcThreshold = 0, altHypothesis = "greater", test = "Wald")
  dfRes <- data.frame("Genes"=res@rownames, "pValue"=res@listData$pvalue, "logFC"=res@listData$log2FoldChange, "padj"= res@listData$padj)
  return(dfRes)
  #return(dfRes$Genes[which(dfRes$padj<0.05)])
}

#' function for Differentiall expression output with names
#' It identifies on differentially expressed TFs
deBoolOutTfs <- function(data1,data2){
  data1<-filter_tfs(data1)
  data2<-filter_tfs(data2)
  options(warn = -1)
  v=DEG_pvalue_singlecell(data1,data2)
  #row.names(v)=row.names(data1)
  #v[1]=row.names(v)
  sub_v=subset(v, adj_p <0.1& (Log_Fold_change > 1 | Log_Fold_change < -1))
  bool=sub_v[,c("Gene","Log_Fold_change")]
  bool$Log_Fold_change[bool$Log_Fold_change>0]<-1
  bool$Log_Fold_change[bool$Log_Fold_change<0]<--1
  names(bool)=c("Gene","Boolean_DE")
  #row.names(bool)<-NULL
  return(bool)
  options(warn = getOption("warn"))
}

#' Filter only TFs
#'
#' @param x single cell RNA seq matrix (first column header must be "Gene") | genes in row, cells in column
#' @param org organism - human (default) or mouse
#'
#' @return filtered matrix (only TFs)
#'
filter_tfs <- function(x,Gene = F){
  tfs <- read.table("Mus_musculus_TF.txt",sep="\t",quote="",header=T)
  colnames(tfs)[2] <- "Gene"
  tfs <- unique(tfs[2])
  x=x[x$Gene %in% as.vector(tfs$Gene),]
  return(x)
}

#' Filter only TFs
#'
#' @param x single cell RNA seq matrix (first column header must be "Gene") | genes in row, cells in column
#' @param org organism - human (default) or mouse
#'
#' @return filtered matrix (only TFs)
#'
filter_tfs_r <- function(x){
  tfs <- read.table("Mus_musculus_TF.txt",sep="\t",quote="",header=T)
  colnames(tfs)[2] <- "Gene"
  tfs <- unique(tfs[2])
  x=x[rownames(x) %in% as.vector(tfs$Gene),]
  return(x)
}

#' function to run sighotspotter
#' @param phenoData1,phenoData2 the two input phenotypes
#' @param cutoff expression cutoff (default 20 for scRNAseq and 100 for bulk)
#' @param percentile is the percentile cutoff (default 90, reduce to 70 if no TFs found)
#' @param FAS T is Feedback arc set is to be removed
#' @param acc is the DE data or ACC data, if NULL t-test is used to calculate DEG
runSigHotSpotter <- function(phenoData1,phenoData2,cutoff,percentile, FAS=NULL,acc=NULL){
  #converting the input file into suitable format
  if (colnames(phenoData1)[1]=="Gene"){
    phenoData1=phenoData1
    phenoData2=phenoData2
  } else {
  phenoData1=as.data.frame(setDT(as.data.frame(phenoData1),keep.rownames = "Gene"))
  phenoData2=as.data.frame(setDT(as.data.frame(phenoData2),keep.rownames = "Gene"))
  }
  if (is.null(acc))
  {
    deTfs=deBoolOutTfs(phenoData1,phenoData2)
  } else {
    deTfs=acc
  }
  # Running SigHotSpotter
  results_1 <<- SigHotSpotter_standalone ("MOUSE", phenoData1, cutoff, deTfs, percentile, invert_DE = FALSE, FAS)
  results_2 <<- SigHotSpotter_standalone ("MOUSE", phenoData2, cutoff, deTfs, percentile, invert_DE = TRUE, FAS)
  final_results <<- list(results_1,results_2)
  return(final_results)
}

#' function for running sighotspotter by considering the accessibility information
#' @param phenoData1 data1
#' @tfExpData..... 
#' @param acc is booleanized differential accessability
runSigHotSpotter_zsc <- function(phenoData1,tfExpData,zscCutoff,cutoff,percentile){
  #converting the input file into suitable format
  phenoData1=as.data.frame(setDT(as.data.frame(phenoData1),keep.rownames = "Gene"))
  deTfs=DeTfsZscore(tfExpData,zscCutoff)
  # Running SigHotSpotter
  final_results <<- SigHotSpotter_standalone ("MOUSE", phenoData1, cutoff, deTfs, percentile, invert_DE = FALSE)
  return(final_results)
}

SigHotSpotter_standalone <- function(species,idata,cutoff,DE_Genes,percentile,invert_DE = FALSE,FAS=NULL){
  ## Choose correct dataset according to species
  if(species == "MOUSE"){
    load("C:/Users/olive/OneDrive/Desktop/MYCELIUM_BIO/mycelial_scaffold/models/signalingfactorscd34positive-master/MOUSE_Background_Network_omnipath_reactome_metacore_01042019.RData")
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
  
  subg=Data_preprocessing(idata,cutoff,species,FAS)
  
  ## Calculate stationary distribution of the MC
  Steady_state_true=Markov_chain_stationary_distribution(subg)
  
  ## Retrieves high probability intermediates
  int=high_probability_intermediates(Steady_state_true, intermediates, percentile)
  ## Retrieves high probability intTFs
    tflist<-as.data.frame(tflist$Symbol)
    names(tflist)<-"Gene"
    intTF<-high_probability_intermediates(Steady_state_true, tflist, percentile)
    #gintg=integrate_sig_TF(subg,Steady_state_true,DE_Genes, non_interface_TFs, TF_TF_interactions )
    #nTF=nonterminal_DE_TFs(gintg,DE_Genes,non_interface_TFs)
  #' integrate Signalinh and TF to cover more nonintDETFs
  gintg_old=integrate_sig_TFGeneral(subg,Steady_state_true,DE_Genes, non_interface_TFs, TF_TF_interactions) # this is just to include the other missing TFs in the network
  trng=buildTrn(DE_Genes,idata,cutoff) # this is to find paths via genes or TFs like Foxos
  gintg=sigTrnIntegratedNet(subg,trng,Steady_state_true,DE_Genes) # this is to find paths via genes or TFs like Foxos i.e. more than just downstream TF of interface TFs
  #' to find more nTFs
  non_interface_DE_TFs <- nonIntDeTfs(subg,DE_Genes)
  nTF=intersect(V(gintg)$name,as.vector(t(non_interface_DE_TFs$Gene)))
  nTF=RmNonreachNtFs(gintg,int,nTF)
  # use only for treating Foxos as nTF
  #nTF=c("Foxo1","Foxo3","Foxo4")
  #comp_score for each TF individually
  #scoreCompatability=lapply(nTF,comp_score_tf, int, gintg)
  #score=lapply(nTF,comp_score_tf, int, gintg)
  scorearray=sapply(nTF,comp_score_tf, int, gintg,simplify = "array")
  score=apply(scorearray,1,simplify2array)
  score_m_means=(sapply(score,calActProb))
  pVals<-sapply(score,calSig,score)
  #out=as.data.frame(cbind(int,sapply(score,calActProb)))
  #colnames(out)<-c("Gene","Activation_probability")
  #final_score=join(out,Steady_state_true,by=c("Gene"),type="left")
  
  #score=lapply(scoreCompatability, function(xl) xl$weight)
  #converting the nested list into a matrix whose row sum will give probability of each intermediate
  #score_m=(matrix(unlist(score), ncol=length(score), byrow=F))
  #score_m_means=as.list(rowMeans(score_m))
  final_score=compatability_score(score_m_means,Steady_state_true,pVals,int)
  
  ## Computing networks for visualization
  trimmed_score_I=.trimResults(final_score,F)
  trimmed_score_A=.trimResults(final_score,T)
  toiintI=c(as.matrix(trimmed_score_I$`Inactive signaling hotspots`))
  toiintA=c(as.matrix(trimmed_score_A$`Active signaling hotspots`))
  twoi=c(toiintI,toiintA)
  
  #pruning the integrated networks
  gintg.p=prun.int.g(gintg)
  #pruning and reversing the edges
  #gintg.p=revertEdgesPrun(gintg,DE_Genes)
  
  #building networks for all intermediates for active signaling hotspots
  #sp_int_A <- lapply(toiintA,to_sp_net_int,gintg.p,nTF,DE_Genes,non_interface_TFs)
  sp_int_A <- lapply(toiintA,to_sp_net_int,gintg.p,nTF,non_interface_DE_TFs,non_interface_DE_TFs[1])
  
  #building networks for inactive signaling hotspots
  #sp_int_I <- lapply(toiintI,to_sp_net_int,gintg.p,nTF,DE_Genes,non_interface_TFs)
  sp_int_I <- lapply(toiintI,to_sp_net_int,gintg.p,nTF,non_interface_DE_TFs,non_interface_DE_TFs[1])
  
  #converting to visNetwork
  vis_net_A <- lapply(sp_int_A,toVisNetworkData)
  vis_net_I <- lapply(sp_int_I,toVisNetworkData)
  
  #for edge color
  vis_net_A <- lapply(vis_net_A,vis.edge.color)
  vis_net_I <- lapply(vis_net_I,vis.edge.color)
  ret_value = list(subg=subg,trng=trng,gintg=gintg,gintg_old=gintg_old,int=int,intTF=intTF,nTF=nTF,subg=subg,SS=Steady_state_true,final_score=final_score,
                   trimmed_score_A=trimmed_score_A,
                   trimmed_score_I=trimmed_score_I,
                   vis_net_A=vis_net_A,
                   vis_net_I=vis_net_I )
  return(ret_value)
}

#' function to run sighotspotter for Foxo
#' @param phenoData1,phenoData2 the two input phenotypes
#' @param cutoff expression cutoff (default 20 for scRNAseq and 100 for bulk)
#' @param percentile is the percentile cutoff (default 90, reduce to 70 if no TFs found)
#' @param FAS T is Feedback arc set is to be removed
#' @param acc is the DE data or ACC data, if NULL t-test is used to calculate DEG
runSigHotSpotterFoxo <- function(phenoData1,phenoData2,cutoff,percentile, FAS=NULL,acc=NULL){
  #converting the input file into suitable format
  if (colnames(phenoData1)[1]=="Gene"){
    phenoData1=phenoData1
    phenoData2=phenoData2
  } else {
    phenoData1=as.data.frame(setDT(as.data.frame(phenoData1),keep.rownames = "Gene"))
    phenoData2=as.data.frame(setDT(as.data.frame(phenoData2),keep.rownames = "Gene"))
  }
  if (is.null(acc))
  {
    deTfs=deBoolOutTfs(phenoData1,phenoData2)
  } else {
    deTfs=acc
  }
  # Running SigHotSpotter
  results_1 <<- SigHotSpotter_standalone_foxo ("MOUSE", phenoData1, cutoff, deTfs, percentile, invert_DE = FALSE, FAS)
  results_2 <<- SigHotSpotter_standalone_foxo ("MOUSE", phenoData2, cutoff, deTfs, percentile, invert_DE = TRUE, FAS)
  final_results <<- list(results_1,results_2)
  return(final_results)
}

SigHotSpotter_standalone_foxo <- function(species,idata,cutoff,DE_Genes,percentile,invert_DE = FALSE,FAS=NULL){
  ## Choose correct dataset according to species
  if(species == "MOUSE"){
    load("C:/Users/olive/OneDrive/Desktop/MYCELIUM_BIO/mycelial_scaffold/models/signalingfactorscd34positive-master/MOUSE_Background_Network_omnipath_reactome_metacore_01042019.RData")
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
  
  subg=Data_preprocessing(idata,cutoff,species,FAS)
  
  ## Calculate stationary distribution of the MC
  Steady_state_true=Markov_chain_stationary_distribution(subg)
  
  ## Retrieves high probability intermediates
  int=high_probability_intermediates(Steady_state_true, intermediates, percentile)
  ## Retrieves high probability intTFs
  tflist<-as.data.frame(tflist$Symbol)
  names(tflist)<-"Gene"
  intTF<-high_probability_intermediates(Steady_state_true, tflist, percentile)
  #gintg=integrate_sig_TF(subg,Steady_state_true,DE_Genes, non_interface_TFs, TF_TF_interactions )
  #nTF=nonterminal_DE_TFs(gintg,DE_Genes,non_interface_TFs)
  #' integrate Signalinh and TF to cover more nonintDETFs
  gintg_old=integrate_sig_TFGeneral(subg,Steady_state_true,DE_Genes, non_interface_TFs, TF_TF_interactions) # this is just to include the other missing TFs in the network
  trng=buildTrn(DE_Genes,idata,cutoff) # this is to find paths via genes or TFs like Foxos
  gintg=sigTrnIntegratedNet(subg,trng,Steady_state_true,DE_Genes) # this is to find paths via genes or TFs like Foxos i.e. more than just downstream TF of interface TFs
  #' to find more nTFs
  #non_interface_DE_TFs <- nonIntDeTfs(subg,DE_Genes)
  non_interface_DE_TFs <- DE_Genes
  nTF=intersect(V(gintg)$name,as.vector(t(non_interface_DE_TFs$Gene)))
  nTF=RmNonreachNtFs(gintg,int,nTF)
  # use only for treating Foxos as nTF
  nTF=c("Foxo1","Foxo3","Foxo4")
  #comp_score for each TF individually
  #scoreCompatability=lapply(nTF,comp_score_tf, int, gintg)
  #score=lapply(nTF,comp_score_tf, int, gintg)
  scorearray=sapply(nTF,comp_score_tf, int, gintg,simplify = "array")
  score=apply(scorearray,1,simplify2array)
  score_m_means=(sapply(score,calActProb))

  pVals<-sapply(score,calSig,score)
  #out=as.data.frame(cbind(int,sapply(score,calActProb)))
  #colnames(out)<-c("Gene","Activation_probability")
  #final_score=join(out,Steady_state_true,by=c("Gene"),type="left")
  
  #score=lapply(scoreCompatability, function(xl) xl$weight)
  #converting the nested list into a matrix whose row sum will give probability of each intermediate
  #score_m=(matrix(unlist(score), ncol=length(score), byrow=F))
  #score_m_means=as.list(rowMeans(score_m))
  final_score=compatability_score(score_m_means,Steady_state_true,pVals,int)
    if(invert_DE)
     {
      final_score[,2]=1-(final_score[,2])
      final_score=final_score[order(final_score$Activation_probability,decreasing = TRUE),]
      }

  ## Computing networks for visualization
  trimmed_score_I=.trimResults(final_score,F)
  trimmed_score_A=.trimResults(final_score,T)
#  if(invert_DE)
 # {
  #  trimmed_score_I[,2] = 1-trimmed_score_I[,2]
   # trimmed_score_A[,2] = 1-trimmed_score_A[,2]
  #}

  toiintI=c(as.matrix(trimmed_score_I$`Inactive signaling hotspots`))
  toiintA=c(as.matrix(trimmed_score_A$`Active signaling hotspots`))
  twoi=c(toiintI,toiintA)
  
  #pruning the integrated networks
  gintg.p=prun.int.g(gintg)
  #pruning and reversing the edges
  #gintg.p=revertEdgesPrun(gintg,DE_Genes)
  
  #building networks for all intermediates for active signaling hotspots
  #sp_int_A <- lapply(toiintA,to_sp_net_int,gintg.p,nTF,DE_Genes,non_interface_TFs)
  sp_int_A <- lapply(toiintA,to_sp_net_int,gintg.p,nTF,non_interface_DE_TFs,non_interface_DE_TFs[1])
  
  #building networks for inactive signaling hotspots
  #sp_int_I <- lapply(toiintI,to_sp_net_int,gintg.p,nTF,DE_Genes,non_interface_TFs)
  sp_int_I <- lapply(toiintI,to_sp_net_int,gintg.p,nTF,non_interface_DE_TFs,non_interface_DE_TFs[1])
  
  #converting to visNetwork
  vis_net_A <- lapply(sp_int_A,toVisNetworkData)
  vis_net_I <- lapply(sp_int_I,toVisNetworkData)
  
  #for edge color
  vis_net_A <- lapply(vis_net_A,vis.edge.color)
  vis_net_I <- lapply(vis_net_I,vis.edge.color)
  ret_value = list(subg=subg,trng=trng,gintg=gintg,gintg_old=gintg_old,int=int,intTF=intTF,nTF=nTF,subg=subg,SS=Steady_state_true,final_score=final_score,
                   trimmed_score_A=trimmed_score_A,
                   trimmed_score_I=trimmed_score_I,
                   vis_net_A=vis_net_A,
                   vis_net_I=vis_net_I )
  return(ret_value)
}


#' Function for extracting Up or down regulated TF for each subpopulation based on z-score
#' @param:x matrix of TF expression levels 
#' @param:y zscore cutoff for defining TFs to be up or downregulated, use 1 as default
DeTfsZscore <- function(x,y){
  m <- rowMeans(x)
  m <- c(subset(m,m<(-y)),subset(m,m>y))
  m[m>0]<-1
  m[m<0]<--1
  m=as.data.frame(m)
  m=cbind(row.names(m),m)
  names(m) <- c("Gene","ZscoreDE")
  return(m)
}

#' Function for reading all sheets in a excel WB
#' @param filename: file to be imported
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#' Function to extract and booleanize the gene expression data
#' @param x,y: expression data, celltype 
#' @return booleanized expression data for input to sighotspotter
cd34deBool<-function(x){
  up<-x[x$log2FoldChange>0,1]
  down<-x[x$log2FoldChange<0,1]
  up<-as.data.frame(cbind(up,1))
  names(up)<-c("Gene","State")
  down<-as.data.frame(cbind(down,-1))
  names(down)<-c("Gene","State")
  deBoolgenes<-data.frame(rbind(up,down))
  deBoolgenes$State=as.numeric(as.character(deBoolgenes$State))
  return(deBoolgenes)
}

#' function to exctract updownFs based on Refbool discritization
#' 
#' 
refboolDisTfs <- function (x){
  x=filter_tfs(x)
  xmean=cbind(x[1],rowMeans(x[2:ncol(x)]))
  xmean[2][xmean[2]>0.5]=1
  xmean[2][xmean[2]<(-0.5)]=-1
  xUD<-subset(xmean,xmean[2]==1|xmean[2]==-1)
  return(xUD)
}


#' function to compare the perturbation effect
#' @param r receptor/intermediate to perturb
#' @param g subg of unperturbed
#' @param SS real SS
perturbCompare<-function(r,g,SS) {
rP=perturb_receptor(r,g)
ssReal=merge(interface_TFs,SS)
ssRp=merge(interface_TFs,rP$Steady_state_true)
merge_ss=merge(ssReal,ssRp,by="Gene")
merge_ss=cbind(merge_ss,merge_ss[3]/merge_ss[2])
return(merge_ss)
}

#' Function to revert back the original edge signs for the interfaceTF-TF interaction
#' @param g: original gintg
#' @param deg: differentially expressed TFs
revertEdgesPrun <- function(g,deg){
el=as_edgelist(g)
graph_markov=as.data.frame(cbind(el,E(g)$Effect))
colnames(graph_markov)=c("Source","Target","Effect")
colnames(deg)=c("Gene","Bool_DEG")
non_interface_DE_TFs=join(deg,non_interface_TFs,by=c("Gene"),type="inner")
#Get the phenotype from the non_interface_DE_TFs and map it to the TF-TF interaction network
colnames(non_interface_DE_TFs)[1]="Target"
graph_markov=join(graph_markov,non_interface_DE_TFs,by=c("Target"),type="left")
graph_markov[is.na(graph_markov)]<-1
ab=graph_markov
#ab$Effect=ab$Effect*ab$Phenotype1 #second Phenotype will be opposite of this #must inpust second data for second phenotype
names(ab)<-NULL
ab=as.matrix(ab)
ab[,3]=as.numeric(ab[,3])*as.numeric(ab[,4])
# adding the original graph weight back to the edge reverted graph
ab=as.data.frame(ab)[1:3]
ab=as.data.frame(cbind(ab,E(g)$weight))
names(ab)=c("Source","Target","Effect","weight")
g=graph.data.frame(ab)
E(g)$weight=1-(abs(E(g)$weight))
#function to delet edges from tf to dummy un the network
del=incident(g, "Dummy", mode = c("in"))
g <- delete.edges(g,del)
g <- set.vertex.attribute(g,"name","Dummy","NICHE")
return(g)
}
