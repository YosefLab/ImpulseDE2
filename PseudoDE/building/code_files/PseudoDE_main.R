source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/ImpulseDE2_main.R")
library(compiler)
library(parallel)
library(DESeq2)
library(BiocParallel)

library(scone)

source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/srcImpulseDE2_CostFunctionsFit.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/srcImpulseDE2_runDESeq2.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/clusterCellsInPseudotime.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/formatDataClusters.R")

runPseudoDE <- function(matCounts,vecPseudotime ){
  
  # 1. Data preprocessing
  print("1. Data preprocessing:")
  tm_preproc <- system.time({
    # Take out NA cells
    vecPseudotime <- vecPseudotime[!is.na(vecPseudotime)]
    # Adjust ording of cells in objects
    matCounts <- matCounts[,names(vecPseudotime)]
    # Remove all zero genes: Mu is initialised
    # as log(sum counts)
    matCounts <- matCounts[apply(matCounts,1,function(gene){any(gene!=0)}),]
  })
  print("DONE")
  print(paste("Consumed time: ",round(tm_preproc["elapsed"]/60,2),
    " min",sep=""))
  print("##################################################################")
  
  # 2. Cluster cells in pseudo-time
  print("2. Clustering:")
  tm_clustering <- system.time({
    lsResultsClustering <- clusterCellsInPseudotime(vecPseudotime=vecPseudotime)
    dfAnnotationClusters <- formatDataClusters(matCounts=matCounts,
      vecPseudotime=vecPseudotime,lsResultsClustering=lsResultsClustering)
  })
  print("DONE")
  print(paste("Consumed time: ",round(tm_clustering["elapsed"]/60,2),
    " min",sep=""))
  print("##################################################################")
  
  # 3. Fit mixture model
  print("3. Fit mixture model:")
  tm_fitmm <- system.time({
    # Fit zero-inflated negative binomial model to each cluster
    MAXITER <- 10
    # Imputed count data matrix
    matCountsImputed <- array(0,c(dim(matCounts)[1],dim(matCounts)[2]))
    # Genes which are all 0 in one cluster recive pi=1, i.e. they
    # drop out of the cost function during fitting as their likelihood
    # is defined by the hyperparameter pi=1 independently of the
    # mean model.
    matDropout <- array(1,c(dim(matCounts)[1],dim(matCounts)[2]))
    # Negative binomial over-dispersion factor: 1 per gene per cluster
    matDispersion <- array(NA,c(dim(matCounts)[1],lsResultsClustering$K))
    for(k in 1:lsResultsClustering$K){
      # Remove all zero genes: Mu is initialised
      # as log(sum counts)
      vecidxCluster <- lsResultsClustering$Assignments==k
      matCountsCluster <- matCounts[,vecidxCluster]
      vecboolNonzeroGenes <- apply(matCountsCluster,1,
        function(gene){any(gene!=0)})
      matCountsCluster <- matCountsCluster[vecboolNonzeroGenes,]
      # Fit zinb model
      lsZINBparam <- estimate_zinb(
        Y = matCountsCluster[vecboolNonzeroGenes,], 
        maxiter = 10, 
        verbose = TRUE)
      # Record parameters
      matDropout[vecboolNonzeroGenes,vecidxCluster] <- lsZINBparam$pi
      matDispersion[vecboolNonzeroGenes,k] <- lsZINBparam$theta
      # Impute count data matrix
      w <- lsZINBparam$p_z
      matCountsClusterImputed <- matCountsCluster * w + lsZINBparam$mu * (1 - w)
      matCountsImputed[vecboolNonzeroGenes,vecidxCluster] <- matCountsClusterImputed
    }
  })
  print("DONE")
  print(paste("Consumed time: ",round(tm_fitmm["elapsed"]/60,2),
    " min",sep=""))
  print("##################################################################")
  
  # 4. Differential expression analysis on clusters
  print("3. Differential expression analysis on cluster:")
  tm_deanalysis <- system.time({  
    lsDEresults <- runImpulseDE2(
      matCountData = matCounts, 
      dfAnnotationFull = dfAnnotationClusters,
      strCaseName = "case", 
      strControlName = NULL, 
      strMode = "batch", 
      nProc = 3, 
      Q_value=0.01,
      boolPlotting=TRUE)
  })
  print("DONE")
  print(paste("Consumed time: ",round(tm_deanalysis["elapsed"]/60,2),
    " min",sep=""))
  print("##################################################################")
  
  print("Completed PseudoDE.")
  return(lsDEresults)
}