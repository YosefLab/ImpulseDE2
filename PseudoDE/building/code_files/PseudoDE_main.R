source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/ImpulseDE2_main.R")
library(compiler)
library(parallel)
library(DESeq2)
library(BiocParallel)
library(MASS)

library(scone)

source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/srcImpulseDE2_CostFunctionsFit.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/srcImpulseDE2_runDESeq2.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/clusterCellsInPseudotime.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/formatDataClusters.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/srcSCONE_zinb.R")

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
    for(k in 1:lsResultsClustering$K){
      # Remove all zero genes: Mu is initialised
      # as log(sum counts)
      matCountsCluster <- matCounts[1:500,lsResultsClustering$Assignments==k]
      matCountsCluster <- matCountsCluster[apply(matCountsCluster,1,
        function(gene){any(gene!=0)}),]
      source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/srcSCONE_zinb.R")
      lsZINBparam <- estimate_zinb_copy(
        Y = matCountsCluster, 
        maxiter = 10, 
        verbose = TRUE)
      lsSCONEout <- scone(matCountsCluster, 
        imputation=identity, 
        scaling=identity, 
        k_ruv=0, 
        k_qc=0,
        evaluate=TRUE, 
        run=TRUE)
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