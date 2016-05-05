source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/ImpulseDE2_main.R")
library(compiler)
library(parallel)
library(DESeq2)
library(BiocParallel)
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/srcImpulseDE2_CostFunctionsFit.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/srcImpulseDE2_runDESeq2.R")
source("/Users/davidsebastianfischer/MasterThesis/code/PseudoDE/building/code_files/clusterCellsInPseudotime.R")
source("/Users/davidsebastianfischer/MasterThesis/code/PseudoDE/building/code_files/formatDataClusters.R")
source("/Users/davidsebastianfischer/MasterThesis/code/PseudoDE/building/code_files/fitHurdleModel.R")

evalLogLikHurdleNB_comp <- cmpfun(evalLogLikHurdleNB)
evalLogLikHurdleDrop_comp <- cmpfun(evalLogLikHurdleDrop)
evalLogLikHurdle_comp <- cmpfun(evalLogLikHurdle)

runPseudoDE <- function(matCounts,vecPseudotime){
  
  # 1. Data preprocessing
  print("1. Data preprocessing:")
  tm_preproc <- system.time({
    # Take out NA cells
    vecPseudotime <- vecPseudotime[!is.na(vecPseudotime)]
    # Adjust ording of cells in objects
    matCounts <- matCounts[,names(vecPseudotime)]
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
    lsHurdleParamters <- fitHurdleModel(matCounts=matCounts, 
      lsResultsClustering=lsResultsClustering, dfAnnotationClusters=dfAnnotationClusters)
  })
  print("DONE")
  print(paste("Consumed time: ",round(tm_fitmm["elapsed"]/60,2),
    " min",sep=""))
  print("##################################################################")
  
  # 4. Differential expression analysis on clusters
  print("3. Differential expression analysis on cluster:")
  tm_deanalysis <- system.time({  
    lsDEresults <- runImpulseDE2(matCountData=matCounts, dfAnnotationFull=dfAnnotationImpulseDE2,
      strCaseName = "case", strControlName=NULL, nProc=3, Q_value=0.01)
  })
  print("DONE")
  print(paste("Consumed time: ",round(tm_deanalysis["elapsed"]/60,2),
    " min",sep=""))
  print("##################################################################")
  
  print("Completed PseudoDE.")
  return(lsDEresults)
}