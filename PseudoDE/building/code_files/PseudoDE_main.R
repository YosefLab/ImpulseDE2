source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/ImpulseDE2_main.R")
source("/Users/davidsebastianfischer/MasterThesis/code/PseudoDE/building/code_files/clusterCellsInPseudotime.R")

runPseudoDE <- function(matCounts,vecPseudotime){
  
  # 1. Data preprocessing
  print("1. Data preprocessing:")
  tm_preproc <- system.time({
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
    lsResultsClustering <- clusterCellsInPseudotime(matCounts=matCounts, 
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