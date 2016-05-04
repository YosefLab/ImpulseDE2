formatDataClusters <- function(matCounts,vecPseudotime,lsResultsClustering){
  
  # Create annotataion table treating a cluster as a sample
  # and its members as replicates
  vecTimepoints <- lsResultsClustering$Centroids[lsResultsClustering$Assignments]
  vecSampleNames <- paste0(rep("cluster_"),seq(1:lsResultsClustering$K))
  vecSamples <- vecSampleNames[lsResultsClustering$Assignments]
  dfAnnotationImpulseDE2 <- as.data.frame(cbind(
    names(vecPseudotime),
    vecSamples,
    "case",
    vecTimepoints ))
  colnames(dfAnnotationImpulseDE2) <- c("Replicate","Sample","Condition","Time")
  dfAnnotationImpulseDE2$Time <- as.numeric(as.vector( dfAnnotationImpulseDE2$Time ))
  
  return(dfAnnotationImpulseDE2)
}