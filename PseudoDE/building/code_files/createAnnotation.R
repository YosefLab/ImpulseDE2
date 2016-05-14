createAnnotationByCluster <- function(matCounts,vecPseudotime,lsResultsClustering){
  
  # Create annotataion table treating a cluster as a sample
  # and its members as replicates
  vecTimepoints <- lsResultsClustering$Centroids[lsResultsClustering$Assignments]
  vecSamples <- ( paste0(rep("Cluster_"),seq(1:lsResultsClustering$K)) )[
    lsResultsClustering$Assignments]
  dfAnnotationImpulseDE2 <- as.data.frame(cbind(
    Replicate=names(vecPseudotime),
    Sample=vecSamples,
    Condition="case",
    Time=vecTimepoints ))
  dfAnnotationImpulseDE2$Time <- as.numeric(as.vector( dfAnnotationImpulseDE2$Time ))
  
  return(dfAnnotationImpulseDE2)
}

createAnnotationByCell <- function(matCounts, vecPseudotime){
  
  # Create annotataion table treating a cell as a sample
  # with only one replicate (itself)
  vecTimepoints <- as.vector(vecPseudotime)
  vecSamples <- paste0(rep("Cell_"),seq(1:length(vecPseudotime)))
  dfAnnotationImpulseDE2 <- as.data.frame(cbind(
    Replicate=names(vecPseudotime),
    Sample=vecSamples,
    Condition="case",
    Time=vecTimepoints ))
  dfAnnotationImpulseDE2$Time <- as.numeric(as.vector( dfAnnotationImpulseDE2$Time ))
  
  return(dfAnnotationImpulseDE2)
}