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
  # with only one replicate (itself), unless multiple cells
  # map to exactly the same pseudotime point.
  vecTimepointsUnique <- sort(unique( as.vector(vecPseudotime) ))
  vecSamplesUnique <- paste0(rep("Sample_",length(vecTimepointsUnique)),vecTimepointsUnique)
  vecReplicates <- names(vecPseudotime)
  vecindTimepoints <- match(vecPseudotime,vecTimepointsUnique)
  vecTimepoints <- vecTimepointsUnique[vecindTimepoints]
  vecSamples <- vecSamplesUnique[vecindTimepoints]
  dfAnnotationImpulseDE2 <- as.data.frame(cbind(
    Replicate=vecReplicates,
    Sample=vecSamples,
    Condition="case",
    Time=vecTimepoints ))
  dfAnnotationImpulseDE2$Time <- as.numeric(as.vector( dfAnnotationImpulseDE2$Time ))
  
  return(dfAnnotationImpulseDE2)
}