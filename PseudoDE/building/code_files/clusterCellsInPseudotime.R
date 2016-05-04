clusterCellsInPseudotime <- function(vecPseudotime){
  
  # Range of K for K-means is number of unique cells
  vecPseudotimeUnique <- unique(vecPseudotime)
  scaK <- length(vecPseudotimeUnique)
  vecWithinSS <- array(NA,scaK)
  lsCentroids <- list()
  matAssignments <- array(NA,c(scaK,length(vecPseudotime)))
  # Loop over range of K
  for(k in 1:scaK){
    if(k==1){
      lsKmeansResults <- kmeans(x=vecPseudotime,centers=1)
    } else {
      lsKmeansResults <- kmeans(x=vecPseudotime,centers=vecPseudotimeUnique[
        round(seq(1,length(vecPseudotimeUnique),by=(length(vecPseudotimeUnique)-1)/(k-1)))]) 
    }
    vecWithinSS[k] <- lsKmeansResults$tot.withinss
    lsCentroids[[k]] <- lsKmeansResults$centers
    matAssignments[k,] <- lsKmeansResults$cluster
  }
  # Select K:
  # Chose model with largest K for which all models with larger
  # K have a smaller within cluster sum of squares.
  # This is a stability measure: If models with larger degrees
  # of freedom perform worse, this is likely due to model instability.
  # -1 to get the last TRUE before the the first FALSE
  scaKhat <- min(which( sapply( 1:scaK, function(k){
    all(vecWithinSS[c(1:scaK)>k] < vecWithinSS[k])
  } ) == FALSE )) - 1
  
  # Summarise output:
  lsResultsKhat <- list()
  lsResultsKhat[[1]] <- matAssignments[scaKhat,]
  lsResultsKhat[[2]] <- sort( lsCentroids[[scaKhat]] )
  lsResultsKhat[[3]] <- scaKhat
  names(lsResultsKhat) <- c("Assignments","Centroids","K")
  return(lsResultsKhat) 
}