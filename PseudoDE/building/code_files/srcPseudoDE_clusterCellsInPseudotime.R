#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++     Cluster cells in pseudotime    ++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cluster cells in pseudotime
#' 
#' This function clusters cells in the 1D space pseudotime with K-means
#' and selects the number of clusters to represent the sample based
#' on the gap statistic.
#' 
#' @seealso Called by \code{runPseudoDE}.
#' 
#' @param vecPseudotime: (numerical vector length number of cells)
#'    Pseudotime coordinates (1D) of cells: One scalar per cell.
#' 
#' @return (list {"Assignments","Centroids","K"})
#'    \itemize{
#'      \item   Assignments: (integer vector length number of
#'        cells) Index of cluster assigned to each cell.
#'      \item   Centroids: 1D Coordinates of cluster centroids,
#'        one scalar per centroid.
#'      \item   K: (scalar) Number of clusters selected.
#'      }
#' @export

clusterCellsInPseudotime <- function(vecPseudotime){
  
  # Range of K for K-means is number of unique cells
  vecPseudotimeUnique <- unique(vecPseudotime)
  K <- max(round( log(length(vecPseudotimeUnique)/2 )),10)
  vecWithinSS <- array(NA,K)
  lsCentroids <- list()
  matAssignments <- array(NA,c(K,length(vecPseudotime)))
  # Loop over range of K
  for(k in 1:K){
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
  # I.
  # Chose model with largest K for which all models with larger
  # K have a smaller within cluster sum of squares.
  # This is a stability measure: If models with larger degrees
  # of freedom perform worse, this is likely due to model instability.
  # -1 to get the last TRUE before the the first FALSE
  #Khat <- min(which( sapply( 1:K, function(k){
  #  all(vecWithinSS[c(1:K)>k] < vecWithinSS[k])
  #} ) == FALSE )) - 1
  # II.
  # Gap statistic (Tibshirani)
  # Bootstraps from a uniform null distribution of cells in pseudo-time.
  B <- 100
  N <- length(vecPseudotime)
  lsPseudotimeBoot <- lapply(seq(1,B),function(b){runif(N, min=min(vecPseudotime), max=max(vecPseudotime))})
  vecTmax <- unlist(lapply( lsPseudotimeBoot, max ))
  vecTmin  <- unlist(lapply( lsPseudotimeBoot, min ))
  matLogWithinSSBoot <- array(NA,c(B,K))
  # Loop over range of K
  for(b in seq(1,B)){
    lsCenters <- lapply(seq(1,k),function(k){
      if(k==1){ return(1) }
      else{ return(vecTmin[b] + (seq(0,(k-1))+0.5)/k*(vecTmax[b]-vecTmin[b])) }
    })
    vecWithinSSBoot <- sapply(seq(1,K), function(k){
      kmeans( x=lsPseudotimeBoot[[b]], centers=lsCenters[[k]] )$tot.withinss
    })
    matLogWithinSSBoot[b,] <- log(vecWithinSSBoot/seq(1:K))
  }
  vecLogWithinSSBootExpect <- 1/B*apply(matLogWithinSSBoot,2,sum)
  vecGapStat <- vecLogWithinSSBootExpect - log(vecWithinSS/seq(1:K))
  vecSigmaK <- sapply(seq(1,K), function(k){
    sqrt(1/B*sum(sapply( seq(1,B), function(b){(matLogWithinSSBoot[b,k]-vecLogWithinSSBootExpect[k])^2} )))
  })
  Khat <- min(which(sapply( seq(1:(K-1)), function(k){
    vecGapStat[k] >= vecGapStat[k+1] - vecSigmaK[k+1]
  })))
  
  # Summarise output:
  lsResultsKhat <- list()
  lsResultsKhat[[1]] <- matAssignments[Khat,]
  lsResultsKhat[[2]] <- sort( lsCentroids[[Khat]] )
  lsResultsKhat[[3]] <- Khat
  names(lsResultsKhat) <- c("Assignments","Centroids","K")
  return(lsResultsKhat) 
}