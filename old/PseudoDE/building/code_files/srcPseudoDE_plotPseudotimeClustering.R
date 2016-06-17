#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++     Plot pseudotime clustering of cells    +++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Plot pseudotime clustering of cells
#' 
#' Summarises the 1D clustering of cells in pseudotime by two plots:
#' The empirical distribution function of cells in pseudotime and
#' the empricial cumulative distribution function of cells in 
#' pseudotime with cluster borders indicated in each.
#' 
#' @seealso Called by \code{runPseudoDE}.
#' 
#' @param vecPseudotime: (numerical vector length number of cells)
#'    Pseudotime coordinates (1D) of cells: One scalar per cell.
#'    Has to be named: Names of elements are cell names.
#' @param lsResultsClustering (list {"Assignments","Centroids","K"})
#'    \itemize{
#'      \item   Assignments: (integer vector length number of
#'        cells) Index of cluster assigned to each cell.
#'      \item   Centroids: 1D Coordinates of cluster centroids,
#'        one scalar per centroid.
#'      \item   K: (scalar) Number of clusters selected.
#'      }
#' @param strPDFname: (str) Name of .pdf with plots.
#'          
#' @return NULL
#' @export

plotPseudotimeClustering <- function(vecPseudotime, 
  lsResultsClustering,
  strPDFname="PseudoDE_ZINBfitsPseudotimeClustering.pdf"){
  
  dfClusterBorders <- data.frame( borders=sapply( seq(1,length(lsResultsClustering$Centroids)-1), 
    function(centroid){(lsResultsClustering$Centroids[centroid]+lsResultsClustering$Centroids[centroid+1])/2} ))
  dfPseudotime <- data.frame( pseudotime=as.vector(vecPseudotime) )
  
  plotEDF <- ggplot() +
    geom_density(data=dfPseudotime, aes(x=pseudotime), colour="black", bw=1) +
    geom_vline(data=dfClusterBorders, aes(xintercept=borders), colour="green", linetype = "longdash") +
    labs(title="Density estimation of cells in pseudotime") +
    xlab("pseudotime") +
    ylab("empirical probability density")
  
  plotECDF <- ggplot() +
    stat_ecdf(data=dfPseudotime, aes(x=pseudotime), colour="red") +
    geom_vline(data=dfClusterBorders, aes(xintercept=borders), colour="green", linetype = "longdash") +
    labs(title="Empirical cumulative density of cells in pseudotime") +
    xlab("pseudotime") +
    ylab("empirical cumulative probability density")
  
  pdf(strPDFname)
  print(plotEDF)
  print(plotECDF)
  dev.off()
  
  return(NULL)
}