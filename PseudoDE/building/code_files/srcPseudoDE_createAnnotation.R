#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++     Create annotation tables    +++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Create annotation table for clustered single cells for 
#' differential expression analysis 
#' 
#' Create annotation table for clustered single cells for differential expression 
#' analysis based on specifications required by ImpulseDE2. Cells within 
#' a cluster are considered to be biological replicates and the cluster
#' centroids (which have pseudotime coordinates) are the time points assigned
#' to each set of replicates. All cells are labelled to originate from
#' condition "case".
#' 
#' @seealso Called by \code{runPseudoDE}. \code{createAnnotationByCell} creates
#' annotation table for each cell as a separate sample, not summarising cells
#' into clusters.
#' 
#' @param matCounts: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param lsResultsClustering (list {"Assignments","Centroids","K"})
#'    \itemize{
#'      \item   Assignments: (integer vector length number of
#'        cells) Index of cluster assigned to each cell.
#'      \item   Centroids: 1D Coordinates of cluster centroids,
#'        one scalar per centroid.
#'      \item   K: (scalar) Number of clusters selected.
#'      }
#' @param vecPseudotime: (numerical vector length number of cells)
#'    Pseudotime coordinates (1D) of cells: One scalar per cell.
#'    Has to be named: Names of elements are cell names.
#' 
#' @return dfAnnotationImpulseDE2: (Table)
#'    Annotation table. Lists covariates of samples: 
#'    Sample, Condition, Time. Time must be numeric. No longitudinal
#'    series given as scRNA-seq data are not longitudinal.
#' @export
 
createAnnotationByCluster <- function(matCounts,
  vecPseudotime,
  lsResultsClustering){
  
  vecTimepoints <- lsResultsClustering$Centroids[lsResultsClustering$Assignments]
  dfAnnotationImpulseDE2 <- as.data.frame(cbind(
    Sample=names(vecPseudotime),
    Condition="case",
    Time=vecTimepoints ))
  dfAnnotationImpulseDE2$Time <- as.numeric(as.vector( dfAnnotationImpulseDE2$Time ))
  
  return(dfAnnotationImpulseDE2)
}

#' Create annotation table for clustered single cells for 
#' differential expression analysis 
#' 
#' Create annotation table for individual single cells for differential expression 
#' analysis based on specifications required by ImpulseDE2. Cells within a cluster 
#' are not considered to be biological replicates. All cells are labelled to 
#' originate from condition "case".
#' 
#' @seealso Called by \code{runPseudoDE}. \code{createAnnotationByCluster} creates
#' annotation table for clusters of cells.
#' 
#' @param matCounts: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param vecPseudotime: (numerical vector length number of cells)
#'    Pseudotime coordinates (1D) of cells: One scalar per cell.
#'    Has to be named: Names of elements are cell names.
#' 
#' @return dfAnnotationImpulseDE2: (Table)
#'    Annotation table. Lists covariates of samples: 
#'    Sample, Condition, Time. Time must be numeric. No longitudinal
#'    series given as scRNA-seq data are not longitudinal.
#' @export

createAnnotationByCell <- function(matCounts, 
  vecPseudotime){
  
  dfAnnotationImpulseDE2 <- as.data.frame(cbind(
    Sample=names(vecPseudotime),
    Condition="case",
    Time=as.vector(vecPseudotime) ))
  dfAnnotationImpulseDE2$Time <- as.numeric(as.vector( dfAnnotationImpulseDE2$Time ))
  
  return(dfAnnotationImpulseDE2)
}