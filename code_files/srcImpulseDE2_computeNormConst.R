#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++     Compute normalisation constants    ++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute normalisation constant for each replicate
#' 
#' The normalisation constant is the median of the ratio of gene counts versus
#' the geomtric gene count mean. There is one normalisation constant per 
#' replicate. An intuitive alternative would be the sequencing depth, the median
#' ratio is however less sensitive to highly differentially expressed genes
#' with high counts (ref. DESeq).
#' The normalisation constants are used to scale the mean of the
#' negative binomial model inferred during fitting to the sequencing depth 
#' of the given sample. The normalisation constants therefore replace
#' normalisation at the count data level, which is not supposed to be done 
#' in the framework of ImpulseDE2.
#' 
#' @seealso Called by \code{runImpulseDE2}.
#' 
#' @param arr2DCountData (2D array genes x replicates) Count data: Reduced 
#'    version of \code{matCountData}. For internal use.
#' @param dfAnnotationFull (Table) Lists co-variables of individual replicates:
#'    Replicate, Sample, Condition, Time. Time must be numeric.
#' 
#' @return vecSizeFactors: (numeric vector number of replicates) 
#'    Normalisation constants for counts for each replicate. 
#' @export

computeNormConst <- function(arr2DCountData, dfAnnotationFull){
  
  # Compute geometric count mean over replicates
  # for each gene: Set zero counts to one
  arr2DCountDataNoZeros <- arr2DCountData
  arr2DCountDataNoZeros[arr2DCountDataNoZeros==0] <- 1
  vecGeomMean <- apply(arr2DCountDataNoZeros, 1, 
    function(gene){
      ( prod(gene, na.rm=TRUE) )^( 1/(sum(!is.na(gene))) )
    })
  matGeomMeans <- matrix(vecGeomMean, 
    nrow=dim(arr2DCountData)[1], 
    ncol=dim(arr2DCountData)[2], 
    byrow=FALSE)
  
  # Compute ratio of each observation to geometric
  # mean.
  matSizeRatios <- arr2DCountData / matGeomMeans
  
  # Chose median of ratios over genes as size factor
  vecSizeFactors <- apply(matSizeRatios, 2,
    function(replicate){
      median(replicate, na.rm=TRUE)
    })
  names(vecSizeFactors) <- colnames(arr2DCountData)
  
  return(vecSizeFactors)
}