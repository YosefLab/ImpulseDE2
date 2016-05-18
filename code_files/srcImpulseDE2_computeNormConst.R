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
#' @param matProbNB: (probability vector genes x replicates) 
#'    Probability of observations to come from negative binomial 
#'    component of mixture model.
#' @param strMode: (str) [Default "batch"] {"batch","timecourses","singlecell"}
#'    Mode of model fitting.
#' 
#' @return vecSizeFactors: (numeric vector number of replicates) 
#'    Normalisation constants for counts for each replicate. 
#' @export

computeNormConst <- function(arr2DCountData, dfAnnotationFull,
  matProbNB, strMode="batch"){
  
  # Compute geometric count mean over replicates
  # for each gene: Set zero counts to one
  # In the case of strMode=singlecell, this becomes
  # the weighted geometric count mean, weighted by
  # by the probability of each observation to come
  # from the negative binomial distribution.
  arr2DCountDataNoZeros <- arr2DCountData
  arr2DCountDataNoZeros[arr2DCountDataNoZeros==0] <- 0.1
  boolObserved <- !is.na(arr2DCountData)
  if(strMode=="batch" | strMode=="timecourses"){
    # Take geometric mean
    vecGeomMean <- sapply(c(1:dim(arr2DCountDataNoZeros)[1]), 
      function(gene){
        ( prod(arr2DCountDataNoZeros[gene,boolObserved[gene,]]) )^
          ( 1/sum(boolObserved[gene,]) )
      })
  } else if(strMode=="singlecell"){
    # Take weighted geometric mean
    vecGeomMean <- sapply(c(1:dim(arr2DCountDataNoZeros)[1]), 
      function(gene){
        ( prod(arr2DCountDataNoZeros[gene,boolObserved[gene,]]^matProbNB[gene,boolObserved[gene,]], na.rm=TRUE) )^
          ( 1/sum(matProbNB[gene,boolObserved[gene,]], na.rm=TRUE) )
      })
  }
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
  
  if(any(vecSizeFactors==0)){
    warning("WARNING: Found size factors==0, setting these to 1.")
    vecSizeFactors[vecSizeFactors==0] <- 1
  }
  
  return(vecSizeFactors)
}