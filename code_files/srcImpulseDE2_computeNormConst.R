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
#' @param arr3DCountData (3D array genes x samples x replicates)
#'    Count data: \code{arr2DCountData} reshaped into a 3D array. 
#'    For internal use.
#' 
#' @return matNormConst: (matrix samples x replicates) Normalisation
#'    constants for each replicate. Missing samples are set NA.
#' @export

computeNormConst <- function(arr3DCountData){
  
  # Compute geometric count mean over replicates
  # for each gene:
  vecGeomMean <- apply(arr3DCountData, 1, 
    function(gene){
      ( prod(gene[gene!=0], na.rm=TRUE) )^( 1/(sum(!is.na(gene[gene!=0]))) )
    })
  arr3DGeomMeans <- array(NA, dim(arr3DCountData))
  for(gene in 1:dim(arr3DGeomMeans)[1]){ arr3DGeomMeans[gene,,] <- vecGeomMean[gene] }
  
  # Compute ratio of each observation to geometric
  # mean.
  matSizeRatios <- arr3DCountData / arr3DGeomMeans
  
  # Chose median of ratios over genes as size factor
  matSizeFactors <- apply(matSizeRatios, c(2,3),
        function(rep){
          median(rep, na.rm=TRUE)
        })
  rownames(matSizeFactors) <- colnames(arr3DCountData)
  colnames(matSizeFactors) <- dimnames(arr3DCountData)[[3]]
  
  return(matSizeFactors)
}