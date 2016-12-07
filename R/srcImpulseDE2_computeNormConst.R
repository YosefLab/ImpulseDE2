#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++     Compute normalisation constants    ++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute size factors for a dataset
#' 
#' This function computes size factors for each sample
#' in the dataset and expands them to a matrix of the size
#' of the dataset.
#' Size factors scale the negative binomial likelihood
#' model of a gene to the sequencing depth of each sample.
#' Note that size factors on bulk and single-cell data are 
#' computed differently: Median ratio of data to geometric mean
#' for bul data and normalised relative sequencing depth for
#' single-cell data.
#' 
#' @seealso Called by \code{computeNormConst}.
#' 
#' @param matCountDataProc: (matrix genes x samples)
#'    Count data: Reduced version of \code{matCountData}. 
#'    For internal use.
#' @param scaSmallRun: (integer) [Default NULL] Number of rows
#'    on which ImpulseDE2 is supposed to be run, the full
#'    data set is only used for size factor estimation.
#' @param strMode: (str) [Default "singelbatch"] 
#'    {"singelbatch","batcheffects"}
#'    Batch model.
#' 
#' @return vecSizeFactors: (numeric vectors number of samples) 
#'    Model scaling factors for each sample which take
#'    sequencing depth into account (size factors).
#' @export

computeSizeFactors <- function(matCountDataProc,
                               scaSmallRun=NULL){
  
  # Compute geometric count mean over replicates
  # for genes without zero observations: Samples
  # with more than half zero observations receive 
  # size factor =1 otherwise.
  vecboolZeroObs <- apply(matCountDataProc,1,function(gene){!any(gene==0)})
  # Take geometric mean
  vecGeomMean <- apply(matCountDataProc[vecboolZeroObs,], 1, function(gene){
    ( prod(gene[!is.na(gene)]) )^( 1/sum(!is.na(gene)) )
  })
  
  # Chose median of ratios over genes as size factor
  vecSizeFactors <- apply(matCountDataProc[vecboolZeroObs,], 2, function(sample){
    median(sample/vecGeomMean, na.rm=TRUE)
  })
  
  if(any(vecSizeFactors==0)){
    print("WARNING: Found size factors==0, setting these to 1.")
    vecSizeFactors[vecSizeFactors==0] <- 1
  }
  names(vecSizeFactors) <- colnames(matCountDataProc)
  
  return(vecSizeFactors)
}

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
#' There is the option to supply size factors to this function to override
#' its size factor choice.
#' 
#' @seealso Called by \code{runImpulseDE2}. 
#' Calls \code{computeSizeFactors}.
#' 
#' @param matCountDataProcFull: (matrix genes x samples)
#'    Count data: Reduced version of \code{matCountData}. 
#'    For internal use. This is the entire data set.
#' @param vecSizeFactorsExternal: (numeric vector number of samples) 
#'    Model scaling factors for each sample which take
#'    sequencing depth into account (size factors). One size
#'    factor per sample. These are supplied to ImpulseDE, if this variable
#'    is not set, size factors are computed in this function.
#'     
#' @return (vecSizeFactors: (numeric vector number of samples) 
#'        Model scaling factors for each sample which take
#'        sequencing depth into account (size factors).
#'      } 
#' @export

computeNormConst <- function(matCountDataProc,
                             vecSizeFactorsExternal=NULL){
  
  # Compute size factors
  # Size factors account for differential sequencing depth.
  if(is.null(vecSizeFactorsExternal)){
    # Compute size factors if not supplied to ImpulseDE2.
    vecSizeFactors <- computeSizeFactors(matCountDataProc=matCountDataProc)
  } else {
    # Chose externally supplied size factors if supplied.
    print("Using externally supplied size factors.")
    vecSizeFactors <- vecSizeFactorsExternal[colnames(matCountDataProc)]
  }
  
  return(vecSizeFactors)
}