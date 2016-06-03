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
#' @param matCountDataProc: (matrix genes x samples)
#'    Count data: Reduced version of \code{matCountData}. 
#'    For internal use.
#' @param matProbNB: (probability vector genes x samples) 
#'    Probability of observations to come from negative binomial 
#'    component of mixture model.
#' @param strMode: (str) [Default "batch"] 
#'    {"batch","longitudinal","singlecell"}
#'    Mode of model fitting.
#' 
#' @return matNormConst: (numeric matrix genes x samples) 
#'    Model scaling factors for each observation: Take
#'    sequencing depth and longitudinal time series mean
#'    within a gene into account (size and translation
#'    factors).
#' @export

computeNormConst <- function(matCountDataProc,
  matProbNB=NULL,
  dfAnnotationProc,
  strMode="batch"){
  
  # 1. Compute size factors
  # Size factors account for differential sequencing depth.
  
  # Compute geometric count mean over replicates
  # for each gene: Set zero counts to one
  # In the case of strMode=singlecell, this becomes
  # the weighted geometric count mean, weighted by
  # by the probability of each observation to come
  # from the negative binomial distribution.
  matCountDataProcNoZeros <- matCountDataProc
  matCountDataProcNoZeros[matCountDataProcNoZeros==0] <- 0.1
  boolObserved <- !is.na(matCountDataProc)
  if(strMode=="batch" | strMode=="longitudinal"){
    # Take geometric mean
    vecGeomMean <- sapply(c(1:dim(matCountDataProcNoZeros)[1]), 
      function(gene){
        ( prod(matCountDataProcNoZeros[gene,boolObserved[gene,]]) )^
          ( 1/sum(boolObserved[gene,]) )
      })
  } else if(strMode=="singlecell"){
    # Take weighted geometric mean
    vecGeomMean <- sapply(c(1:dim(matCountDataProcNoZeros)[1]), 
      function(gene){
        ( prod(matCountDataProcNoZeros[gene,boolObserved[gene,]]^matProbNB[gene,boolObserved[gene,]], na.rm=TRUE) )^
          ( 1/sum(matProbNB[gene,boolObserved[gene,]], na.rm=TRUE) )
      })
  }
  matGeomMeans <- matrix(vecGeomMean, 
    nrow=dim(matCountDataProc)[1], 
    ncol=dim(matCountDataProc)[2], 
    byrow=FALSE)
  
  # Compute ratio of each observation to geometric
  # mean.
  matSizeRatios <- matCountDataProc / matGeomMeans
  
  # Chose median of ratios over genes as size factor
  vecSizeFactors <- apply(matSizeRatios, 2,
    function(replicate){
      median(replicate, na.rm=TRUE)
    })
  
  if(any(vecSizeFactors==0)){
    warning("WARNING: Found size factors==0, setting these to 1.")
    vecSizeFactors[vecSizeFactors==0] <- 1
  }
  
  matSizeFactors <- matrix(vecSizeFactors,
    nrow=dim(matCountDataProc)[1],
    ncol=dim(matCountDataProc)[2],
    byrow=TRUE)
  colnames(matSizeFactors) <- colnames(matCountDataProc)
  rownames(matSizeFactors) <- rownames(matCountDataProc)
  
  # 2. Compute translation factors:
  # Translation factors account for different mean expression levels of a
  # gene between longitudinal sample series.
  if(strMode=="longitudinal"){
    vecLongitudinalSeriesAssign <- dfAnnotationProc[match(
      colnames(matCountDataProc),
      dfAnnotationProc$Sample),]$LongitudinalSeries
    names(vecLongitudinalSeriesAssign) <- colnames(matCountDataProc)
    vecLongitudinalSeries <- unique( vecLongitudinalSeriesAssign )
    matCountDataProcNorm <- matCountDataProc/matSizeFactors
    # Fit mean to normalised counts. Count normalisation 
    # corresponds to model normalisation in impulse model 
    # fitting.
    vecMu <- apply(matCountDataProcNorm,1,function(gene){mean(gene, na.rm=TRUE)})
    matMuLongitudinal <- t(apply(matCountDataProcNorm, 1, function(gene){
      sapply(vecLongitudinalSeries, function(longser){
        mean(gene[vecLongitudinalSeriesAssign==longser], na.rm=TRUE)
      })
    }))
    colnames(matMuLongitudinal) <- vecLongitudinalSeries
    
    # Compute translation factors: Normalisation factor
    # to scale impulse model to a longitudinal series.
    # Note: Translation factor is the same for all samples
    # of a gene in a longitudinal series
    
    # Scaled mean ratio per sample
    matMu <- matrix(vecMu,
      nrow=dim(matCountDataProc)[1],
      ncol=dim(matCountDataProc)[2],
      byrow=FALSE)
    matTranslationFactors <- matMuLongitudinal[,as.vector(vecLongitudinalSeriesAssign)]/matMu
    # Note: If all samples are zero, scaMu is zero and
    # vecTranslationFactors are NA. Genes only with zero counts
    # are removed in processData(). The input data to this function
    # ma still consist entirely of zeros if this function is called
    # on the subset of case or control data (if control condition
    # is given). Those subsets may be only zeros. This exception
    # is caught here and vecTranslationFactors set to 1 from NA.
    # This removes the effect of vecTranslationFactors on fitting.
    # The negative binomial density can still be
    # evaluated as it is initialised with values from an impulse model
    # padded with zero counts.
    matTranslationFactors[matTranslationFactors==0] <- 1
    
    matNormConst <- matSizeFactors*matTranslationFactors
  } else {
    matNormConst <- matSizeFactors
  }
  colnames(matNormConst) <- colnames(matCountDataProc)
  rownames(matNormConst) <- rownames(matCountDataProc)
  
  return(list(matNormConst=matNormConst,matSizeFactors=matSizeFactors))
}