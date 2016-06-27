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
#' @param strMode: (str) [Default "batch"] 
#'    {"batch","longitudinal","singlecell"}
#'    Mode of model fitting.
#' 
#' @return matSizeFactors: (numeric matrix genes x samples) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per sample - rows of this matrix are equal.
#' @export

computeSizeFactors <- function(matCountDataProc,
  scaSmallRun=NULL,
  strMode){
  
  if(strMode=="batch" | strMode=="longitudinal"){
    # Compute geometric count mean over replicates
    # for genes without zero observations: Samples
    # with more than half zero observations receive 
    # size factor =1 otherwise.
    # In the case of strMode=singlecell, this becomes
    # the weighted geometric count mean, weighted by
    # by the probability of each observation to come
    # from the negative binomial distribution.
    matCountDataProcNoZeros <- matCountDataProc
    vecboolZeroObs <- apply(matCountDataProcNoZeros,1,function(gene){!any(gene==0)})
    matCountDataProcNoZeros <- matCountDataProcNoZeros[vecboolZeroObs,]
    matboolObserved <- !is.na(matCountDataProcNoZeros) 
    # Take geometric mean
    vecGeomMean <- sapply(seq(1,dim(matCountDataProcNoZeros)[1]), 
      function(gene){
        ( prod(matCountDataProcNoZeros[gene,matboolObserved[gene,]]) )^
          ( 1/sum(matboolObserved[gene,]) )
      })
    
    matGeomMeans <- matrix(vecGeomMean, 
      nrow=dim(matCountDataProcNoZeros)[1], 
      ncol=dim(matCountDataProcNoZeros)[2], 
      byrow=FALSE)
    
    # Compute ratio of each observation to geometric
    # mean.
    matSizeRatios <- matCountDataProcNoZeros / matGeomMeans
    
    # Chose median of ratios over genes as size factor
    vecSizeFactors <- apply(matSizeRatios, 2,
      function(sample){
        median(sample, na.rm=TRUE)
      })
  } else if(strMode=="singlecell"){
    # Size factors directly represent sequencing depth:
    # Normalised relative sequencing depth.
    vecSeqDepth <- apply(matCountDataProc, 2,
      function(cell){ sum(cell, na.rm=TRUE) })
    vecSizeFactors <- vecSeqDepth/sum(vecSeqDepth)*length(vecSeqDepth)
  }
  
  if(any(vecSizeFactors==0)){
    print("WARNING: Found size factors==0, setting these to 1.")
    vecSizeFactors[vecSizeFactors==0] <- 1
  }
  
  # Replicate vector to matrix
  # Matrix size reduced according to scaSmallRun
  # if ImpulseDE2 is not run on entire data set.
  if(!is.null(scaSmallRun)){ scaNRows <- min(scaSmallRun, dim(matCountDataProc)[1])
  } else {scaNRows <- dim(matCountDataProc)[1] }
  matSizeFactors <- matrix(vecSizeFactors,
    nrow=scaNRows,
    ncol=dim(matCountDataProc)[2],
    byrow=TRUE)
  colnames(matSizeFactors) <- colnames(matCountDataProc)
  rownames(matSizeFactors) <- (rownames(matCountDataProc))[1:scaNRows]
  
  return(matSizeFactors)
}

#' Compute translation factors for a dataset
#' 
#' This function computes translation factors for each sample
#' per gene in the dataset. Note that the translation factors
#' are computed with respect to the overall mean, and with 
#' respect to the condition-wise means if control samples
#' are given. The rescaling of translation factors for case and
#' control conditions is not necessary for differential 
#' expression analysis but gives impulse models which directly
#' represent the data of their conditions, without scaling.
#' Translation factors scale the negative binomial likelihood
#' model of a gene to the mean of a longitudinal series of 
#' samples of this gene. This function is only callde if 
#' strMode="longitudinal".
#' 
#' @seealso Called by \code{computeNormConst}.
#' 
#' @param matCountDataProc: (matrix genes x samples)
#'    Count data: Reduced version of \code{matCountData}. 
#'    For internal use.
#' @param matSizeFactors: (numeric matrix genes x samples) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors).
#' @param scaSmallRun: (integer) [Default NULL] Number of rows
#'    on which ImpulseDE2 is supposed to be run, the full
#'    data set is only used for size factor estimation.
#' @param dfAnnotationProc: (Table) Processed annotation table. 
#'    Lists co-variables of samples: 
#'    Sample, Condition, Time (numeric), TimeCateg (str)
#'    (and LongitudinalSeries). For internal use.
#' @param strCaseName (str) Name of the case condition in \code{dfAnnotationRedFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#'    \code{dfAnnotationRedFull}.
#' 
#' @return matTranslationFactors: (numeric matrix genes x samples) 
#'    Model scaling factors for each observation which take
#'    longitudinal time series mean within a gene into account 
#'    (translation factors). Computed based based on all samples.
#' @export

computeTranslationFactors <- function(matCountDataProc,
  matSizeFactors,
  vecDispersions,
  scaSmallRun=NULL,
  dfAnnotationProc,
  strCaseName,
  strControlName=NULL){
  
  vecLongitudinalSeriesAssign <- dfAnnotationProc[match(
    colnames(matCountDataProc),
    dfAnnotationProc$Sample),]$LongitudinalSeries
  names(vecLongitudinalSeriesAssign) <- colnames(matCountDataProc)
  vecLongitudinalSeries <- unique( vecLongitudinalSeriesAssign )
  
  # Fit mean to normalised counts. Count normalisation 
  # corresponds to model normalisation in impulse model 
  # fitting.
  vecMu <- sapply(seq(1,dim(matCountDataProc)[1]),function(i){
    mean(matCountDataProc[i,]/matSizeFactors[i,], na.rm=TRUE)
    #fitNBMean(vecCounts=matCountDataProc[i,],
    #  scaDispEst=vecDispersions[i],
    #  vecNormConst=matSizeFactors[i,])
  })
  matMu <- matrix(vecMu,
    nrow=dim(matCountDataProc)[1],
    ncol=dim(matCountDataProc)[2],
    byrow=FALSE)
  matMuLongitudinal <- array(NA,c(dim(matCountDataProc)[1],length(vecLongitudinalSeries)))
  rownames(matMuLongitudinal) <- rownames(matCountDataProc)
  colnames(matMuLongitudinal) <- vecLongitudinalSeries
  for(longser in vecLongitudinalSeries){
    vecboolCurrentSer <- vecLongitudinalSeriesAssign==longser
    matMuLongitudinal[,longser] <- sapply(seq(1,dim(matCountDataProc)[1]), function(i){
      mean(matCountDataProc[i,vecLongitudinalSeriesAssign==longser]/
          matSizeFactors[i,vecLongitudinalSeriesAssign==longser], na.rm=TRUE)
      #fitNBMean(vecCounts=matCountDataProc[i,vecLongitudinalSeriesAssign==longser],
      #  scaDispEst=vecDispersions[i],
      #  vecNormConst=matSizeFactors[i,vecLongitudinalSeriesAssign==longser])
    })
  }
  
  matTranslationFactors <- matMuLongitudinal[,as.vector(vecLongitudinalSeriesAssign)]/matMu
  matTranslationFactors[matTranslationFactors==0 | matMu==0] <- 1
  
  return(matTranslationFactors)
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
#' Calls \code{computeTranslationFactors} and
#' \code{computeSizeFactors}.
#' 
#' @param matCountDataProcFull: (matrix genes x samples)
#'    Count data: Reduced version of \code{matCountData}. 
#'    For internal use. This is the entire data set.
#' @param matCountDataProc: (matrix genes x samples)
#'    Count data: Reduced version of \code{matCountData}. 
#'    For internal use. This is the data set reduced to the
#'    genes which are supposed to be analysed.
#' @param scaSmallRun: (integer) [Default NULL] Number of rows
#'    on which ImpulseDE2 is supposed to be run, the full
#'    data set is only used for size factor estimation.
#' @param dfAnnotationProc: (Table) Processed annotation table. 
#'    Lists co-variables of samples: 
#'    Sample, Condition, Time (numeric), TimeCateg (str)
#'    (and LongitudinalSeries). For internal use.
#' @param strCaseName (str) Name of the case condition in \code{dfAnnotationRedFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#'    \code{dfAnnotationRedFull}.
#' @param strMode: (str) [Default "batch"] 
#'    {"batch","longitudinal","singlecell"}
#'    Mode of model fitting.
#' @param vecSizeFactorsExternal: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell. These are supplied to ImpulseDE, if this variable
#'    is not set, size factors are computed in this function.
#' 
#' @return (list {lsMatTranslationFactors,matSizeFactors})
#'    \itemize{
#'      \item matTranslationFactors: (numeric matrix genes x samples) 
#'      Model scaling factors for each observation which take
#'      longitudinal time series mean within a gene into account 
#'      (translation factors). Computed based based on all samples.
#'      \item matSizeFactors: (numeric matrix genes x samples) 
#'        Model scaling factors for each observation which take
#'        sequencing depth into account (size factors). One size
#'        factor per sample - rows of this matrix are equal.
#'      } 
#' @export

computeNormConst <- function(matCountDataProcFull,
  matCountDataProc,
  vecDispersions,
  scaSmallRun=NULL,
  dfAnnotationProc,
  strCaseName,
  strControlName=NULL,
  strMode="batch",
  vecSizeFactorsExternal=NULL){
  
  # 1. Compute size factors
  # Size factors account for differential sequencing depth.
  if(is.null(vecSizeFactorsExternal)){
    # Compute size factors if not supplied to ImpulseDE2.
    matSizeFactors <- computeSizeFactors(matCountDataProc=matCountDataProcFull,
      scaSmallRun=scaSmallRun,
      strMode=strMode)
  } else {
    # Chose externally supplied size factors if supplied.
    # Matrix size reduced according to scaSmallRun
    # if ImpulseDE2 is not run on entire data set.
    print("Using externally supplied size factors.")
    if(!is.null(scaSmallRun)){ scaRows <- min(scaSmallRun,dim(matCountDataProc)[1])
    } else {scaRows <- dim(matCountDataProc)[1] }
    matSizeFactors <- matrix(vecSizeFactorsExternal,
      nrow=scaRows,
      ncol=dim(matCountDataProc)[2],
      byrow=TRUE)
    colnames(matSizeFactors) <- colnames(matCountDataProc)
    rownames(matSizeFactors) <- rownames(matCountDataProc)
  }
  
  # 2. Compute translation factors:
  # Translation factors account for different mean expression levels of a
  # gene between longitudinal sample series.
  if(strMode=="longitudinal"){
    matTranslationFactors <- computeTranslationFactors(matCountDataProc=matCountDataProc,
      matSizeFactors=matSizeFactors,
      vecDispersions=vecDispersions,
      scaSmallRun=scaSmallRun,
      dfAnnotationProc=dfAnnotationProc,
      strCaseName=strCaseName,
      strControlName=strControlName)
  } else {
    matTranslationFactors <- NULL
  }
  
  return(list(matTranslationFactors=matTranslationFactors,
    matSizeFactors=matSizeFactors))
}