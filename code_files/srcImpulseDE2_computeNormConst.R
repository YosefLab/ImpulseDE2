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
#' 
#' @seealso Called by \code{computeNormConst}.
#' 
#' @param matCountDataProc: (matrix genes x samples)
#'    Count data: Reduced version of \code{matCountData}. 
#'    For internal use.
#' @param strMode: (str) [Default "batch"] 
#'    {"batch","longitudinal","singlecell"}
#'    Mode of model fitting.
#' 
#' @return matSizeFactors: (numeric matrix genes x samples) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors).
#' @export

computeSizeFactors <- function(matCountDataProc,
  strMode){
  
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
#' @param dfAnnotationProc: (Table) Processed annotation table. 
#'    Lists co-variables of samples: 
#'    Sample, Condition, Time (numeric), TimeCateg (str)
#'    (and LongitudinalSeries). For internal use.
#' @param strCaseName (str) Name of the case condition in \code{dfAnnotationRedFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#'    \code{dfAnnotationRedFull}.
#' 
#' @return (list {case} or {case, ctrl, combined}) List of
#'    translation factor matrices.
#'    \itemize{
#'      \item matTranslationFactorsAll: (numeric matrix genes x samples) 
#'    Model scaling factors for each observation which take
#'    longitudinal time series mean within a gene into account 
#'    (translation factors). Computed based based on all samples.
#'      \item matTranslationFactorsCase: As matTranslationFactorsAll
#'      for case samples only, non-case samples set NA.
#'      \item matTranslationFactorsCtrl: As matTranslationFactorsAll
#'      for control samples only, non-control samples set NA.
#'      }
#' @export

computeTranslationFactors <- function(matCountDataProc,matSizeFactors,
  dfAnnotationProc,strCaseName,strControlName=NULL){
  
  # Compute size factor normalised count data
  matCountDataProcNorm <- matCountDataProc/matSizeFactors
  colnames(matCountDataProcNorm) <- colnames(matCountDataProc)
  
  vecLongitudinalSeriesAssign <- dfAnnotationProc[match(
    colnames(matCountDataProc),
    dfAnnotationProc$Sample),]$LongitudinalSeries
  names(vecLongitudinalSeriesAssign) <- colnames(matCountDataProc)
  vecLongitudinalSeries <- unique( vecLongitudinalSeriesAssign )
  
  # 1. All data: This is case with a single condition and 
  # combined with control condition.
  
  # Fit mean to normalised counts. Count normalisation 
  # corresponds to model normalisation in impulse model 
  # fitting.
  vecMuAll <- apply(matCountDataProcNorm,1,function(gene){mean(gene, na.rm=TRUE)})
  matMuAll <- matrix(vecMuAll,
    nrow=dim(matCountDataProc)[1],
    ncol=dim(matCountDataProc)[2],
    byrow=FALSE)
  matMuLongitudinalAll <- t(apply(matCountDataProcNorm, 1, function(gene){
    sapply(vecLongitudinalSeries, function(longser){
      mean(gene[vecLongitudinalSeriesAssign==longser], na.rm=TRUE)
    })
  }))
  colnames(matMuLongitudinalAll) <- vecLongitudinalSeries
  
  matTranslationFactorsAll <- matMuLongitudinalAll[,as.vector(vecLongitudinalSeriesAssign)]/matMuAll
  matTranslationFactorsAll[matTranslationFactorsAll==0 | matMuAll==0] <- 1
  
  if(!is.null(strControlName)){
    # 2. Case
    vecboolindColsCase <- c(colnames(matCountDataProc) %in% dfAnnotationProc[dfAnnotationProc$Condition==strCaseName,]$Sample)
    vecMuCase <- apply(matCountDataProcNorm[,vecboolindColsCase],1,function(gene){mean(gene, na.rm=TRUE)})
    matMuLongitudinalCase <- matrix(NA,nrow=dim(matCountDataProc)[1],ncol=length(vecLongitudinalSeries))
    colnames(matMuLongitudinalCase) <- vecLongitudinalSeries
    for(longser in vecLongitudinalSeries){
      matMuLongitudinalCase[,longser] <- apply(matCountDataProcNorm[,vecboolindColsCase & vecLongitudinalSeriesAssign==longser], 1, 
        function(gene){mean(gene, na.rm=TRUE)
      })
    }
    colnames(matMuLongitudinalCase) <- vecLongitudinalSeries
    
    # Scaled mean ratio per sample
    matMuCase <- matrix(vecMuCase,
      nrow=dim(matCountDataProc)[1],
      ncol=dim(matCountDataProc)[2],
      byrow=FALSE)
    matTranslationFactorsCase <- matMuLongitudinalCase[,as.vector(vecLongitudinalSeriesAssign)]/matMuCase
    matTranslationFactorsCase[matTranslationFactorsCase==0 | matMuCase==0] <- 1
    # Take out samples which do not belong to case
    matTranslationFactorsCase[,!vecboolindColsCase] <- NA
    
    # 3. Control
    vecboolindColsCtrl <- c(colnames(matCountDataProc) %in% dfAnnotationProc[dfAnnotationProc$Condition==strControlName,]$Sample)
    vecMuCtrl <- apply(matCountDataProcNorm[,vecboolindColsCtrl],1,function(gene){mean(gene, na.rm=TRUE)})
    matMuLongitudinalCtrl <- matrix(NA,nrow=dim(matCountDataProc)[1],ncol=length(vecLongitudinalSeries))
    colnames(matMuLongitudinalCtrl) <- vecLongitudinalSeries
    for(longser in vecLongitudinalSeries){
      matMuLongitudinalCtrl[,longser] <- apply(matCountDataProcNorm[,vecboolindColsCtrl & vecLongitudinalSeriesAssign==longser], 1, 
        function(gene){mean(gene, na.rm=TRUE)
      })
    }
    colnames(matMuLongitudinalCtrl) <- vecLongitudinalSeries
    
    # Scaled mean ratio per sample
    matMuCtrl <- matrix(vecMuCtrl,
      nrow=dim(matCountDataProc)[1],
      ncol=dim(matCountDataProc)[2],
      byrow=FALSE)
    matTranslationFactorsCtrl <- matMuLongitudinalCtrl[,as.vector(vecLongitudinalSeriesAssign)]/matMuCtrl
    matTranslationFactorsCtrl[matTranslationFactorsCtrl==0 | matMuCtrl==0] <- 1
    # Take out samples which do not belong to control
    matTranslationFactorsCtrl[,!vecboolindColsCtrl] <- NA
    
    lsMatTranslationFactors <- list(matTranslationFactorsAll,
      matTranslationFactorsCase, matTranslationFactorsCtrl)
    names(lsMatTranslationFactors) <- c("combined","case","control")
    return(lsMatTranslationFactors)
  } else {
    lsMatTranslationFactors <- list(matTranslationFactorsAll)
    names(lsMatTranslationFactors) <- c("case")
    return(lsMatTranslationFactors)
  }
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
#' 
#' @seealso Called by \code{runImpulseDE2}. 
#' Calls \code{computeTranslationFactors}.
#' 
#' @param matCountDataProc: (matrix genes x samples)
#'    Count data: Reduced version of \code{matCountData}. 
#'    For internal use.
#' @param matProbNB: (probability vector genes x samples) 
#'    Probability of observations to come from negative binomial 
#'    component of mixture model.
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
#' 
#' @return (list {lsMatTranslationFactors,matSizeFactors})
#'    \itemize{
#'      \item lsMatTranslationFactors 
#'      (list {case} or {case, ctrl, combined}) List of
#'      translation factor matrices. NULL if strMode not
#'      "longitudinal".
#'      \itemize{
#'        \item matTranslationFactorsAll: (numeric matrix genes x samples) 
#'      Model scaling factors for each observation which take
#'      longitudinal time series mean within a gene into account 
#'      (translation factors). Computed based based on all samples.
#'        \item matTranslationFactorsCase: As matTranslationFactorsAll
#'        for case samples only, non-case samples set NA.
#'        \item matTranslationFactorsCtrl: As matTranslationFactorsAll
#'        for control samples only, non-control samples set NA.
#'      }
#'      \item matSizeFactors: (numeric matrix genes x samples) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors).
#'      } 
#' @export

computeNormConst <- function(matCountDataProc,
  matProbNB=NULL,
  dfAnnotationProc,
  strCaseName,
  strControlName=NULL,
  strMode="batch"){
  
  # 1. Compute size factors
  # Size factors account for differential sequencing depth.
  matSizeFactors <- computeSizeFactors(matCountDataProc,
    strMode)
  
  # 2. Compute translation factors:
  # Translation factors account for different mean expression levels of a
  # gene between longitudinal sample series.
  if(strMode=="longitudinal"){
    lsMatTranslationFactors <- computeTranslationFactors(matCountDataProc=matCountDataProc,
      matSizeFactors=matSizeFactors,
      dfAnnotationProc=dfAnnotationProc,
      strCaseName=strCaseName,
      strControlName=strControlName)
  } else {
    lsMatTranslationFactors <- NULL
  }
  
  return(list(lsMatTranslationFactors=lsMatTranslationFactors,matSizeFactors=matSizeFactors))
}