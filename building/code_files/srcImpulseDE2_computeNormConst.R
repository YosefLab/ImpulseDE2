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

computeTranslationFactors <- function(matCountDataProc,
  matSizeFactors,
  vecDispersions,
  scaSmallRun=NULL,
  dfAnnotationProc,
  strCaseName,
  strControlName=NULL){
  
  # Compute size factor normalised count data
  #matCountDataProcNorm <- matCountDataProc/matSizeFactors
  #colnames(matCountDataProcNorm) <- colnames(matCountDataProc)
  
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
  #vecMuAll <- apply(matCountDataProcNorm,1,function(gene){mean(gene, na.rm=TRUE)})
  vecMuAll <- sapply(seq(1,dim(matCountDataProc)[1]),function(i){
    exp(unlist(optim(
      par=log(mean(matCountDataProc[i,]/matSizeFactors[i,], na.rm=TRUE)+1),
      fn=evalLogLikNBMean_comp,
      vecCounts=matCountDataProc[i,],
      scaDispEst=vecDispersions[i], 
      vecNormConst=matSizeFactors[i,],
      vecboolObserved=!is.na(matCountDataProc[i,]),
      method="BFGS",
      control=list(maxit=1000,fnscale=-1)
    )["par"]))
  })
  matMuAll <- matrix(vecMuAll,
    nrow=dim(matCountDataProc)[1],
    ncol=dim(matCountDataProc)[2],
    byrow=FALSE)
  #matMuLongitudinalAll <- t(apply(matCountDataProcNorm, 1, function(gene){
  #  sapply(vecLongitudinalSeries, function(longser){
  #    mean(gene[vecLongitudinalSeriesAssign==longser], na.rm=TRUE)
  #  })
  #}))
  matMuLongitudinalAll <- array(NA,c(dim(matCountDataProc)[1],length(vecLongitudinalSeries)))
  rownames(matMuLongitudinalAll) <- rownames(matCountDataProc)
  colnames(matMuLongitudinalAll) <- vecLongitudinalSeries
  for(longser in vecLongitudinalSeries){
    vecboolCurrentSer <- vecLongitudinalSeriesAssign==longser
    matMuLongitudinalAll[,longser] <- sapply(seq(1,dim(matCountDataProc)[1]), function(i){
      exp(unlist(optim(
        par=log(mean(matCountDataProc[i,vecboolCurrentSer]/matSizeFactors[i,vecboolCurrentSer], na.rm=TRUE)+1),
        fn=evalLogLikNBMean_comp,
        vecCounts=matCountDataProc[i,vecboolCurrentSer],
        scaDispEst=vecDispersions[i], 
        vecNormConst=matSizeFactors[i,vecboolCurrentSer],
        vecboolObserved=!is.na(matCountDataProc[i,vecboolCurrentSer]),
        method="BFGS",
        control=list(maxit=1000,fnscale=-1)
      )["par"]))
    })
  }
  
  matTranslationFactorsAll <- matMuLongitudinalAll[,as.vector(vecLongitudinalSeriesAssign)]/matMuAll
  matTranslationFactorsAll[matTranslationFactorsAll==0 | matMuAll==0] <- 1
  
  if(!is.null(strControlName)){
    warning("not tested: control translation factors computation")
    # 2. Case
    vecboolindColsCase <- c(colnames(matCountDataProc) %in% dfAnnotationProc[dfAnnotationProc$Condition==strCaseName,]$Sample)
    vecLongitudinalSeriesCase <- unique(dfAnnotationProc[dfAnnotationProc$Condition==strCaseName,]$LongitudinalSeries)
    #vecMuCase <- apply(matCountDataProcNorm[,vecboolindColsCase],1,function(gene){mean(gene, na.rm=TRUE)})
    vecMuCase <- sapply(seq(1,dim(matCountDataProc)[1]),function(i){
      exp(unlist(optim(
        par=log(mean(matCountDataProc[i,vecboolindColsCase]/matSizeFactors[i,vecboolindColsCase], na.rm=TRUE)+1),
        fn=evalLogLikNBMean_comp,
        vecCounts=matCountDataProc[i,vecboolindColsCase],
        scaDispEst=vecDispersions[i], 
        vecNormConst=matSizeFactors[i,vecboolindColsCase],
        vecboolObserved=!is.na(matCountDataProc[i,vecboolindColsCase]),
        method="BFGS",
        control=list(maxit=1000,fnscale=-1)
      )["par"]))
    })
    matMuLongitudinalCase <- matrix(NA,nrow=dim(matCountDataProc)[1],ncol=length(vecLongitudinalSeries))
    rownames(matMuLongitudinalCase) <- rownames(matCountDataProc)
    colnames(matMuLongitudinalCase) <- vecLongitudinalSeries
    #for(longser in vecLongitudinalSeries){
    #  matMuLongitudinalCase[,longser] <- apply(matCountDataProcNorm[,vecboolindColsCase & vecLongitudinalSeriesAssign==longser], 1, 
    #    function(gene){mean(gene, na.rm=TRUE)
    #    })
    #}
    for(longser in vecLongitudinalSeriesCase){
      vecboolCurrentSer <- vecboolindColsCase & vecLongitudinalSeriesAssign==longser
      matMuLongitudinalCase[,longser] <- sapply(seq(1,dim(matCountDataProc)[1]), function(i){
        exp(unlist(optim(
          par=log(mean(matCountDataProc[i,vecboolCurrentSer]/matSizeFactors[i,vecboolCurrentSer], na.rm=TRUE)+1),
          fn=evalLogLikNBMean_comp,
          vecCounts=matCountDataProc[i,vecboolCurrentSer],
          scaDispEst=vecDispersions[i], 
          vecNormConst=matSizeFactors[i,vecboolCurrentSer],
          vecboolObserved=!is.na(matCountDataProc[i,vecboolCurrentSer]),
          method="BFGS",
          control=list(maxit=1000,fnscale=-1)
        )["par"]))
      })
    }
    
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
    vecLongitudinalSeriesCtrl <- unique(dfAnnotationProc[dfAnnotationProc$Condition==strControlName,]$LongitudinalSeries)
    #vecMuCtrl <- apply(matCountDataProcNorm[,vecboolindColsCtrl],1,function(gene){mean(gene, na.rm=TRUE)})
    vecMuCtrl <- sapply(seq(1,dim(matCountDataProc)[1]),function(i){
      exp(unlist(optim(
        par=log(mean(matCountDataProc[i,vecboolindColsCtrl]/matSizeFactors[i,vecboolindColsCtrl], na.rm=TRUE)+1),
        fn=evalLogLikNBMean_comp,
        vecCounts=matCountDataProc[i,vecboolindColsCtrl],
        scaDispEst=vecDispersions[i], 
        vecNormConst=matSizeFactors[i,vecboolindColsCtrl],
        vecboolObserved=!is.na(matCountDataProc[i,vecboolindColsCtrl]),
        method="BFGS",
        control=list(maxit=1000,fnscale=-1)
      )["par"]))
    })
    matMuLongitudinalCtrl <- matrix(NA,nrow=dim(matCountDataProc)[1],ncol=length(vecLongitudinalSeries))
    rownames(matMuLongitudinalCtrl) <- rownames(matCountDataProc)
    colnames(matMuLongitudinalCtrl) <- vecLongitudinalSeries
    #for(longser in vecLongitudinalSeries){
    #  matMuLongitudinalCtrl[,longser] <- apply(matCountDataProcNorm[,vecboolindColsCtrl & vecLongitudinalSeriesAssign==longser], 1, 
    #    function(gene){mean(gene, na.rm=TRUE)
    #    })
    #}
    for(longser in vecLongitudinalSeriesCtrl){
      vecboolCurrentSer <- vecboolindColsCtrl & vecLongitudinalSeriesAssign==longser
      matMuLongitudinalCtrl[,longser] <- sapply(seq(1,dim(matCountDataProc)[1]), function(i){
        exp(unlist(optim(
          par=log(mean(matCountDataProc[i,vecboolCurrentSer]/matSizeFactors[i,vecboolCurrentSer], na.rm=TRUE)+1),
          fn=evalLogLikNBMean_comp,
          vecCounts=matCountDataProc[i,vecboolCurrentSer],
          scaDispEst=vecDispersions[i], 
          vecNormConst=matSizeFactors[i,vecboolCurrentSer],
          vecboolObserved=!is.na(matCountDataProc[i,vecboolCurrentSer]),
          method="BFGS",
          control=list(maxit=1000,fnscale=-1)
        )["par"]))
      })
    }
    
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
    lsMatTranslationFactors <- computeTranslationFactors(matCountDataProc=matCountDataProc,
      matSizeFactors=matSizeFactors,
      vecDispersions=vecDispersions,
      scaSmallRun=scaSmallRun,
      dfAnnotationProc=dfAnnotationProc,
      strCaseName=strCaseName,
      strControlName=strControlName)
  } else {
    lsMatTranslationFactors <- NULL
  }
  
  return(list(lsMatTranslationFactors=lsMatTranslationFactors,matSizeFactors=matSizeFactors))
}