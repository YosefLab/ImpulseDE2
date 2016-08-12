#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Impulse model fit    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Fit and compute likelihood of null model for a single gene
#' 
#' Fits an impulse model and a mean model to a single gene.
#' The optimisation method is set within this function (\code{optim}:
#' BFGS). This method is divided into four parts: (I) Prepare data,
#' (II) Fit mean model, (III) Fit impulse model, (IV) Process Fits.
#' Internal function of \code{fitImpulse}.
#' 
#' @seealso Called by \code{fitImpulse_gene}.
#' 
#' @param vecCounts: (count vector number of samples) Count data.
#' @param scaDispersionEstimate: (scalar) Negative binomial 
#'    dispersion coefficient for given gene.
#' @param vecDropoutRate: (probability vector number of samples) 
#'    [Default NULL] Dropout rate/mixing probability of zero inflated 
#'    negative binomial mixturemodel for each gene and cell.
#' @param vecProbNB: (probability vector number of samples) [Default NULL] 
#'    Probability of observations to come from negative binomial 
#'    component of mixture model.
#' @param vecNormConst: (numeric vector number of samples) 
#'    Normalisation constants for each sample.
#' @param vecLongitudinalSeries: (string vector number of longitudinal
#'    sample series) Longitudinal sample series.
#' @param vecLongitudinalSeriesAssign (numeric vector number samples) 
#'    Longitudinal sample series assigned to samples. 
#'    NULL if not operating in strMode="longitudinal".
#' @param vecboolObserved: (bool vector number of samples)
#'    Whether sample is observed (not NA).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' @param vecboolNotZeroObserved: (bool vector number of samples)
#'    Whether sample is not zero and observed (not NA).
#' @param strMode: (str) [Default "batch"] 
#'    {"batch","longitudinal","singlecell"}
#'    Mode of model fitting.
#' @param scaWindowRadius: (integer) [Default NULL]
#'    Smoothing interval radius.
#'    
#' @return (list length 3)
#'    \itemize{
#'      \item scaMu: (scalar) Null model in batch or
#'        singlecell mode: One overall mean.
#'      \item vecMuLongitudinalSeries: (numerical vector
#'        length number of longitudinal series) Null model
#'        in longitudinal mode: One mean per longitudinal
#'        series. Null if mode is not longitudinal.
#'      \item scaLogLikNull: (scalar) Log likelihood of null
#'        model.
#'      }
#' @export

fitNullModel <- function(vecCounts, 
  scaDispersionEstimate,
  vecDropoutRate=NULL, 
  vecProbNB=NULL,
  vecNormConst,
  vecLongitudinalSeries=NULL, 
  vecLongitudinalSeriesAssign=NULL,
  vecboolObserved, 
  vecboolZero=NULL,
  vecboolNotZeroObserved=NULL,
  strMode,
  scaWindowRadius=NULL){
  
  vecMuLongitudinalSeries <- NULL
  
  if(strMode=="batch"){
    # Fit null model:
    scaMu <- fitMuNB(vecCounts=vecCounts,
      scaDispEst=scaDispersionEstimate,
      vecNormConst=vecNormConst)
    
    # Evaluate likelihood of null model:
    scaLogLikNull <- sum(dnbinom(
      vecCounts[vecboolObserved], 
      mu=scaMu*vecNormConst[vecboolObserved], 
      size=scaDispersionEstimate, 
      log=TRUE))
  } else if(strMode=="longitudinal"){
    # Fit null model: 
    scaMu <- fitMuNB(vecCounts=vecCounts,
      scaDispEst=scaDispersionEstimate,
      vecNormConst=vecNormConst)
    vecMuLongitudinalSeries <- sapply(vecLongitudinalSeries,function(longser){ 
      fitMuNB(vecCounts=vecCounts[vecLongitudinalSeriesAssign==longser],
        scaDispEst=scaDispersionEstimate,
        vecNormConst=vecNormConst[vecLongitudinalSeriesAssign==longser])
    })
    names(vecMuLongitudinalSeries) <- vecLongitudinalSeries
    
    # Evaluate likelihood of null model
    scaLogLikNull <- sum(dnbinom(
      vecCounts[vecboolObserved], 
      mu=(vecMuLongitudinalSeries[as.vector(vecLongitudinalSeriesAssign)])[vecboolObserved]*vecNormConst[vecboolObserved], 
      size=scaDispersionEstimate, 
      log=TRUE))
  } else if(strMode=="singlecell"){
    # Fit null model:
    scaMu <- fitMuZINB( vecCounts=vecCounts,
      scaDispEst=scaDispersionEstimate,
      vecNormConst=vecNormConst,
      vecDropoutRateEst=vecDropoutRate,
      vecProbNB=vecProbNB,
      scaWindowRadius=scaWindowRadius )
    if(is.null(scaWindowRadius)){      
      scaLogLikNull <- evalLogLikZINB_comp(vecY=vecCounts,
        vecMu=scaMu*vecNormConst,
        vecDispEst=rep(scaDispersionEstimate,length(vecCounts)),
        vecDropoutRateEst=vecDropoutRate, 
        vecboolNotZeroObserved=vecboolNotZeroObserved, 
        vecboolZero=vecboolZero )
    } else {
      # Evaluate likelihood on window of cells for each impulse value at a cell.
      # This is a local smoothing penalty added to the cost function.
      scaLogLikNull <- evalLogLikSmoothZINB_comp(vecY=vecCounts,
        vecMu=rep(scaMu,length(vecCounts)),
        vecSizeFactors=vecNormConst,
        vecDispEst=rep(scaDispersionEstimate,length(vecCounts)), 
        vecDropoutRateEst=vecDropoutRate, 
        vecboolNotZeroObserved=vecboolNotZeroObserved, 
        vecboolZero=vecboolZero,
        scaWindowRadius=scaWindowRadius )
    }
  } else {
    stop(paste0("ERROR: Unrecognised strMode in fitImpulse::fitNullModel(): ",strMode))
  }
  return(list( scaMu=scaMu, 
    vecMuLongitudinalSeries=vecMuLongitudinalSeries,
    scaLogLikNull=scaLogLikNull ))
}

#' Estimate impulse model parameter initialisations
#' 
#' The initialisations reflect intuitive parameter choices corresponding
#' to a peak and to a valley model.
#' 
#' @seealso Called by \code{fitImpulse_gene}.
#' 
#' @param vecTimepoints: (numeric vector number of timepoints) 
#'    Time-points at which gene was sampled.
#' @param vecCounts: (count vector number of samples) Count data.
#' #' @param vecDropoutRate: (probability vector number of samples) 
#'    [Default NULL] Dropout rate/mixing probability of zero inflated 
#'    negative binomial mixturemodel for each gene and cell.
#' @param vecProbNB: (probability vector number of samples) [Default NULL]
#'    Probability of observations to come from negative binomial 
#'    component of mixture model.
#' @param vecTimepointAssign: (numeric vector number samples) 
#'    Timepoints assigned to samples.
#' @param vecNormConst: (numeric vector number of samples) 
#'    Normalisation constants for each sample.
#' @param strMode: (str) [Default "batch"] 
#'    {"batch","longitudinal","singlecell"}
#'    Mode of model fitting.
#'    
#' @return vecParamGuessPeak: (numeric vector number of impulse
#'    model parameters) Impulse model parameter initialisation 
#'    corresponding to a peak.
#' @export

estimateImpulseParam <- function(vecTimepoints, 
  vecCounts,
  vecDropoutRate=NULL,
  vecProbNB=NULL,
  strSCMode="clustered",
  vecTimepointAssign, 
  vecNormConst,
  strMode){
  # Compute general statistics for initialisation:
  # Expression means by timepoint
  if(strMode=="batch" | strMode=="longitudinal"){
    vecExpressionMeans <- sapply(vecTimepoints,
      function(tp){mean(vecCounts[vecTimepointAssign==tp], na.rm=TRUE)})
  } else if(strMode=="singlecell"){
    if(strSCMode=="clustered"){
      vecExpressionMeans <- unlist(sapply( vecTimepoints, function(tp){
        sum((vecCounts/vecNormConst*vecProbNB)[vecTimepointAssign==tp], na.rm=TRUE)/
          sum(vecProbNB[vecTimepointAssign==tp], na.rm=TRUE)
      } ))
      # Catch exception: sum(vecProbNB[vecTimepointAssign==tp])==0
      vecExpressionMeans[is.na(vecExpressionMeans)] <- 0
    } else if(strSCMode=="continuous"){
      vecTimepointAssignSort <- sort(vecTimepointAssign, index.return=TRUE)
      vecCountsSort <- vecCounts[vecTimepointAssignSort$ix]
      vecProbNBSort <- vecProbNB[vecTimepointAssignSort$ix]
      scaWindows <- 10
      scaCellsPerClus <- round(length(vecCountsSort)/scaWindows)
      vecExpressionMeans <- array(NA, scaWindows)
      vecTimepoints <- array(NA, scaWindows)
      scaidxNew <- 0
      for(k in seq(1,scaWindows)){
        # Define clusters as groups of cells of uniform size
        scaidxLast <- scaidxNew + 1
        scaidxNew <- scaidxLast + scaCellsPerClus
        # Pick up remaining cells in last cluster
        if(k==scaWindows){scaidxNew=length(vecCountsSort)}
        vecidxK <- seq(scaidxLast, scaidxNew)
        # Infer negative binomial mean parameter
        # Mean estimation here has to be replaced if size factors
        # are used.
        vecExpressionMeans[k] <- sum((vecCountsSort*vecProbNBSort)[vecidxK], na.rm=TRUE)/
          sum(vecProbNBSort[vecidxK], na.rm=TRUE)
        vecTimepoints[k] <- mean((vecTimepointAssignSort$x)[vecidxK], na.rm=TRUE)
      }
      # Catch exception: sum(vecProbNB[vecTimepointAssign==tp])==0
      vecExpressionMeans[is.na(vecExpressionMeans)] <- 0
    } else {
      stop(paste0("ERROR: Unrecognised strSCMode in fitImpulse(): ",strSCMode))
    }
  } else {
    stop(paste0("ERROR: Unrecognised strMode in fitImpulse(): ",strMode))
  }
  nTimepts <- length(vecTimepoints)
  scaMaxMiddleMean <- max(vecExpressionMeans[2:(nTimepts-1)], na.rm=TRUE)
  scaMinMiddleMean <- min(vecExpressionMeans[2:(nTimepts-1)], na.rm=TRUE)
  # +1 to push indicices up from middle stretch to entire window (first is omitted here)
  indMaxMiddleMean <- match(scaMaxMiddleMean,vecExpressionMeans[2:(nTimepts-1)]) + 1
  indMinMiddleMean <- match(scaMinMiddleMean,vecExpressionMeans[2:(nTimepts-1)]) + 1
  # Gradients between neighbouring points
  vecGradients <- unlist( lapply(c(1:(nTimepts-1)),function(x){
    (vecExpressionMeans[x+1]-vecExpressionMeans[x])/(vecTimepoints[x+1]-vecTimepoints[x])}) )
  vecGradients[is.na(vecGradients) | !is.finite(vecGradients)] <- 0
  
  # Compute peak initialisation
  # Beta: Has to be negative, Theta1: Low, Theta2: High, Theta3: Low
  # t1: Around first observed inflexion point, t2: Around second observed inflexion point
  indLowerInflexionPoint <- match(max(vecGradients[1:(indMaxMiddleMean-1)], na.rm=TRUE), vecGradients[1:(indMaxMiddleMean-1)])
  indUpperInflexionPoint <- indMaxMiddleMean - 1 + match(min(vecGradients[indMaxMiddleMean:length(vecGradients)], na.rm=TRUE), vecGradients[indMaxMiddleMean:length(vecGradients)])
  vecParamGuessPeak <- c(1,log(vecExpressionMeans[1]+1),log(scaMaxMiddleMean+1),log(vecExpressionMeans[nTimepts]+1),
    (vecTimepoints[indLowerInflexionPoint]+vecTimepoints[indLowerInflexionPoint+1])/2,
    (vecTimepoints[indUpperInflexionPoint]+vecTimepoints[indUpperInflexionPoint+1])/2)
  
  # Compute valley initialisation
  # Beta: Has to be negative, Theta1: High, Theta2: Low, Theta3: High
  # t1: Around first observed inflexion point, t2: Around second observed inflexion point
  indLowerInflexionPoint <- match(min(vecGradients[1:(indMinMiddleMean-1)], na.rm=TRUE), vecGradients[1:(indMinMiddleMean-1)])
  indUpperInflexionPoint <- indMinMiddleMean - 1 + match(max(vecGradients[indMinMiddleMean:(nTimepts-1)], na.rm=TRUE), vecGradients[indMinMiddleMean:(nTimepts-1)])
  vecParamGuessValley <- c(1,log(vecExpressionMeans[1]+1),log(scaMinMiddleMean+1),log(vecExpressionMeans[nTimepts]+1),
    (vecTimepoints[indLowerInflexionPoint]+vecTimepoints[indLowerInflexionPoint+1])/2,
    (vecTimepoints[indUpperInflexionPoint]+vecTimepoints[indUpperInflexionPoint+1])/2 )
  
  lsParamGuesses <- list(peak=vecParamGuessPeak, valley=vecParamGuessValley)
  return(lsParamGuesses)
}

#' Fit an impulse model to data of a gene
#' 
#' Given a parameter initialisation, this function
#' performs numerical optimisation using BFGS of the 
#' likelihood function given the impulse model and returns
#' the fitted (maximum likelihood) model.
#' 
#' @seealso Called by \code{fitImpulse_gene}. This function
#' calls the cost function \code{evalLogLikImpulseBatch_comp},
#' \code{evalLogLikImpulseTC_comp} or 
#' \code{evalLogLikImpulseSC_comp} depending on the mode.
#' 
#' @param vecParamGuessPeak (vector number of parameters [6]) 
#'    Impulse model parameters.
#' @param vecTimepoints (numeric vector number of timepoints) 
#'    Time-points at which gene was sampled.
#' @param vecCounts (count vector numer of samples) 
#     Observed expression values for  given gene.
#' @param scaDispersionEstimate: (scalar) 
#'    Dispersion estimate for given gene.
#' @param vecDropoutRate: (probability vector number of samples)
#'    [Default NULL]
#'    Dropout rate estimate for each cell for given gene.
#' @param vecNormConst: (numeric vector number of samples) 
#'    Normalisation constants for each sample.
#' @param vecindTimepointAssign (numeric vector number samples) 
#'    Index of time point assigned to sample in list of sorted
#'    time points (vecX).
#' @param vecboolObserved: (bool vector number of samples)
#'    Whether sample is observed (not NA).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' @param vecboolNotZeroObserved: (bool vector number of samples)
#'    Whether sample is not zero and observed (not NA).
#' @param scaWindowRadius: (integer) [Default NULL]
#'    Smoothing interval radius.
#' @param strMode: (str) [Default "batch"] 
#'    {"batch","longitudinal","singlecell"}
#'    Mode of model fitting.
#' @param MAXIT: (scalar) [Default 100] Number of iterations, which are performed 
#'    to fit the impulse model to the clusters.
#'   
#' @return vecFit: (vector [beta, h0, h1, h2, t1, t2, logL_H1, 
#'    converge_H1) Impulse model parameters, likelihood of
#'    data under fitted impulse model with given initialiation
#'    and convergence status of numerical optimisation.
#' @export

optimiseImpulseModelFit <- function(vecParamGuess,
  vecTimepoints, 
  vecCounts,
  scaDispersionEstimate,
  vecDropoutRate=NULL,
  matLinModelPi=matLinModelPi,
  vecNormConst,
  vecindTimepointAssign,
  vecboolObserved,
  vecboolZero=NULL,
  vecboolNotZeroObserved=NULL,
  scaWindowRadius=NULL,
  strMode="batch", 
  MAXIT=100){
  
  # Chose the cost function for optimisation according
  # to the mode strMode.
  if(strMode=="batch" | strMode=="longitudinal"){
    vecFit <- tryCatch({
      unlist( optim(
        par=vecParamGuess,
        fn=evalLogLikImpulseBatch_comp,
        vecX=vecTimepoints,
        vecY=vecCounts,
        scaDispEst=scaDispersionEstimate,
        vecNormConst=vecNormConst,
        vecindTimepointAssign=vecindTimepointAssign,
        vecboolObserved=vecboolObserved,
        method="BFGS",
        control=list(maxit=MAXIT,fnscale=-1)
      )[c("par","value","convergence")] )
    }, error=function(strErrorMsg){
      print(paste0("ERROR: Fitting impulse model: optimiseImpulseModelFit().",
        " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
      print(paste0("vecParamGuess ", paste(vecParamGuess,collapse=" ")))
      print(paste0("vecTimepoints ", paste(vecTimepoints,collapse=" ")))
      print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
      print(paste0("scaDispersionEstimate ", paste(scaDispersionEstimate,collapse=" ")))
      print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
      print(paste0("vecindTimepointAssign ", paste(vecindTimepointAssign,collapse=" ")))
      print(paste0("vecboolObserved ", paste(vecboolObserved,collapse=" ")))
      print(paste0("strMode ", strMode))
      print(paste0("MAXIT ", MAXIT))
      lsErrorCausingGene <- list(vecParamGuess, vecTimepoints, vecCounts, 
        scaDispersionEstimate, vecNormConst, vecindTimepointAssign, vecboolObserved,
        strMode, MAXIT)
      names(lsErrorCausingGene) <- c("vecParamGuess","vecTimepoints","vecCounts", 
        "scaDispersionEstimate", "vecNormConst", "vecindTimepointAssign", "vecboolObserved",
        "strMode", "MAXIT")
      save(lsErrorCausingGene,file=file.path(getwd(),"ImpulseDE2_lsErrorCausingGene.RData"))
      stop(strErrorMsg)
    })
  }else if(strMode=="singlecell"){
    vecFit <- tryCatch({
      unlist( optim(
        par=vecParamGuess, 
        fn=evalLogLikImpulseSC_comp, 
        vecX=vecTimepoints,
        vecY=vecCounts, 
        scaDispEst=scaDispersionEstimate,
        vecDropoutRateEst=vecDropoutRate,
        matLinModelPi=matLinModelPi,
        vecNormConst=vecNormConst,
        vecindTimepointAssign=vecindTimepointAssign,
        vecboolNotZeroObserved=vecboolNotZeroObserved, 
        vecboolZero=vecboolZero,
        scaWindowRadius=scaWindowRadius,
        method="BFGS", 
        control=list(maxit=MAXIT,fnscale=-1)
      )[c("par","value","convergence")] )
    }, error=function(strErrorMsg){
      print(paste0("ERROR: Fitting impulse model: optimiseImpulseModelFit().",
        " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
      print(paste0("vecParamGuess ", paste(vecParamGuess,collapse=" ")))
      print(paste0("vecTimepoints ", paste(vecTimepoints,collapse=" ")))
      print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
      print(paste0("scaDispersionEstimate ", paste(scaDispersionEstimate,collapse=" ")))
      print(paste0("vecDropoutRate ", paste(vecDropoutRate,collapse=" ")))
      print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
      print(paste0("vecindTimepointAssign ", paste(vecindTimepointAssign,collapse=" ")))
      print(paste0("vecboolObserved ", paste(vecboolObserved,collapse=" ")))
      print(paste0("vecboolNotZeroObserved ", paste(vecboolNotZeroObserved,collapse=" ")))
      print(paste0("scaWindowRadius", paste(scaWindowRadius,collapse=" ")))
      print(paste0("strMode ", strMode))
      print(paste0("MAXIT ", MAXIT))
      lsErrorCausingGene <- list(vecParamGuess, vecTimepoints, vecCounts, 
        scaDispersionEstimate, vecDropoutRate, vecNormConst, vecindTimepointAssign, 
        vecboolObserved, scaWindowRadius, strMode, MAXIT)
      names(lsErrorCausingGene) <- c("vecParamGuess","vecTimepoints","vecCounts", 
        "scaDispersionEstimate", "vecDropoutRate", "vecNormConst", "vecindTimepointAssign", 
        "vecboolObserved", "scaWindowRadius", "strMode", "MAXIT")
      save(lsErrorCausingGene,file=file.path(getwd(),"ImpulseDE2_lsErrorCausingGene.RData"))
      stop(strErrorMsg)
    })
  } else {
    stop(paste0("ERROR: Unrecognised strMode in fitImpulse(): ",strMode))
  }
  return(vecFit)
}

#' Fit an impulse model to a single gene
#' 
#' Fits an impulse model and a mean model to a single gene.
#' The optimisation method is set within \code{optimiseImpulseModelFit} 
#' (\code{optim}: BFGS). This function is divided into four parts: 
#' (I) Prepare data,
#' (II) Fit mean model, 
#' (III) Fit impulse model, 
#' (IV) Process Fits.
#' The body of this function is broken up into 4 helper functions, 
#' which are exclusively called by this function.
#' 
#' @seealso Called by \code{fitImpulse_matrix}. This function
#' calls \code{fitNullModel}, \code{estimateImpulseParam} and
#' \code{optimiseImpulseModelFit}.
#' 
#' @param vecCounts: (count vector number of samples) Count data.
#' @param scaDispersionEstimate: (scalar) Inverse of negative binomial 
#'    dispersion coefficients computed by DESeq2 for given gene.
#' @param vecDropoutRate: (probability vector number of samples) 
#'    [Default NULL] Dropout rate/mixing probability of zero inflated 
#'    negative binomial mixturemodel.
#' @param vecProbNB: (probability vector number of samples) [Default NULL]
#'    Probability of observations to come from negative binomial 
#'    component of mixture model.
#' @param vecNormConst: (numeric vector number of samples) 
#'    Model scaling factors for each observation: Take
#'    sequencing depth and longitudinal time series mean
#'    within a gene into account (size and translation
#'    factors).
#' @param vecTimepointAssign: (numeric vector number samples) 
#'    Timepoints assigned to samples.
#' @param vecLongitudinalSeriesAssign (numeric vector number samples) 
#'    Longitudinal series assigned to samples. 
#'    NULL if not operating in strMode="longitudinal".
#' @param dfAnnotationProc: (Table) Processed annotation table. 
#'    Lists co-variables of samples: 
#'    Sample, Condition, Time (numeric), TimeCateg (str)
#'    (and LongitudinalSeries). For internal use.
#' @param strMode: (str) [Default "batch"] 
#'    {"batch","longitudinal","singlecell"}
#'    Mode of model fitting.
#' @param strSCMode: (str) {"clustered", "continuous"}
#'    Mode in which singlecell data are fit: as clusters of
#'    cells to pseudotime centroids or in continuous pseudotime
#'    coordinates.
#' @param scaWindowRadius: (integer) [Default NULL]
#'    Smoothing interval radius.
#' @param NPARAM: (scalar) [Default 6] Number of impulse model parameters
#' @param MAXIT: (scalar) [Default 100] Number of iterations, which are performed 
#'    to fit the impulse model to the clusters.
#' 
#' @return vecBestFitSummary: (vector [beta, h0, h1, h2, t1, t2, logL_H1, 
#'    converge_H1, mu, logL_H0]) Beta to t2 are parameters of the impulse
#'    model, mu is the single parameter of the mean model, logL are
#'    log likelihoods of full (H1) and reduced model (H0) respectively, converge
#'    is convergence status of numerical optimisation of model fitting by
#'    \code{optim} from \code{stats}.
#' @export

fitImpulse_gene <- function(vecCounts, 
  scaDispersionEstimate,
  vecDropoutRate=NULL,
  matLinModelPi=NULL,
  vecProbNB=NULL,
  vecNormConst,
  vecTimepointAssign, 
  vecLongitudinalSeriesAssign, 
  dfAnnotationProc,
  strMode="batch",
  strSCMode="clustered",
  scaWindowRadius=NULL,
  NPARAM=6, 
  MAXIT=1000){
  
  # (I) Process data
  # Get boolean observation vectors:
  vecboolObserved <- !is.na(vecCounts)
  vecboolZero <- vecCounts == 0
  vecboolNotZeroObserved <- !is.na(vecCounts) & vecCounts > 0
  
  # Compute time point specifc parameters
  vecTimepoints <- sort(unique( vecTimepointAssign ))
  # Get vector of numeric time point assignment indices:
  vecindTimepointAssign <- match(vecTimepointAssign, vecTimepoints)
  
  # Compute time course specifc parameters
  if(strMode=="longitudinal"){
    vecLongitudinalSeries <- unique( vecLongitudinalSeriesAssign )
  }
  
  # (II) Fit null model and compute likelihood of null model
  lsNullModel <- fitNullModel(
    vecCounts=vecCounts, 
    scaDispersionEstimate=scaDispersionEstimate,
    vecDropoutRate=vecDropoutRate, 
    vecProbNB=vecProbNB,
    vecNormConst=vecNormConst,
    vecLongitudinalSeries=vecLongitudinalSeries, 
    vecLongitudinalSeriesAssign=vecLongitudinalSeriesAssign,
    vecboolObserved=vecboolObserved, 
    vecboolZero=vecboolZero,
    vecboolNotZeroObserved=vecboolNotZeroObserved,
    strMode=strMode,
    scaWindowRadius=scaWindowRadius )
  
  scaMu <- lsNullModel$scaMu
  vecMuLongitudinalSeries <- lsNullModel$vecMuLongitudinalSeries
  scaLogLikNull <- lsNullModel$scaLogLikNull
  
  # (III) Fit Impulse model
  # Compute initialisations
  lsParamGuesses <- estimateImpulseParam(
    vecTimepoints=vecTimepoints, 
    vecCounts=vecCounts,
    vecDropoutRate=vecDropoutRate, 
    vecProbNB=vecProbNB,
    vecTimepointAssign=vecTimepointAssign, 
    vecNormConst=vecNormConst,
    strMode=strMode,
    strSCMode=strSCMode )
  vecParamGuessPeak <- lsParamGuesses$peak
  vecParamGuessValley <- lsParamGuesses$valley
  
  # 1. Initialisation: Peak
  vecFitPeak <- optimiseImpulseModelFit(
    vecParamGuess=vecParamGuessPeak,
    vecTimepoints=vecTimepoints, 
    vecCounts=vecCounts,
    scaDispersionEstimate=scaDispersionEstimate,
    vecDropoutRate=vecDropoutRate,
    matLinModelPi=matLinModelPi,
    vecNormConst=vecNormConst,
    vecindTimepointAssign=vecindTimepointAssign,
    vecboolObserved=vecboolObserved, 
    vecboolNotZeroObserved=vecboolNotZeroObserved,
    scaWindowRadius=scaWindowRadius,
    strMode=strMode, 
    MAXIT=MAXIT)
  # 2. Initialisation: Valley
  vecFitValley <- optimiseImpulseModelFit(
    vecParamGuess=vecParamGuessValley,
    vecTimepoints=vecTimepoints, 
    vecCounts=vecCounts,
    scaDispersionEstimate=scaDispersionEstimate,
    vecDropoutRate=vecDropoutRate,
    matLinModelPi=matLinModelPi,
    vecNormConst=vecNormConst,
    vecindTimepointAssign=vecindTimepointAssign,
    vecboolObserved=vecboolObserved, 
    vecboolNotZeroObserved=vecboolNotZeroObserved,
    scaWindowRadius=scaWindowRadius,
    strMode=strMode, 
    MAXIT=MAXIT)
  
  # (IV) Process fits
  dfFitsByInitialisation <- cbind(vecFitPeak, vecFitValley)
  
  # Select best fit and report fit type
  # Report mean fit objective value as null hypothesis, too.
  # match() selects first hit if maximum occurs multiple times
  indBestFit <- match(max(dfFitsByInitialisation["value",]),dfFitsByInitialisation["value",])
  if(strMode=="batch" | strMode=="singlecell"){
    vecBestFitSummary <- c(dfFitsByInitialisation[,indBestFit],scaMu,scaLogLikNull)
    names(vecBestFitSummary) <- c("beta","h0","h1","h2","t1","t2",
      "logL_H1","converge_H1","mu","logL_H0")
  } else if(strMode=="longitudinal") {
    vecBestFitSummary <- c(dfFitsByInitialisation[,indBestFit],
      scaMu,
      scaLogLikNull,
      vecMuLongitudinalSeries)
    vecColnamesMubyLongitudinalSeries <- paste0(rep("mu_",length(vecMuLongitudinalSeries)), names(vecMuLongitudinalSeries))
    names(vecBestFitSummary) <- c("beta","h0","h1","h2","t1","t2",
      "logL_H1","converge_H1", "mu", "logL_H0", vecColnamesMubyLongitudinalSeries)
  }  else {
    stop(paste0("ERROR: Unrecognised strMode in fitImpulse(): ",strMode))
  }
  # Remove log scale from count parameters of impulse model
  vecBestFitSummary[c("h0","h1","h2")] <- exp(vecBestFitSummary[c("h0","h1","h2")])
  
  return(vecBestFitSummary)
}

#' Fits impulse models to all genes of a dataset
#' 
#' Fits impulse models to all genes of a dataset. Performs parallelisation.
#' Internal function of \code{fitImpulse}.
#' 
#' @details Maximum number of iterations for optimisation is set within this 
#'    function as \code{MAXIT}. In the case of single condition differential 
#'    expression over time, this function is called once for the case condition.
#'    In the case of case and control condition, this function is called three 
#'    times: data from case, control and combined conditions.
#'    Calls \code{fitImpulse_gene}.
#' 
#' @seealso Called by \code{fitImpulse}. Calls \code{fitImpulse_gene}.
#'   
#' @param matCountDataProc: (count matrix genes x samples) Count data.
#' @param vecDispersions: (vector number of genes) Gene-wise 
#'    negative binomial dispersion coefficients.
#' @param matDropoutRate: (probability matrix genes x samples) Dropout 
#'    rate/mixing probability of zero inflated negative binomial mixture 
#'    model for each gene and cell.
#' @param matProbNB: (probability vector genes x samples) 
#'    Probability of observations to come from negative binomial 
#'    component of mixture model.
#' @param matNormConst: (numeric matrix genes x samples) 
#'    Model scaling factors for each observation: Take
#'    sequencing depth and longitudinal time series mean
#'    within a gene into account (size and translation
#'    factors).
#' @param vecTimepointAssign: (numeric vector number samples) Timepoints 
#'    assigned to samples.
#' @param vecLongitudinalSeriesAssign: (numeric vector number samples) Time courses 
#'    assigned to samples. NULL if not operating in strMode="longitudinal".
#' @param dfAnnotationProc: (Table) Processed annotation table. 
#'    Lists co-variables of samples: 
#'    Sample, Condition, Time (numeric), TimeCateg (str)
#'    (and LongitudinalSeries). For internal use.
#' @param strCaseName: (str) Name of the case condition in \code{dfAnnotationRedFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#'    \code{dfAnnotationRedFull}.
#' @param strMode: (str) [Default "batch"] 
#'    {"batch","longitudinal","singlecell"}
#'    Mode of model fitting.
#' @param strSCMode: (str) {"clustered", "continuous"}
#'    Mode in which singlecell data are fit: as clusters of
#'    cells to pseudotime centroids or in continuous pseudotime
#'    coordinates.
#' @param scaWindowRadius: (integer) [Default NULL]
#'    Smoothing interval radius.
#' @param nProc: (scalar) [Default 3] Number of processes for parallelisation.
#' @param NPARAM: (scalar) [Default 6] Number of parameters of impulse model.
#' 
#' @return lsFitResults_matrix (list length 2) List of two matrices which
#'    contain parameter fits and model values for given time course for the
#'    given condition.
#' \itemize{
#'    \item \code{parameter}(genes x [beta, h0, h1, h2, t1, t2, logL_H1, 
#'        converge_H1, mu, logL_H0, converge_H0]) 
#'        Beta to t2 are parameters of the impulse model, mu is the single 
#'        parameter of the mean model, logL are log likelihoods of full 
#'        (H1) and reduced model (H0) respectively, converge is convergence 
#'        status of numerical optimisation of model fitting by
#'        \code{optim} from \code{stats} of either model.
#'    \item \code{values} (genes x time points) Contains the counts 
#'        predicted by the impulse model at the observed time points.
#' }
#' @export

fitImpulse_matrix <- function(matCountDataProcCondition, 
  vecDispersions, 
  matDropoutRate=NULL, 
  matProbNB=NULL,
  matNormConst, 
  vecTimepointAssign, 
  vecLongitudinalSeriesAssign, 
  dfAnnotationProc,
  strCaseName, 
  strControlName=NULL, 
  strMode="batch",
  strSCMode="clustered",
  scaWindowRadius=NULL,
  nProc=1, 
  NPARAM=6){
  
  # Maximum number of iterations for numerical optimisation of
  # likelihood function in MLE fitting of impulse model:
  MAXIT <- 1000
  
  if(nrow(matCountDataProcCondition) > max(2*nProc,10)){
    # Use parallelisation if number of genes/centroids to fit is large
    
    # Define partitioning of genes onto nodes: lsGeneIndexByCore
    lsGeneIndexByCore = list()
    bord = floor(nrow(matCountDataProcCondition)/nProc)
    for (i in 1:nProc){
      lsGeneIndexByCore[[i]] <- ((i-1)*bord+1):(i*bord)
    }
    if(nProc*bord < nrow(matCountDataProcCondition)){
      # Add remaining genes to last node
      lsGeneIndexByCore[[nProc]] <-  c(lsGeneIndexByCore[[nProc]],(nProc*bord+1):nrow(matCountDataProcCondition))
    }
    
    cl <- makeCluster(nProc, outfile = "ImpulseDE2_ClusterOut.txt")
    my.env <- new.env()
    # Variables
    assign("matCountDataProcCondition", matCountDataProcCondition, envir = my.env)
    assign("vecDispersions", vecDispersions, envir = my.env)
    assign("matDropoutRate", matDropoutRate, envir = my.env)
    assign("matProbNB", matProbNB, envir = my.env)
    assign("matNormConst", matNormConst, envir = my.env)
    assign("vecTimepointAssign", vecTimepointAssign, envir = my.env)
    assign("vecLongitudinalSeriesAssign", vecLongitudinalSeriesAssign, envir = my.env)
    assign("dfAnnotationProc", dfAnnotationProc, envir = my.env)
    assign("strMode", strMode, envir = my.env)
    assign("strSCMode", strMode, envir = my.env)
    assign("scaWindowRadius", scaWindowRadius, envir = my.env)
    assign("NPARAM", NPARAM, envir = my.env)
    assign("MAXIT", MAXIT, envir = my.env)
    # Functions
    assign("fitImpulse_gene",fitImpulse_gene, envir = my.env)
    assign("calcImpulse_comp", calcImpulse_comp, envir = my.env)
    assign("fitMuNB", fitMuNB, envir = my.env)
    assign("evalLogLikMuNB_comp", evalLogLikMuNB_comp, envir = my.env)
    assign("evalLogLikImpulseBatch_comp", evalLogLikImpulseBatch_comp, envir = my.env)
    assign("evalLogLikZINB_comp", evalLogLikZINB_comp, envir = my.env)
    assign("evalLogLikSmoothZINB_comp", evalLogLikSmoothZINB_comp, envir = my.env)
    assign("fitMuZINB", fitMuZINB, envir = my.env)
    assign("evalLogLikMuZINB_comp", evalLogLikMuZINB_comp, envir = my.env)
    assign("evalLogLikImpulseSC_comp", evalLogLikImpulseSC_comp, envir = my.env)
    assign("fitNullModel", fitNullModel, envir = my.env)
    assign("estimateImpulseParam", estimateImpulseParam, envir = my.env)
    assign("optimiseImpulseModelFit", optimiseImpulseModelFit, envir = my.env)
    
    clusterExport(cl=cl, varlist=c(
      "matCountDataProcCondition",
      "vecDispersions",
      "matDropoutRate",
      "matProbNB",
      "matNormConst",
      "vecTimepointAssign",
      "vecLongitudinalSeriesAssign",
      "dfAnnotationProc",
      "strMode",
      "strSCMode",
      "scaWindowRadius",
      "MAXIT",
      "NPARAM",
      "calcImpulse_comp",
      "fitMuNB",
      "evalLogLikMuNB_comp",
      "evalLogLikImpulseBatch_comp",
      "evalLogLikZINB_comp",
      "evalLogLikSmoothZINB_comp",
      "fitMuZINB",
      "evalLogLikMuZINB_comp",
      "evalLogLikImpulseSC_comp", 
      "fitImpulse_gene",
      "fitNullModel",
      "estimateImpulseParam",
      "optimiseImpulseModelFit"
    ), envir = my.env)
    
    # Set default length of output vector for each gene,
    # used upon encountering error in optimisation for gracefull exit.
    # 6 parameteres, 4 meta, one mean for each longitudinal series
    scaSizeReport <- 6 + 4
    if(strMode=="longitudinal"){
      scaNumberLongitudinalSeries <- length(unique( vecLongitudinalSeriesAssign ))
      scaSizeReport <- scaSizeReport + scaNumberLongitudinalSeries
    }
    # Fit impulse model to each gene of matrix and get impulse parameters:
    # clusterApply runs the function impulse_fit_gene_wise
    # The input data are distributed to nodes by lsGeneIndexByCore partitioning
    if(strMode=="singlecell"){
      # Call fitting with single cell parameters
      lsmatFits <- clusterApply(cl, 1:length(lsGeneIndexByCore),
        function(z){ t(sapply( lsGeneIndexByCore[[z]],
          function(x){
            impulsefitGene <- NULL
            tryCatch({
              impulsefitGene <- fitImpulse_gene(
                vecCounts=matCountDataProcCondition[x,],
                scaDispersionEstimate=vecDispersions[x],
                vecDropoutRate=matDropoutRate[x,],
                vecProbNB=matProbNB[x,],
                vecNormConst=matNormConst[x,],
                vecTimepointAssign=vecTimepointAssign,
                vecLongitudinalSeriesAssign=vecLongitudinalSeriesAssign,
                dfAnnotationProc=dfAnnotationProc,
                strMode=strMode,
                strSCMode=strSCMode,
                scaWindowRadius=scaWindowRadius,
                NPARAM=NPARAM,
                MAXIT=MAXIT )
            }, error=function(strErrorMsg){
              print("ERROR: Impulse fitting failed in one instance with error message:")
              print(strErrorMsg)
              warning(paste0("Fitting of impuluse model failed for one gene.",
                " Consider contacting developers if you care about all genes.",
                " david.seb.fischer@gmail.com"))
            }, finally={
              if(is.null(impulsefitGene)){
                impulsefitGene <- array(NA, scaSizeReport)
              }
            })
            return(impulsefitGene)
          }))
        })
    } else {
      # Call fitting without single cell parameters
      lsmatFits <- clusterApply(cl, 1:length(lsGeneIndexByCore),
        function(z){ t(sapply( lsGeneIndexByCore[[z]],
          function(x){
            impulsefitGene <- NULL
            tryCatch({
              impulsefitGene <- fitImpulse_gene(
                vecCounts=matCountDataProcCondition[x,],
                scaDispersionEstimate=vecDispersions[x],
                vecNormConst=matNormConst[x,],
                vecTimepointAssign=vecTimepointAssign,
                vecLongitudinalSeriesAssign=vecLongitudinalSeriesAssign,
                dfAnnotationProc=dfAnnotationProc,
                strMode=strMode,
                NPARAM=NPARAM,
                MAXIT=MAXIT )
            }, error=function(strErrorMsg){
              print("ERROR: Impulse fitting failed in one instance with error message:")
              print(strErrorMsg)
              warning(paste0("Fitting of impuluse model failed for one gene.",
                " Consider contacting developers if you care about all genes.",
                " david.seb.fischer@gmail.com"))
            }, finally={
              if(is.null(impulsefitGene)){
                impulsefitGene <- array(NA, scaSizeReport)
              }
            })
            return(impulsefitGene)
          }))
        })      
    }
    # Give output rownames again, which are lost above
    for(i in 1:length(lsGeneIndexByCore)){
      rownames(lsmatFits[[i]]) <- rownames(matCountDataProcCondition[lsGeneIndexByCore[[i]],])
    }
    # Concatenate the output objects of each node
    matFits <- do.call(rbind,lsmatFits)
    # Parallelisation ends here
    stopCluster(cl)
    
  } else {
    # Do not use parallelisation if number of genes to fit is small
    
    # Fit impulse model to each gene of matrix and get impulse parameters
    lsGeneIndexByCore_all <- 1:nrow(matCountDataProcCondition)
    if(strMode=="singlecell"){
      # Call fitting with single cell parameters
      matFits <- lapply(lsGeneIndexByCore_all,function(x){
        fitImpulse_gene(
          vecCounts=matCountDataProcCondition[x,],
          scaDispersionEstimate=vecDispersions[x],
          vecDropoutRate=matDropoutRate[x,],
          vecProbNB=matProbNB[x,],
          vecNormConst=matNormConst[x,],
          vecTimepointAssign=vecTimepointAssign,
          vecLongitudinalSeriesAssign=vecLongitudinalSeriesAssign,
          dfAnnotationProc=dfAnnotationProc,
          strMode=strMode,
          strSCMode=strSCMode,
          scaWindowRadius=scaWindowRadius,
          NPARAM=NPARAM,
          MAXIT=MAXIT )})
    } else {
      # Call fitting without single cell parameters
      matFits <- lapply(lsGeneIndexByCore_all,function(x){
        fitImpulse_gene(
          vecCounts=matCountDataProcCondition[x,],
          scaDispersionEstimate=vecDispersions[x],
          vecNormConst=matNormConst[x,],
          vecTimepointAssign=vecTimepointAssign,
          vecLongitudinalSeriesAssign=vecLongitudinalSeriesAssign,
          dfAnnotationProc=dfAnnotationProc,
          strMode=strMode,
          NPARAM=NPARAM,
          MAXIT=MAXIT )})
    }
    # Call fitting without single cell parameters
    matFits <- do.call(rbind,matFits)
    rownames(matFits) <- rownames(matCountDataProcCondition)
  }
  
  # Use obtained impulse parameters to calculate impulse fit values
  matTheta <- matFits[,1:NPARAM]
  matTheta[,2:4] <- log(matTheta[,2:4])
  if(nrow(matFits) == 1){
    # If matrix contains only one gene
    matImpulseValues <- as.data.frame(t(calcImpulse_comp(matTheta,
      unique(sort(vecTimepointAssign)))),row.names=rownames(matFits))
  } else {                    
    # If matrix contains > 1 genes
    matImpulseValues <- t(apply(matTheta,1,function(x){calcImpulse_comp(x,
      unique(sort(vecTimepointAssign)))}))
  } 
  colnames(matImpulseValues) = unique(sort(vecTimepointAssign))
  rownames(matImpulseValues) <- rownames(matCountDataProcCondition)
  
  # Report results.
  lsFitResults_matrix <- list()
  lsFitResults_matrix[[1]] <- matFits
  lsFitResults_matrix[[2]] <- matImpulseValues
  names(lsFitResults_matrix) <- c("parameters","values")
  return(lsFitResults_matrix)
}

#' Fits impulse model to a timecourse dataset
#' 
#' This function processes the input matrix and
#' coordinates impulse model fitting through
#' \code{impulse_fit_matrix}.
#' [Helper \code{fitImpulse_matrix}] Fit impulse model to matrix of genes.
#' Calls \code{fitImpulse_gene}. \code{fitImpulse_matrix} fits impulse
#' models to a matrix of samples from one condition.
#'  
#' @seealso Called by \code{runImpulseDE2}. Calls \code{fitImpulse_matrix}.
#'  
#' @param matCountDataProc: (matrix genes x samples)
#'    Count data: Reduced version of \code{matCountData}. 
#'    For internal use.
#' @param dfAnnotationProc: (Table) Processed annotation table. 
#'    Lists co-variables of samples: 
#'    Sample, Condition, Time (numeric), TimeCateg (str)
#'    (and LongitudinalSeries). For internal use.
#' @param matTranslationFactors: (numeric matrix genes x samples) 
#'    Model scaling factors for each observation which take
#'    longitudinal time series mean within a gene into account 
#'    (translation factors). Computed based based on all samples.
#' @param matSizeFactors: (numeric matrix genes x samples) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per sample - rows of this matrix are equal.
#' @param vecDispersions: (vector number of genes) Gene-wise 
#'    negative binomial dispersion coefficients.
#' @param matDropoutRate: (probability matrix genes x samples) Dropout 
#'    rate/mixing probability of zero inflated negative binomial mixture 
#'    model for each gene and cell.
#' @param matProbNB: (probability matrix genes x samples) 
#'    Probability of observations to come from negative binomial 
#'    component of mixture model.
#' @param strCaseName (str) Name of the case condition in \code{dfAnnotationRedFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#'    \code{dfAnnotationRedFull}.
#' @param strMode: (str) [Default "batch"] 
#'    {"batch","longitudinal","singlecell"}
#'    Mode of model fitting.
#' @param strSCMode: (str) {"clustered", "continuous"}
#'    Mode in which singlecell data are fit: as clusters of
#'    cells to pseudotime centroids or in continuous pseudotime
#'    coordinates.
#' @param scaWindowRadius: (integer) [Default NULL]
#'    Smoothing interval radius.
#' @param nProc: (scalar) [Default 1] Number of processes for parallelisation.
#' @param NPARAM: (scalar) [Default 6] Number of parameters of impulse model.
#' 
#' @return lsFitResults_all (list length 2 or 6) List of matrices which
#'    contain parameter fits and model values for given time course for the
#'    case condition (and control and combined if control is present).
#'    Each parameter matrix is called parameter_'condition' and has the form
#'    (genes x [beta, h0, h1, h2, t1, t2, logL_H1, converge_H1, mu, logL_H0, 
#'    converge_H0]) where beta to t2 are parameters of the impulse
#'    model, mu is the single parameter of the mean model, logL are
#'    log likelihoods of full (H1) and reduced model (H0) respectively, converge
#'    is convergence status of numerical optimisation of model fitting by
#'    \code{optim} from \code{stats} of either model. Each value matrix is called
#'    value_'condition' and has the form (genes x time points) and contains the
#'    counts predicted by the impulse model at the observed time points.
#' @export

fitImpulse <- function(matCountDataProc, 
  dfAnnotationProc, 
  matTranslationFactors, 
  matSizeFactors,
  vecDispersions, 
  matDropoutRate, 
  matProbNB,
  vecClusterAssignments,
  strCaseName, 
  strControlName=NULL, 
  strMode="batch",
  strSCMode="clustered",
  scaWindowRadius=NULL,
  nProc=1, 
  NPARAM=6){
  
  #g_names = rownames(arr3DCountData)
  lsFitResults_all = list()
  
  # Condition labels to be used in runs:
  # These labels are used to group samples
  if(!is.null(strControlName)){
    lsLabels <- c("combined","case","control")
    names(lsLabels) <- c("combined","case","control")
  } else {
    lsLabels <- c("case")
    names(lsLabels) <- c("case")
  }
  
  # Create lists of samples to be used per run
  if(!is.null(strControlName)){    
    lsSamplesByCond <- list(
      colnames(matCountDataProc),
      (colnames(matCountDataProc))[dfAnnotationProc[match(
        colnames(matCountDataProc),
        dfAnnotationProc$Sample),
        ]$Condition %in% strCaseName],
      (colnames(matCountDataProc))[dfAnnotationProc[match(
        colnames(matCountDataProc),
        dfAnnotationProc$Sample),
        ]$Condition %in% strControlName] )
    names(lsSamplesByCond) <- c("combined","case","control")   
  } else {
    lsSamplesByCond <- list(
      colnames(matCountDataProc) )
    names(lsSamplesByCond) <- c("case")
  }
  # Get time course assignments of samples
  if(strMode=="longitudinal"){
    vecLongitudinalSeriesAssign <- dfAnnotationProc[match(
      colnames(matCountDataProc),
      dfAnnotationProc$Sample),]$LongitudinalSeries
    names(vecLongitudinalSeriesAssign) <- colnames(matCountDataProc)
  } else {
    vecLongitudinalSeriesAssign <- NULL
  }
  # Get time point assignments of samples
  vecTimepointAssign <- dfAnnotationProc[match(
    colnames(matCountDataProc),
    dfAnnotationProc$Sample),]$Time
  names(vecTimepointAssign) <- colnames(matCountDataProc)
  if(strMode=="longitudinal"){
    matNormConst <- matSizeFactors*matTranslationFactors
  } else {
    matNormConst <- matSizeFactors
  }
  
  # Fitting for different runs
  for (label in lsLabels){
    if(strMode=="singlecell"){
      # Call fitting with single cell parameters
      lsFitResults_run <- fitImpulse_matrix(
        matCountDataProcCondition=matCountDataProc[,lsSamplesByCond[[label]]],
        matNormConst=matNormConst[,lsSamplesByCond[[label]]],
        vecDispersions=vecDispersions,
        matDropoutRate=matDropoutRate[,lsSamplesByCond[[label]]],
        matProbNB=matProbNB[,lsSamplesByCond[[label]]],
        vecTimepointAssign=vecTimepointAssign[lsSamplesByCond[[label]]],
        vecLongitudinalSeriesAssign=vecLongitudinalSeriesAssign[lsSamplesByCond[[label]]],
        dfAnnotationProc=dfAnnotationProc,
        strCaseName=strCaseName,
        strControlName=strControlName,
        strMode=strMode,
        strSCMode=strSCMode,
        scaWindowRadius=scaWindowRadius,
        nProc=nProc,
        NPARAM=NPARAM )
    } else {
      # Call fitting without singlecell parameters
      lsFitResults_run <- fitImpulse_matrix(
        matCountDataProcCondition=matCountDataProc[,lsSamplesByCond[[label]]],
        matNormConst=matNormConst[,lsSamplesByCond[[label]]],
        vecDispersions=vecDispersions,
        vecTimepointAssign=vecTimepointAssign[lsSamplesByCond[[label]]],
        vecLongitudinalSeriesAssign=vecLongitudinalSeriesAssign[lsSamplesByCond[[label]]],
        dfAnnotationProc=dfAnnotationProc,
        strCaseName=strCaseName,
        strControlName=strControlName,
        strMode=strMode,
        nProc=nProc,
        NPARAM=NPARAM )
    }
    
    names(lsFitResults_run) <- paste(c("parameters","values"),lsLabels[label],sep="_")
    lsFitResults_all[[names(lsFitResults_run)[1]]] <- lsFitResults_run[[1]]
    lsFitResults_all[[names(lsFitResults_run)[2]]] <- lsFitResults_run[[2]]
  }
  
  return(lsFitResults_all)
}