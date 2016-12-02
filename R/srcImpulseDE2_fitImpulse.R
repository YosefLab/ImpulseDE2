#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Impulse model fit    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Estimate impulse model parameter initialisations
#' 
#' The initialisations reflect intuitive parameter choices corresponding
#' to a peak and to a valley model.
#' 
#' @seealso Called by \code{fitImpulse_gene}.
#' 
#' @param vecCounts: (count vector number of samples) Count data.
#' @param vecTimepoints: (numeric vector number samples) 
#'    Timepoints assigned to samples.
#' @param vecSizeFactors: (numeric vector number of samples) 
#'    Normalisation constants for each sample.
#' @param strMode: (str) [Default "batch"] 
#'    {"batch","longitudinal","singlecell"}
#'    Mode of model fitting.
#'    
#' @return vecParamGuessPeak: (numeric vector number of impulse
#'    model parameters) Impulse model parameter initialisation 
#'    corresponding to a peak.
#' @export

estimateImpulseParam <- function(vecCounts,
                                 vecTimepoints,
                                 vecSizeFactors,
                                 vecindBatchAssign,
                                 strMode){
  # Compute general statistics for initialisation:
  # Expression means by timepoint
  vecCountsSFcorrected <- vecCounts/vecSizeFactors
  if(strMode=="batcheffects"){
    # Estimate batch factors
    vecBatchFactors <- sapply(unique(vecindBatchAssign), function(batch){
      mean(vecCountsSFcorrected[vecindBatchAssign==batch]/mean(vecCounts, na.rm=TRUE), na.rm=TRUE)
    })
    # Catch exception that all observations of a batch are zero:
    vecBatchFactors[is.na(vecBatchFactors) | vecBatchFactors==0] <- 1
    vecCountsSFBatchcorrected <- vecCountsSFcorrected/vecBatchFactors[vecindBatchAssign]
    vecExpressionMeans <- sapply(vecTimepoints, function(tp){
      mean(vecCountsSFBatchcorrected[vecTimepoints==tp], na.rm=TRUE)
    })
  } else {
    vecExpressionMeans <- sapply(vecTimepoints, function(tp){
      mean(vecCountsSFcorrected[vecTimepoints==tp], na.rm=TRUE)
    })
  }
  scaNTimepoints <- length(vecTimepoints)
  scaMaxMiddleMean <- max(vecExpressionMeans[2:(scaNTimepoints-1)], na.rm=TRUE)
  scaMinMiddleMean <- min(vecExpressionMeans[2:(scaNTimepoints-1)], na.rm=TRUE)
  # +1 to push indicices up from middle stretch to entire window (first is omitted here)
  indMaxMiddleMean <- match(scaMaxMiddleMean,vecExpressionMeans[2:(scaNTimepoints-1)]) + 1
  indMinMiddleMean <- match(scaMinMiddleMean,vecExpressionMeans[2:(scaNTimepoints-1)]) + 1
  # Gradients between neighbouring points
  vecGradients <- unlist( lapply(c(1:(scaNTimepoints-1)),function(x){
    (vecExpressionMeans[x+1]-vecExpressionMeans[x])/(vecTimepoints[x+1]-vecTimepoints[x])}) )
  vecGradients[is.na(vecGradients) | !is.finite(vecGradients)] <- 0
  
  # Compute peak initialisation
  # Beta: Has to be negative, Theta1: Low, Theta2: High, Theta3: Low
  # t1: Around first observed inflexion point, t2: Around second observed inflexion point
  indLowerInflexionPoint <- match(
    max(vecGradients[1:(indMaxMiddleMean-1)], na.rm=TRUE), 
    vecGradients[1:(indMaxMiddleMean-1)])
  indUpperInflexionPoint <- indMaxMiddleMean - 1 + match(
    min(vecGradients[indMaxMiddleMean:length(vecGradients)], na.rm=TRUE), 
    vecGradients[indMaxMiddleMean:length(vecGradients)])
  vecParamGuessPeak <- c(1,log(vecExpressionMeans[1]+1),
                         log(scaMaxMiddleMean+1),log(vecExpressionMeans[scaNTimepoints]+1),
                         (vecTimepoints[indLowerInflexionPoint]+vecTimepoints[indLowerInflexionPoint+1])/2,
                         (vecTimepoints[indUpperInflexionPoint]+vecTimepoints[indUpperInflexionPoint+1])/2)
  
  # Compute valley initialisation
  # Beta: Has to be negative, Theta1: High, Theta2: Low, Theta3: High
  # t1: Around first observed inflexion point, t2: Around second observed inflexion point
  indLowerInflexionPoint <- match(
    min(vecGradients[1:(indMinMiddleMean-1)], na.rm=TRUE), 
    vecGradients[1:(indMinMiddleMean-1)])
  indUpperInflexionPoint <- indMinMiddleMean - 1 + match(
    max(vecGradients[indMinMiddleMean:(scaNTimepoints-1)], na.rm=TRUE), 
    vecGradients[indMinMiddleMean:(scaNTimepoints-1)])
  vecParamGuessValley <- c(1,
                           log(vecExpressionMeans[1]+1),
                           log(scaMinMiddleMean+1),
                           log(vecExpressionMeans[scaNTimepoints]+1),
                           (vecTimepoints[indLowerInflexionPoint]+vecTimepoints[indLowerInflexionPoint+1])/2,
                           (vecTimepoints[indUpperInflexionPoint]+vecTimepoints[indUpperInflexionPoint+1])/2 )
  
  lsParamGuesses <- list(peak=vecParamGuessPeak, valley=vecParamGuessValley)
  return(lsParamGuesses)
}

#' Fit negative binomial mean parameter
#' 
#' Numerical optimisation of negative binomial mean model fit.
#' Note that the closed form solution of the maximum likelihood 
#' estimator of the negative binomial mean parameter 
#' (the weighted average) only holds if all normalisation
#' factors are 1. Catches numerical errors. Fit in log space
#' to guarantee positive mean paramter.
#' 
#' @seealso Called by \code{computeTranslationFactors()}
#' and \code{computeLogLikNull()}.
#' 
#' @param vecCounts (count vector samples) 
#'    Observed expression values for given gene.
#' @param scaDisp: (scalar) Dispersion estimate for given gene.
#' @param vecSizeFactors: (numeric vector number of samples) 
#'    Normalisation constants for each sample.
#' @param vecWeights: (probability vector number of samples) 
#'    Weights for inference on mixture models.
#' @param vecboolObserved: (bool vector number of samples)
#'    Stores bool of sample being not NA (observed).
#'     
#' @return scaMu: (scalar) Maximum likelihood estimator of
#'    negative binomial mean parameter.
#' @export

fitConstModel <- function(vecCounts,
                          scaDisp,
                          vecSizeFactors,
                          vecBatchesUnique,
                          vecindBatchAssign,
                          strMode,
                          MAXIT=100,
                          RELTOL=10^(-8),
                          trace=0,
                          REPORT=10){
  
  if(all(vecSizeFactors==1) & strMode=="singlebatch"){
    # Closed form maximum likelihood estimator
    scaMu <- log(mean(vecCounts, na.rm=TRUE))
  } else {
    # Numerical maximum likelihood estimator
    if(strMode=="batcheffects"){
      # One mean parameter and #batches-1 batch factors
      # initialised to 1 (0 in log space)
      vecParamGuess <- c(log(mean(vecCounts, na.rm=TRUE)), 
                         rep(0, length(vecBatchesUnique)-1) )
    } else { vecParamGuess <- log(mean(vecCounts, na.rm=TRUE)) }
  }
  
  lsFit <- tryCatch({
    optim(
      par=vecParamGuess,
      evalLogLikMu_comp,
      vecCounts=vecCounts,
      scaDisp=scaDisp, 
      vecSizeFactors=vecSizeFactors,
      vecindBatch=vecindBatchAssign,
      strMode=strMode,
      vecboolObserved=!is.na(vecCounts),
      method="BFGS",
      control=list(maxit=MAXIT,
                   reltol=RELTOL,
                   fnscale=-1)
    )[c("par","value","convergence")]
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting null model: fitrMu().",
                 " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
    print(paste0("vecParamGuess ", paste(vecParamGuess,collapse=" ")))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("scaDisp ", paste(scaDisp,collapse=" ")))
    print(paste0("vecSizeFactors ", paste(vecSizeFactors,collapse=" ")))
    print(paste0("strMode ", strMode))
    print(paste0("MAXIT ", MAXIT))
    print(strErrorMsg)
    stop(strErrorMsg)
  })
  
  # Extract parameter estimates
  scaMu <- exp(lsFit$par[1])
  # Catch boundary of likelihood domain on mu space:
  if(scaMu < 10^(-10)){scaMu <- 10^(-10)}
  if(strMode=="batcheffects"){
    vecBatchFactors <- c(1, exp(lsFit$par[2:length(lsFit$par)]))
    # Catch boundary of likelihood domain on batch factor space:
    vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
    vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
    names(vecBatchFactors) <- vecBatchesUnique
  } else { vecBatchFactors <- NULL }
  
  return(list( scaMu=scaMu,
               vecBatchFactors=vecBatchFactors,
               scaLL=lsFit$value,
               scaConvergence=lsFit$convergence))
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
#' @param scaDisp: (scalar) 
#'    Dispersion estimate for given gene.
#' @param vecDropoutRate: (probability vector number of samples)
#'    [Default NULL]
#'    Dropout rate estimate for each cell for given gene.
#' @param vecSizeFactors: (numeric vector number of samples) 
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

fitImpulseModel <- function(vecImpulseParamGuess,
                            vecCounts,
                            scaDisp,
                            vecSizeFactors,
                            vecBatchesUnique,
                            vecindBatchAssign,
                            vecTimepointsUnique,
                            vecindTimepointAssign,
                            strMode, 
                            MAXIT=100,
                            RELTOL=10^(-8),
                            trace=0,
                            REPORT=10 ){
  
  
  if(strMode=="batcheffects"){
    # Impulse parameter guesses (6) and #batches-1 batch factor parameters
    # initialised to 1 (0 in log space)
    vecParamGuess <- c(vecImpulseParamGuess, rep(0, length(vecBatchesUnique)-1))
  } else { vecParamGuess <- vecImpulseParamGuess }
  
  lsFit <- tryCatch({
    optim(
      par=vecParamGuess,
      fn=evalLogLikImpulse_comp,
      vecCounts=vecCounts,
      scaDisp=scaDisp,
      vecSizeFactors=vecSizeFactors,
      vecTimepointsUnique=vecTimepointsUnique,
      vecindTimepointAssign=vecindTimepointAssign,
      vecindBatch=vecindBatchAssign,
      strMode=strMode,
      vecboolObserved=!is.na(vecCounts),
      method="BFGS",
      control=list(maxit=MAXIT,
                   reltol=RELTOL,
                   fnscale=-1)
    )[c("par","value","convergence")]
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting impulse model: optimiseImpulseModelFit().",
                 " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
    print(paste0("vecParamGuess ", paste(vecParamGuess,collapse=" ")))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("scaDisp ", paste(scaDisp,collapse=" ")))
    print(paste0("vecSizeFactors ", paste(vecSizeFactors,collapse=" ")))
    print(paste0("vecTimepointsUnique ", paste(vecTimepointsUnique,collapse=" ")))
    print(paste0("vecindTimepointAssign ", paste(vecindTimepointAssign,collapse=" ")))
    print(paste0("strMode ", strMode))
    print(paste0("MAXIT ", MAXIT))
    print(strErrorMsg)
    stop(strErrorMsg)
  })
  
  # Extract parameter estimates
  vecImpulseParam <- lsFit$par[1:6]
  vecImpulseParam[2:4] <- exp(vecImpulseParam[2:4])
  names(vecImpulseParam) <- c("beta", "h0", "h1", "h2", "t1", "t2")
  vecImpulseValue <- calcImpulse_comp(vecImpulseParam=vecImpulseParam,
                                      vecTimepoints=vecTimepointsUnique)[vecindTimepointAssign]
  names(vecImpulseValue) <- names(vecCounts)
  if(strMode=="batcheffects"){
    vecBatchFactors <- c(1, exp(lsFit$par[7:length(lsFit$par)]))
    # Catch boundary of likelihood domain on batch factor space:
    vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
    vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
    names(vecBatchFactors) <- vecBatchesUnique
  } else { vecBatchFactors <- NULL }
  
  return(list( vecImpulseParam=vecImpulseParam,
               vecImpulseValue=vecImpulseValue,
               vecBatchFactors=vecBatchFactors,
               scaLL=lsFit$value,
               scaConvergence=lsFit$convergence))
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
#' @param scaDisp: (scalar) Inverse of negative binomial 
#'    dispersion coefficients computed by DESeq2 for given gene.
#' @param vecDropoutRate: (probability vector number of samples) 
#'    [Default NULL] Dropout rate/mixing probability of zero inflated 
#'    negative binomial mixturemodel.
#' @param vecProbNB: (probability vector number of samples) [Default NULL]
#'    Probability of observations to come from negative binomial 
#'    component of mixture model.
#' @param vecSizeFactors: (numeric vector number of samples) 
#'    Model scaling factors for each observation: Take
#'    sequencing depth and longitudinal time series mean
#'    within a gene into account (size and translation
#'    factors).
#' @param vecTimepoints: (numeric vector number samples) 
#'    Timepoints assigned to samples.
#' @param vecBatch (numeric vector number samples) 
#'    Longitudinal series assigned to samples. 
#'    NULL if not operating in strMode="longitudinal".
#' @param dfAnnotationProc: (Table) Processed annotation table. 
#'    Lists co-variables of samples: 
#'    Sample, Condition, Time (numeric), TimeCateg (str)
#'    (and Batch). For internal use.
#' @param strMode: (str) [Default "singelbatch"] 
#'    {"singelbatch","batcheffects"}
#'    Batch model.
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

fitConstImpulseGene <- function(vecCounts, 
                                scaDisp,
                                vecSizeFactors,
                                vecTimepoints, 
                                vecBatches=vecBatches,
                                dfAnnotationProc,
                                strMode,
                                boolFitConst,
                                MAXIT=1000){
  
  # (I) Process data
  # Check whether all counts are zero: this
  # Compute time point specifc parameters
  vecTimepointsUnique <- sort(unique( vecTimepoints ))
  vecindTimepointUniqueAssign <- match(vecTimepoints, vecTimepointsUnique)
  # Get batch assignments
  vecBatchesUnique <- unique(vecBatches)
  vecindBatchAssign <- match(vecBatches, vecBatchesUnique)
  
  if(boolFitConst){
    # (II) Fit null model and compute loglikelihood of null model
    # Only necessary if running case only DE analysis as
    # otherwise two impulse models are compared.
    lsConstFit <- fitConstModel(
      vecCounts=vecCounts, 
      scaDisp=scaDisp,
      vecSizeFactors=vecSizeFactors,
      vecBatchesUnique=vecBatchesUnique,
      vecindBatchAssign=vecindBatchAssign,
      strMode=strMode)
  } else {
    lsConstFit <- NULL
  }
  
  # (III) Fit Impulse model
  # 1. Compute initialisations
  lsParamGuesses <- estimateImpulseParam(
    vecCounts=vecCounts,
    vecTimepoints=vecTimepointsUnique, 
    vecindBatchAssign=vecindBatchAssign,
    vecSizeFactors=vecSizeFactors,
    strMode=strMode )
  vecParamGuessPeak <- lsParamGuesses$peak
  vecParamGuessValley <- lsParamGuesses$valley
  
  # 2. Initialisation: Peak
  lsFitPeak <- fitImpulseModel(vecImpulseParamGuess=vecParamGuessPeak,
                               vecCounts=vecCounts,
                               scaDisp=scaDisp,
                               vecSizeFactors=vecSizeFactors,
                               vecTimepointsUnique=vecTimepointsUnique, 
                               vecindTimepointAssign=vecindTimepointUniqueAssign,
                               vecBatchesUnique=vecBatchesUnique,
                               vecindBatchAssign=vecindBatchAssign,
                               strMode=strMode, 
                               MAXIT=MAXIT)
  # 3. Initialisation: Valley
  lsFitValley <- fitImpulseModel(vecImpulseParamGuess=vecParamGuessValley,
                                 vecCounts=vecCounts,
                                 scaDisp=scaDisp,
                                 vecSizeFactors=vecSizeFactors,
                                 vecTimepointsUnique=vecTimepointsUnique, 
                                 vecindTimepointAssign=vecindTimepointUniqueAssign,
                                 vecBatchesUnique=vecBatchesUnique,
                                 vecindBatchAssign=vecindBatchAssign,
                                 strMode=strMode, 
                                 MAXIT=MAXIT)
  
  # (IV) Select best fit and report fit type
  if(lsFitValley$scaLL > lsFitPeak$scaLL){
    lsBestImpulseFit <- lsFitValley
  } else {
    lsBestImpulseFit <- lsFitPeak
  }
  
  return(list(
    lsImpulseFit=lsBestImpulseFit,
    lsConstFit=lsConstFit
  ))
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
#' @param vecTimepoints: (numeric vector number samples) Timepoints 
#'    assigned to samples.
#' @param vecindBatch: (numeric vector number samples) Time courses 
#'    assigned to samples. NULL if not operating in strMode="longitudinal".
#' @param dfAnnotationProc: (Table) Processed annotation table. 
#'    Lists co-variables of samples: 
#'    Sample, Condition, Time (numeric), TimeCateg (str)
#'    (and Batch). For internal use.
#' @param strCaseName: (str) Name of the case condition in \code{dfAnnotationRedFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#'    \code{dfAnnotationRedFull}.
#' @param strMode: (str) [Default "singelbatch"] 
#'    {"singelbatch","batcheffects"}
#'    Batch model.
#' @param nProc: (scalar) [Default 3] Number of processes for parallelisation.
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

fitConstImpulse <- function(matCountDataProcCondition, 
                            vecDispersions, 
                            vecSizeFactors,
                            vecTimepoints, 
                            vecBatches, 
                            dfAnnotationProc,
                            boolFitConst,
                            strMode){
  
  # Maximum number of iterations for numerical optimisation of
  # likelihood function in MLE fitting of impulse model:
  MAXIT <- 1000
  
  lsFits <- bplapply(rownames(matCountDataProcCondition),function(x){
    fitConstImpulseGene(
      vecCounts=matCountDataProcCondition[x,],
      scaDisp=vecDispersions[x],
      vecSizeFactors=vecSizeFactors,
      vecTimepoints=vecTimepoints,
      vecBatches=vecBatches,
      dfAnnotationProc=dfAnnotationProc,
      strMode=strMode,
      boolFitConst=boolFitConst,
      MAXIT=MAXIT )
  })
  names(lsFits) <- rownames(matCountDataProcCondition)
  
  return(lsFits)
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
#'    (and Batch). For internal use.
#' @param vecSizeFactors: (numeric vector number of samples) 
#'    Model scaling factors for each sample which take
#'    sequencing depth into account (size factors).
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
#' @param strMode: (str) [Default "singelbatch"] 
#'    {"singelbatch","batcheffects"}
#'    Batch model.
#' @param nProc: (scalar) [Default 1] Number of processes for parallelisation.
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

fitModels <- function(matCountDataProc, 
                      dfAnnotationProc, 
                      vecSizeFactors,
                      vecDispersions, 
                      strCaseName, 
                      strControlName=NULL, 
                      strMode="singlebatch",
                      nProc=1 ){
  
  lsFitResults_all = list()
  
  # Condition labels to be used in runs:
  # These labels are used to group samples
  if(!is.null(strControlName)){
    lsLabels <- c("combined","case","control")
  } else {
    lsLabels <- c("case")
  }
  
  # Create lists of samples to be used per run
  if(!is.null(strControlName)){    
    lsSamplesByCond <- list(
      combined=colnames(matCountDataProc),
      case=(colnames(matCountDataProc))[dfAnnotationProc[match(
        colnames(matCountDataProc),
        dfAnnotationProc$Sample),
        ]$Condition %in% strCaseName],
      control=(colnames(matCountDataProc))[dfAnnotationProc[match(
        colnames(matCountDataProc),
        dfAnnotationProc$Sample),
        ]$Condition %in% strControlName] )
  } else {
    lsSamplesByCond <- list(
      case=colnames(matCountDataProc) )
  }
  # Get batch assignments of samples
  if(strMode=="batcheffects"){
    vecBatches <- dfAnnotationProc[match(
      colnames(matCountDataProc),
      dfAnnotationProc$Sample),]$Batch
    names(vecBatches) <- colnames(matCountDataProc)
  } else {
    vecBatches <- NULL
  }
  # Get time point assignments of samples
  vecTimepoints <- dfAnnotationProc[match(
    colnames(matCountDataProc),
    dfAnnotationProc$Sample),]$Time
  names(vecTimepoints) <- colnames(matCountDataProc)
  
  # Fitting for different runs
  # Fit constant model only if doing case-only analysis
  # for which the constant model is the null model or 
  # in the combined case of case-ctrl to derive a inferred
  # mean parameter for reference.
  for (label in lsLabels){
    lsFitResults <- fitConstImpulse(
      matCountDataProcCondition=matCountDataProc[,lsSamplesByCond[[label]]],
      vecSizeFactors=vecSizeFactors[lsSamplesByCond[[label]]],
      vecDispersions=vecDispersions,
      vecTimepoints=vecTimepoints[lsSamplesByCond[[label]]],
      vecBatches=vecBatches[lsSamplesByCond[[label]]],
      dfAnnotationProc=dfAnnotationProc,
      boolFitConst= is.null(strControlName) | label=="combined",
      strMode=strMode )
    
    lsFitResults_all[[label]] <- lsFitResults
  }
  
  return(lsFitResults_all)
}