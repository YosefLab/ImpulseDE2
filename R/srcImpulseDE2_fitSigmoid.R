#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Sigmoidal  model fit    +++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Estimate sigmoidal model parameter initialisations
#' 
#' The initialisations reflect intuitive parameter choices corresponding
#' to a peak and to a valley model.
#' 
#' @seealso Called by \code{fit_gene}.
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

estimateSigmoidParam <- function(vecCounts,
                                 vecTimepoints,
                                 vecSizeFactors,
                                 lsvecindBatch){
  # Compute general statistics for initialisation:
  # Expression means by timepoint
  vecCountsSFcorrected <- vecCounts/vecSizeFactors
  if(!is.null(lsvecindBatch)){
    # Estimate batch factors
    vecBatchFactors <- array(1, length(vecCounts))
    for(vecindBatch in lsvecindBatch){
      vecBatchFactorsConfounder <- sapply(unique(vecindBatch), function(batch){
        mean(vecCountsSFcorrected[vecindBatch==batch]/mean(vecCounts, na.rm=TRUE), na.rm=TRUE)
      })
      # Catch exception that all observations of a batch are zero or all observations are zero:
      vecBatchFactorsConfounder[is.na(vecBatchFactorsConfounder) | vecBatchFactorsConfounder==0] <- 1
      vecBatchFactors <- vecBatchFactors*vecBatchFactorsConfounder[vecindBatch]
    }
    vecCountsSFBatchcorrected <- vecCountsSFcorrected/vecBatchFactors
    vecExpressionMeans <- sapply(vecTimepoints, function(tp){
      mean(vecCountsSFBatchcorrected[vecTimepoints==tp], na.rm=TRUE)
    })
  } else {
    vecExpressionMeans <- sapply(vecTimepoints, function(tp){
      mean(vecCountsSFcorrected[vecTimepoints==tp], na.rm=TRUE)
    })
  }
  scaNTimepoints <- length(vecTimepoints)
  idxMiddleTP <- round(scaNTimepoints/2)
  scaMaxEarlyMean <- max(vecExpressionMeans[1:(idxMiddleTP-1)], na.rm=TRUE)
  scaMinEarlyMean <- min(vecExpressionMeans[1:(idxMiddleTP-1)], na.rm=TRUE)
  scaMaxLateMean <- max(vecExpressionMeans[idxMiddleTP:scaNTimepoints], na.rm=TRUE)
  scaMinLateMean <- min(vecExpressionMeans[idxMiddleTP:scaNTimepoints], na.rm=TRUE)
  
  # Compute peak initialisation
  vecParamGuessPeak <- c(1,
                         log(scaMinEarlyMean+1),
                         log(scaMaxLateMean+1),
                         vecTimepoints[idxMiddleTP])
  
  # Compute valley initialisation
  vecParamGuessValley <- c(1,
                           log(scaMaxEarlyMean+1),
                           log(scaMinLateMean+1),
                           vecTimepoints[idxMiddleTP])
  
  return(list(peak=vecParamGuessPeak, 
              valley=vecParamGuessValley))
}

#' Fit a sigmoidal model to data of a gene
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
#' @param vecindTimepoint (numeric vector number samples) 
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

fitSigmoidModel <- function(vecSigmoidParamGuess,
                            vecCounts,
                            scaDisp,
                            vecSizeFactors,
                            lsvecindBatch,
                            vecTimepointsUnique,
                            vecindTimepoint,
                            MAXIT=100,
                            RELTOL=10^(-8),
                            trace=0,
                            REPORT=10 ){
  
  
  vecParamGuess <- vecSigmoidParamGuess
  if(!is.null(lsvecindBatch)){
    for(vecindConfounder in lsvecindBatch){
      vecParamGuess <- c(vecParamGuess, rep(0, length(unique(vecindConfounder))-1))
    }
  }
  
  lsFit <- tryCatch({
    optim(
      par=vecParamGuess,
      fn=evalLogLikSigmoid_comp,
      vecCounts=vecCounts,
      scaDisp=scaDisp,
      vecSizeFactors=vecSizeFactors,
      vecTimepointsUnique=vecTimepointsUnique,
      vecindTimepoint=vecindTimepoint,
      lsvecindBatch=lsvecindBatch,
      vecboolObserved=!is.na(vecCounts),
      method="BFGS",
      control=list(maxit=MAXIT,
                   reltol=RELTOL,
                   fnscale=-1)
    )[c("par","value","convergence")]
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting sigmoid model: fitSigmoidModel().",
                 " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
    print(paste0("vecParamGuess ", paste(vecParamGuess,collapse=" ")))
    print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
    print(paste0("scaDisp ", paste(scaDisp,collapse=" ")))
    print(paste0("vecSizeFactors ", paste(vecSizeFactors,collapse=" ")))
    print(paste0("vecTimepointsUnique ", paste(vecTimepointsUnique,collapse=" ")))
    print(paste0("vecindTimepoint ", paste(vecindTimepoint,collapse=" ")))
    print(paste0("lsvecindBatch ", paste(lsvecindBatch,collapse=" ")))
    print(paste0("MAXIT ", MAXIT))
    print(strErrorMsg)
    stop(strErrorMsg)
  })
  
  # Extract parameter estimates
  vecSigmoidParam <- lsFit$par[1:6]
  vecSigmoidParam[2:4] <- exp(vecSigmoidParam[2:4])
  names(vecSigmoidParam) <- c("beta", "h0", "h1", "t")
  vecSigmoidValue <- evalSigmoid_comp(vecSigmoidParam=vecSigmoidParam,
                                      vecTimepoints=vecTimepointsUnique)[vecindTimepoint]
  names(vecSigmoidValue) <- names(vecCounts)
  scaNParamUsed <- 6
  if(!is.null(lsvecindBatch)){
    lsvecBatchFactors <- lapply(lsvecindBatch, function(vecindConfounder){
      scaNBatchFactors <- max(vecindConfounder)-1 # Batches are counted from 1
      # Factor of first batch is one (constant), the remaining
      # factors scale based on the first batch.
      vecBatchFactorsConfounder <- c(1, exp(lsFit$par[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
      scaNParamUsed <- scaNParamUsed+scaNBatchFactors
      # Catch boundary of likelihood domain on batch factor space:
      vecBatchFactorsConfounder[vecBatchFactorsConfounder < 10^(-10)] <- 10^(-10)
      vecBatchFactorsConfounder[vecBatchFactorsConfounder > 10^(10)] <- 10^(10)
      return(vecBatchFactorsConfounder)
    })
  }
  
  return(list( vecSigmoidParam=vecSigmoidParam,
               vecSigmoidValue=vecSigmoidValue,
               lsvecBatchFactors=lsvecBatchFactors,
               scaDispParam=scaDisp,
               scaLL=lsFit$value,
               scaConvergence=lsFit$convergence))
}

#' Fit an sigmoidal model to a single gene
#' 
#' Fits an sigmoidal model and a mean model to a single gene.
#' The optimisation method is set within \code{optimiseImpulseModelFit} 
#' (\code{optim}: BFGS). This function is divided into four parts: 
#' (I) Prepare data,
#' (II) Fit sigmoidal models, 
#' (III) Process Fits.
#' The body of this function is broken up into 4 helper functions, 
#' which are exclusively called by this function.
#' 
#' @seealso Called by \code{fitImpulse_matrix}. This function
#' calls \code{fitNullModel}, \code{estimateSigmoidParam} and
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

fitSigmoidGene <- function(vecCounts, 
                           scaDisp,
                           vecSizeFactors,
                           vecTimepoints, 
                           lsvecBatches=lsvecBatches,
                           MAXIT=1000){
  
  # (I) Process data
  # Compute time point specifc parameters
  vecTimepointsUnique <- sort(unique( vecTimepoints ))
  vecindTimepoint <- match(vecTimepoints, vecTimepointsUnique)
  # Get batch assignments
  if(!is.null(lsvecBatches)){
    lsvecindBatch <- list()
    for(batchfactor in lsvecBatches){
      lsvecindBatch[[length(lsvecindBatch)+1]] <- match(batchfactor, unique(batchfactor))
    }
  } else { lsvecindBatch <- NULL }
  
  # (II) Fit sigmoidal model
  # 1. Compute initialisations
  lsParamGuesses <- estimateSigmoidParam(
    vecCounts=vecCounts,
    vecTimepoints=vecTimepointsUnique, 
    lsvecindBatch=lsvecindBatch,
    vecSizeFactors=vecSizeFactors )
  vecParamGuessPeak <- lsParamGuesses$peak
  vecParamGuessValley <- lsParamGuesses$valley
  
  # 2. Initialisation: Peak
  lsFitPeak <- fitSigmoidModel(vecSigmoidParamGuess=vecParamGuessPeak,
                               vecCounts=vecCounts,
                               scaDisp=scaDisp,
                               vecSizeFactors=vecSizeFactors,
                               vecTimepointsUnique=vecTimepointsUnique, 
                               vecindTimepoint=vecindTimepoint,
                               lsvecindBatch=lsvecindBatch,
                               MAXIT=MAXIT)
  # 3. Initialisation: Valley
  lsFitValley <- fitSigmoidModel(vecSigmoidParamGuess=vecParamGuessValley,
                                 vecCounts=vecCounts,
                                 scaDisp=scaDisp,
                                 vecSizeFactors=vecSizeFactors,
                                 vecTimepointsUnique=vecTimepointsUnique, 
                                 vecindTimepoint=vecindTimepoint,
                                 lsvecindBatch=lsvecindBatch,
                                 MAXIT=MAXIT)
  
  # (III) Select best fit and report fit type
  if(lsFitValley$scaLL > lsFitPeak$scaLL){
    lsbestSigmoidFit <- lsFitValley
  } else {
    lsbestSigmoidFit <- lsFitPeak
  }
  
  return(lsbestSigmoidFit=lsbestSigmoidFit)
}

#' Fits sigmoid model to a timecourse dataset
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
#' @export

fitSigmoidModels <- function(matCountDataProc, 
                             dfAnnotationProc, 
                             lsModelFits,
                             vecSizeFactors,
                             vecDispersions,
                             vecConfounders,
                             strCondition){
  
  vecSamplesCond <- dfAnnotationProc[dfAnnotationProc$Condition==strCondition,]$Sample
  
  # Get batch assignments of samples
  if(!is.null(vecConfounders)){
    lsvecBatches <- lapply(vecConfounders, function(confounder){
      vecBatches <- dfAnnotationProc[,confounder]
      names(vecBatches) <- 	colnames(matCountDataProc)
      return(vecBatches)
    })
  } else { lsvecBatches <- NULL }
  
  # Get time point assignments of samples
  vecTimepoints <- dfAnnotationProc$Time
  names(vecTimepoints) <- colnames(matCountDataProc)
  
  # Developmental note: Compared to impulse/constant fitting,
  # this function does not iterate over conditions as this is likely
  # only used for one condition (case). Therefore merge two wrappers
  # used for impulse/const fit into one here.
  
  # Maximum number of iterations for numerical optimisation of
  # likelihood function in MLE fitting of impulse model:
  MAXIT <- 1000
  
  lsSigmoidFits <- bplapply(rownames(matCountDataProc),function(x){
    fitSigmoidGene(
      vecCounts=matCountDataProc[x,vecSamplesCond],
      scaDisp=vecDispersions[x],
      vecSizeFactors=vecSizeFactors[vecSamplesCond],
      vecTimepoints=vecTimepoints[vecSamplesCond],
      lsvecBatches=lapply(lsvecBatches, function(confounder) confounder[vecSamplesCond] ),
      MAXIT=MAXIT )
  })
  names(lsSigmoidFits) <- rownames(matCountDataProc)
  
  # Add sigmoid fits into model fit data structure to preexisting impulse
  # (and constant) fits.
  for(x in rownames(matCountDataProc)){
    lsModelFits[[strCondition]][[x]]$lsSigmoidFit <- lsSigmoidFits[[x]]
  }
  
  return(lsModelFits)
}