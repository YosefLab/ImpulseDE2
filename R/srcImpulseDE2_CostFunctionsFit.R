#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++     Cost Functions    +++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cost function for negative binomial mean parameter fit
#' 
#' Log likelihood cost function for numerical optimisation of mean model fit 
#' based on negative binomial model. Note that the closed form solution
#' of the negative binomial mean parameter only holds if all normalisation
#' factors are 1.
#' In analogy to generalised linear models, a log linker
#' function is used for the mean. The inferred negative
#' binomial mean model is scaled by the factors in vecSizeFactors,
#' which represent size factors (and translation factors),
#' for evaluation of the likelihood on the data. 
#' 
#' @aliases evalLogLikMuNB_comp
#' 
#' @seealso
#' 
#' @param scaTheta (vector number of parameters [6]) 
#'    Negative binomial mean parameter estimate (log).
#' @param vecCounts (count vector samples) 
#'    Observed expression values for given gene.
#' @param scaDispEst: (scalar) Dispersion estimate for given gene.
#' @param vecSizeFactors: (numeric vector number of samples) 
#'    Normalisation constants for each sample.
#' @param vecWeights: (probability vector number of samples) 
#'    Weights for inference on mixture models.
#' @param vecboolObserved: (bool vector number of samples)
#'    Stores bool of sample being not NA (observed).
#'     
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikMu <- function(vecTheta,
  vecCounts,
  scaDisp, 
  vecSizeFactors,
  lsvecindBatch,
  vecboolObserved){
  
  scaMu <- exp(vecTheta[1])
  scaNParamUsed <- 1
  # Prevent mean shrinkage to zero:
  if(scaMu < 10^(-10)){ scaMu <- 10^(-10) }
  
  vecBatchFactors <- array(1, length(vecCounts))
  if(!is.null(lsvecindBatch)){
  	for(vecindConfounder in lsvecindBatch){
  		scaNBatchFactors <- max(vecindConfounder)-1 # Batches are counted from 1
  		# Factor of first batch is one (constant), the remaining
  		# factors scale based on the first batch.
  		vecBatchFactorsConfounder <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))[vecindConfounder]
  		scaNParamUsed <- scaNParamUsed+scaNBatchFactors
  		# Prevent batch factor shrinkage and explosion:
  		vecBatchFactorsConfounder[vecBatchFactorsConfounder < 10^(-10)] <- 10^(-10)
  		vecBatchFactorsConfounder[vecBatchFactorsConfounder > 10^(10)] <- 10^(10)
  		
  		vecBatchFactors <- vecBatchFactors*vecBatchFactorsConfounder
  	}
  }
  
  # Compute log likelihood under constant model by
  # adding log likelihood of model at each timepoint.
  scaLogLik <- sum(dnbinom(
  	vecCounts[vecboolObserved], 
  	mu=scaMu*vecBatchFactors[vecboolObserved]*
  		vecSizeFactors[vecboolObserved], 
  	size=scaDisp, 
  	log=TRUE))
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Cost function impulse model fit
#' 
#' Log likelihood cost function for impulse model fit based on negative 
#' binomial model.
#' In analogy to generalised linear models, a log linker
#' function is used for the count parameters. The inferred negative
#' binomial model is scaled by the factors in vecSizeFactors,
#' which represent size factors (and translation factors),
#' for evaluation of the likelihood on the data.
#' 
#' @aliases evalLogLikImpulse_comp
#' 
#' @seealso Called by \code{fitImpulse}:\code{fitImpulse_gene}.
#' Calls \code{evalImpulse}. \code{evalLogLikImpulseByTC} for
#' dependent residuals (i.e. time course experiments)
#' 
#' @param vecTheta (vector number of parameters [6]) 
#'    Impulse model parameters.
#' @param vecX (numeric vector number of timepoints) 
#'    Time-points at which gene was sampled.
#' @param vecY (count vector samples) 
#'    Observed expression values for given gene.
#' @param scaDispEst: (scalar) Dispersion estimate for given gene.
#' @param vecSizeFactors: (numeric vector number of samples) 
#'    Normalisation constants for each sample.
#' @param vecindTimepoint (numeric vector number samples) 
#'    Index of time point assigned to sample in list of sorted
#'    time points (vecX).
#' @param vecboolObserved: (bool vector number of samples)
#'    Stores bool of sample being not NA (observed).
#'     
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikImpulse <- function(vecTheta,
  vecCounts,
  scaDisp, 
  vecSizeFactors,
  vecTimepointsUnique,
  vecindTimepoint,
  lsvecindBatch,
  vecboolObserved){
  
  # Compute normalised impulse function value:
  vecImpulseParam <- vecTheta[1:6]
  vecImpulseParam[2:4] <- exp(vecImpulseParam[2:4])
  vecImpulseValue <- evalImpulse_comp(vecImpulseParam=vecImpulseParam,
                                      vecTimepoints=vecTimepointsUnique)[vecindTimepoint]
  scaNParamUsed <- 6
  
  vecBatchFactors <- array(1, length(vecCounts))
  if(!is.null(lsvecindBatch)){
  	for(vecindConfounder in lsvecindBatch){
  		scaNBatchFactors <- max(vecindConfounder)-1 # Batches are counted from 1
  		# Factor of first batch is one (constant), the remaining
  		# factors scale based on the first batch.
  		vecBatchFactorsConfounder <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))[vecindConfounder]
  		scaNParamUsed <- scaNParamUsed+scaNBatchFactors
  		# Prevent batch factor shrinkage and explosion:
  		vecBatchFactorsConfounder[vecBatchFactorsConfounder < 10^(-10)] <- 10^(-10)
  		vecBatchFactorsConfounder[vecBatchFactorsConfounder > 10^(10)] <- 10^(10)
  		
  		vecBatchFactors <- vecBatchFactors*vecBatchFactorsConfounder
  	}
  }
  
  # Compute log likelihood under impulse model by
  # adding log likelihood of model at each timepoint.
  scaLogLik <- sum(dnbinom(
  	vecCounts[vecboolObserved], 
  	mu=vecImpulseValue[vecboolObserved]* 
  		vecBatchFactors[vecboolObserved]*
  		vecSizeFactors[vecboolObserved], 
  	size=scaDisp, 
  	log=TRUE))
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Cost function impulse model fit - Batch mode
#' 
#' Log likelihood cost function for impulse model fit based on negative 
#' binomial model. This cost function is called in the modes
#' "batch" and "longitudinal".
#' In analogy to generalised linear models, a log linker
#' function is used for the count parameters. The inferred negative
#' binomial model is scaled by the factors in vecSizeFactors,
#' which represent size factors (and translation factors),
#' for evaluation of the likelihood on the data.
#' 
#' @aliases evalLogLikImpulseBatch_comp
#' 
#' @seealso Called by \code{fitImpulse}:\code{fitImpulse_gene}.
#' Calls \code{evalImpulse}. \code{evalLogLikImpulseByTC} for
#' dependent residuals (i.e. time course experiments)
#' 
#' @param vecTheta (vector number of parameters [6]) 
#'    Impulse model parameters.
#' @param vecX (numeric vector number of timepoints) 
#'    Time-points at which gene was sampled.
#' @param vecY (count vector samples) 
#'    Observed expression values for given gene.
#' @param scaDispEst: (scalar) Dispersion estimate for given gene.
#' @param vecSizeFactors: (numeric vector number of samples) 
#'    Normalisation constants for each sample.
#' @param vecindTimepoint (numeric vector number samples) 
#'    Index of time point assigned to sample in list of sorted
#'    time points (vecX).
#' @param vecboolObserved: (bool vector number of samples)
#'    Stores bool of sample being not NA (observed).
#'     
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikSigmoid <- function(vecTheta,
                              vecCounts,
                              scaDisp, 
                              vecSizeFactors,
                              vecTimepointsUnique,
                              vecindTimepoint,
                              lsvecindBatch,
                              vecboolObserved){
  
  # Compute normalised impulse function value:
  vecSigmoidParam <- vecTheta[1:4]
  vecSigmoidParam[2:3] <- exp(vecSigmoidParam[2:3])
  vecSigmoidValue <- evalImpulse_comp(vecSigmoidParam=vecSigmoidParam,
                                      vecTimepoints=vecTimepointsUnique)[vecindTimepoint]
  scaNParamUsed <- 4
  
  vecBatchFactors <- array(1, length(vecCounts))
  if(!is.null(lsvecindBatch)){
    for(vecindConfounder in lsvecindBatch){
      scaNBatchFactors <- max(vecindConfounder)-1 # Batches are counted from 1
      # Factor of first batch is one (constant), the remaining
      # factors scale based on the first batch.
      vecBatchFactorsConfounder <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))[vecindConfounder]
      scaNParamUsed <- scaNParamUsed+scaNBatchFactors
      # Prevent batch factor shrinkage and explosion:
      vecBatchFactorsConfounder[vecBatchFactorsConfounder < 10^(-10)] <- 10^(-10)
      vecBatchFactorsConfounder[vecBatchFactorsConfounder > 10^(10)] <- 10^(10)
      
      vecBatchFactors <- vecBatchFactors*vecBatchFactorsConfounder
    }
  }
  
  # Compute log likelihood under impulse model by
  # adding log likelihood of model at each timepoint.
  scaLogLik <- sum(dnbinom(
    vecCounts[vecboolObserved], 
    mu=vecSigmoidValue[vecboolObserved]* 
      vecBatchFactors[vecboolObserved]*
      vecSizeFactors[vecboolObserved], 
    size=scaDisp, 
    log=TRUE))
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}