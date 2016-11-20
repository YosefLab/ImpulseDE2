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
  scaDispEst, 
  vecSizeFactors,
  vecindBatch,
  strMode,
  vecboolObserved){
  
  scaMu <- exp(vecTheta[1])
  # Prevent means shrinkage to zero:
  if(scaMu < 10^(-10)){ scaMu <- 10^(-10) }
  
  if(strMode=="singlebatch"){
    # Compute log likelihood under constant model by
    # adding log likelihood of model at each timepoint.
    scaLogLik <- sum(dnbinom(
      vecCounts[vecboolObserved], 
      mu=scaMu * vecSizeFactors[vecboolObserved], 
      size=scaDispEst, 
      log=TRUE))
  } else {
    # Factor of first batch is one (constant), the remaining
    # factors scale based on the first batch.
    vecBatchFactors <- c(1, exp(vecTheta[2:length(vecTheta)]))
    # Prevent batch factor shrinkage and explosion:
    vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
    vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
    
    # Compute log likelihood under constant model by
    # adding log likelihood of model at each timepoint.
    scaLogLik <- sum(dnbinom(
      vecCounts[vecboolObserved], 
      mu=scaMu*vecBatchFactors[vecindBatch[vecboolObserved]]*
        vecSizeFactors[vecboolObserved], 
      size=scaDispEst, 
      log=TRUE))
  }
  
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
#' Calls \code{calcImpulse}. \code{evalLogLikImpulseByTC} for
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
#' @param vecindTimepointAssign (numeric vector number samples) 
#'    Index of time point assigned to sample in list of sorted
#'    time points (vecX).
#' @param vecboolObserved: (bool vector number of samples)
#'    Stores bool of sample being not NA (observed).
#'     
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikImpulse <- function(vecTheta,
  vecCounts,
  scaDispEst, 
  vecSizeFactors,
  vecTimepointsUnique,
  vecindTimepointAssign,
  vecindBatch,
  strMode,
  vecboolObserved){
  
  # Compute normalised impulse function value:
  vecImpulseParam <- vecTheta[1:6]
  vecImpulseParam[2:4] <- exp(vecImpulseParam[2:4])
  vecImpulseValue <- calcImpulse_comp(vecImpulseParam=vecImpulseParam,
                                      vecTimepoints=vecTimepointsUnique)[vecindTimepointAssign]
  
  if(strMode=="singlebatch"){
    # Compute log likelihood under impulse model by
    # adding log likelihood of model at each timepoint.
    scaLogLik <- sum(dnbinom(
      vecCounts[vecboolObserved], 
      mu=vecImpulseValue[vecboolObserved]*
        vecSizeFactors[vecboolObserved], 
      size=scaDispEst, 
      log=TRUE))
  } else {
    # Factor of first batch is one (constant), the remaining
    # factors scale based on the first batch.
    vecBatchFactors <- c(1, exp(vecTheta[7:length(vecTheta)]))
    # Prevent batch factor shrinkage and explosion:
    vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
    vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
    
    # Compute log likelihood under impulse model by
    # adding log likelihood of model at each timepoint.
    scaLogLik <- sum(dnbinom(
      vecCounts[vecboolObserved], 
      mu=vecImpulseValue[vecboolObserved]* 
        vecBatchFactors[vecindBatch[vecboolObserved]]*
        vecSizeFactors[vecboolObserved], 
      size=scaDispEst, 
      log=TRUE))
  }
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}