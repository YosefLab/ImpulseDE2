#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++     Cost Functions    +++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cost function impulse model fit - Batch mode
#' 
#' Log likelihood cost function for impulse model fit based on negative 
#' binomial model. This cost function is called in the modes
#' "batch" and "longitudinal".
#' In analogy to generalised linear models, a log linker
#' function is used for the count parameters. The inferred negative
#' binomial model is scaled by the factors in vecNormConst,
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
#' @param vecNormConst: (numeric vector number of samples) 
#'    Normalisation constants for each sample.
#' @param vecindTimepointAssign (numeric vector number samples) 
#'    Index of time point assigned to sample in list of sorted
#'    time points (vecX).
#' @param vecboolObserved: (bool vector number of samples)
#'    Stores bool of sample being not NA (observed).
#'     
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikImpulseBatch <- function(vecTheta, 
  vecX, 
  vecY,
  scaDispEst, 
  vecNormConst, 
  vecindTimepointAssign, 
  vecboolObserved){

  # Compute normalised impulse function value: 
  # Mean of negative binomial density at each time point,
  # scaled by normalisation factor of each sample.
  vecImpulseValue <- calcImpulse_comp(vecTheta,vecX)[vecindTimepointAssign]*
    vecNormConst
  
  # Compute log likelihood under impulse model by
  # adding log likelihood of model at each timepoint.
  scaLogLik <- sum(dnbinom(
    vecY[vecboolObserved], 
    mu=vecImpulseValue[vecboolObserved], 
    size=scaDispEst, 
    log=TRUE))
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute log likelihood of zero-inflated negative binomial model
#' 
#' This liklihood function is appropriate for sequencing data with high drop 
#' out rate, commonly observed in single cell data (e.g. scRNA-seq).
#' 
#' @aliases evalLogLikZINB_comp
#' 
#' @seealso Called by \code{fitImpulse}::\code{fitImpulse_matrix}::
#' \code{fitImpulse_gene}::\code{computeLogLikNull} and
#' \code{evalLogLikImpulseSC}.
#'
#' @param vecY (count vector number of amples)
#'    Observed expression values for  given gene.
#' @param vecMu (vector number of samples) Negative binomial
#'    mean parameter for each sample.
#' @param scaDispEst: (scalar) Negative binomial dispersion 
#'    parameter for given gene.
#' @param vecDropoutRateEst: (probability vector number of samples) 
#'    Dropout rate estimate for each cell for given gene.
#' @param vecboolNotZeroObserved: (bool vector number of samples)
#'    Whether sample is not zero and observed (not NA).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#'    
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikZINB <- function(vecY,
  vecMu,
  scaDispEst, 
  vecDropoutRateEst, 
  vecboolNotZeroObserved, 
  vecboolZero){  

  # Note on handling very low probabilities: vecLikZeros
  # typically does not have zero elements as it has the 
  # the summand drop-out rate. Also the log cannot be
  # broken up over the sum to dnbinom. In contrast to that,
  # the log is taken in dnbinom for vecLikNonzeros to avoid 
  # zero probabilities. Zero probabilities are handled
  # through substitution of the minimal probability under
  # machine precision.
  
  # Likelihood of zero counts:
  vecLikZeros <- (1-vecDropoutRateEst[vecboolZero])*
    dnbinom(
      vecY[vecboolZero], 
      mu=vecMu[vecboolZero], 
      size=scaDispEst, 
      log=FALSE) +
    vecDropoutRateEst[vecboolZero]
  # Replace zero likelihood observation with machine precision
  # for taking log.
  scaLogLikZeros <- sum( log(vecLikZeros[vecLikZeros!=0]) +
      sum(vecLikZeros==0)*log(.Machine$double.eps) )
  # Likelihood of non-zero counts:
  vecLikNonzeros <- log(1-vecDropoutRateEst[vecboolNotZeroObserved]) +
    dnbinom(
      vecY[vecboolNotZeroObserved], 
      mu=vecMu[vecboolNotZeroObserved], 
      size=scaDispEst, 
      log=TRUE)
  # Replace zero likelihood observation with machine precision
  # for taking log.
  scaLogLikNonzeros <- sum( vecLikNonzeros[is.finite(vecLikNonzeros)]) +
      sum(!is.finite(vecLikNonzeros))*log(.Machine$double.eps)
  # Compute likelihood of all data:
  scaLogLik <- scaLogLikZeros + scaLogLikNonzeros
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Cost function impulse model fit - Single cell mode
#' 
#' Log likelihood cost function for impulse model fit based on zero inflated  
#' negative binomial model ("hurdle model"). This cost function is appropriate
#' for sequencing data with high drop out rate, commonly observed in single
#' cell data (e.g. scRNA-seq).
#' In analogy to generalised linear models, a log linker
#' function is used for the count parameters. The impulse model values
#' (negative binomial mean parameters) are normalised by vecNormConst
#' for evaluation of the likelihood on the data.
#' 
#' @aliases evalLogLikImpulseSC_comp
#' 
#' @seealso Called by \code{fitImpulse}::\code{fitImpulse_matrix}::
#' \code{fitImpulse_gene}::\code{optimiseImpulseModelFit}.
#' Calls \code{calcImpulse} and \code{evalLogLikZINB_comp}.
#' 
#' @param vecTheta (vector number of parameters [6]) Impulse model parameters.
#' @param vecX (numeric vector number of timepoints) 
#'    Time-points at which gene was sampled.
#' @param vecY (count vector samples) Observed expression values for 
#'    given gene.
#' @param scaDispEst: (scalar) Negative binomial dispersion 
#'    parameter for given gene.
#' @param vecDropoutRateEst: (probability vector number of samples) 
#'    Dropout rate estimate for each cell for given gene.
#' @param vecNormConst: (numeric vector number of samples) 
#'    Normalisation constants for each sample.
#' @param vecindTimepointAssign (numeric vector number samples) 
#'    Index of time point assigned to sample in list of sorted
#'    time points (vecX).
#' @param vecboolNotZeroObserved: (bool vector number of samples)
#'    Whether sample is not zero and observed (not NA).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#'    
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikImpulseSC <- function(vecTheta,
  vecX,
  vecY,
  scaDispEst, 
  vecDropoutRateEst, 
  vecNormConst,
  vecindTimepointAssign, 
  vecboolNotZeroObserved, 
  vecboolZero){  

  # Generate negative binomial mean parameters from impulse model:
  # This links the parameters vecTheta which are changed in optimisation
  # to the cost function.
  vecImpulseValue <- calcImpulse_comp(vecTheta,vecX)[vecindTimepointAssign]*
    vecNormConst
  
  # Evaluate likelihood (this is the cost function):  
  scaLogLik <- evalLogLikZINB_comp(vecY=vecY,
    vecMu=vecImpulseValue,
    scaDispEst=scaDispEst, 
    vecDropoutRateEst=vecDropoutRateEst, 
    vecboolNotZeroObserved=vecboolNotZeroObserved, 
    vecboolZero=vecboolZero)
  
  return(scaLogLik)
}