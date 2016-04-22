#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++     Cost Functions    +++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Log likelihood cost function for impulse model fit based on negative 
#' binomial model.
#' 
#' @aliases evalLogLikImpulse_comp
#' 
#' @seealso Called by \code{fitImpulse}:\code{fitImpulse_gene}.
#' Calls \code{calcImpulse}.
#' 
#' @param vecTheta (vector number of parameters [6]) Impulse model parameters.
#' @param vecX (vector number of timepoints) Time-points at which gene was sampled
#' @param matY (2D array timepoints x replicates) Observed expression values for 
#     given gene.
#' @param scaDispEst: (scalar) Dispersion estimate for given gene.
#' 
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikImpulse <- function(vecTheta,vecX,matY,scaDispEst){  
  # Compute function value - mean
  vecImpulseValue = calcImpulse_comp(vecTheta,vecX)
  
  # Compute log likelihood under impulse model by
  # adding log likelihood of model at each timepoint.
  # Only consider timepoints with not all entries NA.
  scaLogLik <- sum( sapply(c(1:length(vecX)), function(tp){
    sum(dnbinom(matY[tp,!is.na(matY[tp,])], mu=vecImpulseValue[tp], size=scaDispEst, log=TRUE))}) )
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Log likelihood cost function for mean model fit based on negative binomial
#' model.
#' 
#' @aliases evalLogLikMean_comp
#' 
#' @seealso Called by \code{fitImpulse}:\code{fitImpulse_gene}.
#' Calls \code{calcImpulse}.
#' 
#' @param scaMuEst (scalar) Mean for single negative binomial model.
#' @param matY (2D array timepoints x replicates) Observed expression values for 
#     given gene.
#' @param scaDispEst: (scalar) Dispersion estimate for given gene.
#' 
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikMean <- function(scaMuEst,matY,scaDispEst){
  # Compute log likelihood assuming constant mean of negative binomial
  scaLogLik <- sum (dnbinom(matY[!is.na(matY)], mu=exp(scaMuEst), size=scaDispEst, log=TRUE))
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}
