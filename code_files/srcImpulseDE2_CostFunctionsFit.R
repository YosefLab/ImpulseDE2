#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++     Cost Functions    +++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cost function impulse model fit - Batch mode
#' 
#' Log likelihood cost function for impulse model fit based on negative 
#' binomial model. Batch mode means that residuals are assumed to be 
#' independent. 
#' 
#' @aliases evalLogLikImpulseBatch_comp
#' 
#' @seealso Called by \code{fitImpulse}:\code{fitImpulse_gene}.
#' Calls \code{calcImpulse}. \code{evalLogLikImpulseByTC} for
#' dependent residuals (i.e. time course experiments)
#' 
#' @param vecTheta (vector number of parameters [6]) Impulse model parameters.
#' @param vecX (vector number of timepoints) Time-points at which gene was sampled
#' @param matY (2D array timepoints x replicates) Observed expression values for 
#     given gene.
#' @param scaDispEst: (scalar) Dispersion estimate for given gene.
#' 
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikImpulseBatch <- function(vecTheta,vecX,matY,scaDispEst){  
  # Compute impulse function value: Mean of negative binomial
  # likelihood function for each time point.
  vecImpulseValue = calcImpulse_comp(vecTheta,vecX)
  
  # Compute log likelihood under impulse model by
  # adding log likelihood of model at each timepoint.
  scaLogLik <- sum( sapply(c(1:length(vecX)), function(tp){
    sum(dnbinom(matY[tp,!is.na(matY[tp,])], mu=vecImpulseValue[tp], size=scaDispEst, log=TRUE))}) )
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Cost function impulse model fit - Time course mode
#' 
#' Log likelihood cost function for impulse model fit based on negative 
#' binomial model. Time course means that replicates of samples can be
#' grouped into time course experiments, i.e. residuals are dependent
#' within a time course. The impulse model is scaled to the
#' mean of each time course.
#' 
#' @aliases evalLogLikImpulseByTC_comp
#' 
#' @seealso Called by \code{fitImpulse}:\code{fitImpulse_gene}.
#' Calls \code{calcImpulse}. \code{evalLogLikImpulseBatch} for
#' independent residuals.
#' 
#' @param vecTheta (vector number of parameters [6]) Impulse model parameters.
#' @param vecX (vector number of timepoints) Time-points at which gene was sampled
#' @param matY (2D array timepoints x replicates) Observed expression values for 
#     given gene.
#' @param scaDispEst: (scalar) Dispersion estimate for given gene.
#' @param scaMu: (scalar) MLE of overall mean of negative binomial
#'    constant model on all data.
#' @param vecMuTimecourses: (vector time courses) MLEs of meana of negative binomial
#'    constant models by time course.
#'    
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikImpulseByTC <- function(vecTheta,vecX,matY,scaDispEst,
  scaMu,vecMuTimecourses){  
  # Compute impulse function value: Mean of negative binomial
  # likelihood function for each time point.
  vecImpulseValue = calcImpulse_comp(vecTheta,vecX)
  
  # Compute log likelihood under impulse model by
  # adding log likelihood of model at each timepoint.
  vecTranslationFactors <- vecMuTimecourses/scaMu
  scaLogLik <- sum(sapply( c(1:length(vecX)), function(tp){
    sum(dnbinom(matY[tp,!is.na(matY[tp,])], 
      mu=vecImpulseValue[tp] * vecTranslationFactors[!is.na(matY[tp,])], 
      size=scaDispEst, log=TRUE))
    }) )
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cost function impulse model fit - Single cell mode
#' 
#' Log likelihood cost function for impulse model fit based on zero inflated  
#' negative binomial model ("hurdle model"). This cost function is appropriate
#' for sequencing data with high drop out rate, commonly observed in single
#' cell data (e.g. scRNA-seq).
#' 
#' @aliases evalLogLikImpulseSC_comp
#' 
#' @seealso Called by \code{fitImpulse}:\code{fitImpulse_gene}.
#' Calls \code{calcImpulse}.
#' 
#' @param vecTheta (vector number of parameters [6]) Impulse model parameters.
#' @param vecX (vector number of timepoints) Time-points at which gene was sampled
#' @param matY (2D array timepoints x replicates) Observed expression values for 
#     given gene.
#' @param scaDispEst: (scalar) Dispersion estimate for given gene.
#' @param scaDropoutEst: (scalar) Dropout rate estimate for given cell.
#' 
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

### DAVID developmental note: think about taking this from 2D into 1D, only need 2D if work in clusters
evalLogLikImpulseSC <- function(vecTheta,vecX,matY,
  scaDispEst,scaDropoutEst){  
  # Compute impulse function value: Mean of negative binomial
  # likelihood function for each time point.
  vecImpulseValue = calcImpulse_comp(vecTheta,vecX)
  
  # Compute log likelihood under impulse model by
  # adding log likelihood of model at each timepoint.
  # Likelihood function of hurdle modle differs between
  # zero and non-zero counts: Add likelihood of both together.
  # Note that the log is taken over the sum of mixture model 
  # components of the likelihood model and accordingly log=FALSE
  # inside dnbinom.
  # Likelihood of zero counts:
  scaLogLikZeros <- sum( sapply(c(1:length(vecX)), function(tp){ sum(log(
    (1-scaDropoutEst) * dnbinom(matY[tp,matY[tp,]==0], mu=vecImpulseValue[tp], size=scaDispEst, log=FALSE) +
      scaDropoutEst
  ))}) )
  # Likelihood of non-zero counts:
  scaLogLikNonzeros <- sum( sapply(c(1:length(vecX)), function(tp){ sum(log(
    (1-scaDropoutEst) * dnbinom(matY[tp,matY[tp,]!=0 & !is.na(matY[tp,])], mu=vecImpulseValue[tp], size=scaDispEst, log=FALSE)
  ))}) )
  # Compute likelihood of all data:
  scaLogLik <- scaLogLikZeros + scaLogLikNonzeros
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}
