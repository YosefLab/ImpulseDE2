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
#' binomial mean model is scaled by the factors in vecNormConst,
#' which represent size factors (and translation factors),
#' for evaluation of the likelihood on the data. 
#' 
#' @aliases evalLogLikNBMean_comp
#' 
#' @seealso
#' 
#' @param scaTheta (vector number of parameters [6]) 
#'    Negative binomial mean parameter.
#' @param vecCounts (count vector samples) 
#'    Observed expression values for given gene.
#' @param scaDispEst: (scalar) Dispersion estimate for given gene.
#' @param vecNormConst: (numeric vector number of samples) 
#'    Normalisation constants for each sample.
#' @param vecWeights: (probability vector number of samples) 
#'    Weights for inference on mixture models.
#' @param vecboolObserved: (bool vector number of samples)
#'    Stores bool of sample being not NA (observed).
#'     
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikNBMean <- function(scaTheta,
  vecCounts,
  scaDispEst, 
  vecNormConst,
  vecboolObserved){
  
  scaNBMean <- exp(scaTheta)
  # Compute log likelihood under impulse model by
  # adding log likelihood of model at each timepoint.
  scaLogLik <- sum(dnbinom(
    vecCounts[vecboolObserved], 
    mu=scaNBMean*vecNormConst[vecboolObserved], 
    size=scaDispEst, 
    log=TRUE))
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Fit negative binomial mean parameter
#' 
#' Numerical optimisation of negative binomial mean model fit.
#' Note that the closed form solution of the maximum likelihood 
#' estimator of the negative binomial mean parameter 
#' (the weighted average) only holds if all normalisation
#' factors are 1. Catches numerical errors.
#' 
#' @seealso Called by \code{computeTranslationFactors()}
#' and \code{computeLogLikNull()}.
#' 
#' @param vecCounts (count vector samples) 
#'    Observed expression values for given gene.
#' @param scaDispEst: (scalar) Dispersion estimate for given gene.
#' @param vecNormConst: (numeric vector number of samples) 
#'    Normalisation constants for each sample.
#' @param vecWeights: (probability vector number of samples) 
#'    Weights for inference on mixture models.
#' @param vecboolObserved: (bool vector number of samples)
#'    Stores bool of sample being not NA (observed).
#'     
#' @return scaMu: (scalar) Maximum likelihood estimator of
#'    negative binomial mean parameter.
#' @export

fitNBMean <- function(vecCounts,
  scaDispEst,
  vecNormConst){
  
  scaMu <- tryCatch({
    exp(unlist(optimise(
      evalLogLikNBMean_comp,
      vecCounts=vecCounts,
      scaDispEst=scaDispEst, 
      vecNormConst=vecNormConst,
      vecboolObserved=!is.na(vecCounts),
      lower = log(.Machine$double.eps),
      upper = log(max(vecNormConst*vecCounts, na.rm=TRUE)+1),
      maximum = TRUE)["maximum"]))
  }, error=function(strErrorMsg){
    print(paste0("ERROR: Fitting negative binomial mean parameter: fitNBMean().",
      " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
    print(paste0("vecCounts ", paste(vecCounts,sep=" ")))
    print(paste0("scaDispEst ", paste(scaDispEst,sep=" ")))
    print(paste0("vecNormConst ", paste(vecNormConst,sep=" ")))
    lsErrorCausingGene <- list(vecCounts, scaDispEst, vecNormConst)
    names(lsErrorCausingGene) <- c("vecCounts", "scaDispEst", "vecNormConst")
    save(lsErrorCausingGene,file=file.path(getwd(),"ImpulseDE2_lsErrorCausingGene.RData"))
    stop(strErrorMsg)
  })
  
  return(scaMu)
}

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
  
  # Catch cases in which entire sample is zero observation
  # and impulse value is pushed to close to zero to be evaluated
  # as a negative binomial mean without numerical error.
  vecImpulseValue[vecImpulseValue < .Machine$double.eps] <- .Machine$double.eps
  
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
  #vecLikZeros <- (1-vecDropoutRateEst[vecboolZero])*
  #  dnbinom(
  #    vecY[vecboolZero], 
  #    mu=vecMu[vecboolZero], 
  #    size=scaDispEst, 
  #    log=FALSE) +
  #  vecDropoutRateEst[vecboolZero]
  # Use closed form solution to negative binomial likelihood at zero here.
  vecLikZeros <- (1-vecDropoutRateEst[vecboolZero])*
    (scaDispEst/(scaDispEst + vecMu[vecboolZero]))^scaDispEst +
    vecDropoutRateEst[vecboolZero]
  # Replace zero likelihood observation with machine precision
  # for taking log.
  scaLogLikZeros <- sum( log(vecLikZeros[vecLikZeros > .Machine$double.eps]) +
      sum(vecLikZeros <= .Machine$double.eps)*log(.Machine$double.eps) )
  #scaLogLikZeros <- sum( log(vecLikZeros[vecLikZeros!=0]) +
  #    sum(vecLikZeros==0)*log(.Machine$double.eps) )
  # Likelihood of non-zero counts:
  vecLogLikNonzeros <- log(1-vecDropoutRateEst[vecboolNotZeroObserved]) +
    dnbinom(
      vecY[vecboolNotZeroObserved], 
      mu=vecMu[vecboolNotZeroObserved], 
      size=scaDispEst, 
      log=TRUE)
  # Replace zero likelihood observation with machine precision
  # for taking log.
  #scaLogLikNonzeros <- sum( vecLogLikNonzeros[is.finite(vecLogLikNonzeros)]) +
  #    sum(!is.finite(vecLogLikNonzeros))*log(.Machine$double.eps)
  scaLogLikNonzeros <- sum( vecLogLikNonzeros[vecLogLikNonzeros > log(.Machine$double.eps)] ) +
    sum(vecLogLikNonzeros <= log(.Machine$double.eps))*log(.Machine$double.eps)
  # Compute likelihood of all data:
  scaLogLik <- scaLogLikZeros + scaLogLikNonzeros
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Cost function impulse model fit - Single cell mode
#' 
#' Log likelihood cost function for impulse model fit based on zero inflated  
#' negative binomial mode. This cost function is appropriate
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
  vecboolZero,
  scaWindowRadius=NULL){  

  # Generate negative binomial mean parameters from impulse model:
  # This links the parameters vecTheta which are changed in optimisation
  # to the cost function.
  vecImpulseValue <- calcImpulse_comp(vecTheta,vecX)[vecindTimepointAssign]*
    vecNormConst
  
  # Catch cases in which entire sample is zero observation
  # and impulse value is pushed to close to zero to be evaluated
  # as a negative binomial mean without numerical error.
  vecImpulseValue[vecImpulseValue < .Machine$double.eps] <- .Machine$double.eps
  
  # Evaluate likelihood (this is the cost function): 
  if(is.null(scaWindowRadius)){
    # Evaluate likelihood on individual cells.  
    scaLogLik <- evalLogLikZINB_comp(vecY=vecY,
      vecMu=vecImpulseValue,
      scaDispEst=scaDispEst, 
      vecDropoutRateEst=vecDropoutRateEst, 
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero)
    
  } else {
    # Evaluate likelihood on window of cells for each impulse value at a cell.
    # This is a local smoothing penalty added to the cost function.
    scaLogLik <- 0
    for(indcell in seq(1,length(vecY))){
      scaWindowStart <- max(0,indcell-scaWindowRadius)
      scaWindowEnd <- min(length(vecY),indcell+scaWindowRadius)
      scaLogLik <- scaLogLik + evalLogLikZINB_comp(
        vecY=vecY[scaWindowStart:scaWindowEnd],
        vecMu=vecImpulseValue[scaWindowStart:scaWindowEnd],
        scaDispEst=scaDispEst, 
        vecDropoutRateEst=vecDropoutRateEst[scaWindowStart:scaWindowEnd], 
        vecboolNotZeroObserved=vecboolNotZeroObserved[scaWindowStart:scaWindowEnd], 
        vecboolZero=vecboolZero[scaWindowStart:scaWindowEnd])
    }
  }
  
  return(scaLogLik)
}