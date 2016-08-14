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
#' @aliases evalLogLikMuNB_comp
#' 
#' @seealso
#' 
#' @param scaTheta (vector number of parameters [6]) 
#'    Negative binomial mean parameter estimate (log).
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

evalLogLikMuNB <- function(scaTheta,
  vecCounts,
  scaDispEst, 
  vecNormConst,
  vecboolObserved){
  
  scaMu <- exp(scaTheta)
  
  # Prevent means estimate from shrinking to zero
  # to avoid numerical errors:
  if(scaMu < .Machine$double.eps){ scaMu <- .Machine$double.eps }
  
  # Compute log likelihood under impulse model by
  # adding log likelihood of model at each timepoint.
  scaLogLik <- sum(dnbinom(
    vecCounts[vecboolObserved], 
    mu=scaMu*vecNormConst[vecboolObserved], 
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
#' factors are 1. Catches numerical errors. Fit in log space
#' to guarantee positive mean paramter.
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

fitMuNB <- function(vecCounts,
  scaDispEst,
  vecNormConst){
  
  if(all(vecNormConst==1)){
    # Closed form maximum likelihood estimator
    scaMu <- mean(vecCounts, na.rm=TRUE)
  } else {
    # Numerical maximum likelihood estimator
    scaMu <- tryCatch({
      exp(unlist(optimise(
        evalLogLikMuNB_comp,
        vecCounts=vecCounts,
        scaDispEst=scaDispEst, 
        vecNormConst=vecNormConst,
        vecboolObserved=!is.na(vecCounts),
        lower = log(.Machine$double.eps),
        upper = log(max(vecNormConst*vecCounts, na.rm=TRUE)+1),
        maximum = TRUE)["maximum"]))
    }, error=function(strErrorMsg){
      print(paste0("ERROR: Fitting negative binomial mean parameter: fitMuNB().",
        " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
      print(paste0("vecCounts ", paste(vecCounts,sep=" ")))
      print(paste0("scaDispEst ", paste(scaDispEst,sep=" ")))
      print(paste0("vecNormConst ", paste(vecNormConst,sep=" ")))
      lsErrorCausingGene <- list(vecCounts, scaDispEst, vecNormConst)
      names(lsErrorCausingGene) <- c("vecCounts", "scaDispEst", "vecNormConst")
      save(lsErrorCausingGene,file=file.path(getwd(),"ImpulseDE2_lsErrorCausingGene.RData"))
      stop(strErrorMsg)
    })
  }
  
  # Catch boundary of likelihood domain on mu space
  if(scaMu < .Machine$double.eps){scaMu <- .Machine$double.eps}
  
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
#' @param vecMu (vector number of cell) Negative binomial
#'    mean parameter for each cell.
#' @param vecDispEst: (numerical vector number of samples) 
#'    Negative binomial dispersion parameter estimates for given gene.
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
  vecDispEst, 
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
  scaLogPrecLim <- -323*log(10)

  # Likelihood of zero counts:
  vecLogLikZeros <- log((1-vecDropoutRateEst[vecboolZero])*
    (vecDispEst[vecboolZero]/(vecDispEst[vecboolZero] + 
        vecMu[vecboolZero]))^vecDispEst[vecboolZero] +
    vecDropoutRateEst[vecboolZero])
  vecLogLikZeros[vecLogLikZeros < scaLogPrecLim] <- scaLogPrecLim
  scaLogLikZeros <- sum(vecLogLikZeros)
  # Likelihood of non-zero counts:
  vecLogLikNonzerosDropout <- log(1-vecDropoutRateEst[vecboolNotZeroObserved])
  vecLogLikNonzerosNB <- dnbinom(
      vecY[vecboolNotZeroObserved], 
      mu=vecMu[vecboolNotZeroObserved], 
      size=vecDispEst[vecboolNotZeroObserved], 
      log=TRUE)
  vecLogLikNonzerosDropout[vecLogLikNonzerosDropout < scaLogPrecLim] <- scaLogPrecLim
  vecLogLikNonzerosNB[vecLogLikNonzerosNB < scaLogPrecLim] <- scaLogPrecLim
  # Replace zero likelihood observation with machine precision
  # for taking log.
  scaLogLikNonzeros <- sum( vecLogLikNonzerosDropout+vecLogLikNonzerosNB )
  # Compute likelihood of all data:
  scaLogLik <- scaLogLikZeros + scaLogLikNonzeros
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Compute smoothed log likelihood of zero-inflated negative binomial model for one gene
#' 
#' This liklihood function is appropriate for sequencing data with high drop 
#' out rate, commonly observed in single cell data (e.g. scRNA-seq). It includes
#' a smoothin penalty on the means.
#' 
#' @aliases evalLogLikSmoothZINB_comp
#' 
#' @seealso Called by \code{fitZINB} and
#' \code{runModelFreeDEAnalysis}.
#'
#' @param vecY (count vector number of amples)
#'    Observed expression values for  given gene.
#' @param vecMu (vector number of samples) Negative binomial
#'    mean parameter for each sample.
#' @param vecSizeFactors: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecDispEst: (scalar vector number of samples) 
#'    Negative binomial dispersion  parameter for given 
#'    gene and observations.
#' @param vecDropoutRateEst: (probability vector number of samples) 
#'    Dropout rate estimate for each cell for given gene.
#' @param vecboolNotZeroObserved: (bool vector number of samples)
#'    Whether sample is not zero and observed (not NA).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval length.
#'    
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikSmoothZINB <- function(vecY,
  vecMu,
  vecSizeFactors,
  vecDispEst, 
  vecDropoutRateEst, 
  vecboolNotZeroObserved, 
  vecboolZero,
  scaWindowRadius=NULL ){
  
  scaNumCells <- length(vecY)
  scaLogLik <- sum(sapply(seq(1,scaNumCells), 
    function(j){
      scaindIntervalStart <- max(1,j-scaWindowRadius)
      scaindIntervalEnd <- min(scaNumCells,j+scaWindowRadius)
      vecInterval <- seq(scaindIntervalStart,scaindIntervalEnd)
      scaLogLikCell <- evalLogLikZINB_comp(vecY=vecY[vecInterval],
        vecMu=vecMu[j]*vecSizeFactors[vecInterval],
        vecDispEst=rep(vecDispEst[j], length(vecInterval)), 
        vecDropoutRateEst=vecDropoutRateEst[vecInterval], 
        vecboolNotZeroObserved[vecInterval], 
        vecboolZero[vecInterval])
      return(scaLogLikCell)
    }))
  return(scaLogLik)
}

#' Cost function zero-inflated negative binomial model for mean fitting
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of negative binomial mean paramater on single gene given
#' the drop-out rate and negative binomial dispersion parameter. The
#' mean parameter is fit in log space and is therefore fit
#' as a positive scalar. The cost function is insensitive to the
#' dispersion factor shrinking beyond a numerical threshold to zero
#' to avoid shrinkage of the dispersion factor to zero which 
#' may cause numerical errors. Accordingly, growth above a numerical
#' threshold to infinity (this correponds to Poissonian noise) is 
#' also guarded against.
#' Smoothing is applied if scaWindowRadius is handed as not NULL.
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param scaTheta: (scalar) Log of mean parameter estimate.
#' @param vecY: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecDispEst: (scalar vector number of samples) 
#'    Negative binomial dispersion  parameter for given 
#'    gene and observations.
#' @param vecSizeFactors: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecDropoutRateEst: (vector number of cells) Dropout estimate of cell.
#' @param vecboolObserved: (bool vector number of samples)
#'    Whether sample is not NA (observed).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval length.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

evalLogLikMuZINB <- function(scaTheta,
  vecY,
  vecDispEst,
  vecNormConst,
  vecDropoutRateEst,
  vecboolNotZeroObserved,
  vecboolZero,
  scaWindowRadius=NULL ){ 
  
  # Log linker function to fit positive means
  scaMu <- exp(scaTheta)
  
  # Prevent means estimate from shrinking to zero
  # to avoid numerical errors:
  if(scaMu < .Machine$double.eps){ scaMu <- .Machine$double.eps }
  
  if(is.null(scaWindowRadius)){
    scaLogLik <- evalLogLikZINB_comp( vecY=vecY,
      vecMu=scaMu*vecNormConst,
      vecDispEst=vecDispEst, 
      vecDropoutRateEst=vecDropoutRateEst,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero )
  } else {
    scaLogLik <- evalLogLikSmoothZINB_comp( vecY=vecY,
      vecMu=rep(scaMu, length(vecY)),
      vecSizeFactors=vecNormConst,
      vecDispEst=vecDispEst, 
      vecDropoutRateEst=vecDropoutRateEst,
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero,
      scaWindowRadius=scaWindowRadius)
  }
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Cost function zero-inflated negative binomial model for mean fitting
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of negative binomial mean paramater on single gene given
#' the drop-out rate and negative binomial dispersion parameter. The
#' mean parameter is fit in log space and is therefore fit
#' as a positive scalar. The cost function is insensitive to the
#' dispersion factor shrinking beyond a numerical threshold to zero
#' to avoid shrinkage of the dispersion factor to zero which 
#' may cause numerical errors. Accordingly, growth above a numerical
#' threshold to infinity (this correponds to Poissonian noise) is 
#' also guarded against.
#' 
#' @seealso Called by \code{computeLogLikNull}.
#' 
#' @param vecCounts: (vector number of cells) Observed expression values 
#'    of gene in cells.
#' @param scaDispEst: (vector number of cells) Negative binomial
#'    dispersion parameter estimate.
#' @param vecNormConst: (numeric vector number of cells) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per cell.
#' @param vecDropoutRateEst: (vector number of cells) Dropout estimate of cell.
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval length.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

fitMuZINB <- function(vecCounts,
  scaDispEst,
  vecNormConst,
  vecDropoutRateEst,
  vecProbNB,
  scaWindowRadius){ 
  
  if(all(vecNormConst==1) & is.null(scaWindowRadius)){
    # Closed form maximum likelihood estimator
    scaMu <- sum(vecCounts*vecProbNB, na.rm=TRUE)/sum(vecProbNB, na.rm=TRUE)
  } else {
    # Numerical maximum likelihood estimator
    scaMu <- tryCatch({
      exp(unlist(optimise(
        evalLogLikMuZINB_comp,
        vecCounts=vecCounts,
        vecDispEst=rep(scaDispEst, length(vecCounts)),
        vecDropoutRateEst=vecDropoutRateEst,
        vecNormConst=vecNormConst,
        vecboolNotZeroObserved=!is.na(vecCounts) & vecCounts>0,
        vecboolObserved=!is.na(vecCounts),
        scaWindowRadius=scaWindowRadius,
        lower = log(.Machine$double.eps),
        upper = log(max(vecNormConst*vecCounts, na.rm=TRUE)+1),
        maximum = TRUE)["maximum"]))
    }, error=function(strErrorMsg){
      print(paste0("ERROR: Fitting zero-inflated negative binomial mean parameter: fitMuZINB().",
        " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
      print(paste0("vecCounts ", paste(vecCounts,collapse=" ")))
      print(paste0("scaDispEst ", paste(scaDispEst,collapse=" ")))
      print(paste0("vecDropoutRateEst ", paste(vecDropoutRateEst,collapse=" ")))
      print(paste0("vecNormConst ", paste(vecNormConst,collapse=" ")))
      lsErrorCausingGene <- list(vecCounts, scaDispEst, vecDropoutRateEst, vecNormConst)
      names(lsErrorCausingGene) <- c("vecCounts", "scaDispEst", "vecDropoutRateEst","vecNormConst")
      save(lsErrorCausingGene,file=file.path(getwd(),"ImpulseDE2_lsErrorCausingGene.RData"))
      stop(strErrorMsg)
    })
  }
  
  # Catch boundary of likelihood domain on mu space
  if(scaMu < .Machine$double.eps){scaMu <- .Machine$double.eps}
  
  return(scaMu)
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
#' @param scaWindowRadius: (integer) 
#'    Smoothing interval length.
#'    
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikImpulseSC <- function(vecTheta,
  vecX,
  vecY,
  scaDispEst, 
  vecDropoutRateEst,
  matLinModelPi=NULL,
  vecNormConst,
  vecindTimepointAssign, 
  vecboolNotZeroObserved, 
  vecboolZero,
  scaWindowRadius=NULL){  
  
  # Generate negative binomial mean parameters from impulse model:
  # This links the parameters vecTheta which are changed in optimisation
  # to the cost function.
  vecImpulseValue <- calcImpulse_comp(vecTheta,vecX)[vecindTimepointAssign]
  
  # Catch cases in which entire sample is zero observation
  # and impulse value is pushed to close to zero to be evaluated
  # as a negative binomial mean without numerical error.
  vecImpulseValue[vecImpulseValue < .Machine$double.eps] <- .Machine$double.eps
  
  if(!is.null(matLinModelPi)){
    vecLinModelOut <- sapply(seq(1,length(vecImpulseValue)), function(cell){
      sum(matLinModelPi[cell,] * c(1,log(vecImpulseValue[cell])))
    })
    vecDropoutRateEst <- 1/(1+exp(-vecLinModelOut))
  }
  
  # Evaluate likelihood (this is the cost function): 
  if(is.null(scaWindowRadius)){
    # Evaluate likelihood on individual cells.  
    scaLogLik <- evalLogLikZINB_comp(vecY=vecY,
      vecMu=vecImpulseValue*vecNormConst,
      vecDispEst=rep(scaDispEst, length(vecY)), 
      vecDropoutRateEst=vecDropoutRateEst, 
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero)  
  } else {
    # Evaluate likelihood on individual cells.  
    scaLogLik <- evalLogLikSmoothZINB_comp(vecY=vecY,
      vecMu=vecImpulseValue,
      vecSizeFactors=vecNormConst,
      vecDispEst=rep(scaDispEst, length(vecY)), 
      vecDropoutRateEst=vecDropoutRateEst, 
      vecboolNotZeroObserved=vecboolNotZeroObserved, 
      vecboolZero=vecboolZero,
      scaWindowRadius=scaWindowRadius )
  }
  
  return(scaLogLik)
}