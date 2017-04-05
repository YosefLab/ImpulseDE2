### Cost Functions

#' Cost function for constant model
#' 
#' Log likelihood cost function for numerical optimisation of
#' constant model.
#' Implements log linker function for the constant mean parameter
#' and the batch correction factors.
#' Implements lower sensitivity bound of likelihood with respect
#' to constant mean parameter.
#' Implements upper and lower sensitivity bound of likelihood with respect
#' to batch correction factors.
#' 
#' @seealso Compiled version: \link{evalLogLikMu_comp}
#' 
#' @param vecTheta (numeric vector number of parameters to be estimated) 
#' Constant model parameter and batch correction factor estimates.
#' @param vecCounts (numeric vector number of samples)
#' Read count data.
#' @param scaDisp (scalar) Gene-wise 
#' negative binomial dispersion hyper-parameter.
#' @param vecSizeFactors (numeric vector number of samples) 
#' Model scaling factors for each sample which take
#' sequencing depth into account (size factors).
#' @param lsvecidxBatch (list length number of confounding variables)
#' List of index vectors. 
#' One vector per confounding variable.
#' Each vector has one entry per sample with the index batch
#' within the given confounding variable of the given sample.
#' Batches are enumerated from 1 to number of batches.
#' @param vecboolObserved (bool vector number of samples)
#' Whether sample is observed (finite and not NA).
#'  
#' @return scaLogLik (scalar) Value of cost function 
#' (loglikelihood) for given gene.
#' 
#' @author David Sebastian Fischer
evalLogLikMu <- function(
    vecTheta, vecCounts, scaDisp, vecSizeFactors, lsvecidxBatch, 
    vecboolObserved) {
    
    scaMu <- exp(vecTheta[1])
    scaNParamUsed <- 1
    # Prevent mean shrinkage to zero:
    if (scaMu < 10^(-10)) {
        scaMu <- 10^(-10)
    }
    
    vecBatchFactors <- array(1, length(vecCounts))
    if (!is.null(lsvecidxBatch)) {
        for (vecidxConfounder in lsvecidxBatch) {
            scaNBatchFactors <- max(vecidxConfounder) - 1  
            # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining 
            # factors scale based on the first batch.
            vecBatchFacConf <- c(
                1, exp(vecTheta[(scaNParamUsed + 1):
                                    (scaNParamUsed + scaNBatchFactors)])
            )[vecidxConfounder]
            scaNParamUsed <- scaNParamUsed + scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchFacConf[vecBatchFacConf < 10^(-10)] <- 10^(-10)
            vecBatchFacConf[vecBatchFacConf > 10^(10)] <- 10^(10)
            
            vecBatchFactors <- vecBatchFactors * vecBatchFacConf
        }
    }
    
    # Compute log likelihood under constant model by adding log likelihood
    # of model at each timepoint.
    scaLogLik <- sum(dnbinom(
        vecCounts[vecboolObserved], 
        mu = scaMu * vecBatchFactors[vecboolObserved] * 
            vecSizeFactors[vecboolObserved], 
        size = scaDisp, log = TRUE))
    
    # Maximise log likelihood: Return likelihood as value to optimisation
    # routine
    return(scaLogLik)
}

#' Compiled function: evalLogLikMu
#' 
#' Pre-compile heavily used functions.
#' Refer to \link{evalLogLikMu}.
#' 
#' @param vecTheta (numeric vector number of parameters to be estimated) 
#' Constant model parameter and batch correction factor estimates.
#' @param vecCounts (numeric vector number of samples)
#' Read count data.
#' @param scaDisp (scalar) Gene-wise 
#' negative binomial dispersion hyper-parameter.
#' @param vecSizeFactors (numeric vector number of samples) 
#' Model scaling factors for each sample which take
#' sequencing depth into account (size factors).
#' @param lsvecidxBatch (list length number of confounding variables)
#' List of index vectors. 
#' One vector per confounding variable.
#' Each vector has one entry per sample with the index batch
#' within the given confounding variable of the given sample.
#' Batches are enumerated from 1 to number of batches.
#' @param vecboolObserved (bool vector number of samples)
#' Whether sample is observed (finite and not NA).
#'  
#' @return scaLogLik (scalar) Value of cost function 
#' (loglikelihood) for given gene.
#' 
#' @author David Sebastian Fischer
evalLogLikMu_comp <- cmpfun(evalLogLikMu)

#' Cost function for impulse model
#' 
#' Log likelihood cost function for numerical optimisation of
#' impulse model.
#' Implements log linker function for the amplitude parameters
#' and the batch correction factors.
#' Implements upper and lower sensitivity bound of likelihood with respect
#' to batch correction factors and lower bound for amplitude paramters.
#' 
#' @seealso Compiled version: \link{evalLogLikImpulse_comp}
#' 
#' @param vecTheta (numeric vector number of parameters to be estimated) 
#' Impulse model parameter and batch correction factor estimates.
#' @param vecCounts (numeric vector number of samples)
#' Read count data.
#' @param scaDisp (scalar) Gene-wise 
#' negative binomial dispersion hyper-parameter.
#' @param vecSizeFactors (numeric vector number of samples) 
#' Model scaling factors for each sample which take
#' sequencing depth into account (size factors).
#' @param vecTimepointsUnique 
#' (numeric vector length number of unique time points)
#' Unique time points of set of time points of given samples.
#' @param vecidxTimepoint (index vector length number of samples)
#' Index of of time point assigned to each sample in vector
#' vecTimepointsUnique.
#' @param lsvecidxBatch (list length number of confounding variables)
#' List of index vectors. 
#' One vector per confounding variable.
#' Each vector has one entry per sample with the index batch
#' within the given confounding variable of the given sample.
#' Batches are enumerated from 1 to number of batches.
#' @param vecboolObserved (bool vector number of samples)
#' Whether sample is observed (finite and not NA).
#'  
#' @return scaLogLik (scalar) Value of cost function 
#' (loglikelihood) for given gene.
#' 
#' @author David Sebastian Fischer
evalLogLikImpulse <- function(
    vecTheta, vecCounts, scaDisp, vecSizeFactors, 
    vecTimepointsUnique, vecidxTimepoint, lsvecidxBatch, vecboolObserved) {
    
    # Compute normalised impulse function value:
    vecImpulseParam <- vecTheta[1:6]
    vecImpulseParam[2:4] <- exp(vecImpulseParam[2:4])
    vecImpulseValue <- evalImpulse_comp(
        vecImpulseParam = vecImpulseParam, 
        vecTimepoints = vecTimepointsUnique)[vecidxTimepoint]
    scaNParamUsed <- 6
    
    vecBatchFactors <- array(1, length(vecCounts))
    if (!is.null(lsvecidxBatch)) {
        for (vecidxConfounder in lsvecidxBatch) {
            scaNBatchFactors <- max(vecidxConfounder) - 1  
            # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining 
            # factors scale based on the first batch.
            vecBatchFacConf <- c(
                1, exp(vecTheta[(scaNParamUsed +1):
                                    (scaNParamUsed + scaNBatchFactors)])
            )[vecidxConfounder]
            scaNParamUsed <- scaNParamUsed + scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchFacConf[vecBatchFacConf < 10^(-10)] <- 10^(-10)
            vecBatchFacConf[vecBatchFacConf > 10^(10)] <- 10^(10)
            
            vecBatchFactors <- vecBatchFactors * vecBatchFacConf
        }
    }
    
    # Compute log likelihood under impulse model by adding log likelihood of
    # model at each timepoint.
    scaLogLik <- sum(dnbinom(
        vecCounts[vecboolObserved], 
        mu = vecImpulseValue[vecboolObserved] * vecBatchFactors[vecboolObserved] * 
            vecSizeFactors[vecboolObserved], 
        size = scaDisp, log = TRUE))
    
    # Maximise log likelihood: Return likelihood as value to optimisation
    # routine
    return(scaLogLik)
}

#' Compiled function: evalLogLikImpulse
#' 
#' Pre-compile heavily used functions.
#' Refer to \link{evalLogLikImpulse}.
#' 
#' @param vecTheta (numeric vector number of parameters to be estimated) 
#' Impulse model parameter and batch correction factor estimates.
#' @param vecCounts (numeric vector number of samples)
#' Read count data.
#' @param scaDisp (scalar) Gene-wise 
#' negative binomial dispersion hyper-parameter.
#' @param vecSizeFactors (numeric vector number of samples) 
#' Model scaling factors for each sample which take
#' sequencing depth into account (size factors).
#' @param vecTimepointsUnique
#' (numeric vector length number of unique time points)
#' Unique time points of set of time points of given samples.
#' @param vecidxTimepoint (index vector length number of samples)
#' Index of of time point assigned to each sample in vector
#' vecTimepointsUnique.
#' @param lsvecidxBatch (list length number of confounding variables)
#' List of index vectors. 
#' One vector per confounding variable.
#' Each vector has one entry per sample with the index batch
#' within the given confounding variable of the given sample.
#' Batches are enumerated from 1 to number of batches.
#' @param vecboolObserved (bool vector number of samples)
#' Whether sample is observed (finite and not NA).
#'  
#' @return scaLogLik (scalar) Value of cost function 
#' (loglikelihood) for given gene.
#' 
#' @author David Sebastian Fischer
evalLogLikImpulse_comp <- cmpfun(evalLogLikImpulse)

#' Cost function for sigmoidal model
#' 
#' Log likelihood cost function for numerical optimisation of
#' sigmoidal model.
#' Implements log linker function for the amplitude parameters
#' and the batch correction factors.
#' Implements upper and lower sensitivity bound of likelihood with respect
#' to batch correction factors and lower bound for amplitude paramters.
#' 
#' @seealso Compiled version: \link{evalLogLikSigmoid_comp}
#' 
#' @param vecTheta (numeric vector number of parameters to be estimated) 
#' Sigmoid model parameter and batch correction factor estimates.
#' @param vecCounts (numeric vector number of samples)
#' Read count data.
#' @param scaDisp (scalar) Gene-wise 
#' negative binomial dispersion hyper-parameter.
#' @param vecSizeFactors (numeric vector number of samples) 
#' Model scaling factors for each sample which take
#' sequencing depth into account (size factors).
#' @param vecTimepointsUnique
#' (numeric vector length number of unique time points)
#' Unique time points of set of time points of given samples.
#' @param vecidxTimepoint (index vector length number of samples)
#' Index of of time point assigned to each sample in vector
#' vecTimepointsUnique.
#' @param lsvecidxBatch (list length number of confounding variables)
#' List of index vectors. 
#' One vector per confounding variable.
#' Each vector has one entry per sample with the index batch
#' within the given confounding variable of the given sample.
#' Batches are enumerated from 1 to number of batches.
#' @param vecboolObserved (bool vector number of samples)
#' Whether sample is observed (finite and not NA).
#'  
#' @return scaLogLik (scalar) Value of cost function 
#' (loglikelihood) for given gene.
#' 
#' @author David Sebastian Fischer
evalLogLikSigmoid <- function(
    vecTheta, vecCounts, scaDisp, vecSizeFactors, 
    vecTimepointsUnique, vecidxTimepoint, lsvecidxBatch, vecboolObserved) {
    
    # Compute normalised impulse function value:
    vecSigmoidParam <- vecTheta[1:4]
    vecSigmoidParam[2:3] <- exp(vecSigmoidParam[2:3])
    vecSigmoidValue <- evalSigmoid_comp(
        vecSigmoidParam = vecSigmoidParam, 
        vecTimepoints = vecTimepointsUnique)[vecidxTimepoint]
    scaNParamUsed <- 4
    
    vecBatchFactors <- array(1, length(vecCounts))
    if (!is.null(lsvecidxBatch)) {
        for (vecidxConfounder in lsvecidxBatch) {
            scaNBatchFactors <- max(vecidxConfounder) - 1  
            # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining 
            # factors scale based on the first batch.
            vecBatchFacConf <- c(
                1, exp(vecTheta[(scaNParamUsed +1):
                                    (scaNParamUsed + scaNBatchFactors)])
            )[vecidxConfounder]
            scaNParamUsed <- scaNParamUsed + scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchFacConf[vecBatchFacConf < 10^(-10)] <- 10^(-10)
            vecBatchFacConf[vecBatchFacConf > 10^(10)] <- 10^(10)
            
            vecBatchFactors <- vecBatchFactors * vecBatchFacConf
        }
    }
    
    # Compute log likelihood under impulse model by adding log likelihood of
    # model at each timepoint.
    scaLogLik <- sum(dnbinom(
        vecCounts[vecboolObserved], 
        mu = vecSigmoidValue[vecboolObserved] *
            vecBatchFactors[vecboolObserved] * 
            vecSizeFactors[vecboolObserved], 
        size = scaDisp, log = TRUE))
    
    # Maximise log likelihood: Return likelihood as value to optimisation
    # routine
    return(scaLogLik)
}

#' Compiled function: evalLogLikSigmoid
#' 
#' Pre-compile heavily used functions.
#' Refer to \link{evalLogLikSigmoid}.
#' 
#' @param vecTheta (numeric vector number of parameters to be estimated) 
#' Sigmoid model parameter and batch correction factor estimates.
#' @param vecCounts (numeric vector number of samples)
#' Read count data.
#' @param scaDisp (scalar) Gene-wise 
#' negative binomial dispersion hyper-parameter.
#' @param vecSizeFactors (numeric vector number of samples) 
#' Model scaling factors for each sample which take
#' sequencing depth into account (size factors).
#' @param vecTimepointsUnique
#' (numeric vector length number of unique time points)
#' Unique time points of set of time points of given samples.
#' @param vecidxTimepoint (index vector length number of samples)
#' Index of of time point assigned to each sample in vector
#' vecTimepointsUnique.
#' @param lsvecidxBatch (list length number of confounding variables)
#' List of index vectors. 
#' One vector per confounding variable.
#' Each vector has one entry per sample with the index batch
#' within the given confounding variable of the given sample.
#' Batches are enumerated from 1 to number of batches.
#' @param vecboolObserved (bool vector number of samples)
#' Whether sample is observed (finite and not NA).
#'  
#' @return scaLogLik (scalar) Value of cost function 
#' (loglikelihood) for given gene.
#' 
#' @author David Sebastian Fischer
evalLogLikSigmoid_comp <- cmpfun(evalLogLikSigmoid)
