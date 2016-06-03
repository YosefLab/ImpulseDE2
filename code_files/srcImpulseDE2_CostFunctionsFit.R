#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++     Cost Functions    +++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cost function impulse model fit - Batch mode
#' 
#' Log likelihood cost function for impulse model fit based on negative 
#' binomial model. Batch mode means that residuals are assumed to be 
#' independent.
#' In analogy to generalised linear models, a log linker
#' function is used for the count parameters. The inferred negative
#' binomial model is normalised by the factors in vecNormConst
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

# to be deprecated DAVID
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cost function impulse model fit - Time course mode
#' 
#' Log likelihood cost function for impulse model fit based on negative 
#' binomial model. Time course means that samples of samples can be
#' grouped into time course experiments, i.e. residuals are dependent
#' within a time course. The impulse model is scaled to the
#' mean of each time course.
#' In analogy to generalised linear models, a log linker
#' function is used for the count parameters. The inferred negative
#' binomial model is normalised by the factors in vecNormConst
#' for evaluation of the likelihood on the data.
#' 
#' @aliases evalLogLikImpulseByTC_comp
#' 
#' @seealso Called by \code{fitImpulse}:\code{fitImpulse_gene}.
#' Calls \code{calcImpulse}. \code{evalLogLikImpulseBatch} for
#' independent residuals.
#' 
#' @param vecTheta (vector number of parameters [6]) Impulse model parameters.
#' @param vecX (numeric vector number of timepoints) 
#'    Time-points at which gene was sampled.
#' @param vecY (numeric vector samples) Observed expression values for 
#     given gene.
#' @param scaDispEst: (scalar) Dispersion estimate for given gene.
#' @param vecNormConst: (count vector number of samples) 
#'    Normalisation constants for each sample.
#' @param vecTranslationFactors: (numeric vector number of samples)
#'    Scaling factors for impulse model between different timecourses:
#'    Mean timecourse/overall mean, each scaled by size factors.
#' @param vecindTimepointAssign (numeric vector number samples) 
#'    Index of time point assigned to sample in list of sorted
#'    time points (vecX).
#' @param vecboolObserved: (bool vector number of samples)
#'    Stores bool of sample being not NA (observed).
#'    
#' @return scaLogLik: (scalar) Value of cost function (likelihood) for given gene.
#' @export

evalLogLikImpulseByTC <- function(vecTheta,
  vecX,
  vecY,
  scaDispEst, 
  vecNormConst,
  vecTranslationFactors, 
  vecindTimepointAssign, 
  vecboolObserved){  
  # Compute normalised impulse function value: 
  # Mean of negative binomial density at each time point,
  # scaled by normalisation factor of each sample
  # and scaled by translation factor (one for each
  # longitudinal series).
  vecImpulseValue <- calcImpulse_comp(vecTheta,vecX)[vecindTimepointAssign]*
    vecNormConst*
    vecTranslationFactors
  
  # Compute log likelihood under impulse model by
  # adding log likelihood of model at each timepoint.
  scaLogLik <- sum(dnbinom(
    vecY[vecboolObserved], 
    mu=vecImpulseValue[vecboolObserved],   
    size=scaDispEst, 
    log=TRUE))
  
  #DAVID to do: take vecboolObserved out of these functions and reduce data at input to cost function
  
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
#' In analogy to generalised linear models, a log linker
#' function is used for the count parameters. The inferred negative
#' binomial model is normalised by the factors in vecNormConst
#' for evaluation of the likelihood on the data.
#' 
#' @aliases evalLogLikImpulseSC_comp
#' 
#' @seealso Called by \code{fitImpulse}:\code{fitImpulse_gene}.
#' Calls \code{calcImpulse}.
#' 
#' @param vecTheta (vector number of parameters [6]) Impulse model parameters.
#' @param vecX (numeric vector number of timepoints) 
#'    Time-points at which gene was sampled.
#' @param vecY (count vector samples) Observed expression values for 
#     given gene.
#' @param scaDispEst: (scalar) Dispersion estimate for given gene.
#' @param vecDropoutRateEst: (probability vector number of samples) 
#'    Dropout rate estimate for each cell for given gene.
#' @param vecNormConst: (numeric vector number of samples) 
#'    Normalisation constants for each sample.
#' @param vecindTimepointAssign (numeric vector number samples) 
#'    Index of time point assigned to sample in list of sorted
#'    time points (vecX).
#' @param vecboolObserved: (bool vector number of samples)
#'    Stores bool of sample being not NA (observed).
#' @param vecboolZero: (bool vector number of samples)
#'    Stores bool of sample having a count of 0.
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
  # Compute normalised impulse function value: 
  # Mean of negative binomial density at each time point,
  # scaled by normalisation factor of each sample.
  vecImpulseValue <- calcImpulse_comp(vecTheta,vecX)[vecindTimepointAssign]*
    vecNormConst
  
  # Compute log likelihood under impulse model by
  # adding log likelihood of model at each timepoint.
  # Likelihood function of hurdle modle differs between
  # zero and non-zero counts: Add likelihood of both together.
  # Note that the log is taken over the sum of mixture model 
  # components of the likelihood model and accordingly log=FALSE
  # inside dnbinom.
  # Likelihood of zero counts:
  vecLikZeros <- (1-vecDropoutRateEst[vecboolZero])*
    dnbinom(
      vecY[vecboolZero], 
      mu=vecImpulseValue[vecboolZero], 
      size=scaDispEst, 
      log=FALSE) +
    vecDropoutRateEst[vecboolZero]
  # Replace zero likelihood observation with machine precision
  # for taking log.
  scaLogLikZeros <- sum( log(vecLikZeros[vecLikZeros!=0]) +
      sum(vecLikZeros==0)*log(.Machine$double.eps) )
  # Likelihood of non-zero counts:
  vecLikNonzeros <- (1-vecDropoutRateEst[vecboolNotZeroObserved])*
    dnbinom(
      vecY[vecboolNotZeroObserved], 
      mu=vecImpulseValue[vecboolNotZeroObserved], 
      size=scaDispEst, 
      log=FALSE)
  # Replace zero likelihood observation with machine precision
  # for taking log.
  scaLogLikNonzeros <- sum( log(vecLikNonzeros[vecLikNonzeros!=0]) +
      sum(vecLikNonzeros==0)*log(.Machine$double.eps) )
  # Compute likelihood of all data:
  scaLogLik <- scaLogLikZeros + scaLogLikNonzeros
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#DAVID to be deprecated from here
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cost function hurdle model fit
#' 
#' Log likelihood cost function for hurdle model fit (zero inflated  negative 
#' binomial model). This function is designed to allow numerical optimisation
#' of negative binomial mean and drop-out rate on the entire data set given
#' the gene-wise overdispersion estimates.
#' 
#' @aliases evalLogLikHurdle_com
#' 
#' @seealso Called by \code{fitHurdleModel}.
#' 
#' @param lsThetaHurdle (list {vecDropoutRateEst, matMeanEst}) Parmeters
#'    of hurdle model to be estimated:
#'    vecDropoutRate: (vector cells) Dropout rate estimate for all cells.
#'    matMuEst: (matrix genes x clusters) Negative binomial 
#'        component mean estimates for each gene and each cluster.
#' @param vecY (matrix genes x cells) Observed expression values of all genes
#'    in all cells.
#' @param vecDispEst: (vector genes) Dispersion estimates for all genes. 
#' @param lsResultsClustering:
#' 
#' @return scaLogLik: (scalar) Value of cost function (likelihood)
#'    of hurdle model for dataset.
#' @export

evalLogLikHurdle <- function(lsThetaHurdle,vecY,
  vecDispEst,lsResultsClustering){ 
  
  vecDropoutRateEst <- lsThetaHurdle$vecDropoutRateEst
  matMeanEst <- lsThetaHurdle$matMeanEst
  # Likelihood function of hurdle modle differs between
  # zero and non-zero counts: Add likelihood of both together.
  # Note that the log is taken over the sum of mixture model 
  # components of the likelihood model and accordingly log=FALSE
  # inside dnbinom.
  
  # Create matrices of zero and non-zero counts with all elements
  # outside the target interval ({0} or {x>0 and finite}) set
  # to NA. NA will be ignored in likelihood computation below.
  vecYZeros <- vecY
  vecYZeros[vecYZeros!=0] <- NA
  vecYNonzeros <- vecY
  vecYNonzeros[vecYNonzeros <= 0 | is.na(vecYNonzeros) | !is.finite(vecYNonzeros)] <- NA
  
  # Loop over clusters
  scaLogLik <- sum(sapply( seq(1:lsResultsClustering$K), function(cluster){ 
    # Get indices of cells in current cluster
    vecCells <- which(lsResultsClustering$Assignments==cluster)
    
    # Evaluate on all zero count genes in cells
    scaLogLikZerosCluster <- sum(log(
      (1-vecDropoutRateEst[vecCells]) * dnbinom(vecYZeros[,vecCells], 
        mu=matMeanEst[,cluster], size=vecDispEst, log=FALSE) +
        vecDropoutRateEst[vecCells]
    ), na.rm=TRUE)
    # Evaluate on all nonzero count genes in cells
    scaLogLikNonzerosCluster <- sum(log(
      (1-vecDropoutRateEst[vecCells]) * dnbinom(vecYNonzeros[,vecCells], 
        mu=matMeanEst[,cluster], size=vecDispEst, log=FALSE)
    ), na.rm=TRUE)
    
    # Get likelihood of all data in cluster
    scaLogLikCluster <- scaLogLikZerosCluster + scaLogLikNonzerosCluster
    
    return(scaLogLikCluster)
  } ))
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Cost function hurdle model fit - Negative binomial
#' 
#' Log likelihood cost function for hurdle model fit (zero inflated  negative 
#' binomial model). This function is designed to allow numerical optimisation
#' of negative binomial mean and overdispersion on single gene given
#' the drop-out rate.
#' Note: During optimisation, this function does not make use of the
#' fact that there exists a closed-form MLE for the mean of a negative
#' binomial because the negative binomial is estimated within the hurdle model
#' here.
#' This function is used for iterative hurdle model fitting (separate fitting
#' of negative binomial and drop-out rate). This separation allows 
#' parallelisation of the individual fitting steps over genes/cells and
#' reduces the dimensionality of the individual optimisation problems
#' drastically.
#' 
#' @aliases evalLogLikHurdleNB_com
#' 
#' @seealso Called by \code{fitHurdleModel}.
#' 
#' @param vecTheta (vector parameters) Parmeters of negative
#'    binomial in hurdle model of one gene to be estimated:
#'    Overdispersion and one mean per cluster.
#' @param vecY (vector cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecDropoutRateEst: (vector cells) Dropout estimate of cell. 
#' @param lsResultsClustering:
#' 
#' @return scaLogLik: (scalar) Value of cost function (likelihood)
#'    of hurdle model for gene
#' @export

evalLogLikHurdleNB <- function(vecTheta,vecY,
  vecDropoutRateEst,lsResultsClustering){ 
  
  scaDispEst <- vecTheta[1]
  vecMeanEst <- vecTheta[2:length(vecTheta)]
  # Likelihood function of hurdle modle differs between
  # zero and non-zero counts: Add likelihood of both together.
  # Note that the log is taken over the sum of mixture model 
  # components of the likelihood model and accordingly log=FALSE
  # inside dnbinom.
  
  # Get indices of zero and non-zero counts.
  indZeros <- which(vecY==0)
  indNonzeros <- which(vecY > 0 & !is.na(vecY) & is.finite(vecY))
  
  # Evaluate on all zero count genes in cells
  scaLogLikZeros <- sum(log(
    (1-vecDropoutRateEst[indZeros]) * dnbinom(vecY[indZeros],
      mu=exp(vecMeanEst[lsResultsClustering$Assignments[indZeros]]), 
      size=scaDispEst[indZeros],
      log=FALSE) + vecDropoutRateEst[indZeros]
  ), na.rm=TRUE)
  # Evaluate on all nonzero count genes in cells
  scaLogLikNonzeros <- sum(log(
    (1-vecDropoutRateEst[indNonzeros]) * dnbinom(vecY[indNonzeros],
      mu=exp(vecMeanEst[lsResultsClustering$Assignments[indNonzeros]]),
      size=scaDispEst[indNonzeros],
      log=FALSE)
  ), na.rm=TRUE)
  
  # Get likelihood of all data in cluster
  scaLogLik <- scaLogLikZeros + scaLogLikNonzeros
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

#' Cost function hurdle model fit - Dropout
#' 
#' Log likelihood cost function for hurdle model fit (zero inflated  negative 
#' binomial model). This function is designed to allow numerical optimisation
#' of drop-out rate of a cell given the parameterisation of the negative 
#' binomial distribution of each gene (mean and overdispersion).
#' This function is used for iterative hurdle model fitting (separate fitting
#' of negative binomial and drop-out rate). This separation allows 
#' parallelisation of the individual fitting steps over genes/cells and
#' reduces the dimensionality of the individual optimisation problems
#' drastically.
#' 
#' @aliases evalLogLikHurdleDrop_com
#' 
#' @seealso Called by \code{fitHurdleModel}.
#' 
#' @param scaTheta (scalar) Estimate of dropout rate.
#' @param vecY (vector genes) Observed expression values 
#'    of all genes in cell.
#' @param vecMeanEst: (vector genes) Negative binomial 
#'    mean estimates of genes.
#' @param vecDispEst: (vector genes) Negative binomial 
#'    overdispersion estimates of all genes.
#' 
#' @return scaLogLik: (scalar) Value of cost function (likelihood)
#'    of hurdle model for gene
#' @export

evalLogLikHurdleDrop <- function(scaTheta, vecY,
  vecMeanEst, vecDispEst,
  vecboolNotZeroObserved, vecboolZero){ 
  
  scaDropoutRateEst <- scaTheta
  # Likelihood function of hurdle modle differs between
  # zero and non-zero counts: Add likelihood of both together.
  # Note that the log is taken over the sum of mixture model 
  # components of the likelihood model and accordingly log=FALSE
  # inside dnbinom.
  
  # Likelihood of zero counts:
  vecLikZeros <- (1-scaDropoutRateEst)*
    dnbinom(
      vecY[vecboolZero], 
      mu=vecMeanEst[vecboolZero], 
      size=vecDispEst[vecboolZero], 
      log=FALSE) +
    scaDropoutRateEst
  # Replace zero likelihood observation with machine precision
  # for taking log.
  scaLogLikZeros <- sum( log(vecLikZeros[vecLikZeros!=0]) +
      sum(vecLikZeros==0)*log(.Machine$double.eps) )
  # Likelihood of non-zero counts:
  vecLikNonzeros <- (1-scaDropoutRateEst)*
    dnbinom(
      vecY[vecboolNotZeroObserved], 
      mu=vecMeanEst[vecboolNotZeroObserved], 
      size=vecDispEst[vecboolNotZeroObserved], 
      log=FALSE)
  # Replace zero likelihood observation with machine precision
  # for taking log.
  scaLogLikNonzeros <- sum( log(vecLikNonzeros[vecLikNonzeros!=0]) +
      sum(vecLikNonzeros==0)*log(.Machine$double.eps) )
  # Compute likelihood of all data:
  scaLogLik <- scaLogLikZeros + scaLogLikNonzeros
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}
