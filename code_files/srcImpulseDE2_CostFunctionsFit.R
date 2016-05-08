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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

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
#' @param lsThetaHurdle (list {vecDropoutEst, matMeanEst}) Parmeters
#'    of hurdle model to be estimated:
#'    vecDropoutRate: (vector cells) Dropout rate estimate for all cells.
#'    matMuEst: (matrix genes x clusters) Negative binomial 
#'        component mean estimates for each gene and each cluster.
#' @param matY (matrix genes x cells) Observed expression values of all genes
#'    in all cells.
#' @param vecDispEst: (vector genes) Dispersion estimates for all genes. 
#' @param lsResultsClustering:
#' 
#' @return scaLogLik: (scalar) Value of cost function (likelihood)
#'    of hurdle model for dataset.
#' @export

evalLogLikHurdle <- function(lsThetaHurdle,matY,
  vecDispEst,lsResultsClustering){ 
  
  vecDropoutEst <- lsThetaHurdle$vecDropoutEst
  matMeanEst <- lsThetaHurdle$matMeanEst
  # Likelihood function of hurdle modle differs between
  # zero and non-zero counts: Add likelihood of both together.
  # Note that the log is taken over the sum of mixture model 
  # components of the likelihood model and accordingly log=FALSE
  # inside dnbinom.
  
  # Create matrices of zero and non-zero counts with all elements
  # outside the target interval ({0} or {x>0 and finite}) set
  # to NA. NA will be ignored in likelihood computation below.
  matYZeros <- matY
  matYZeros[matYZeros!=0] <- NA
  matYNonzeros <- matY
  matYNonzeros[matYNonzeros <= 0 | is.na(matYNonzeros) | !is.finite(matYNonzeros)] <- NA
  
  # Loop over clusters
  scaLogLik <- sum(sapply( c(1:lsResultsClustering$K), function(cluster){ 
    # Get indices of cells in current cluster
    vecCells <- which(lsResultsClustering$Assignments==cluster)
    
    # Evaluate on all zero count genes in cells
    scaLogLikZerosCluster <- sum(log(
      (1-vecDropoutEst[vecCells]) * dnbinom(matYZeros[,vecCells], 
        mu=matMeanEst[,cluster], size=vecDispEst, log=FALSE) +
        vecDropoutEst[vecCells]
    ), na.rm=TRUE)
    # Evaluate on all nonzero count genes in cells
    scaLogLikNonzerosCluster <- sum(log(
      (1-vecDropoutEst[vecCells]) * dnbinom(matYNonzeros[,vecCells], 
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
#' @param vecDropoutEst: (vector cells) Dropout estimate of cell. 
#' @param lsResultsClustering:
#' 
#' @return scaLogLik: (scalar) Value of cost function (likelihood)
#'    of hurdle model for gene
#' @export

evalLogLikHurdleNB <- function(vecTheta,vecY,
  vecDropoutEst,lsResultsClustering){ 
  
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
    (1-vecDropoutEst[indZeros]) * dnbinom(vecY[indZeros],
      mu=exp(vecMeanEst[lsResultsClustering$Assignments[indZeros]]), 
      size=scaDispEst[indZeros],
      log=FALSE) + vecDropoutEst[indZeros]
  ), na.rm=TRUE)
  # Evaluate on all nonzero count genes in cells
  scaLogLikNonzeros <- sum(log(
    (1-vecDropoutEst[indNonzeros]) * dnbinom(vecY[indNonzeros],
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

evalLogLikHurdleDrop <- function(scaTheta,vecY,
  vecMeanEst,vecDispEst){ 
  
  scaDropoutEst <- scaTheta
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
    (1-scaDropoutEst) * dnbinom(vecY[indZeros],
      mu=exp(vecMeanEst[indZeros]), 
      size=vecDispEst[indZeros],
      log=FALSE) + scaDropoutEst
  ), na.rm=TRUE)
  # Evaluate on all nonzero count genes in cells
  scaLogLikNonzeros <- sum(log(
    (1-scaDropoutEst) * dnbinom(vecY[indNonzeros], 
      mu=exp(vecMeanEst[indNonzeros]),
      size=vecDispEst[indNonzeros],
      log=FALSE)
  ), na.rm=TRUE)
  
  # Get likelihood of all data in cluster
  scaLogLik <- scaLogLikZeros + scaLogLikNonzeros
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}
