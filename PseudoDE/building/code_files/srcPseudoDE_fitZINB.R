#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++     Fit ZINB model    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Cost function zero-inflated negative binomial model for dispersion fitting
#' 
#' Log likelihood of zero inflated  negative binomial model. 
#' This function is designed to allow numerical optimisation
#' of negative binomial overdispersion on single gene given
#' the drop-out rate and negative binomial mean parameter.
#' 
#' @aliases evalLogLikHurdleNB_com
#' 
#' @seealso Called by \code{fitZINB}.
#' 
#' @param scaTheta: (scalar) Log of dispersion estimate.
#' @param vecY: (vector number of cells) Observed expression values 
#'    of gene in cells in cluster.
#' @param vecMuEst: (vector number of cells) Negative binomial
#'    mean parameter estimate of clusters to which cells
#'    belong.
#' @param vecDropoutRateEst: (vector number of cells) Dropout estimate of cell. 
#' @param vecboolObserved: (bool vector number of samples)
#'    Whether sample is not NA (observed).
#' @param vecboolZero: (bool vector number of samples)
#'    Whether sample has zero count.
#' 
#' @return scaLogLik: (scalar) Value of cost function:
#'    zero-inflated negative binomial likelihood.
#' @export

evalLogLikDispNB <- function(scaTheta,
  vecY,
  vecMuEst,
  vecDropoutRateEst,
  vecboolNotZeroObserved, 
  vecboolZero){ 
  
  scaDispEst <- exp(scaTheta)
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
      mu=vecMuEst[vecboolZero], 
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
      mu=vecMuEst[vecboolNotZeroObserved], 
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

#' Fit zero-inflated negative binomial model to data
#' 
#' Fit zero-inflated negative binomial model to data: One mean per cluster 
#' and either one dispersion parameter across all observations in all cluster
#' for a gene or one dispersion parameter per cluster per gene. Dropout rate
#' and dispersion factor inferred here are used as hyperparamters in the
#' impulse fitting stage of PseudoDE.
#' To parallelise this code, replace lapply by bplapply in dispersion
#' factor and drop-out rate estimation and uncomment the BiocParallel
#' register() command at the beginning of this function.
#' 
#' @seealso Called by \code{runPseudoDE}.
#' 
#' @param matCountsProc: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param lsResultsClustering (list {"Assignments","Centroids","K"})
#'    \itemize{
#'      \item   Assignments: (integer vector length number of
#'        cells) Index of cluster assigned to each cell.
#'      \item   Centroids: 1D Coordinates of cluster centroids,
#'        one scalar per centroid.
#'      \item   K: (scalar) Number of clusters selected.
#'      }
#' @param boolOneDispPerGene: (bool) [Default TRUE]
#'    Whether one negative binomial dispersion factor is fitted
#'    per gene or per gene for each cluster.
#' @param nProc: (scalar) [Default 1] Number of processes for 
#'    parallelisation.
#' @param scaMaxiterEM: (scalar) Maximum number of EM-iterations to
#'    be performed in ZINB model fitting.
#' @param verbose: (bool) Whether to follow EM-algorithm
#'    convergence.
#' @param boolSuperVerbose: (bool) Whether to follow EM-algorithm
#'    progress in high detail with local convergence flags. 
#' 
#' @return (list length 6)
#'    \itemize{
#'      \item vecDispersions: (numeric matrix genes x clusters)
#'        Inferred negative binomial dispersions.
#'      \item matDropout: (numeric matrix genes x cells)
#'        Inferred zero inflated negative binomial drop out rates.
#'      \item matProbNB: (numeric matrix genes x cells)
#'        Inferred probabilities of every observation to be generated
#'        from the negative binomial component of the zero-inflated
#'        negative binomial mixture model.
#'      \item matCountsProcImputed: (numeric matrix genes x cells)
#'        Data predicted with inferred zero-inflated negative binomial
#'        model.
#'      \item matMuCluster: (numeric matrix genes x clusters)
#'        Inferred negative binomial cluster means.
#'      \item boolConvergence: (bool) Convergence of EM algorithm.
#'    }
#' @export

fitZINB <- function(matCountsProc, 
  lsResultsClustering,
  nProc=1,
  strDropoutTrainingSet="PoissonVar",
  vecHousekeepingGenes=NULL,
  vecSpikeInGenes=NULL,
  boolOneDispPerGene=TRUE,
  scaMaxiterEM=20,
  verbose=FALSE,
  boolSuperVerbose=FALSE ){
  
  # Parameters
  # Minimim fractional liklihood increment necessary to
  # continue EM-iterations.
  scaPrecEM <- 1-10^(-4)
  
  # Set number of processes to be used for parallelisation
  # This function is currently not parallelised to reduce memory usage.
  # Read function description for further information.
  # register(MulticoreParam(nProc))
  
  vecClusterAssign <- paste0(rep("cluster_",length(lsResultsClustering$Assignments)),lsResultsClustering$Assignments)
  vecClusters <- unique(vecClusterAssign)
  vecindClusterAssign <- match(vecClusterAssign, vecClusters)
  scaNumGenes <- dim(matCountsProc)[1]
  scaNumCells <- dim(matCountsProc)[2]  
  
  # (I) Initialisation
  # The initial model is estimated under the assumption that all zero-counts
  # are drop-outs.
  matMuCluster <- array(NA,c(scaNumGenes,lsResultsClustering$K))
  lsPiCoefs <- NULL
  
  # E-step:
  # Posterior of dropout: matZ
  if(boolSuperVerbose){print("Initialisation E-step: Estimtate posterior of dropout")}
  matZ <- t(apply(matCountsProc==0, 1, as.numeric))
  
  # (II) EM itertion
  # Allow two iteration to get reference likelihood
  scaIter <- 1
  scaLogLikNew <- 0
  scaLogLikOld <- 0
  while(scaIter == 1 | scaIter == 2 | (scaLogLikNew > scaLogLikOld*scaPrecEM & scaIter <= scaMaxiterEM)){
    # M-step:
    # a) Negative binomial mean parameter
    # Use MLE of mean parameter of negative binomial distribution: weighted average
    if(boolSuperVerbose){print("M-step: Estimtate negative binomial mean parameters")}
    for(k in seq(1,max(vecindClusterAssign))){
      matMuCluster[,k] <- rowSums(
        matCountsProc[,vecindClusterAssign==k]*(1-matZ)[,vecindClusterAssign==k],
        na.rm=TRUE) / rowSums((1-matZ)[,vecindClusterAssign==k])
    }
    # Add pseudocounts to zeros
    matMuCluster[matMuCluster==0 | is.na(matMuCluster)] <- 1/scaNumCells
    matMu <- matMuCluster[,vecindClusterAssign]
    
    # b) Dropout rate
    # Fit dropout rate with GLM
    if(boolSuperVerbose){print("M-step: Estimtate dropout rate")}
    vecPiFit <- lapply(seq(1,scaNumCells), function(j) {
      glm( matZ[,j] ~ log(matMu[,j]),
        family=binomial(link=logit),
        control=list(maxit=1000)
      )[c("converged","fitted.values")]
    })
    vecboolConvergedGLMpi <- sapply(vecPiFit, function(x) x[[1]])
    matDropout <- sapply(vecPiFit, function(x) x[[2]])
    if(boolSuperVerbose){
      print(paste0("GLM to estimate drop-out rate for cells did not converge in ", 
        sum(!vecboolConvergedGLMpi), " cases."))
    }
    
    # c) Negative binomial dispersion parameter
    # Use MLE of dispersion factor: numeric optimisation of likelihood
    if(boolSuperVerbose){print("M-step: Estimtate negative binomial dispersion parameters")}
    matboolNotZeroObserved <- matCountsProc > 0 & !is.na(matCountsProc) & is.finite(matCountsProc)
    matboolZero <- matCountsProc == 0
    scaDispGues <- 1
    if(boolOneDispPerGene){  
      vecDispFit <- lapply(seq(1,scaNumGenes), function(i){
        unlist(optim(
          par=log(scaDispGues),
          fn=evalLogLikDispNB,
          vecY=matCountsProc[i,],
          vecMuEst=matMu[i,],
          vecDropoutRateEst=matDropout[i,],
          vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
          vecboolZero=matboolZero[i,],
          method="BFGS",
          control=list(maxit=1000, fnscale=-1)
        )[c("par","convergence")] )
      })
    } else {
      #  not coded
    }
    vecboolConvergedGLMdisp <- sapply(vecDispFit, function(x) x["convergence"])
    vecDispersions <- sapply(vecDispFit, function(fit){ exp(fit["par"]) })
    matDispersions <- matrix(vecDispersions, nrow=length(vecDispersions), ncol=scaNumCells, byrow=FALSE)
    matDispersions[matDispersions<=0] <- min(matDispersions[matDispersions>0])
    if(boolSuperVerbose){
      print(paste0("GLM to estimate drop-out rate for cells did not converge in ", 
        sum(vecboolConvergedGLMdisp), " cases."))
    }
    
    # E-step:
    if(boolSuperVerbose){print("E-step: Estimtate posterior of dropout")}
    matZ <- matDropout/(matDropout + (1-matDropout)*
        dnbinom(0, mu = matMu, size = matDispersions) )
    matZ[matCountsProc > 0] <- 0
    
    # Evaluate Likelihood
    scaLogLikOld <- scaLogLikNew
    matLikNew <- matDropout*(matCountsProc==0) + (1-matDropout)*
      dnbinom(matCountsProc, mu = matMu, size = matDispersions)
    scaLogLikNew <- sum( log(matLikNew[matLikNew!=0]) +
        sum(matLikNew==0)*log(.Machine$double.eps) )
    if(verbose){print(paste0("Completed Iteration ", scaIter, " with data log likelihood of ", scaLogLikNew))}
    scaIter <- scaIter+1
  }
  # Evaluate convergence
  if(all(vecboolConvergedGLMdisp) & all(vecboolConvergedGLMpi) &
      scaLogLikNew < scaLogLikOld*scaPrecEM & scaLogLikNew > scaLogLikOld){
    boolConvergence <- TRUE
  } else {
    boolConvergence <- FALSE
  }
  
  # Compute mixture probabilities and imputed counts
  matProbNB <- 1 - matZ
  matCountsProcImputed <- matDropout * (1 - matProbNB) + matMu * matProbNB
  
  # Name rows and columns of output
  rownames(matDropout) <- rownames(matCountsProc)
  colnames(matDropout) <- colnames(matCountsProc)
  rownames(matMu) <- rownames(matCountsProc)
  colnames(matMu) <- colnames(matCountsProc)
  rownames(matMuCluster) <- rownames(matCountsProc)
  names(vecDispersions) <- rownames(matCountsProc)
  rownames(matDispersions) <- rownames(matCountsProc)
  colnames(matDispersions) <- colnames(matCountsProc)
  rownames(matProbNB) <- rownames(matCountsProc)
  colnames(matProbNB) <- colnames(matCountsProc)
  rownames(matCountsProcImputed) <- rownames(matCountsProc)
  colnames(matCountsProcImputed) <- colnames(matCountsProc)
  
  # Check dispersions
  if(any(is.na(matDispersions) | !is.finite(matDispersions))){
    matDispersions[is.na(matDispersions) | !is.finite(matDispersions)] <- 1
    print("WARNING: Found NA/inf dispersions")
    print(matDispersions)
  }
  
  return(list( vecDispersions=vecDispersions,
    matDropout=matDropout, 
    matProbNB=matProbNB, 
    matCountsProcImputed=matCountsProcImputed, 
    matMuCluster=matMuCluster,
    boolConvergence=boolConvergence ))
}