#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++    Model-free DE analysis  ++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Model-free differential expression analysis
#' 
#' Performs model-free differential expression analysis. This functiion is
#' broken up into:
#' (I) Fit null model.
#' (II) Compute likelihood of full and reduced model.
#' (II) Differential expression analysis.
#' 
#' @seealso Called by \code{runPseudoDE}.
#' 
#' @param matCountsProc: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param vecPseudotime: (numerical vector length number of cells)
#'    Pseudotime coordinates (1D) of cells: One scalar per cell.
#'    Has to be named: Names of elements are cell names.
#' @param vecDispersions: (vector number of genes) Gene-wise 
#'    negative binomial dispersion coefficients.
#' @param vecDispersions: (numeric matrix genes x clusters)
#'    Inferred negative binomial dispersions.
#' @param matMuCluster: (numeric matrix genes x clusters)
#'    Inferred negative binomial cluster means.
#' @param matDropout: (numeric matrix genes x cells)
#'    Inferred zero inflated negative binomial drop out rates.
#' @param boolConvergenceZINBH1: (bool) Convergence of EM algorithm for
#'    full model.
#' @param nProc: (scalar) [Default 1] Number of processes for 
#'    parallelisation.
#' @param scaMaxiterEM: (scalar) Maximum number of EM-iterations to
#'    be performed in ZINB model fitting.
#' @param verbose: (bool) Whether to follow EM-algorithm
#'    convergence.
#' 
#' @return dfModelFreeDEAnalysis: (data frame genes) 
#'    Summary of model-free differential expression analysis.
#' @export

runModelFreeDEAnalysis <- function(matCountsProc,
  vecPseudotime,
  lsResultsClustering,
  vecDispersions,
  matMuCluster,
  matDropout,
  boolConvergenceZINBH1,
  boolOneDispPerGene=TRUE,
  nProc=1,
  scaMaxiterEM = 20,
  verbose=TRUE ){
  
  # (I) Fit null model
  # Null model is a single cluster: Create Clustering.
  lsNullClustering <- list()
  lsNullClustering[[1]] <- rep(1,dim(matCountsProc)[2])
  lsNullClustering[[2]] <- max(vecPseudotime)-min(vecPseudotime)
  lsNullClustering[[3]] <- 1
  names(lsNullClustering) <- c("Assignments","Centroids","K")
  
  # Fit zero-inflated negative binomial null model
  lsResZINBFitsNULL <- fitZINB( matCountsProc=matCountsProc, 
    lsResultsClustering=lsNullClustering,
    vecSpikeInGenes=NULL,
    boolOneDispPerGene=boolOneDispPerGene,
    nProc=nProc,
    scaMaxiterEM=scaMaxiterEM,
    verbose=verbose )
  vecDispersionsNULL <- lsResZINBFitsNULL$vecDispersions
  matDropoutNULL <- lsResZINBFitsNULL$matDropout
  matProbNBNULL  <- lsResZINBFitsNULL$matProbNB
  matCountsProcImputedNULL  <- lsResZINBFitsNULL$matCountsProcImputed
  vecMuClusterNULL  <- lsResZINBFitsNULL$matMuCluster
  boolConvergenceZINBH0 <- lsResZINBFitsNULL$boolConvergence
  matMuNULL <- matrix(vecMuClusterNULL, nrow=dim(matCountsProc)[1], 
    ncol=dim(matCountsProc)[2], byrow=FALSE)
  matDispersionsNULL <- matrix(vecDispersionsNULL, nrow=dim(matCountsProc)[1], 
    ncol=dim(matCountsProc)[2], byrow=FALSE)
  
  matDispersions <- matrix(vecDispersions, nrow=dim(matCountsProc)[1], 
    ncol=dim(matCountsProc)[2], byrow=FALSE)
  
  # (II) Compute log likelihoods
  matboolZeros <- matCountsProc==0
  matMu <- matMuCluster[,lsResultsClustering$Assignments]
  vecLogLikFull <- array(NA,dim(matCountsProc)[1])
  vecLogLikRed <- array(NA,dim(matCountsProc)[1])
  for(i in dim(matCountsProc)[1]){
    vecCounts <- matCountsProc[i,]
    vecboolZero <- vecCounts==0
    vecboolNotZeroObserved <- vecCounts>0 & !is.na(vecCounts)
    # Evaluate likelihood of null model
    # Likelihood of zero counts:
    vecLikH0Zeros <- (1-matDropoutNULL[i,vecboolZero])*
      dnbinom( vecCounts[vecboolZero], 
        mu=matMuNULL[i,vecboolZero], 
        size=matDispersionsNULL[i,vecboolZero], 
        log=FALSE) +
      matDropoutNULL[i,vecboolZero]
    # Replace zero likelihood observation with machine precision
    # for taking log.
    scaLogLikH0Zeros <- sum( log(vecLikH0Zeros[vecLikH0Zeros!=0]) +
        sum(vecLikH0Zeros==0)*log(.Machine$double.eps) )
    # Likelihood of non-zero counts:
    vecLikH0Nonzeros <- log(1-matDropoutNULL[i,vecboolNotZeroObserved]) + 
      dnbinom(vecCounts[vecboolNotZeroObserved], 
        mu=matMuNULL[i,vecboolNotZeroObserved], 
        size=matDispersionsNULL[i,vecboolNotZeroObserved], 
        log=TRUE)
    # Replace zero likelihood observation with machine precision
    # for taking log.
    scaLogLikH0Nonzeros <- sum( vecLikH0Nonzeros[is.finite(vecLikH0Nonzeros)]) +
      sum(!is.finite(vecLikH0Nonzeros))*log(.Machine$double.eps)
    # Compute likelihood of all data:
    vecLogLikRed[i] <- scaLogLikH0Zeros + scaLogLikH0Nonzeros
    
    # Evaluate likelihood of alternative model
    # Likelihood of zero counts:
    vecLikH1Zeros <- (1-matDropout[i,vecboolZero])*
      dnbinom( vecCounts[vecboolZero], 
        mu=matMu[i,vecboolZero], 
        size=matDispersions[i,vecboolZero], 
        log=FALSE) +
      matDropout[i,vecboolZero]
    # Replace zero likelihood observation with machine precision
    # for taking log.
    scaLogLikH1Zeros <- sum( log(vecLikH1Zeros[vecLikH1Zeros!=0]) +
        sum(vecLikH1Zeros==0)*log(.Machine$double.eps) )
    # Likelihood of non-zero counts:
    vecLikH1Nonzeros <- log(1-matDropout[i,vecboolNotZeroObserved]) + 
      dnbinom(vecCounts[vecboolNotZeroObserved], 
        mu=matMu[i,vecboolNotZeroObserved], 
        size=matDispersions[i,vecboolNotZeroObserved], 
        log=TRUE)
    # Replace zero likelihood observation with machine precision
    # for taking log.
    scaLogLikH1Nonzeros <- sum( vecLikH1Nonzeros[is.finite(vecLikH1Nonzeros)]) +
      sum(!is.finite(vecLikH1Nonzeros))*log(.Machine$double.eps)    
    # Compute likelihood of all data:
    vecLogLikFull[i] <- scaLogLikH1Zeros + scaLogLikH1Nonzeros
  }
  
  # (III) Differential expression analysis
  # scaK: Number of clusters used in full model
  scaK <- lsResultsClustering$K
  # K x mean and 1 or K x dispersion parameters
  if(boolOneDispPerGene){ scaDegFreedomFull <- scaK + 1
  }else{ scaDegFreedomFull <- scaK*2 }
  # One dispersion estimate and one overall mean estimate
  scaDegFreedomRed <- 2
  # Compute difference in degrees of freedom between null model and alternative model.
  scaDeltaDegFreedom <- scaDegFreedomFull - scaDegFreedomRed
  # Compute test statistic: Deviance
  vecDeviance <- 2*(vecLogLikFull - vecLogLikRed)
  # Get p-values from Chi-square distribution (assumption about null model)
  vecPvalue <- pchisq(vecDeviance,df=scaDeltaDegFreedom,lower.tail=FALSE)
  # Multiple testing correction (Benjamini-Hochberg)
  vecPvalueBH = p.adjust(vecPvalue, method = "BH")
  
  dfModelFreeDEAnalysis =   as.data.frame(cbind(
    "Gene" = row.names(matCountsProc),
    "p"=as.numeric(vecPvalue),
    "adj.p"=as.numeric(vecPvalueBH),
    "loglik_full"=vecLogLikFull,
    "loglik_red"=vecLogLikRed,
    "deviance"=vecDeviance,
    "mean"=vecMuClusterNULL,
    "dispersion_null"=vecDispersionsNULL,
    "converged_null"=rep(boolConvergenceZINBH0!=0,dim(matCountsProc)[1]),
    "converged_full"=rep(all(boolConvergenceZINBH1),dim(matCountsProc)[1]),
    stringsAsFactors = FALSE))
  
  # Order data frame by adjusted p-value
  dfModelFreeDEAnalysis$adj.p <- as.numeric(as.character(dfModelFreeDEAnalysis$adj.p))
  dfModelFreeDEAnalysis = dfModelFreeDEAnalysis[order(dfModelFreeDEAnalysis$adj.p),]
  
  return(dfModelFreeDEAnalysis)
}