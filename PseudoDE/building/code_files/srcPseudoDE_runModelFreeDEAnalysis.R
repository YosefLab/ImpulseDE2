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
  lsResultsClusteringH1,
  vecDispersionsH1,
  matMuClusterH1,
  matDropoutH1,
  boolConvergenceZINBH1,
  boolOneDispPerGene=TRUE,
  nProc=1,
  scaMaxiterEM = 20,
  verbose=TRUE ){
  
  # (I) Fit null model
  # Null model is a single cluster: Create Clustering.
  lsH0Clustering <- list()
  lsH0Clustering[[1]] <- rep(1,dim(matCountsProc)[2])
  lsH0Clustering[[2]] <- max(vecPseudotime)-min(vecPseudotime)
  lsH0Clustering[[3]] <- 1
  names(lsH0Clustering) <- c("Assignments","Centroids","K")
  
  # Fit zero-inflated negative binomial null model
  lsResZINBFitsH0 <- fitZINB( matCountsProc=matCountsProc, 
    lsResultsClustering=lsH0Clustering,
    vecSpikeInGenes=NULL,
    boolOneDispPerGene=boolOneDispPerGene,
    nProc=nProc,
    scaMaxiterEM=scaMaxiterEM,
    verbose=verbose )
  vecDispersionsH0 <- lsResZINBFitsH0$vecDispersions
  matDropoutH0 <- lsResZINBFitsH0$matDropout
  matProbNBH0  <- lsResZINBFitsH0$matProbNB
  vecMuClusterH0  <- lsResZINBFitsH0$matMuCluster
  boolConvergenceZINBH0 <- lsResZINBFitsH0$boolConvergence
  matMuH0 <- matrix(vecMuClusterH0, nrow=dim(matCountsProc)[1], 
    ncol=dim(matCountsProc)[2], byrow=FALSE)
  matDispersionsH0 <- matrix(vecDispersionsH0, nrow=dim(matCountsProc)[1], 
    ncol=dim(matCountsProc)[2], byrow=FALSE)
  
  matMuH1 <- matMuClusterH1[,lsResultsClusteringH1$Assignments]
  matDispersionsH1 <- matrix(vecDispersionsH1, nrow=dim(matCountsProc)[1], 
    ncol=dim(matCountsProc)[2], byrow=FALSE)
  
  # (II) Compute log likelihoods
  matboolNotZeroObserved <- matCountsProc >0 & !is.na(matCountsProc)
  matboolZero <- matCountsProc==0
  vecLogLikFull <- sapply( seq(1,dim(matCountsProc)[1]), function(i){
    evalLogLikZINB_PseudoDE_comp(vecY=matCountsProc[i,],
      vecMu=matMuH1[i,],
      vecDispEst=matDispersionsH1[i,], 
      vecDropoutRateEst=matDropoutH1[i,],
      vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
      vecboolZero=matboolZero[i,])
  })
  vecLogLikRed <- sapply( seq(1,dim(matCountsProc)[1]), function(i){
    evalLogLikZINB_PseudoDE_comp(vecY=matCountsProc[i,],
      vecMu=matMuH0[i,],
      vecDispEst=matDispersionsH0[i,], 
      vecDropoutRateEst=matDropoutH0[i,],
      vecboolNotZeroObserved=matboolNotZeroObserved[i,], 
      vecboolZero=matboolZero[i,])
  })
  
  # (III) Differential expression analysis
  # scaK: Number of clusters used in full model
  scaK <- lsResultsClusteringH1$K
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
    "mean"=vecMuClusterH0,
    "dispersion_H0"=vecDispersionsH0,
    "converged_H0"=rep(boolConvergenceZINBH0!=0,dim(matCountsProc)[1]),
    "converged_full"=rep(all(boolConvergenceZINBH1),dim(matCountsProc)[1]),
    stringsAsFactors = FALSE))
  
  # Order data frame by adjusted p-value
  dfModelFreeDEAnalysis$adj.p <- as.numeric(as.character(dfModelFreeDEAnalysis$adj.p))
  dfModelFreeDEAnalysis = dfModelFreeDEAnalysis[order(dfModelFreeDEAnalysis$adj.p),]
  
  return(dfModelFreeDEAnalysis)
}