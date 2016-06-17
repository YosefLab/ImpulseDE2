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
  matMu <- matMuCluster[,lsResultsClustering$Assignments]
  matLikFull <- matDropout*(matCountsProc==0) + (1-matDropout)*
      dnbinom(matCountsProc, mu = matMu, size = matDispersions)
  vecLogLikFull <- apply(matLikFull, 1, function(gene){
    sum( log(gene[gene!=0]) +
        sum(gene==0)*log(.Machine$double.eps) )
  })
  matLikRed <- matDropoutNULL*(matCountsProc==0) + (1-matDropoutNULL)*
      dnbinom(matCountsProc, mu = matMuNULL, size = matDispersionsNULL)
  vecLogLikRed <- apply(matLikRed, 1, function(gene){
    sum( log(gene[gene!=0]) +
        sum(gene==0)*log(.Machine$double.eps) )
  })
  
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