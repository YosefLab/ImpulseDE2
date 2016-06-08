#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++    Model-free DE analysis  ++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Model-free differential expression analysis
#' 
#' Performs model-free differential expression analysis. This functiion is
#' broken up into:
#' (I) Fit null model
#' (II) Differential expression analysis
#' 
#' @seealso Called by \code{runPseudoDE}.
#' 
#' @param matCounts: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param vecDispersions: (vector number of genes) Gene-wise 
#'    negative binomial dispersion coefficients.
#' @param matClusterMeans:
#' @param vecLogLikFull: (numerical vector numer of genes)
#'    Log likelihood of full model which is a ZINB model for each cluster.
#' @param vecConvergenceFull: (numeric vector number of clusters)
#'    Convergence value for each zero-inflated negative binomial
#'    model fitting process (one for each cluster).
#' @param Q_value: (scalar) Q-value threshold for calling differentially
#'    expressed genes.
#' @param MAXITER: (scalar) Maximum number of EM-iterations to
#'    be performed in ZINB model fitting.
#' 
#' @return dfModelFreeDEAnalysis: (data frame genes) 
#'    Summary of model-free differential expression analysis.
#' @export

runModelFreeDEAnalysis <- function(matCounts,
  vecDispersions,
  matClusterMeans,
  scaLogLikFull,
  vecConvergenceFull,
  Q_value = 0.01,
  MAXITER = 20 ){
  
  # (I) Fit null model
  # Fit once to all data, one mean
  
  # Negative binomial over-dispersion factor: 1 per gene per cluster
  matDispersion <- array(NA,c(dim(matCounts)[1],1))
  
  # this doesnt work for 20 genes:
  vecboolNonzeroGenes <- apply(matCountsCluster,1,
    function(gene){any(gene!=0)})
  # this works for 20 genes:
  #vecboolNonzeroGenes <- apply(matCounts,1,
  #  function(gene){mean(gene)>10})
  matCountsClean <- matCounts[vecboolNonzeroGenes,]
  
  # Fit zinb model
  lsZINBparam <- estimate_zinb(
    Y = matCountsCluster, 
    maxiter = MAXITER, 
    verbose = TRUE)
  
  # Record parameters
  vecDispersions <- lsZINBparam$theta
  vecMeans <- lsZINBparam$mu
  vecLogLikRed <- lsZINBparam$loglikelihood
  scaConvergedRed <- lsZINBparam$converged
  
  # (II) Differential expression analysis
  # scaK: Number of clusters used in full model
  scaK <- dim(matClusterMeans)[2]
  # K mean and dispersion parameters
  scaDegFreedomFull <- K*2
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
    "Gene" = row.names(matCounts),
    "p"=as.numeric(vecPvalue),
    "adj.p"=as.numeric(vecPvalueBH),
    "loglik_full"=vecLogLikFull,
    "loglik_red"=vecLogLikRed,
    "deviance"=vecDeviance,
    "mean"=vecMeans,
    "dispersion_null"=vecDispersions,
    "converged_null"=rep(scaConvergedRed!=0,dim(matCounts)[1]),
    "converged_full"=rep(all(vecConvergenceFull),dim(matCounts)[1]),
    stringsAsFactors = FALSE))
  
  # Order data frame by adjusted p-value
  dfModelFreeDEAnalysis$adj.p <- as.numeric(as.character(dfDEAnalysis$adj.p))
  dfModelFreeDEAnalysis = dfModelFreeDEAnalysis[order(dfDEAnalysis$adj.p),]
  
  return(dfModelFreeDEAnalysis)
}