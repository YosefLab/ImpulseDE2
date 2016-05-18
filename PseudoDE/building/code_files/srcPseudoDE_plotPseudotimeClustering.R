#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++     Plot ZINB fits    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Plots the zero inflated negative binomial fits and data to pdf
#' 
#' Plots the zero inflated negative binomial fits and data to pdf
#' 
#' @seealso Called by \code{runPseudoDE}.
#' 
#' @param lsGeneIDs (character vector) of gene names to be plotted; must be
#     trownames of arr3DCountData.
#' @param arr2DCountData (2D array genes x replicates) Count data: Reduced 
#'    version of \code{matCountData}. For internal use.
#' @param vecNormConst: (numeric vector number of replicates) 
#'    Normalisation constants for each replicate.
#' @param dfAnnotation (Table) Lists co-variables of individual replicates:
#'    Replicate, Sample, Condition, Time. Time must be numeric.
#' @param lsImpulseFits (list length 2 or 6) List of matrices which
#'    contain parameter fits and model values for given time course for the
#'    case condition (and control and combined if control is present).
#'    Each parameter matrix is called parameter_'condition' and has the form
#'    (genes x [beta, h0, h1, h2, t1, t2, logL_H1, converge_H1, mu, logL_H0, 
#'    ) where beta to t2 are parameters of the impulse
#'    model, mu is the single parameter of the mean model, logL are
#'    log likelihoods of full (H1) and reduced model (H0) respectively, converge
#'    is convergence status of numerical optimisation of model fitting by
#'    \code{optim} from \code{stats}. Each value matrix is called
#'    value_'condition' and has the form (genes x time points) and contains the
#'    counts predicted by the impulse model at the observed time points.
#' @param dfDEAnalysis (data frame genes x fitting characteristics) 
#'    Summary of fitting procedure for each gene.
#' @param vecRefPval (vec length genes) Method 2 (DESeq2) adjusted p-values
#' @param strCaseName (str) Name of the case condition in \code{dfAnnotationFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#'    \code{dfAnnotationFull}.
#' @param strMode: (str) [Default "batch"] {"batch","timecourses","singlecell"}
#'    Mode of model fitting.
#' @param NPARAM (scalar) [Default 6] Number of parameters of impulse model.
#' @param strFileNameSuffix (character string) [Default ""] File extention.
#' @param title_string (character string) [Default ""] Title for each plot.
#' @param strPlotSubtitle (character string) [Default ""] Subtitle for each plot.
#' @param NPARAM (scalar) [Default 6] Number of parameters of impulse model.
#' 
#' @return NULL
#' @export

plotPseudotimeClustering <- function(vecPseudotime, lsResultsClustering,
  strPDFname="PseudoDE_ZINBfitsPseudotimeClustering.pdf"){
  
  dfClusterBorders <- data.frame( borders=sapply( seq(1,length(lsResultsClustering$Centroids)-1), 
    function(centroid){(lsResultsClustering$Centroids[centroid]+lsResultsClustering$Centroids[centroid+1])/2} ))
  dfPseudotime <- data.frame( pseudotime=as.vector(vecPseudotime) )
  
  plotEDF <- ggplot() +
    geom_density(data=dfPseudotime, aes(x=pseudotime), colour="black", bw=1) +
    geom_vline(data=dfClusterBorders, aes(xintercept=borders), colour="green", linetype = "longdash") +
    labs(title="Density estimation of cells in pseudotime") +
    xlab("pseudotime") +
    ylab("empirical probability density")
  
  plotECDF <- ggplot() +
    stat_ecdf(data=dfPseudotime, aes(x=pseudotime), colour="red") +
    geom_vline(data=dfClusterBorders, aes(xintercept=borders), colour="green", linetype = "longdash") +
    labs(title="Empirical cumulative density of cells in pseudotime") +
    xlab("pseudotime") +
    ylab("empirical cumulative probability density")
    # Open .pdf
  pdf(strPDFname)
  print(plotEDF)
  print(plotECDF)
  # Close .pdf
  dev.off()
}