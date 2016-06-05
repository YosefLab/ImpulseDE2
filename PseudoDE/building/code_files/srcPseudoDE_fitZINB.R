#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++     Plot ZINB fits    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Plots the zero inflated negative binomial fits and data to pdf.
#' 
#' Plots the zero inflated negative binomial fits and data to pdf.
#' 
#' @seealso Called by \code{runPseudoDE}.
#' 
#' @param matCounts: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param lsResultsClustering (list {"Assignments","Centroids","K"})
#'    \itemize{
#'      \item   Assignments: (integer vector length number of
#'        cells) Index of cluster assigned to each cell.
#'      \item   Centroids: 1D Coordinates of cluster centroids,
#'        one scalar per centroid.
#'      \item   K: (scalar) Number of clusters selected.
#'      }
#' @param nProc: (scalar) [Default 1] Number of processes for 
#'    parallelisation. The specified value is set
#'    to \code{min(detectCores() - 1, nProc)} using the 
#'    \code{detectCores} function from the package \code{parallel} 
#'    to avoid overload.
#' @param MAXITER: (scalar) Maximum number of EM-iterations to
#'    be performed in ZINB model fitting.
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
#'      \item matCountsImputed: (numeric matrix genes x cells)
#'        Data predicted with inferred zero-inflated negative binomial
#'        model.
#'      \item matClusterMeansFitted: (numeric matrix genes x clusters)
#'        Inferred negative binomial cluster means.
#'      \item vecConvergence: (numeric vector number of clusters)
#'        Convergence value for each zero-inflated negative binomial
#'        model fitting process.
#'    }
#' @export

fitZINB <- function(matCounts, 
  lsResultsClustering, 
  nProc=1,
  MAXITER=20){
  
  # Set number of processes to be used in this step:
  # register(MulticoreParam()) controls the number of processes used for 
  # BiocParallel, used in zinb().
  nProcesses <- min(detectCores() - 1, nProc)
  print(paste0("Number of processes: ", nProcesses))
  register(MulticoreParam(nProcesses))
  
  J <- dim(matCounts)[1]
  N <- dim(matCounts)[2]
  
  if(TRUE){
    # Fit to clusters seperately, use original verision of fitting function
    
    matDispersion <- matrix(NA,nrow=J,ncol=lsResultsClustering$K)
    matDropout <- matrix(NA,nrow=J,ncol=N)
    matProbNB <- matrix(NA,nrow=J,ncol=N)
    matDispersion <- matrix(NA,nrow=J,ncol=N)
    matCountsImputed <- matrix(NA,nrow=J,ncol=N)
    vecConvergence <- array(NA,lsResultsClustering$K)
    
    # Name objects
    rownames(matDropout) <- rownames(matCounts)
    colnames(matDropout) <- colnames(matCounts)
    rownames(matProbNB) <- rownames(matCounts)
    colnames(matProbNB) <- colnames(matCounts)
    rownames(matCountsImputed) <- rownames(matCounts)
    colnames(matCountsImputed) <- colnames(matCounts)
    rownames(vecDispersions) <- rownames(matCounts)
    rownames(matClusterMeansFitted) <- rownames(matCounts)
    
    lsZinbOutputByCluster <- list()
    for(k in 1:lsResultsClustering$K){
      print(paste0("Fitting cluster ",k))
      
      ### remove with new version?
      # Remove all zero genes: Mu is initialised
      # as log(sum counts)
      vecidxCluster <- lsResultsClustering$Assignments==k
      matCountsCluster <- matCounts[,vecidxCluster]
      # this doesnt work for 20 genes:
      #vecboolNonzeroGenes <- apply(matCountsCluster,1,
      #  function(gene){any(gene!=0)})
      # this works for 20 genes:
      vecboolNonzeroGenes <- apply(matCountsCluster,1,
        function(gene){mean(gene)>10})
      matCountsCluster <- matCountsCluster[vecboolNonzeroGenes,]
      
      # Fit zinb model
      lsZINBparam <- estimate_zinb(
        Y = matCountsCluster, 
        maxiter = MAXITER, 
        verbose = TRUE)
      lsZinbOutputByCluster[[k]] <- lsZINBparam
      
      # Record inferred parameters
      matClusterMeansFitted[vecboolNonzeroGenes,k] <- lsZINBparam$mu
      matDispersion[vecboolNonzeroGenes,k] <- lsZINBparam$theta
      matDropout[vecboolNonzeroGenes,vecidxCluster] <- lsZINBparam$pi
      matProbNB[vecboolNonzeroGenes,vecidxCluster] <- lsZINBparam$p_z
      # Impute count data matrix
      matCountsImputed[vecboolNonzeroGenes,vecidxCluster] <- 
        matCountsCluster * lsZINBparam$p_z + 
        lsZINBparam$mu * (1 - lsZINBparam$p_z)
      vecConvergence[k] <- lsZINBparam$converged
    }
  }
  
  if(FALSE){
    # Old code on rebuild estimate_zinb
    # this doesnt work for 20 genes:
    vecboolNonzeroGenes <- apply(matCounts,1,
      function(gene){ all(unlist(sapply( seq(1,lsResultsClustering$K), 
        function(cl){any(gene[lsResultsClustering$Assignments==cl]>10)} 
      ))) } )
    # this works for 20 genes:
    #vecboolNonzeroGenes <- apply(matCounts,1,
    #  function(gene){ all(unlist(sapply( seq(1,lsResultsClustering$K), 
    #    function(cl){any(gene[lsResultsClustering$Assignments==cl]>20)} 
    #  ))) } )
    matCountsClean <- matCounts[vecboolNonzeroGenes,]
    # Fit zinb model
    vecClusterAssign <- paste0(rep("cluster_",length(lsResultsClustering$Assignments)),lsResultsClustering$Assignments)
    lsZINBparam <- estimate_zinb(
      Y = matCountsClean,
      vecClusterAssign = vecClusterAssign,
      maxiter = MAXITER, 
      verbose = TRUE)
    lsZinbOutput <- lsZINBparam
    
    # Record parameters
    matDropout <- lsZINBparam$pi
    matProbNB <- lsZINBparam$p_z
    vecDispersions <- lsZINBparam$theta
    
    # Impute count data matrix
    vecClusters <- unique(vecClusterAssign)
    vecindClusterAssing <- match(vecClusterAssign, vecClusters)
    matCountsImputed <- 
      matCountsClean * lsZINBparam$p_z + 
      lsZINBparam$mu[,vecindClusterAssing] * (1 - lsZINBparam$p_z)
    
    # Name objects
    rownames(matDropout) <- rownames(matCountsClean)
    colnames(matDropout) <- colnames(matCountsClean)
    rownames(matProbNB) <- rownames(matCountsClean)
    colnames(matProbNB) <- colnames(matCountsClean)
    rownames(matCountsImputed) <- rownames(matCountsClean)
    colnames(matCountsImputed) <- colnames(matCountsClean)
    matCountsImputed <- round(matCountsImputed)
    names(vecDispersions) <- rownames(matCountsClean)
    matClusterMeansFitted <- (lsZINBparam$mu)[,match(seq(1,length(vecClusters)),lsResultsClustering$Assignments)]
    rownames(matClusterMeansFitted) <- rownames(matCountsClean)
  }
  
  # Check dispersions
  if(any(is.na(vecDispersions) | !is.finite(vecDispersions))){
    vecDispersions[is.na(vecDispersions) | !is.finite(vecDispersions)] <- 1
    print("WARNING: Found NA/inf dispersions")
    print(vecDispersions)
  }
  
  return(list( vecDispersions=vecDispersions,
    matDropout=matDropout, 
    matProbNB=matProbNB, 
    matCountsImputed=matCountsImputed, 
    matClusterMeansFitted=matClusterMeansFitted,
    vecConvergence=vecConvergence ))
}