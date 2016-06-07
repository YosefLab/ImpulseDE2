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
  strDropoutTraining="PoissonVar",
  vecHousekeepingGenes=NULL,
  vecSpikeInGenes=NULL,
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
    matDispersion <- matrix(NA,nrow=J,ncol=lsResultsClustering$K)
    matDropout <- matrix(NA,nrow=J,ncol=N)
    matProbNB <- matrix(NA,nrow=J,ncol=N)
    matDispersion <- matrix(NA,nrow=J,ncol=N)
    matCountsImputed <- matrix(NA,nrow=J,ncol=N)
    vecConvergence <- array(NA,lsResultsClustering$K)
    
    # 1. Identify dropout rates as Bernoulli model:
    # Dropout rates are modelled as a sample-specific
    # linear model of the gene-wise mean expression level.
    # This step is performed on the entire data set.

    # a) Define the set of genes on which dropout Bernoulli
    # model is trained for each cell.
    if(strDropoutTraining=="PoissonVar"){
      vecTargetGenes <- 
    } else if(strDropoutTraining=="Housekeeping"){
      vecTargetGenes <- vecHousekeepingGenes
    } else if(strDropoutTraining=="SpikeIns"){
      vecTargetGenes <- vecSpikeInGenes
    } else if(strDropoutTraining=="All"){
      vecTargetGenes <- as.vector(rownames(matCounts))
    } else {
      stop(paste0("ERROR: Did not recognise strDropoutTraining=",strDropoutTraining," in fitZINB()."))
    }

    # b) Identify zero-inflated Bernoulli model
    # using an EM-algorithm implementation from SCONE.
    lsZIBERparam <- estimate_ziber(
        x = matCounts,
        fp_thres = 0,
        bulk_model = TRUE,
        gfeatM = NULL,
        pos_controls = vecTargetGenes,
        maxiter = MAXITER, 
        verbose = TRUE)
    matDropout <- 1 - lsZIBERparam$p_nodrop
    vecConvergence <- lsZIBERparam$convergence

    # 2. Estimate dispersion factors:
    # This step is performed on each cluster separately:
    # Either as a GLM with cluster index as predictor,
    # cluster-wise mean parameters and gene-wise dispersion
    # parameter or K separate GLMs with one dispersion factor each.

    # a) Initialisation
    vecClusterAssign <- paste0(rep("cluster_",length(lsResultsClustering$Assignments)),lsResultsClustering$Assignments)
    vecindClusterAssign <- match(vecClusterAssign, unique(vecClusterAssign))
    scaNumGenes <- dim(matCounts)[1]
    # Initialise mean parameters
    lsMu0 <- lapply(vecClusters ,function(cluster){apply(matCounts[,vecClusterAssign==cluster],1,
      function(gene){mean(gene[gene>0])}
    )})
    matMu0 <- do.call(cbind, lsMu0)
    if(dim(matMu0)[2]>1){
      matMuLinearCoeff0 <- log(matMu0[,2:dim(matMu0)[2]]) - log(matMu0[,1])
      matMuCoeff <- cbind( log(matMu0[,1]), matMuLinearCoeff0 ) # Intersection, linear coefficients
    } else {
      matMuCoeff <- log(matMu0[,1])
    }
    # Initialise dispersion parameters
    vecTheta0 <- array(1,scaNumGenes)

    # b) GLM fitting
    fit_mu <- bplapply(seq(1,scaNumGenes), function(i) {
      fit <- glm.nb(matCounts[i,] ~ vecClusterAssign, weights = (1 - matDropout[i,]), init.theta = vecTheta0[i], start=matMuCoeff[i,])
      return(list(fitted=fit$fitted.values, theta=fit$theta))
    })
    # Extract negative binomial parameters into matrix (genes x clusters)
    matClusterMeansFitted <- do.call( rbind, lapply(fit_mu, function(x) x$fitted) )
    vecDispersions <- unlist( lapply(fit_mu, function(x) x$theta) )
    matDispersions <- matrix(thetahat, nrow=scaNumGenes, ncol=lsResultsClustering$K, byrow=FALSE)
    # Compute mixture probabilities and imputed counts
    matProbNB <- 1 - (matCounts == 0) * matDropout / (matDropout + (1 - matDropout) * (1 + muhat[,vecindClusterAssign] / matDispersions[,vecindClusterAssign])^(-matDispersions[,vecindClusterAssign]))
    matCountsImputed <- 
      matDropout * (1 - matProbNB) + 
      matClusterMeansFitted[,vecindClusterAssign] * matProbNB

    # Name objects
    rownames(matDropout) <- rownames(matCounts)
    colnames(matDropout) <- colnames(matCounts)
    rownames(matProbNB) <- rownames(matCounts)
    colnames(matProbNB) <- colnames(matCounts)
    rownames(matCountsImputed) <- rownames(matCounts)
    colnames(matCountsImputed) <- colnames(matCounts)
    rownames(vecDispersions) <- rownames(matCounts)
    rownames(matClusterMeansFitted) <- rownames(matCounts)
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