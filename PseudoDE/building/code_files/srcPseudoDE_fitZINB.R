#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++     Fit ZINB model    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

selectPoissonGenes <- function(matCounts){
  # Tolerance level: Factor by which variance may be
  # larger than mean to be considered Poisson.
  scaPoissonTol <- 1.1
  
  # Only non-zero counts are tested for whether they are Poisson distributed.
  vecMu <- apply(matCounts, 1, function(gene){mean(gene[gene>0], na.rm=TRUE)})
  vecStdv <- apply(matCounts, 1, function(gene){sd(gene[gene>0], na.rm=TRUE)})
  vecPoissonGenes <- as.vector( (rownames(matCounts))[vecStdv^2 <= scaPoissonTol*vecMu] )
  # Print distribution of gene-wise coefficients of
  # variation to file.
  pdf("Genewise_CV",height=6.0,width=9.0)
  
  dev.off()
  # Handle case of few near Poisson-distributed genes.
  if(length(vecPoissonGenes) < 50){
    warning(paste0("WARNING: Found few (", length(vecPoissonGenes),") Poisson distributed genes in dataset.",
      "Consider identifying drop-out rate through different gene set (strDropoutTrainingSet)."))
  }
  if(length(vecPoissonGenes) < 10){
    stop(paste0("ERROR: Found too few (", length(vecPoissonGenes),") Poisson distributed genes in dataset.",
      "Switch to different gene set for identification of drop-out rates (strDropoutTrainingSet)."))
  }
  return(vecPoissonGenes)
}

#' Fit zero-inflated negative binomial model to data
#' 
#' Fit zero-inflated negative binomial model to data: One mean per cluster 
#' and either one dispersion parameter across all observations in all cluster
#' for a gene or one dispersion parameter per cluster per gene. Dropout rate
#' and dispersion factor inferred here are used as hyperparamters in the
#'impulse fitting stage of PseudoDE.
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

library(scone)
#source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/srcSCONE_zinb.R")

fitZINB <- function(matCounts, 
  lsResultsClustering, 
  nProc=1,
  strDropoutTrainingSet="PoissonVar",
  vecHousekeepingGenes=NULL,
  vecSpikeInGenes=NULL,
  boolOneDispPerGene=TRUE,
  MAXITER=20){
  
  # Set number of processes to be used in this step:
  # register(MulticoreParam()) controls the number of processes used for 
  # BiocParallel, used in zinb().
  nProcesses <- min(detectCores() - 1, nProc)
  print(paste0("Number of processes: ", nProcesses))
  register(MulticoreParam(nProcesses))
  
  # 1. Identify dropout rates as Bernoulli model:
  # Dropout rates are modelled as a sample-specific
  # linear model of the gene-wise mean expression level.
  # This step is performed on the entire data set.
  
  # a) Define the set of genes on which dropout Bernoulli
  # model is trained for each cell.
  if(strDropoutTrainingSet=="PoissonVar"){
    # Training genes are genes with count variance
    # close to Poissonian variance (stdv=mean).
    # Poissonian variance can be interpreted as noise
    # entirely derived from sampling, i.e. no biological noise.
    # Poissonian noise is therefore an identifier for
    # constantly expressed "housekeeping" genes (see below).
    vecTargetGenes <- selectPoissonGenes(matCounts)
  } else if(strDropoutTrainingSet=="Housekeeping"){
    # Training genes are externally supplied "housekeeping"
    # which are thought to have constant expression levels
    # across genes.
    vecTargetGenes <- vecHousekeepingGenes
  } else if(strDropoutTrainingSet=="SpikeIns"){
    # Training genes are the external RNA spike-ins,
    # which were given at constant concentration to each cell.
    vecTargetGenes <- vecSpikeInGenes
  } else if(strDropoutTrainingSet=="All"){
    # All genes are training genes. Note that the estimation
    # of the drop-out rate is impeded by overdispersed 
    # distribution with high density at zero.
    vecTargetGenes <- as.vector(rownames(matCounts))
  } else {
    stop(paste0("ERROR: Did not recognise strDropoutTrainingSet=",strDropoutTrainingSet," in fitZINB()."))
  }
  
  # b) Identify zero-inflated Bernoulli model
  # using an EM-algorithm implementation from SCONE.
  lsZIBERparam <- estimate_ziber(
    x = matCounts,
    fp_thres = 0,
    bulk_model = strDropoutTrainingSet!="All",
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
  vecClusterAssign <- paste0(rep("cluster_",length(lsResultsClustering$Assignments)),lsResultsClustering$Assignments)
  vecClusters <- unique(vecClusterAssign)
  vecindClusterAssign <- match(vecClusterAssign, vecClusters)
  scaNumGenes <- dim(matCounts)[1]
  
  # b) GLM fitting
  # Initialise dispersion parameters
  vecTheta0 <- array(1,scaNumGenes)
  
  if(boolOneDispPerGene){
    # b1) Fit negative binomial with one mean per cluster and one
    # dispersion factor per gene (i.e. global for all clusters).
    
    # Initialise mean parameters as log weighted cluster mean.
    # Note that the coefficients of the linear model matMuCoeff
    # have the following columns: intersection (mean of cluster 1),
    # (k in 2:K elements with) difference of mean of cluster 1 to cluster k.
    # The structure of this linear model is required by the model 
    # formulation in glm.nb: Cluster assignment is given as 
    # a categorial variable.
    matMu0 <- do.call(cbind, lapply(vecClusters, function(cluster){
      apply( (matCounts*(1-matDropout))[,vecClusterAssign==cluster], 1,
        function(gene){ mean(gene[gene>0], na.rm=TRUE) 
        })
    }))
    
    if(dim(matMu0)[2]>1){
      matMuLinearCoeff0 <- log(matMu0[,2:dim(matMu0)[2]]) - log(matMu0[,1])
      matMuCoeff <- cbind( log(matMu0[,1]), matMuLinearCoeff0 ) # Intersection, linear coefficients
    } else {
      matMuCoeff <- log(matMu0[,1])
    }
    
    # Fit GLM
    lsGLMNB <- bplapply(seq(1,scaNumGenes), function(i) {
      fit <- glm.nb( matCounts[i,] ~ vecClusterAssign, 
        weights = (1 - matDropout[i,]), 
        init.theta = vecTheta0[i], 
        start=matMuCoeff[i,], 
        link=log )
      return(list(mu_estimate=fit$fitted.values, theta_estimate=fit$theta))
    })
    # Extract negative binomial parameters into matrix (genes x clusters)
    matMu <- do.call( rbind, lapply(lsGLMNB, function(x) x$mu_estimate) )
    vecDispersions <- unlist( lapply(lsGLMNB, function(x) x$theta_estimate) )
    matDispersions <- matrix(thetahat, nrow=scaNumGenes, ncol=lsResultsClustering$K, byrow=FALSE)
    
  } else {
    # b2) Fit negative binomial with one mean per cluster and one
    # dispersion factor per cluster: The negative binomial models
    # for the clusters of one gene are entirely independent.
    matMu <- matrix(NA, nrow=scaNumGenes, ncol=lsResultsClustering$K)
    matDispersion <- matrix(NA, nrow=scaNumGenes, ncol=lsResultsClustering$K)
    
    # Initialise mean parameters as log weighted cluster mean:
    # Note that coefficients are the weighted cluster means
    # because this model set does not have categorial predictors
    # (cluster assignment) but only intersections. Each column
    # in matMuCoeff represents the intersection of a different linear
    # model, one linear model for each cluster.
    matMu0 <- do.call(cbind, lapply(vecClusters, function(cluster){
      apply( (matCounts*(1-matDropout))[,vecClusterAssign==cluster], 1,
        function(gene){ mean(gene[gene>0], na.rm=TRUE) 
        })
    }))
    matMuCoeff <- log(matMu)
    
    # Fit GLM
    for(cluster in vecClusters){
      vecidxCells <- vecClusterAssign==cluster
      idxCluster <- match(cluster, vecClusters)
      lsGLMNB <- bplapply(seq(1,scaNumGenes), function(i) {
        fit <- glm.nb( matCounts[i,vecidxCells] ~ 1, 
          weights = (1 - matDropout[i,vecidxCells]), 
          init.theta = vecTheta0[i], 
          start=matMuCoeff[i,idxCluster], 
          link=log )
        return(list(mu_estimate=fit$fitted.values, theta_estimate=fit$theta))
      })
      matMu[,idxCluster] <- unlist( lapply(lsGLMNB, function(x) x$mu_estimate) )
      matDispersion[,idxCluster] <- unlist( lapply(lsGLMNB, function(x) x$theta_estimate) )
    }
  }
  # Expand matrices from (genes x clusters)
  # to (genes x cells) by column duplication.
  # Now each element represents one observation.
  matMu <- matMu[,vecindClusterAssign]
  matDispersions <- matDispersions[,vecindClusterAssign]
  
  # Compute mixture probabilities and imputed counts
  matProbNB <- 1 - (matCounts == 0) * matDropout / (matDropout + (1 - matDropout) *
      (1 + matMu / matDispersions)^(-matDispersions))
  matCountsImputed <- matDropout * (1 - matProbNB) + matMu * matProbNB
  
  # Compute the overall log likelihood of the inferred zero-inflated
  # Bernoulli model. This model is the alternative model for
  # model-free differential expression analysis
  matboolZeros <- matCounts==0
  scaLogLik <- sum(log( matDropout*(matCounts==0) + (1-matDropout) * dnbinom(matCounts, mu=matMu, size=matDispersions, log=FALSE) ))
  #scalLogLikZeros <- sum(log( matDropout[matCounts==0] + (1 - matDropout[matboolZeros]) * dnbinom(matCounts[matboolZeros], mu=matMu[matboolZeros], size=matDispersions[matboolZeros], log=FALSE) ))
  #scalLogLikNonZeros <- sum(log( (1 - matDropout[!matboolZeros]) * dnbinom(matCounts[!matboolZeros], mu=matMu[!matboolZeros], size=matDispersions[!matboolZeros], log=FALSE) ))
  #scaLogLik <- scalLogLikZeros + scalLogLikNonZeros
  
  # Rows and columns of output matrix accoring to matCounts
  rownames(matDropout) <- rownames(matCounts)
  colnames(matDropout) <- colnames(matCounts)
  rownames(matMu) <- rownames(matCounts)
  colnames(matMu) <- colnames(matCounts)
  rownames(matDispersions) <- rownames(matCounts)
  colnames(matDispersions) <- colnames(matCounts)
  rownames(matProbNB) <- rownames(matCounts)
  colnames(matProbNB) <- colnames(matCounts)
  rownames(matCountsImputed) <- rownames(matCounts)
  colnames(matCountsImputed) <- colnames(matCounts)
  
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
  if(any(is.na(matDispersions) | !is.finite(matDispersions))){
    matDispersions[is.na(matDispersions) | !is.finite(matDispersions)] <- 1
    print("WARNING: Found NA/inf dispersions")
    print(matDispersions)
  }
  
  return(list( matDispersions=matDispersions,
    matDropout=matDropout, 
    matProbNB=matProbNB, 
    matCountsImputed=matCountsImputed, 
    matClusterMeansFitted=matClusterMeansFitted,
    vecConvergence=vecConvergence ))
}