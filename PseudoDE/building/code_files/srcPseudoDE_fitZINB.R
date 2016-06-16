#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++     Fit ZINB model    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

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
#' @param vecTheta: (vector number of parameters) Parmeters of negative
#'    binomial in hurdle model of one gene to be estimated:
#'    Overdispersion and one mean per cluster.
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
#' @return scaLogLik: (scalar) Value of cost function (likelihood)
#'    of hurdle model for gene
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

selectPoissonGenes <- function(matCountsProc){
  # Tolerance level: Factor by which variance may be
  # larger than mean to be considered Poisson.
  scaPoissonTol <- 10
  
  # Only non-zero counts are tested for whether they are Poisson distributed.
  vecMu <- apply(matCountsProc, 1, function(gene){mean(gene[gene>=1], na.rm=TRUE)})
  vecStdv <- apply(matCountsProc, 1, function(gene){sd(gene[gene>=1], na.rm=TRUE)})
  vecPoissonGenes <- c(vecStdv^2 <= scaPoissonTol*vecMu)
  # Standard deviation is NA if only single count is non-zero:
  vecPoissonGenes[is.na(vecPoissonGenes)] <- FALSE
  
  vecCV <- vecStdv/vecMu
  dfScatter <- data.frame( logmu=log(vecMu)/log(2),
    logcv=log(vecCV)/log(2) )
  pdf("PseudoDE_Genewise_CV.pdf",height=6.0,width=9.0)
  #hist(vecStdv^2/vecMu)
  plotScatter <- ggplot(dfScatter, aes(x = logmu, y = logcv)) + 
    geom_point() +
    labs(title="Mean-variance relationship") +
    xlab("log_2 Mean") +
    ylab("log_2 CV")
  print(plotScatter)
  #scatter(x=log(vecMu)/log(2),y=log(vecStdv^2)/log(vecMu))
  dev.off()
  # Handle case of few near Poisson-distributed genes.
  scaNWarning <- 50
  scaNStop <- 10
  if(sum(vecPoissonGenes) < scaNWarning & sum(vecPoissonGenes) > scaNStop){
    warning(paste0("WARNING: Found few (", sum(vecPoissonGenes),") Poisson distributed genes in dataset.",
      "Consider identifying drop-out rate through different gene set (strDropoutTrainingSet)."))
  } else if(sum(vecPoissonGenes) <= scaNStop){
    stop(paste0("ERROR: Found too few (", sum(vecPoissonGenes),") Poisson distributed genes in dataset.",
      " Switch to different gene set for identification of drop-out rates (strDropoutTrainingSet)."))
  } else {
    print(paste0("Found ", sum(vecPoissonGenes)," Poisson distributed genes in dataset."))
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
#'      \item matCountsProcImputed: (numeric matrix genes x cells)
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

fitZINB <- function(matCountsProc, 
  lsResultsClustering,
  nProc=1,
  strDropoutTrainingSet="PoissonVar",
  vecHousekeepingGenes=NULL,
  vecSpikeInGenes=NULL,
  boolOneDispPerGene=TRUE,
  MAXITER=100,
  verbose=FALSE){
  
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
  
  if(FALSE){
    # a) Define the set of genes on which dropout Bernoulli
    # model is trained for each cell.
    if(strDropoutTrainingSet=="PoissonVar"){
      # Training genes are genes with count variance
      # close to Poissonian variance (stdv=mean).
      # Poissonian variance can be interpreted as noise
      # entirely derived from sampling, i.e. no biological noise.
      # Poissonian noise is therefore an identifier for
      # constantly expressed "housekeeping" genes (see below).
      vecTargetGenes <- selectPoissonGenes(matCountsProc)
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
      vecTargetGenes <- as.vector(rownames(matCountsProc))
    } else {
      stop(paste0("ERROR: Did not recognise strDropoutTrainingSet=",strDropoutTrainingSet," in fitZINB()."))
    }
    
    # b) Identify zero-inflated Bernoulli model
    # using an EM-algorithm implementation from SCONE.
    if(verbose){
      print("4.1 Fitting zero-inflated Bernoulli model for drop-out rates.")
      print("### SCONE::estimate_ziber() output starts here: #####################")
    }
    lsZIBERparam <- estimate_ziber(
      x = matCountsProc,
      fp_tresh = 0,
      bulk_model = strDropoutTrainingSet!="All",
      gfeatM = NULL,
      pos_controls = vecTargetGenes,
      maxiter = MAXITER, 
      verbose = TRUE)
    if(verbose){print("### SCONE::estimate_ziber() output ends here. ######################")}
    matDropout <- 1 - lsZIBERparam$p_nodrop
    vecConvergence <- lsZIBERparam$convergence
  }
  
  vecClusterAssign <- paste0(rep("cluster_",length(lsResultsClustering$Assignments)),lsResultsClustering$Assignments)
  vecClusters <- unique(vecClusterAssign)
  vecindClusterAssign <- match(vecClusterAssign, vecClusters)
  scaNumGenes <- dim(matCountsProc)[1]
  scaNumCells <- dim(matCountsProc)[2]
  
  # b) GLM fitting
  # Initialise dispersion parameters
  #vecTheta0 <- array(1,scaNumGenes)
  #if(verbose){print("4.2 Fitting negative model to drop-out weighted data with glm.nb().")}
  
  if(boolOneDispPerGene){
    # b1) Fit negative binomial with one mean per cluster and one
    # dispersion factor per gene (i.e. global for all clusters).
    
    # This is GLM
    #matMu0 <- do.call(cbind, lapply(vecClusters, function(cluster){
    #  apply( (matCountsProc*matWeights)[,vecClusterAssign==cluster], 1,
    #    function(gene){ mean(gene, na.rm=TRUE) 
    #    })
    #}))
    # Add pseudocounts for taking log
    #matMu0 <- matMu0+1
    #matMu0[is.na(matMu0)] <- 1
    
    #if(dim(matMu0)[2]>1){
    #  matMuLinearCoeff0 <- log(matMu0[,2:dim(matMu0)[2]]) - log(matMu0[,1])
    #  matMuCoeff0 <- cbind( log(matMu0[,1]), matMuLinearCoeff0 ) # Intersection, linear coefficients
    #} else {
    #  matMuCoeff0 <- log(matMu0[,1])
    #}
    
    # Fit GLM
    #lsGLMNB <- bplapply(seq(1,scaNumGenes), function(i) {
    #  fit <- glm.nb( matCountsProc[i,] ~ vecClusterAssign, 
    #    weights = (1-matDropout)[i,], 
    #    init.theta = vecTheta0[i], 
    #    start=matMuCoeff0[i,], 
    #    link=log )
    #  return(list(mu_estimate=fit$fitted.values, theta_estimate=fit$theta))
    #})
    # Extract negative binomial parameters into matrix (genes x clusters)
    #matMu <- do.call( rbind, lapply(lsGLMNB, function(x) x$mu_estimate) )
    #vecDispersions <- unlist( lapply(lsGLMNB, function(x) x$theta_estimate) )
    #matDispersions <- matrix(vecDispersions, nrow=length(vecDispersions), ncol=scaNumCells, byrow=FALSE)
    
    matMuCluster <- array(NA,c(scaNumGenes,lsResultsClustering$K))
    vecDispersions <- array(NA,c(scaNumGenes))
    
    scaIter <- 1
    scaLogLikNew <- 0
    scaLogLikOld <- 0
    scaPrec <- 1-10^(-4)
    scaMaxiterEM <- 20
    superverbose = FALSE
    
    # (I) Initialisation
    # The initial model is estimated under the assumption that all zero-counts
    # are drop-outs.
    
    # E-step:
    # Posterior of dropout: matZ
    if(superverbose){print("Initialisation E-step: Estimtate posterior of dropout")}
    matZ <- t(apply(matCountsProc==0, 1, as.numeric))
    
    # (II) EM itertion
    # Have to allow second iteration independent of likelihood as the first iteration
    # represents a reinitialisation based on a different model.
    while(scaIter == 1 | scaIter == 2 | (scaLogLikNew > scaLogLikOld*scaPrec & scaIter <= scaMaxiterEM)){
      # M-step:
      # a) Negative binomial mean parameter
      # Use MLE of mean parameter of negative binomial distribution: weighted average
      if(superverbose){print("M-step: Estimtate negative binomial mean parameters")}
      for(k in seq(1,max(vecindClusterAssign))){
        matMuCluster[,k] <- rowSums(
          matCountsProc[,vecindClusterAssign==k]*(1-matZ)[,vecindClusterAssign==k],
          na.rm=TRUE) / rowSums((1-matZ)[,vecindClusterAssign==k])
      }
      # Add pseudocounts to zeros
      matMuCluster[matMuCluster==0 | is.na(matMuCluster)] <- 1/scaNumCells
      matMu <- matMuCluster[,vecindClusterAssign]
      
      # b) Dropout rate
      if(superverbose){print("M-step: Estimtate dropout rate")}
      # Fit dropout rate with GLM
      fit_pi <- lapply(seq(1,scaNumCells), function(j) {
        fit <- glm( matZ[,j] ~ log(matMu[,j]), family=binomial(link=logit))
        return(list(coefs=coefficients(fit), fitted=fitted.values(fit)))
      })
      coefs_pi <- sapply(fit_pi, function(x) x$coefs)
      matDropout <- sapply(fit_pi, function(x) x$fitted)

      # c) Negative binomial dispersion parameter
      # Use MLE of dispersion factor: numeric optimisation of likelihood
      if(superverbose){print("M-step: Estimtate negative binomial dispersion parameters")}
      for(i in seq(1,scaNumGenes)){
        if(FALSE){
          print(matCountsProc[i,matCountsProc[i,]>0])
          print(unique(matMu[i,]))
          print(matDropout[i,])
        }
        vecboolNotZeroObserved <- matCountsProc[i,] > 0 & !is.na(matCountsProc[i,]) & is.finite(matCountsProc[i,])
        vecboolZero <- matCountsProc[i,] == 0
        vecFit <- unlist( optim(
          par=log(1),
          fn=evalLogLikDispNB,
          vecY=matCountsProc[i,],
          vecMuEst=matMu[i,],
          vecDropoutRateEst=matDropout[i,],
          vecboolNotZeroObserved=vecboolNotZeroObserved, 
          vecboolZero=vecboolZero,
          method="BFGS",
          control=list(maxit=1000, fnscale=-1)
        )[c("par","convergence")] )
        if(superverbose & vecFit["convergence"]){
          print(paste0(i, " convergence: ", vecFit["convergence"] ))
        }
        vecDispersions[i] <- exp(vecFit["par"])
      }
      matDispersions <- matrix(vecDispersions, nrow=length(vecDispersions), ncol=scaNumCells, byrow=FALSE)
      matDispersions[matDispersions<=0] <- min(matDispersions[matDispersions>0])
      
      # E-step:
      if(superverbose){print("E-step: Estimtate posterior of dropout")}
      matZ <- matDropout/(matDropout + (1-matDropout)*
          dnbinom(0, mu = matMu, size = matDispersions) )
      matZ[matCountsProc > 0] <- 0
      
      # Evaluate Likelihood
      scaLogLikOld <- scaLogLikNew
      matLikNew <- matDropout*(matCountsProc==0) + (1-matDropout)*
          dnbinom(matCountsProc, mu = matMu, size = matDispersions)
      scaLogLikNew <- sum( log(matLikNew[matLikNew!=0]) +
          sum(matLikNew==0)*log(.Machine$double.eps) )
      #print(sum(matLikNew==0))
      if(verbose){print(paste0("Completed Iteration ", scaIter, " with data log likelihood of ", scaLogLikNew))}
      scaIter <- scaIter+1
    }
    
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
      apply( (matCountsProc*(1-matDropout))[,vecClusterAssign==cluster], 1,
        function(gene){ mean(gene[gene>0], na.rm=TRUE) 
        })
    }))
    matMu0[is.na(matMu0) | matMu0==0] <- 0.01
    matMuCoeff0 <- log(matMu0)
    
    # Fit GLM
    for(cluster in vecClusters){
      vecidxCells <- vecClusterAssign==cluster
      idxCluster <- match(cluster, vecClusters)
      lsGLMNB <- bplapply(seq(1,scaNumGenes), function(i) {
        fit <- glm.nb( matCountsProc[i,vecidxCells] ~ 1, 
          weights = (1 - matDropout[i,vecidxCells]), 
          init.theta = vecTheta0[i], 
          start=matMuCoeff0[i,idxCluster], 
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
  #matMu <- matMu[,vecindClusterAssign]
  #matDispersions <- matDispersions[,vecindClusterAssign]
  
  # Compute mixture probabilities and imputed counts
  matProbNB <- 1 - matZ
  matCountsProcImputed <- matDropout * (1 - matProbNB) + matMu * matProbNB
  
  # Compute the overall log likelihood of the inferred zero-inflated
  # Bernoulli model. This model is the alternative model for
  # model-free differential expression analysis
  matboolZeros <- matCountsProc==0
  scaLogLik <- sum(log( matDropout*(matCountsProc==0) + (1-matDropout) * dnbinom(matCountsProc, mu=matMu, size=matDispersions, log=FALSE) ))
  #scalLogLikZeros <- sum(log( matDropout[matCountsProc==0] + (1 - matDropout[matboolZeros]) * dnbinom(matCountsProc[matboolZeros], mu=matMu[matboolZeros], size=matDispersions[matboolZeros], log=FALSE) ))
  #scalLogLikNonZeros <- sum(log( (1 - matDropout[!matboolZeros]) * dnbinom(matCountsProc[!matboolZeros], mu=matMu[!matboolZeros], size=matDispersions[!matboolZeros], log=FALSE) ))
  #scaLogLik <- scalLogLikZeros + scalLogLikNonZeros
  
  # Rows and columns of output matrix accoring to matCountsProc
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
  
  if(FALSE){
    # Old code on rebuild estimate_zinb
    # this doesnt work for 20 genes:
    vecboolNonzeroGenes <- apply(matCountsProc,1,
      function(gene){ all(unlist(sapply( seq(1,lsResultsClustering$K), 
        function(cl){any(gene[lsResultsClustering$Assignments==cl]>10)} 
      ))) } )
    # this works for 20 genes:
    #vecboolNonzeroGenes <- apply(matCountsProc,1,
    #  function(gene){ all(unlist(sapply( seq(1,lsResultsClustering$K), 
    #    function(cl){any(gene[lsResultsClustering$Assignments==cl]>20)} 
    #  ))) } )
    matCountsProcClean <- matCountsProc[vecboolNonzeroGenes,]
    # Fit zinb model
    vecClusterAssign <- paste0(rep("cluster_",length(lsResultsClustering$Assignments)),lsResultsClustering$Assignments)
    lsZINBparam <- estimate_zinb(
      Y = matCountsProcClean,
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
    matCountsProcImputed <- 
      matCountsProcClean * lsZINBparam$p_z + 
      lsZINBparam$mu[,vecindClusterAssing] * (1 - lsZINBparam$p_z)
    
    # Name objects
    rownames(matDropout) <- rownames(matCountsProcClean)
    colnames(matDropout) <- colnames(matCountsProcClean)
    rownames(matProbNB) <- rownames(matCountsProcClean)
    colnames(matProbNB) <- colnames(matCountsProcClean)
    rownames(matCountsProcImputed) <- rownames(matCountsProcClean)
    colnames(matCountsProcImputed) <- colnames(matCountsProcClean)
    matCountsProcImputed <- round(matCountsProcImputed)
    names(vecDispersions) <- rownames(matCountsProcClean)
    matClusterMeansFitted <- (lsZINBparam$mu)[,match(seq(1,length(vecClusters)),lsResultsClustering$Assignments)]
    rownames(matClusterMeansFitted) <- rownames(matCountsProcClean)
  }
  
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
    matMuCluster=matMuCluster))
    #vecConvergence=vecConvergence ))
}