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

#---------
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
#vecTheta0 <- array(1,scaNumGenes)
#if(verbose){print("4.2 Fitting negative model to drop-out weighted data with glm.nb().")}
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

#=====

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

# Expand matrices from (genes x clusters)
# to (genes x cells) by column duplication.
# Now each element represents one observation.
#matMu <- matMu[,vecindClusterAssign]
#matDispersions <- matDispersions[,vecindClusterAssign]

#====
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
# === 
for(i in seq(1,scaNumGenes)){
  vecboolNotZeroObserved <- matCountsProc[i,] > 0 & !is.na(matCountsProc[i,]) & is.finite(matCountsProc[i,])
  vecboolZero <- matCountsProc[i,] == 0
  vecFit <- bplapunlist( optim(
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
  if(boolSuperVerbose & vecFit["convergence"]){
    print(paste0(i, " convergence: ", vecFit["convergence"] ))
  }
  vecDispersions[i] <- exp(vecFit["par"])
}