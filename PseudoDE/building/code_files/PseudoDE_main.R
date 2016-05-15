library(compiler)
library(parallel)
library(DESeq2)
library(BiocParallel)
library(scone)

source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/ImpulseDE2_main.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/srcImpulseDE2_CostFunctionsFit.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/srcImpulseDE2_runDESeq2.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/clusterCellsInPseudotime.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/createAnnotation.R")

runPseudoDE <- function(matCounts, vecPseudotime,
  nProc=1){
  
  MAXITER <- 3
  
  # 1. Data preprocessing
  print("1. Data preprocessing:")
  tm_preproc <- system.time({
    # Check data
    # 1. matCounts
    if(any(matCounts %% 1 != 0)){
      stop("ERROR: matCounts contains non-integer elements. Requires count data.")
    }
    if(any(is.na(matCounts) | !is.finite(matCounts) | matCounts<0 )){
      stop("ERROR: Non positive counts (NA/Inf/-Inf/negative).")
    }
    # 2. vecPseudotime
    if(!all(names(vecPseudotime) %in% colnames(matCounts))){
      stop("ERROR: Cell names in vecPseudotime and matCounts do not match.")
    }
    
    # Convert from data frame to matrix for zinb()
    if(is.list(matCounts)){
      matCounts <- data.matrix(matCounts)
    }
    # Take out NA cells
    vecPseudotime <- vecPseudotime[!is.na(vecPseudotime)]
    # Adjust ording of cells in objects
    matCounts <- matCounts[,names(vecPseudotime)]
    # Remove all zero genes: Mu is initialised
    # as log(sum counts)
    matCounts <- matCounts[apply(matCounts,1,function(gene){any(gene!=0)}),]
    
  })
  
  # 2. Cluster cells in pseudo-time
  print("2. Clustering:")
  tm_clustering <- system.time({
    lsResultsClustering <- clusterCellsInPseudotime(vecPseudotime=vecPseudotime)
    # Fit by cluster
    dfAnnotation <- createAnnotationByCluster(matCounts=matCounts,
      vecPseudotime=vecPseudotime,
      lsResultsClustering=lsResultsClustering)
    # Fit by cell
    #dfAnnotation <- createAnnotationByCell(matCounts=matCounts,
    #  vecPseudotime=vecPseudotime)
  })
  save(lsResultsClustering,file=file.path(getwd(),"PseudoDE_lsResultsClustering.RData"))
  save(dfAnnotation,file=file.path(getwd(),"PseudoDE_dfAnnotation.RData"))
  print(paste("Consumed time: ",round(tm_clustering["elapsed"]/60,2),
    " min",sep=""))
  
  # 3. Fit mixture model
  print("3. Fit mixture model:")
  tm_fitmm <- system.time({
    # Set number of processes to be used in this step:
    # register(MulticoreParam()) controls the number of processes used for 
    # BiocParallel, used in zinb().
    #nProcesses <- min(detectCores() - 1, nProc)
    #register(MulticoreParam(nProcesses))
    
    # Fit zero-inflated negative binomial model to each cluster
    # Imputed count data matrix
    matCountsImputed <- array(0,c(dim(matCounts)[1],dim(matCounts)[2]))
    # Genes which are all 0 in one cluster recive pi=1, i.e. they
    # drop out of the cost function during fitting as their likelihood
    # is defined by the hyperparameter pi=1 independently of the
    # mean model.
    matDropout <- array(1,c(dim(matCounts)[1],dim(matCounts)[2]))
    matProbNB <- array(1,c(dim(matCounts)[1],dim(matCounts)[2]))
    # Negative binomial over-dispersion factor: 1 per gene per cluster
    matDispersion <- array(NA,c(dim(matCounts)[1],lsResultsClustering$K))
    lsZinbOutputByCluster <- list()
    for(k in 1:lsResultsClustering$K){
      print(paste0("Fitting cluster ",k))
      
      # Remove all zero genes: Mu is initialised
      # as log(sum counts)
      vecidxCluster <- lsResultsClustering$Assignments==k
      matCountsCluster <- matCounts[,vecidxCluster]
      #vecboolNonzeroGenes <- apply(matCountsCluster,1,
      #  function(gene){any(gene!=0)})
      vecboolNonzeroGenes <- apply(matCountsCluster,1,
        function(gene){mean(gene)>10})
      matCountsCluster <- matCountsCluster[vecboolNonzeroGenes,]
      
      # Fit zinb model
      lsZINBparam <- estimate_zinb(
        Y = matCountsCluster, 
        maxiter = MAXITER, 
        verbose = TRUE)
      lsZinbOutputByCluster[[k]] <- lsZINBparam
      
      # Record parameters
      matDropout[vecboolNonzeroGenes,vecidxCluster] <- lsZINBparam$pi
      matProbNB[vecboolNonzeroGenes,vecidxCluster] <- lsZINBparam$p_z
      matDispersion[vecboolNonzeroGenes,k] <- lsZINBparam$theta
    
      # Impute count data matrix
      matCountsImputed[vecboolNonzeroGenes,vecidxCluster] <- 
        matCountsCluster * lsZINBparam$p_z + 
        lsZINBparam$mu * (1 - lsZINBparam$p_z)
    }
    # Put objects together
    lsInputToImpulseDE2 <- list(matDropout, matProbNB, matCountsImputed)
    names(lsInputToImpulseDE2) <- c("matDropout", "matProbNB", "matCountsImputed")
  })
  save(lsZinbOutputByCluster,file=file.path(getwd(),"PseudoDE_lsZinbOutputByCluster.RData"))
  save(lsInputToImpulseDE2,file=file.path(getwd(),"PseudoDE_lsInputToImpulseDE2.RData"))
  print(paste("Consumed time: ",round(tm_fitmm["elapsed"]/60,2),
    " min",sep=""))
  
  # 4. Differential expression analysis: Run ImpulseDE2
  print("4. Differential expression analysis: Run ImpulseDE2:")
  print("### Begin ImpulseDE2 output ################################")
  tm_deanalysis <- system.time({
    lsImpulseDE2results <- runImpulseDE2(
      matCountData = matCounts, 
      dfAnnotationFull = dfAnnotation,
      strCaseName = "case", 
      strControlName = NULL, 
      strMode = "singlecell", 
      nProc = nProc, 
      Q_value = 0.01,
      boolPlotting = TRUE,
      lsPseudo = lsInputToImpulseDE2)
  })
  print("### End ImpulseDE2 output ##################################")
  save(lsImpulseDE2results,file=file.path(getwd(),"PseudoDE_lsImpulseDE2results.RData"))
  print(paste("Consumed time: ",round(tm_deanalysis["elapsed"]/60,2),
    " min",sep=""))
  
  print("Completed PseudoDE.")
  return(lsImpulseDE2results)
}