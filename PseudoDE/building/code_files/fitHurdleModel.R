fitHurdleModel <- function(matCounts, lsResultsClustering, 
  dfAnnotationClusters){
  
  MAXIT <- 1000
  
  #########################
  ### large optimisation of everything together
  #########################
  if(FALSE){
    # 1.Identify overdispersion factor of negative binomial
    # by gene:
    # From GLM based on categorial independent variable
    # (cluster index). Dispersion estimation is performed by 
    # DESeq2, which extends GLM fitting of overdispersion by
    # an Empririal Bayes correction.
    matCountsDESeq2 <- matCounts
    matCountsDESeq2[matCountsDESeq2==0] <- NA
    lsDESeq2ResultsClusters <- runDESeq2(
      dfAnnotationFull=dfAnnotationClusters, 
      arr2DCountData=matCounts, 
      strMode="singlecell")
    vecDESeq2Dispersions <- lsDESeq2ResultsClusters[[1]]
    dfDESeq2ResultsClusters <- lsDESeq2ResultsClusters[[2]]
    
    # 2. MLE of mean of negative binomials of each cluster
    # of each gene and MLE of drop-out rate for each cell:
    # Means and dropout rates are fitted to data in the same
    # step through numerical optimisation of the hurdle
    # likelihood function.
    
    # a) Initialisation: Initialise based on the assumption
    # that all zero counts are drop-outs.
    
    # Means are initialised to cluster averages over NON-ZERO
    # counts.
    matMeanEst <- array(NA,c(dim(matCounts)[1],lsResultsClustering$K))
    rownames(matMeanEst) <- rownames(matCounts)
    colnames(matMeanEst) <- c(1:lsResultsClustering$K)
    for(gene in rownames(matCounts)){
      matMeanEst[gene,] <- sapply(c(1:lsResultsClustering$K),
        function(cluster){ matMeanEst[gene,cluster] <- 
            mean( matCounts[gene,lsResultsClustering$Assignments==cluster], na.rm=TRUE )})
    }
    # Dropout rates are initialised to fraction of 
    # ZERO counts in cell.
    vecDropoutGuess <- apply(matCounts,2,
      function(cell){ sum(cell==0)/sum(!is.na(cell) )})
    # Assemble parameter object to be handed to fitting function
    lsThetaHurdleGuess <- list(vecDropoutGuess, matMeanEst)
    names(lsThetaHurdleGuess) <- c("vecDropoutEst","matMeanEst")  
    
    # b) Fitting of hurdle model by numerical optimisation
    # of likelihood function.
    lsThetaHurdleEst <- unlist(optim( 
      par=vecDropoutGuess, 
      fn=evalLogLikImpulseBatch_comp,
      matY=matCounts, 
      vecDispEst=vecDESeq2Dispersions,
      lsResultsClustering=lsResultsClustering,
      method="BFGS", 
      control=list(maxit=MAXIT,fnscale=-1) 
    )["par"])
  }
  #########################
  ### iterative optimisation
  #########################
  print("Initialise")
  # 1. Initialisation: Initialise based on the assumption
  # that all zero counts are drop-outs.
  
  # Means are initialised to cluster averages over NON-ZERO
  # counts.
  matMeanEst <- array(NA,c(dim(matCounts)[1],lsResultsClustering$K))
  rownames(matMeanEst) <- rownames(matCounts)
  colnames(matMeanEst) <- c(1:lsResultsClustering$K)
  lsMeanEst <- lapply( matMeanEst,1, function(gene){
    sapply(c(1:lsResultsClustering$K), function(cluster){ 
      mean( gene[lsResultsClustering$Assignments==cluster], na.rm=TRUE )
    })
  })
  matMeanEst <- do.call(rbind,lsMeanEst)
  
  # Random guess of dispersions: All 1
  vecDispEst <- array(1,dim(matCounts)[1])
  # Dropout rates are initialised to fraction of 
  # ZERO counts in cell.
  vecDropoutEst <- apply(matCounts,2,
    function(cell){ sum(cell==0)/sum(!is.na(cell) )})
  
  # 2. Iterative fitting:
  boolConverved <- FALSE
  cntr <- 0
  while(!boolConverved){
    cntr <- cntr +1
    print(cntr)
    
    # During fitting, parameters are initialised
    # with prior values.
    # Fit negative binomial
    print("NB")
    lsNBEst_new <- lapply( c(1:dim(matCounts)[1]), function(gene){
      unlist(optim(
        par=c(vecDispEst[gene],matMeanEst[gene,]), 
        fn=evalLogLikHurdleNB_comp,
        vecY=matCounts[gene,], 
        vecDropoutEst=vecDropoutEst,
        lsResultsClustering=lsResultsClustering,
        method="BFGS", 
        control=list(maxit=MAXIT,fnscale=-1) 
      )["par"])
    })
    matNBEst_new <- do.call(rbind,lsNBEst_new)
    vecDispEst_new <- matNBEst_new[,1]
    matMeanEst_new <- matNBEst_new[,2:(dim(matNBEst_new)[2])]
    # Fit drop-out rates
    print("DO")
    lsDropoutEst_new <- lapply( c(1:dim(matCounts)[2]), function(cell){
      unlist(optim( 
        par=vecDropoutEst[cell], 
        fn=evalLogLikHurdleDrop_comp,
        vecY=matCounts[,cell], 
        vecMeanEst=matMeanEst_new[,lsResultsClustering$Assignments[cell]],
        vecDispEst=vecDispEst_new,
        method="BFGS", 
        control=list(maxit=MAXIT,fnscale=-1) 
      )["par"])
    })
    vecDropoutEst_new <- unlist(lsDropoutEst_new)
    
    # Test for convergence:
    if(cntr == 3){
      boolConverged <- TRUE
    }
    
    # Update parameter guesses
    matMeanEst <- matMeanEst_new
    vecDispEst <- vecDispEst_new
    vecDropoutEst <- vecDropoutEst_new
  }
  
  #########################
  
  # Assemble output list.
  #lsParamHurdle <- list(vecDESeq2Dispersions,lsThetaHurdleEst$vecDropoutEst)
  lsParamHurdle <- list(vecDispEst,vecDropoutEst)
  names(lsParamHurdle) <- c("Dispersions","DropoutRates")
  
  return(lsParamHurdle)
}