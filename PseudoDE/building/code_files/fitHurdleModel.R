fitHurdleModel <- function(matCounts, lsResultsClustering, 
  dfAnnotationClusters){
  
  MAXIT <- 1000
  nProcesses <- 3
  
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
  # counts. With one pseudo count.
  matMeanEst <- array(NA,c(dim(matCounts)[1],lsResultsClustering$K))
  rownames(matMeanEst) <- rownames(matCounts)
  colnames(matMeanEst) <- c(1:lsResultsClustering$K)
  lsMeanEst <- lapply( c(1:dim(matCounts)[1]), function(gene){
    sapply(c(1:lsResultsClustering$K), function(cluster){ 
      log(mean( matCounts[gene,lsResultsClustering$Assignments==cluster], na.rm=TRUE ) +1)
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
    
    # Define partitioning of genes onto nodes: lsGeneIndexByCore
    lsGeneIndexByCore = list()
    bord = floor(nrow(matCounts)/nProcesses)
    for (i in 1:nProcesses){
      lsGeneIndexByCore[[i]] <- ((i-1)*bord+1):(i*bord)
    }
    if(nProcesses*bord < nrow(matCounts)){
      # Add remaining genes to last node
      lsGeneIndexByCore[[nProcesses]] <-  c(lsGeneIndexByCore[[nProcesses]],(nProcesses*bord+1):nrow(matCounts))
    }
    
    cl <- makeCluster(nProcesses)
    
    my.env <- new.env()
    assign("vecDispEst", vecDispEst, envir = my.env)
    assign("matMeanEst", matMeanEst, envir = my.env)
    assign("vecDropoutEst", vecDropoutEst, envir = my.env)
    assign("matCounts", matCounts, envir = my.env)
    assign("lsResultsClustering", lsResultsClustering, envir = my.env)
    assign("lsGeneIndexByCore", lsGeneIndexByCore, envir = my.env)
    assign("evalLogLikHurdleNB_comp", evalLogLikHurdleNB_comp, envir = my.env)
    assign("evalLogLikHurdleDrop_comp", evalLogLikHurdleDrop_comp, envir = my.env)
    assign("MAXIT", MAXIT, envir = my.env)
    
    clusterExport(cl=cl, varlist=c("vecDispEst","matMeanEst", "vecDropoutEst",
      "matCounts", "lsResultsClustering", "lsGeneIndexByCore",
      "evalLogLikHurdleNB_comp", "evalLogLikHurdleDrop_comp", 
      "MAXIT"), envir = my.env)
    
    lsNBEst_new <- clusterApply(cl, c(1:length(lsGeneIndexByCore)), function(core){
      t(sapply( lsGeneIndexByCore[[core]], function(gene){
        unlist(optim(
          par=c(vecDispEst[gene],matMeanEst[gene,]), 
          fn=evalLogLikHurdleNB_comp,
          vecY=matCounts[gene,], 
          vecDropoutEst=vecDropoutEst,
          lsResultsClustering=lsResultsClustering,
          method="BFGS", 
          control=list(maxit=MAXIT,fnscale=-1) 
        )["par"])
      }))
    })
    # Give output rownames
    for(i in 1:length(lsNBEst_new)){
      rownames(lsNBEst_new[[i]]) <- rownames(matCounts[lsGeneIndexByCore[[i]],])
    }
    matNBEst_new <- do.call(rbind,lsNBEst_new)
    
    stopCluster(cl)
    
    vecDispEst_new <- matNBEst_new[,1]
    matMeanEst_new <- matNBEst_new[,2:(dim(matNBEst_new)[2])]
    
    # Fit drop-out rates
    print("DO")
    
    # Define partitioning of cells onto nodes: lsCellIndexByCore
    lsCellIndexByCore = list()
    bord = floor(ncol(matCounts)/nProcesses)
    for (i in 1:nProcesses){
      lsCellIndexByCore[[i]] <- ((i-1)*bord+1):(i*bord)
    }
    if(nProcesses*bord < ncol(matCounts)){
      # Add remaining genes to last node
      lsCellIndexByCore[[nProcesses]] <-  c(lsCellIndexByCore[[nProcesses]],(nProcesses*bord+1):ncol(matCounts))
    }
    
    cl <- makeCluster(nProcesses)
    
    my.env <- new.env()
    assign("vecDispEst_new", vecDispEst, envir = my.env)
    assign("matMeanEst_new", matMeanEst, envir = my.env)
    assign("vecDropoutEst", vecDropoutEst, envir = my.env)
    assign("matCounts", matCounts, envir = my.env)
    assign("lsResultsClustering", lsResultsClustering, envir = my.env)
    assign("lsCellIndexByCore", lsCellIndexByCore, envir = my.env)
    assign("evalLogLikHurdleNB_comp", evalLogLikHurdleNB_comp, envir = my.env)
    assign("evalLogLikHurdleDrop_comp", evalLogLikHurdleDrop_comp, envir = my.env)
    assign("MAXIT", MAXIT, envir = my.env)
    
    clusterExport(cl=cl, varlist=c("vecDispEst_new","matMeanEst_new", "vecDropoutEst",
      "matCounts", "lsResultsClustering", "lsCellIndexByCore",
      "evalLogLikHurdleNB_comp", "evalLogLikHurdleDrop_comp", 
      "MAXIT"), envir = my.env)
    lsDropoutEst_new <- clusterApply(cl, c(1:length(lsCellIndexByCore)), function(core){
      t(sapply( lsCellIndexByCore[[core]], function(cell){
        unlist(optim( 
          par=vecDropoutEst[cell], 
          fn=evalLogLikHurdleDrop_comp,
          vecY=matCounts[,cell], 
          vecMeanEst=matMeanEst_new[,lsResultsClustering$Assignments[cell]],
          vecDispEst=vecDispEst_new,
          method="L-BFGS-B",
          lower=0, upper=1,
          control=list(maxit=MAXIT,fnscale=-1) 
        )["par"])
      }))
    })
    # Give output elment names
    for(i in 1:length(lsDropoutEst_new)){
      names(lsDropoutEst_new[[i]]) <- colnames(matCounts[,lsCellIndexByCore[[i]]])
    }
    vecDropoutEst_new <- unlist(lsDropoutEst_new)
    
    stopCluster(cl)
    
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