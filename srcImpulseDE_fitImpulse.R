#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Impulse model fit    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Fits impulse model to a timecourse dataset
### Function is split into 3 parts
###   Check data consistency
###   1. [Subfunction impulse_fit_gene_wise] Fit impulse model to a single gene
###   2. [Subfunction impulse_fit_matrix] Fit impulse model to matrix of genes.
###       Calls impulse_fit_gene_wise, calc_impulse_comp
###   3. ["Script-body"] Prepare data and fit model by calling impulse_fit_matrix.

### Developmental note:
### Three-way separation was chosen to use 'apply' instead of for-loops

# INPUT:
#   arrCountData: Two possibilities, depending on cluster or gene-wise fitting:
#       1.  Cluster. cluster_results as below.
#       2.  Gene-wise. (Numeric 3D array genes x samples x replicates).
#           Contains expression values or similar locus-specific read-outs.
#   dfAnnotation: (Table samples x 2[time and condition]) 
#       Co-variables for the samples including condition and time points.
#       Time points must be numeric.
#   n_iter: (scalar) [Defaul 100] Number of iterations, which are performed 
#       to fit the impulse model to the clusters.
#   boolControlTimecourse: (bool) [Default FALSE] Control time timecourse is 
#       part of the data set (TRUE) or not (FALSE).
#   strControlName: (str) name of the control condition in annotation_table.
#   cluster_results: (list ["kmeans_clus","cluster_means", "n_pre_clus",
#       "fine_clus"] x number of runs [combined, case, control])
#       kmeans_clus: (vector number of genes) [1, number of clusters] 
#           Indicates assignment (value of function Kmeans) 
#       cluster_means: (matrix number of clusters x gene dimensions 
#           [number of samples per run, i.e. in case])
#           Centroids of fine (2nd) clustering.
#       n_pre_clus: (scalar) Number of clusters used in pre-clustering
#       n_fine_clusters: (scalar) Number of clusters used in fine-clustering
#       NOTE: This variable is only passed in gene-wise fitting, the results
#       from clustering are passed as arrCountData in cluster fitting. This is
#       done because cluster results serve as data to be fitted in cluster
#       fitting, while they are only auxillary in gene-wise fitting.
#   start_values:
#   fit_backg: (bool) [Defaul FALSE]
#   nProcessesAssigned: (scalar) [Default 4] number of processes, which can be used on 
#       the machine to run the background calculation in parallel
# OUTPUT:
#   results: (list runs*2 ["impulse_parameters","impulse_fits"])
#       impulse_parameters: (matrix genes x (NPARAM+1)*runs) +1 is value of 
#           objective of optimisation (i.e. sum of squares or weighted SS)
#       impulse_fits: (matrix genes x vecTimepoints) model values for gene data

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

fitImpulse <- function(arrCountData, dfAnnotation,
  boolControlTimecourse = FALSE, strControlName = NULL, nProcessesAssigned = 4,
  vecDispersions=NULL, NPARAM=6){
  
  # Check dataset for consistency
  if(ncol(dfAnnotation) != 2 ||
      FALSE %in% (colnames(dfAnnotation) == c("Time","Condition"))){
    stop("Please use function 'annotation_preparation' to prepare
            the annotation table")
  }
  if(boolControlTimecourse == TRUE & is.null(strControlName)){
    stop("Please specify the name of the control in the column 'Condition'")
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~ 1. Fit impulse model to gene ~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  ### Fits an impulse model to a single gene. Fitted are 1. peak, 2. valley,
  ### 3. up, 4. down and 5. mean model, the best fit is kept.
  
  # INPUT:
  #   matCounts: (numeric 2D array vecTimepoints x replicates) 
  #   dfAnnotation: (Table samples x 2[time and condition]) 
  #       Co-variables for the samples including condition and time points.
  #       Time points must be numeric.
  #   NPARAM: (scalar) Number of model parameters
  #   n_iter: (scalar) [Default 100] Number of iterations, which are performed 
  #       to fit the impulse model to the clusters.
  #   nProcessesAssigned: (scalar) [Default 4] number of processes, which can be used on 
  #       the machine to run the background calculation in parallel
  # OUTPUT:
  #   pvec_and_objective: (vector NPARAM+1) +1 is value of objective of optimisation 
  #       (i.e. sum of squares or weighted SS)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  fitImpulse_gene_wise <- function(matCounts, vecTimepoints, 
    NPARAM=6, MAXIT=1000, scaDispersionEstimate=NULL,...){
    
    optim_method <- "optim"
    #optim_method <- "nlminb"
    #optim_method <- c("optim","nlminb")
    
    # If handed a list, i.e. not replicates
    if (is.vector(matCounts)==TRUE){
      # Transform into 2D array
      matCounts <- array(matCounts,c(length(matCounts),1))
      # Do intitial guesses based on centroids
      vecExpressionMeans <- matCounts
    } else {
      # Do initial guesses based on mean over replicates:
      vecExpressionMeans <- apply(matCounts,1,mean)
    }
    
    nTimepts <- length(vecExpressionMeans)
    scaMaxMiddleMean <- max(vecExpressionMeans[2:(nTimepts-1)])
    scaMinMiddleMean <- min(vecExpressionMeans[2:(nTimepts-1)])
    # +1 to push indicices up from middle stretch to entire window (first is omitted here)
    indMaxMiddleMean <- match(scaMaxMiddleMean,vecExpressionMeans[2:(nTimepts-1)]) + 1
    indMinMiddleMean <- match(scaMinMiddleMean,vecExpressionMeans[2:(nTimepts-1)]) + 1
    # Gradients between neighbouring points
    vecGradients <- unlist( lapply(c(1:(nTimepts-1)),function(x){
      (vecExpressionMeans[x+1]-vecExpressionMeans[x])/(vecTimepoints[x+1]-vecTimepoints[x])}) )
    
    # 1. Alternative model: Peak
    # Beta: Has to be negative
    # Theta1: Low
    # Theta2: High
    # Theta3: Low
    # t1: Around first observed inflexion point
    # t2: Around second observed inflexion point
    indLowerInflexionPoint <- match(max(vecGradients[1:(indMaxMiddleMean-1)]), vecGradients[1:(indMaxMiddleMean-1)])
    indUpperInflexionPoint <- indMaxMiddleMean - 1 + match(min(vecGradients[indMaxMiddleMean:length(vecGradients)]), vecGradients[indMaxMiddleMean:length(vecGradients)])
    vecParamGuess <- c(1,log(vecExpressionMeans[1]+1),log(scaMaxMiddleMean+1),log(vecExpressionMeans[nTimepts]+1),
      (vecTimepoints[indLowerInflexionPoint]+vecTimepoints[indLowerInflexionPoint+1])/2,
      (vecTimepoints[indUpperInflexionPoint]+vecTimepoints[indUpperInflexionPoint+1])/2)
    #vecParamGuess <- c(1,1,2,1,vecTimepoints[2],vecTimepoints[3])
    
    if("optim" %in% optim_method){
      lsFitPeak <- unlist( optim(par=vecParamGuess, fn=evalLogLikImpulse_comp, vecX=vecTimepoints,
        matY=matCounts, scaDispEst=scaDispersionEstimate,
        method="BFGS", control=list(maxit=MAXIT,fnscale=-1))[c("par","value","convergence")] )
    }
    if("nlminb" %in% optim_method){
      stop("switched optimisation to maximisation which is not used here")
      lsFitPeak2 <- unlist(nlminb(start=vecParamGuess, objective=evalLogLikImpulse_comp, vecX=vecTimepoints,
        matY=matCounts, scaDispEst=scaDispersionEstimate)[c("par","objective")])
    }
    
    # 2. Alternative model: Valley
    # Beta: Has to be negative
    # Theta1: High
    # Theta2: Low
    # Theta3: High
    # t1: Around first observed inflexion point
    # t2: Around second observed inflexion point
    indLowerInflexionPoint <- match(min(vecGradients[1:(indMinMiddleMean-1)]), vecGradients[1:(indMinMiddleMean-1)])
    indUpperInflexionPoint <- indMinMiddleMean - 1 + match(max(vecGradients[indMinMiddleMean:(nTimepts-1)]), vecGradients[indMinMiddleMean:(nTimepts-1)])
    vecParamGuess <- c(1,log(vecExpressionMeans[1]+1),log(scaMinMiddleMean+1),log(vecExpressionMeans[nTimepts]+1),
      (vecTimepoints[indLowerInflexionPoint]+vecTimepoints[indLowerInflexionPoint+1])/2,
      (vecTimepoints[indUpperInflexionPoint]+vecTimepoints[indUpperInflexionPoint+1])/2 )
    #vecParamGuess <- c(1,1,2,1,vecTimepoints[2],vecTimepoints[3])
    
    if("optim" %in% optim_method){
      lsFitValley <- unlist( optim(par=vecParamGuess, fn=evalLogLikImpulse_comp, vecX=vecTimepoints,
        matY=matCounts, scaDispEst=scaDispersionEstimate,
        method="BFGS", control=list(maxit=MAXIT,fnscale=-1))[c("par","value","convergence")] )
    }
    if("nlminb" %in% optim_method){
      stop("switched optimisation to maximisation which is not used here")
      lsFitValley2 <- unlist(nlminb(start=vecParamGuess, objective=evalLogLikImpulse_comp, vecX=vecTimepoints,
        matY=matCounts, scaDispEst=scaDispersionEstimate)[c("par","objective")])
    }
    
    # 3. Null model: Single mean
    # Parameter estimate: Overall mean
    scaMuGuess <- log(mean(matCounts))
    lsFitMean <- unlist(optim(par=scaMuGuess, fn=evalLogLikMean_comp,
      matY=matCounts, scaDispEst=scaDispersionEstimate,
      method="BFGS", control=list(maxit=MAXIT,fnscale=-1))[c("par","value","convergence")])
    
    # Summarise results
    if(("optim" %in% optim_method) && ("nlminb" %in% optim_method)){
      dfFitsByInitialisation <- cbind(lsFitPeak,lsFitPeak2,lsFitValley,lsFitValley2)
    }
    if(("optim" %in% optim_method) && !("nlminb" %in% optim_method)){
      dfFitsByInitialisation <- cbind(lsFitPeak,lsFitValley)
    }
    if(!("optim" %in% optim_method) && ("nlminb" %in% optim_method)){
      dfFitsByInitialisation <- cbind(lsFitPeak2,lsFitValley2)
    } 
    
    # Name row containing value of objective as value 
    # This name differs between optimisation functions
    # DAVID: take out when take out nlminb
    rownames(dfFitsByInitialisation) <- c(rownames(dfFitsByInitialisation[1:(nrow(dfFitsByInitialisation)-1),]),"value")
    
    # Select best fit and report fit type
    # Report mean fit objective value as null hypothesis, too.
    # match() selects first hit if maximum occurs multiple times
    indBestFit <- match(max(dfFitsByInitialisation["value",]),dfFitsByInitialisation["value",])
    lsBestFitSummary <- c(dfFitsByInitialisation[,indBestFit],lsFitMean)
    names(lsBestFitSummary) <- c("beta","h0","h1","h2","t1","t2",
      "logL_H1","converge_H1","mu","logL_H0","converge_H0")
    
    return(lsBestFitSummary)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~ 2. Fit impulse model to matrix ~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  ### Fits an impulse model to all genes of a dataset
  
  # INPUT:
  #   arrCounts: (Numeric 3D array genes x samples x replicates).
  #       Contains expression values or similar locus-specific read-outs.
  #   vecTimepoints: (Table samples x 2 [time and condition]) 
  #       Co-variables for the samples including condition and time points.
  #       Time points must be numeric.
  #   NPARAM: 
  #   n_it: (scalar) [Defaul 100] Number of iterations, which are performed 
  #       to fit the impulse model to the clusters.
  #   ctrl_tc
  #   ctrl: (bool)
  #   nProcessesAssigned: (scalar) [Default 4] number of processes, which can be used on 
  #       the machine to run the background calculation in parallel
  #   ...
  # OUTPUT:
  #   res: (list 2 ["impulse_parameters","impulse_fits"])
  #       impulse_parameters: (matrix genes x (NPARAM+1)*runs) +1 is value of 
  #           objective of optimisation (i.e. sum of squares or weighted SS)
  #       impulse_fits: (matrix genes x vecTimepoints) model values for gene data
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  fitImpulse_matrix <- function(arrCounts, vecTimepoints,
    ctrl_tc = FALSE, ctrl = NULL, nProcessesAssigned = 3, vecDispersions=NULL,NPARAM=6){
    
    # Maximum number of iterations for numerical optimisation of
    # likelihood function in MLE fitting of Impulse model:
    MAXIT <- 1000
    
    # Set number of processes to number of cores assigned if available
    nProcesses <- min(detectCores() - 1, nProcessesAssigned)
    # Use parallelisation if number of genes/centroids to fit is large
    if(nrow(arrCounts) > max(2*nProcesses,10)){
      # Define partitioning of genes onto nodes: lsGeneIndexByCore
      lsGeneIndexByCore = list()
      bord = floor(nrow(arrCounts)/nProcesses)
      for (i in 1:nProcesses){
        lsGeneIndexByCore[[i]] <- ((i-1)*bord+1):(i*bord)
      }
      if(nProcesses*bord < nrow(arrCounts)){
        # Add remaining genes to last node
        lsGeneIndexByCore[[nProcesses]] <-  c(lsGeneIndexByCore[[nProcesses]],(nProcesses*bord+1):nrow(arrCounts))
      }
      
      cl <- makeCluster(nProcesses, outfile = "clus_out_impulse_fit.txt")
      #DAVDI take fitImpulse out?
      my.env <- new.env()
      assign("arrCounts", arrCounts, envir = my.env)
      assign("calcImpulse_comp", calcImpulse_comp, envir = my.env)
      assign("fitImpulse", fitImpulse, envir = my.env)
      assign("evalLogLikImpulse_comp", evalLogLikImpulse_comp, envir = my.env)
      assign("evalLogLikMean_comp", evalLogLikMean_comp, envir = my.env)
      assign("fitImpulse_gene_wise",fitImpulse_gene_wise, envir = my.env)
      assign("fitImpulse_matrix",fitImpulse_matrix, envir = my.env)
      assign("vecTimepoints", vecTimepoints, envir = my.env)
      assign("NPARAM", NPARAM, envir = my.env)
      assign("MAXIT", MAXIT, envir = my.env)
      assign("vecDispersions", vecDispersions, envir = my.env)
      
      clusterExport(cl=cl, varlist=c("arrCounts","vecDispersions", "vecTimepoints",
        "calcImpulse_comp", "evalLogLikImpulse_comp", "evalLogLikMean_comp",
        "fitImpulse", "fitImpulse_matrix", "fitImpulse_gene_wise",
        "MAXIT", "NPARAM"), envir = my.env)
      
      # Fit impulse model to each gene of matrix and get impulse parameters:
      # clusterApply runs the function impulse_fit_gene_wise
      # The input data are distributed to nodes by lsGeneIndexByCore partitioning
      lsmatFits <- clusterApply(cl, 1:length(lsGeneIndexByCore), function(z){t(sapply(lsGeneIndexByCore[[z]],
        function(x){fitImpulse_gene_wise(matCounts=arrCounts[x,,],
          vecTimepoints=vecTimepoints,NPARAM=NPARAM,MAXIT=MAXIT,
          scaDispersionEstimate=vecDispersions[x])}))})
      # Give output rownames again, which are lost above
      for(i in 1:length(lsGeneIndexByCore)){
        rownames(lsmatFits[[i]]) <- rownames(arrCounts[lsGeneIndexByCore[[i]],,])
      }
      # Concatenate the output objects of each node
      matFits <- do.call(rbind,lsmatFits)
      stopCluster(cl)
      # Parallelisation ends here
      
      # Do not use parallelisation if number of genes to fit is small
    } else {
      # Fit impulse model to each gene of matrix and get impulse parameters
      lsGeneIndexByCore_all <- 1:nrow(arrCounts)
      matFits <- lapply(lsGeneIndexByCore_all,function(x){
        fitImpulse_gene_wise(matCounts=arrCounts[x,,],
          vecTimepoints=vecTimepoints,NPARAM=NPARAM,MAXIT=MAXIT,
          scaDispersionEstimate=vecDispersions[x])})
      matFits_temp <- array(NA,c(length(matFits),length(matFits[[1]])))
      for (i in 1:length(matFits)){
        matFits_temp[i,] <- matFits[[i]]
      }
      matFits <- matFits_temp
    }   
    colnames(matFits) <- c("beta","h0","h1","h2","t1","t2","logL_H1",
      "converge_H1","mu","logL_H0","converge_H0")
    rownames(matFits) <- rownames(arrCounts)
    
    ### DAVID: can reduce this if clause?
    # Use obtained impulse parameters to calculate impulse fit values   
    if(nrow(matFits) == 1){
      # if matrix contains only one gene
      matImpulseValues <- as.data.frame(t(calcImpulse_comp(matFits[,1:NPARAM],
        unique(sort(vecTimepoints)))),row.names=rownames(matFits))
    } else {                    
      # if matrix contains > 1 genes
      matImpulseValues <- t(apply(matFits[,1:NPARAM],1,function(x){calc_impulse_comp(x,
        unique(sort(vecTimepoints)))}))
    } 
    colnames(matImpulseValues) = unique(sort(vecTimepoints))
    rownames(matImpulseValues) <- rownames(arrCounts)
    
    # Report results.
    lsFitResults_matrix <- list()
    lsFitResults_matrix[[1]] <- matFits   
    # ["beta","h0","h1","h2","t1","t2","logL_H1",
    # "converge_H1","mu","logL_H0","converge_H0"] (x3 if with control)
    lsFitResults_matrix[[2]] <- matImpulseValues
    # Values of impulse model on input timepoints by gene
    return(lsFitResults_matrix)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~ 3. Prepare data for impulse model fit ~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  #g_names = rownames(arrCountData)
  lsFitResults_all = list()
  
  if(boolControlTimecourse == TRUE){
    nRuns = 3
    lsLabels = c("combined","case","control")
  } else {
    nRuns = 1
    lsLabels = "case"
  }
  
  # arrCountData is defined in the following as:
  #   1) A list of three data 3D arrays combined, case and control
  #   2) A list of a single data 3D array if only have case
  
  # Perform 3 runs (combined, case and control) if control timecourse is present
  if(boolControlTimecourse == TRUE){    
    # Combined: Count data of both case and control (all columns).
    arrCountDataAll = arrCountData
    # Case: Count data of case.
    arrCountDataCase = arrCountData[,!(dfAnnotation$Condition %in% strControlName),]
    # Control: Count data of control.  
    arrCountDataCtrl = arrCountData[,(dfAnnotation$Condition %in% strControlName),]
    
    # Formatting into 3D breaks if single gene is selected:
    # DAVID rework this, doesn work?
    # If exactly one gene is tested:
    if(dim(arrCountDataAll)[1]==1){ 
      arrCountDataAll_temp <- array(NA,c(1,dim(arrCountData)[2],dim(arrCountData)[3]))
      for(z in 1:dim(arrCountData)[3]){
        arrCountDataAll_temp[1,,z] <- arrCountDataAll[,z]
        arrCountDataCase_temp[1,,z] <- arrCountDataCase[,z]
        arrCountDataCtrl_temp[1,,z] <- arrCountDataCtrl[,z]
      }
      arrCountDataAll <- arrCountDataAll_temp
      arrCountDataCase <- arrCountDataCase_temp
      arrCountDataCtrl <- arrCountDataCtrl_temp
      rownames(arrCountDataAll) <- rownames(arrCountData)
      colnames(arrCountDataAll) <- colnames(arrCountData)
      rownames(arrCountDataCase) <- rownames(arrCountData)
      colnames(arrCountDataCase) <- colnames(arrCountData[,!(dfAnnotation$Condition %in% strControlName),])
      rownames(arrCountDataCtrl) <- rownames(arrCountData)
      colnames(arrCountDataCtrl) <- colnames(arrCountData[,(dfAnnotation$Condition %in% strControlName),])
    }
    arrCountData <- list(arrCountDataAll, arrCountDataCase, arrCountDataCtrl)
    
  } else if(boolControlTimecourse == FALSE){
    arrCountDataAll = arrCountData[,,]
    # Formatting into 3D breaks if only single gene is selected:
    # If exactly one gene is tested:
    if(dim(arrCountDataAll)[1]==1){ 
      arrCountDataAll_temp <- array(NA,c(1,dim(arrCountData)[2],dim(arrCountData)[3]))
      for(z in 1:dim(arrCountData)[3]){
        arrCountDataAll_temp[1,,z] <- arrCountDataAll[,z]
      }
      arrCountDataAll <- arrCountDataAll_temp
      rownames(arrCountDataAll) <- rownames(arrCountData)
      colnames(arrCountDataAll) <- colnames(arrCountData)
    }
    arrCountData <- list(arrCountDataAll)
  }
  
  # Fitting for different runs
  for (iRuns in 1:nRuns){
    lsFitResults_run <- impulse_fit_matrix(arrCounts=arrCountData[[iRuns]],
      vecTimepoints=as.numeric(as.character(
        dfAnnotation[colnames(arrCountData[[iRuns]]),"Time"])), 
      ctrl_tc = boolControlTimecourse,ctrl = strControlName,
      nProcessesAssigned = nProcessesAssigned, vecDispersions=vecDispersions,NPARAM=NPARAM)
    
    # DAVID: do i still need this?
    #imp_res$impulse_parameters <- (imp_res$impulse_parameters)[g_names,]
    #imp_res$impulse_fits <- (imp_res$impulse_fits)[g_names,]
    
    names(lsFitResults_run) <- paste(c("parameters","values"),lsLabels[iRuns],sep="_")
    lsFitResults_all[[names(lsFitResults_run)[1]]] <- lsFitResults_run[[1]]
    lsFitResults_all[[names(lsFitResults_run)[2]]] <- lsFitResults_run[[2]]
  }
  
  return(lsFitResults_all)
}

