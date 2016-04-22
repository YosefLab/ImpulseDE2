#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Impulse model fit    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Fits impulse model to a timecourse dataset with either one condition
### or case and control condition.
### Function is split into 3 parts:
###   1. [Subfunction fitImpulse_gene] Fit impulse model to a single gene
###   2. [Subfunction fitImpulse_matrix] Fit impulse model to matrix of genes.
###       Calls fitImpulse_gene, calcImpulse_comp
###   3. ["Script-body"] Prepare data and fit model by calling impulse_fit_matrix.

# INPUT:
#   arr3DCountData: (3D array genes x samples x replicates) Count data: 
#       \code{matCountData} reshaped into a 3D array. For internal use.
#   dfAnnotationRed:(data frame) Reduced version of 
#       \code{dfAnnotationFull}. For internal use.
#   strCaseName: (str) Name of the case condition in dfAnnotationFull.
#   strControlName: (str) Name of the control condition in dfAnnotationFull.
#   nProcessesAssigned: (scalar) [Default 3] Number of processes for 
#       parallelisation.
#   NPARAM: (scalar) [Default 6] Number of parameters of impulse model.
# OUTPUT:
#   lsFitResults_all: (list) List of matrices which
#       contain parameter fits and model values for given time course for the
#       case condition (and control and combined if control is present).
#       Each parameter matrix is called parameter_'condition' and has the form
#       (genes x \{"beta","h0","h1","h2","t1","t2","logL_H1","converge_H1","mu",
#       "logL_H0","converge_H0"\}) where beta to t2 are parameters of the impulse
#       model, mu is the single parameter of the mean model, logL are
#       log likelihoods of full (H1) and reduced model (H0) respectively, converge
#       is convergence status of numerical optimisation of model fitting by
#       \code{optim} from \code{stats} of either model. Each value matrix is called
#       value_'condition' and has the form (genes x time points) and contains the
#       counts predicted by the impulse model at the observed time points.

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' @param ImpulseDE2_arr3DCountData: (3D array genes x samples x replicates)
#' Count data: \code{arr2DCountData} reshaped into a 3D array. For internal use.
#' @param dfAnnotationRed: (data frame) Reduced version of 
#' \code{dfAnnotationFull}. For internal use.
#' @param strCaseName (str) Name of the case condition in \code{dfAnnotationFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#' \code{dfAnnotationFull}.
#' @param nProc (scalar) [Default 3] Number of processes for parallelisation. The
#' specified value is internally changed to \code{min(detectCores() - 1, nProc)} 
#' using the \code{detectCores} function from the package \code{parallel} to 
#' avoid overload.
#' @param NPARAM (scalar) [Default 6] Number of parameters of impulse model.
#' @return lsFitResults_all (list) List of matrices which
#' contain parameter fits and model values for given time course for the
#' case condition (and control and combined if control is present).
#' Each parameter matrix is called parameter_'condition' and has the form
#' (genes x \{"beta","h0","h1","h2","t1","t2","logL_H1","converge_H1","mu",
#' "logL_H0","converge_H0"\}) where beta to t2 are parameters of the impulse
#' model, mu is the single parameter of the mean model, logL are
#' log likelihoods of full (H1) and reduced model (H0) respectively, converge
#' is convergence status of numerical optimisation of model fitting by
#' \code{optim} from \code{stats} of either model. Each value matrix is called
#' value_'condition' and has the form (genes x time points) and contains the
#' counts predicted by the impulse model at the observed time points.

fitImpulse <- function(arr3DCountData, vecDispersions=NULL, dfAnnotation,
  strCaseName=NULL, strControlName=NULL, nProcessesAssigned=3, NPARAM=6){
  
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
  
  fitImpulse_gene <- function(matCounts, vecTimepoints, 
    NPARAM=6, MAXIT=1000, scaDispersionEstimate=NULL){
    
    optim_method <- "optim"
    #optim_method <- "nlminb"
    #optim_method <- c("optim","nlminb")
    
    # Remove timepoint if any is entirely missing
    vecboolObservedTimepoint <- apply(matCounts,1,function(tp){any(!is.na(tp))})
    if(any( !vecboolObservedTimepoint )){
      vecTimepoints <- vecTimepoints[vecboolObservedTimepoint]
      matCounts <- matCounts[vecboolObservedTimepoint,]
    }
    
    # If handed a list, i.e. not replicates
    if (is.vector(matCounts)==TRUE){
      # Transform into 2D array
      matCounts <- array(matCounts,c(length(matCounts),1))
      # Do intitial guesses based on centroids
      vecExpressionMeans <- matCounts
    } else {
      # Do initial guesses based on mean over replicates:
      vecExpressionMeans <- apply(matCounts,1,function(x){mean(x, na.rm=TRUE)})
    }
    
    nTimepts <- length(vecExpressionMeans)
    scaMaxMiddleMean <- max(vecExpressionMeans[2:(nTimepts-1)], na.rm=TRUE)
    scaMinMiddleMean <- min(vecExpressionMeans[2:(nTimepts-1)], na.rm=TRUE)
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
    indLowerInflexionPoint <- match(max(vecGradients[1:(indMaxMiddleMean-1)], na.rm=TRUE), vecGradients[1:(indMaxMiddleMean-1)])
    indUpperInflexionPoint <- indMaxMiddleMean - 1 + match(min(vecGradients[indMaxMiddleMean:length(vecGradients)], na.rm=TRUE), vecGradients[indMaxMiddleMean:length(vecGradients)])
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
    indLowerInflexionPoint <- match(min(vecGradients[1:(indMinMiddleMean-1)], na.rm=TRUE), vecGradients[1:(indMinMiddleMean-1)])
    indUpperInflexionPoint <- indMinMiddleMean - 1 + match(max(vecGradients[indMinMiddleMean:(nTimepts-1)], na.rm=TRUE), vecGradients[indMinMiddleMean:(nTimepts-1)])
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
    scaMuGuess <- log(mean(matCounts, na.rm=TRUE))
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
  #   arr3DCountDataCondition: (Numeric 3D array genes x samples x replicates).
  #       Contains expression values or similar locus-specific read-outs.
  #   vecTimepoints: (Table samples x 2 [time and condition]) 
  #       Co-variables for the samples including condition and time points.
  #       Time points must be numeric.
  #   NPARAM: 
  #   n_it: (scalar) [Defaul 100] Number of iterations, which are performed 
  #       to fit the impulse model to the clusters.
  #   strControlName
  #   nProcessesAssigned: (scalar) [Default 4] number of processes, which can be used on 
  #       the machine to run the background calculation in parallel
  #   ...
  # OUTPUT:
  #   res: (list 2 ["impulse_parameters","impulse_fits"])
  #       impulse_parameters: (matrix genes x (NPARAM+1)*runs) +1 is value of 
  #           objective of optimisation (i.e. sum of squares or weighted SS)
  #       impulse_fits: (matrix genes x vecTimepoints) model values for gene data
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  fitImpulse_matrix <- function(arr3DCountDataCondition, vecTimepoints, vecDispersions=NULL,
    strCaseName=NULL, strControlName=NULL, nProcessesAssigned=3, NPARAM=6){
    
    # Maximum number of iterations for numerical optimisation of
    # likelihood function in MLE fitting of Impulse model:
    MAXIT <- 1000
    
    # Set number of processes to number of cores assigned if available
    nProcesses <- min(detectCores() - 1, nProcessesAssigned)
    # Use parallelisation if number of genes/centroids to fit is large
    if(nrow(arr3DCountDataCondition) > max(2*nProcesses,10)){
      # Define partitioning of genes onto nodes: lsGeneIndexByCore
      lsGeneIndexByCore = list()
      bord = floor(nrow(arr3DCountDataCondition)/nProcesses)
      for (i in 1:nProcesses){
        lsGeneIndexByCore[[i]] <- ((i-1)*bord+1):(i*bord)
      }
      if(nProcesses*bord < nrow(arr3DCountDataCondition)){
        # Add remaining genes to last node
        lsGeneIndexByCore[[nProcesses]] <-  c(lsGeneIndexByCore[[nProcesses]],(nProcesses*bord+1):nrow(arr3DCountDataCondition))
      }
      
      cl <- makeCluster(nProcesses, outfile = "ImpulseDE2_ClusterOut.txt")
      #DAVID take fitImpulse out?
      my.env <- new.env()
      assign("arr3DCountDataCondition", arr3DCountDataCondition, envir = my.env)
      assign("calcImpulse_comp", calcImpulse_comp, envir = my.env)
      assign("fitImpulse", fitImpulse, envir = my.env)
      assign("evalLogLikImpulse_comp", evalLogLikImpulse_comp, envir = my.env)
      assign("evalLogLikMean_comp", evalLogLikMean_comp, envir = my.env)
      assign("fitImpulse_gene",fitImpulse_gene, envir = my.env)
      assign("fitImpulse_matrix",fitImpulse_matrix, envir = my.env)
      assign("vecTimepoints", vecTimepoints, envir = my.env)
      assign("NPARAM", NPARAM, envir = my.env)
      assign("MAXIT", MAXIT, envir = my.env)
      assign("vecDispersions", vecDispersions, envir = my.env)
      
      clusterExport(cl=cl, varlist=c("arr3DCountDataCondition","vecDispersions", "vecTimepoints",
        "calcImpulse_comp", "evalLogLikImpulse_comp", "evalLogLikMean_comp",
        "fitImpulse", "fitImpulse_matrix", "fitImpulse_gene",
        "MAXIT", "NPARAM"), envir = my.env)
      
      # Fit impulse model to each gene of matrix and get impulse parameters:
      # clusterApply runs the function impulse_fit_gene_wise
      # The input data are distributed to nodes by lsGeneIndexByCore partitioning
      lsmatFits <- clusterApply(cl, 1:length(lsGeneIndexByCore), function(z){t(sapply(lsGeneIndexByCore[[z]],
        function(x){fitImpulse_gene(matCounts=arr3DCountDataCondition[x,,],
          vecTimepoints=vecTimepoints,NPARAM=NPARAM,MAXIT=MAXIT,
          scaDispersionEstimate=vecDispersions[x])}))})
      # Give output rownames again, which are lost above
      for(i in 1:length(lsGeneIndexByCore)){
        rownames(lsmatFits[[i]]) <- rownames(arr3DCountDataCondition[lsGeneIndexByCore[[i]],,])
      }
      # Concatenate the output objects of each node
      matFits <- do.call(rbind,lsmatFits)
      stopCluster(cl)
      # Parallelisation ends here
      
      # Do not use parallelisation if number of genes to fit is small
    } else {
      # Fit impulse model to each gene of matrix and get impulse parameters
      lsGeneIndexByCore_all <- 1:nrow(arr3DCountDataCondition)
      matFits <- lapply(lsGeneIndexByCore_all,function(x){
        fitImpulse_gene(matCounts=arr3DCountDataCondition[x,,],
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
    rownames(matFits) <- rownames(arr3DCountDataCondition)
    
    ### DAVID: can reduce this if clause?
    # Use obtained impulse parameters to calculate impulse fit values   
    if(nrow(matFits) == 1){
      # if matrix contains only one gene
      matImpulseValues <- as.data.frame(t(calcImpulse_comp(matFits[,1:NPARAM],
        unique(sort(vecTimepoints)))),row.names=rownames(matFits))
    } else {                    
      # if matrix contains > 1 genes
      matImpulseValues <- t(apply(matFits[,1:NPARAM],1,function(x){calcImpulse_comp(x,
        unique(sort(vecTimepoints)))}))
    } 
    colnames(matImpulseValues) = unique(sort(vecTimepoints))
    rownames(matImpulseValues) <- rownames(arr3DCountDataCondition)
    
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
  
  #g_names = rownames(arr3DCountData)
  lsFitResults_all = list()
  
  if(!is.null(strControlName)){
    nRuns = 3
    lsLabels = c("combined","case","control")
  } else {
    nRuns = 1
    lsLabels = "case"
  }

  # arr3DCountData is defined in the following as:
  #   1) A list of three data 3D arrays combined, case and control
  #   2) A list of a single data 3D array if only have case
  
  # Perform 3 runs (combined, case and control) if control timecourse is present
  if(!is.null(strControlName)){    
    # (I) Combined: Count data of both case and control
    # Bind case and control as replicates in 3rd dimension of array
    # Find samples corresponding in time first:
    lsTimepointsCase <- as.vector( dfAnnotation[dfAnnotation$Condition %in% strCaseName,"Time"] )
    lsTimepointsCtrl <- as.vector( dfAnnotation[dfAnnotation$Condition %in% strControlName,"Time"] )
    lsTimepointsAll <- unique(c(lsTimepointsCase,lsTimepointsCtrl))
    # Pad with NA sample if one time point not present in one of the two conditions
    # Define new array dimensions
    nGenes <- dim(arr3DCountData[,,])[1]
    nTimepoints <- length(lsTimepointsAll)
    nRep <- dim(arr3DCountData[,,])[3]
    # Leave non-sampled positions in new array (i.e. entire time points of
    # either case or control) as NA.
    arr3DCountDataAll <- array(NA,c(nGenes,nTimepoints,nRep*2))
    arr3DCountDataAll[,lsTimepointsCase %in% lsTimepointsAll,1:nRep] <- arr3DCountData[,(dfAnnotation$Condition %in% strCaseName),]
    arr3DCountDataAll[,lsTimepointsCtrl %in% lsTimepointsAll,(nRep+1):(nRep*2)] <- arr3DCountData[,(dfAnnotation$Condition %in% strControlName),]
    
    # (II) Case: Count data of case
    arr3DCountDataCase = arr3DCountData[,(dfAnnotation$Condition %in% strCaseName),]
    
    # (III) Control: Count data of control
    arr3DCountDataCtrl = arr3DCountData[,(dfAnnotation$Condition %in% strControlName),]
    
    # Formatting into 3D breaks if single gene is selected:
    # DAVID rework this, doesn work?
    # If exactly one gene is tested:
    if(nGenes==1){ 
      arr3DCountDataAll_temp <- array(NA,c(1,dim(arr3DCountData)[2],dim(arr3DCountData)[3]))
      for(z in 1:dim(arr3DCountData)[3]){
        arr3DCountDataAll_temp[1,,z] <- arr3DCountDataAll[,z]
        arr3DCountDataCase_temp[1,,z] <- arr3DCountDataCase[,z]
        arr3DCountDataCtrl_temp[1,,z] <- arr3DCountDataCtrl[,z]
      }
      arr3DCountDataAll <- arr3DCountDataAll_temp
      arr3DCountDataCase <- arr3DCountDataCase_temp
      arr3DCountDataCtrl <- arr3DCountDataCtrl_temp
      rownames(arr3DCountDataAll) <- rownames(arr3DCountData)
      colnames(arr3DCountDataAll) <- colnames(arr3DCountData[,,])
      rownames(arr3DCountDataCase) <- rownames(arr3DCountData)
      colnames(arr3DCountDataCase) <- colnames(arr3DCountData[,(dfAnnotation$Condition %in% strCaseName),])
      rownames(arr3DCountDataCtrl) <- rownames(arr3DCountData)
      colnames(arr3DCountDataCtrl) <- colnames(arr3DCountData[,(dfAnnotation$Condition %in% strControlName),])
    }
    lsarr3DCountData <- list(arr3DCountDataAll, arr3DCountDataCase, arr3DCountDataCtrl)
    
  } else {
    arr3DCountDataAll = arr3DCountData[,(dfAnnotation$Condition %in% strCaseName),]
    # Formatting into 3D breaks if only single gene is selected:
    # If exactly one gene is tested:
    if(dim(arr3DCountDataAll)[1]==1){ 
      arr3DCountDataAll_temp <- array(NA,c(1,dim(arr3DCountData)[2],dim(arr3DCountData)[3]))
      for(z in 1:dim(arr3DCountData)[3]){
        arr3DCountDataAll_temp[1,,z] <- arr3DCountDataAll[,z]
      }
      arr3DCountDataAll <- arr3DCountDataAll_temp
      rownames(arr3DCountDataAll) <- rownames(arr3DCountData)
      colnames(arr3DCountDataAll) <- colnames(arr3DCountData[,(dfAnnotation$Condition %in% strCaseName),])
    }
    lsarr3DCountData <- list(arr3DCountDataAll)
  }
  
  # Fitting for different runs
  for (iRuns in 1:nRuns){
    lsFitResults_run <- fitImpulse_matrix(arr3DCountDataCondition=lsarr3DCountData[[iRuns]],
      vecTimepoints=as.numeric(as.character(
        dfAnnotation[dfAnnotation$Sample %in% colnames(lsarr3DCountData[[iRuns]]),"Time"])),
      vecDispersions=vecDispersions,
      strCaseName=strCaseName, strControlName=strControlName,
      nProcessesAssigned = nProcessesAssigned, NPARAM=NPARAM)
    
    # DAVID: do i still need this?
    #imp_res$impulse_parameters <- (imp_res$impulse_parameters)[g_names,]
    #imp_res$impulse_fits <- (imp_res$impulse_fits)[g_names,]
    
    names(lsFitResults_run) <- paste(c("parameters","values"),lsLabels[iRuns],sep="_")
    lsFitResults_all[[names(lsFitResults_run)[1]]] <- lsFitResults_run[[1]]
    lsFitResults_all[[names(lsFitResults_run)[2]]] <- lsFitResults_run[[2]]
  }
  
  return(lsFitResults_all)
}

