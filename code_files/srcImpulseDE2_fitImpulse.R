#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Impulse model fit    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Fit an impulse model to a single gene
#' 
#' Fits an impulse model and a mean model to a single gene.
#' The optimisation method is set within this function (\code{optim}:
#' BFGS). This method is divided into four parts: (I) Prepare data,
#' (II) Fit mean model, (III) Fit impulse model, (IV) Process Fits.
#' Internal function of \code{fitImpulse}.
#' 
#' @seealso Called by \code{fitImpulse_matrix}. Calls \code{evalLogLikImpulse},
#'    and \code{calcImpulse_comp}.
#' 
#' @param matCounts (numeric 2D array vecTimepoints x replicates) Count data.
#' @param scaDispersionEstimate (scalar) Inverse of negative binomial 
#'    dispersion coefficients computed by DESeq2 for given gene.
#' @param matNormConst: (matrix samples x replicates) Normalisation
#'    constants for each replicate. Missing samples are set NA.
#' @param vecTimepoints (vector number timepoints) Observed timepoints, numeric.
#' @param NPARAM (scalar) [Default 6] Number of impulse model parameters
#' @param MAXIT (scalar) [Default 100] Number of iterations, which are performed 
#'    to fit the impulse model to the clusters.
#' 
#' @return lsBestFitSummary (list [beta, h0, h1, h2, t1, t2, logL_H1, 
#'    converge_H1, mu, logL_H0]) Beta to t2 are parameters of the impulse
#'    model, mu is the single parameter of the mean model, logL are
#'    log likelihoods of full (H1) and reduced model (H0) respectively, converge
#'    is convergence status of numerical optimisation of model fitting by
#'    \code{optim} from \code{stats}. Each value matrix is called
#'    value_'condition' and has the form (genes x time points) and contains the
#'    counts predicted by the impulse model at the observed time points.)
#' @export

fitImpulse_gene <- function(matCounts, scaDispersionEstimate, 
  matNormConst, vecTimepoints,
  strMode="batch", NPARAM=6, MAXIT=1000){
  
  optim_method <- "optim"
  #optim_method <- "nlminb"
  #optim_method <- c("optim","nlminb")
  
  # (I) Prepare data
  # Remove timepoint if any is entirely missing
  vecboolObservedTimepoint <- apply(matCounts,1,function(tp){any(!is.na(tp))})
  if(any( !vecboolObservedTimepoint )){
    vecTimepoints <- vecTimepoints[vecboolObservedTimepoint]
    matCounts <- matCounts[vecboolObservedTimepoint,]
  }
  # Remove time course if is entirely missing
  vecboolObservedTimecourse <- apply(matCounts,2,function(tp){any(!is.na(tp))})
  if(any( !vecboolObservedTimecourse )){
    matCounts <- matCounts[vecboolObservedTimecourse,]
  }
  # Transform into 2D array if handed a list, i.e. one replicate
  if (is.vector(matCounts)==TRUE){
    matCounts <- array(matCounts,c(length(matCounts),1))
  }
  # Compute statistics for initialisation
  vecExpressionMeans <- apply(matCounts,1,function(x){mean(x, na.rm=TRUE)})
  nTimepts <- length(vecExpressionMeans)
  scaMaxMiddleMean <- max(vecExpressionMeans[2:(nTimepts-1)], na.rm=TRUE)
  scaMinMiddleMean <- min(vecExpressionMeans[2:(nTimepts-1)], na.rm=TRUE)
  # +1 to push indicices up from middle stretch to entire window (first is omitted here)
  indMaxMiddleMean <- match(scaMaxMiddleMean,vecExpressionMeans[2:(nTimepts-1)]) + 1
  indMinMiddleMean <- match(scaMinMiddleMean,vecExpressionMeans[2:(nTimepts-1)]) + 1
  # Gradients between neighbouring points
  vecGradients <- unlist( lapply(c(1:(nTimepts-1)),function(x){
    (vecExpressionMeans[x+1]-vecExpressionMeans[x])/(vecTimepoints[x+1]-vecTimepoints[x])}) )
  
  # Precompute matrices used in cost functions:
  matboolObserved <- !is.na(matCounts)
  matboolZero <- matCounts==0
  
  # DAVID: todo export this into extra function and handle entire matrix at once
  # (II) Fit mean model
  # Fit mean to normalised counts. Count normalisation 
  # corresponds to model normalisation in impulse model 
  # fitting.
  scaMu <- mean(matCounts, na.rm=TRUE)
  if(strMode=="batch" | strMode=="singlecell"){
    # Null model (timecourses): One mean for all points
    # Fitting is done above: scaMu is used in different parts
    # of code, too. 
    # Evaluate likelihood of null model
    scaLoglikNull <- sum(dnbinom(
      matCounts[matboolObserved], 
      mu=scaMu, 
      size=scaDispersionEstimate, 
      log=TRUE))
  }else if(strMode=="timecourses"){
    # Null model (timecourses): Fit one mean for each time course
    nTimecourses <- dim(matCounts)[2]
    # Fit null model.
    vecMuTimecourses <- apply(matCounts,2,function(vecTC){mean(vecTC, na.rm=TRUE)})
    matMuTimecourses <- matrix(vecMuTimecourses, nrow=dim(matCounts)[1], ncol=dim(matCounts)[2], byrow=TRUE)
    
    # Evaluate likelihood of null model
    scaLoglikNull <- sum(dnbinom(
      matCounts[matboolObserved], 
      mu=matMuTimecourses[matboolObserved], 
      size=scaDispersionEstimate, 
      log=TRUE))
    
    # Compute translation factors: Normalisation factor
    # to scale impulse model to indivindual time courses.
    # Note: Translation factor is the same for all replicates
    # in a time course.
    # Reference: Overall scaled mean
    scaMuScaled <- mean(matCounts/matNormConst, na.rm=TRUE)
    # Scaled means by timecourse
    vecMuTimecoursesScaled <- apply(matCounts/matNormConst,2,function(vecTC){mean(vecTC, na.rm=TRUE)})
    matMuScaledTimecourses <- matrix(vecMuTimecoursesScaled, nrow=dim(matCounts)[1], ncol=dim(matCounts)[2], byrow=TRUE)
    # Scaled mean ratio
    matTranslationFactors <- matMuScaledTimecourses/scaMuScaled
  } else {
    stop(paste0("ERROR: Unrecognised strMode in fitImpulse(): ",strMode))
  }
  
  # (III) Fit Impulse model
  # 1. Initialisation: Peak
  # Beta: Has to be negative, Theta1: Low, Theta2: High, Theta3: Low
  # t1: Around first observed inflexion point, t2: Around second observed inflexion point
  indLowerInflexionPoint <- match(max(vecGradients[1:(indMaxMiddleMean-1)], na.rm=TRUE), vecGradients[1:(indMaxMiddleMean-1)])
  indUpperInflexionPoint <- indMaxMiddleMean - 1 + match(min(vecGradients[indMaxMiddleMean:length(vecGradients)], na.rm=TRUE), vecGradients[indMaxMiddleMean:length(vecGradients)])
  vecParamGuess <- c(1,log(vecExpressionMeans[1]+1),log(scaMaxMiddleMean+1),log(vecExpressionMeans[nTimepts]+1),
    (vecTimepoints[indLowerInflexionPoint]+vecTimepoints[indLowerInflexionPoint+1])/2,
    (vecTimepoints[indUpperInflexionPoint]+vecTimepoints[indUpperInflexionPoint+1])/2)
  
  if("optim" %in% optim_method){
    if(strMode=="batch"){
      lsFitPeak <- unlist( optim(
        par=vecParamGuess,
        fn=evalLogLikImpulseBatch_comp,
        vecX=vecTimepoints,
        matY=matCounts,
        scaDispEst=scaDispersionEstimate,
        matNormConst=matNormConst,
        matboolObserved=matboolObserved,
        method="BFGS",
        control=list(maxit=MAXIT,fnscale=-1)
      )[c("par","value","convergence")] )
    }else if(strMode=="timecourses"){
      lsFitPeak <- unlist( optim(
        par=vecParamGuess, 
        fn=evalLogLikImpulseByTC_comp, 
        vecX=vecTimepoints,
        matY=matCounts, 
        scaDispEst=scaDispersionEstimate,
        matNormConst=matNormConst,
        matTranslationFactors=matTranslationFactors,
        matboolObserved=matboolObserved,
        method="BFGS", 
        control=list(maxit=MAXIT,fnscale=-1)
      )[c("par","value","convergence")] )
    }else if(strMode=="singlecell"){
      lsFitPeak <- unlist( optim(par=vecParamGuess, 
        fn=evalLogLikImpulseSC_comp, 
        vecX=vecTimepoints,
        matY=matCounts, 
        scaDispEst=scaDispersionEstimate,
        scaDropoutEst=scaDropoutEstimate,
        matNormConst=matNormConst,
        matboolObserved=matboolObserved, 
        matboolZero=matboolZero,
        method="BFGS", 
        control=list(maxit=MAXIT,fnscale=-1)
      )[c("par","value","convergence")] )
    } else {
      stop(paste0("ERROR: Unrecognised strMode in fitImpulse(): ",strMode))
    }
  }
  if("nlminb" %in% optim_method){
    stop("switched optimisation to maximisation which is not used here")
    lsFitPeak2 <- unlist(nlminb(start=vecParamGuess, objective=evalLogLikImpulseByTC_comp, vecX=vecTimepoints,
      matY=matCounts, scaDispEst=scaDispersionEstimate)[c("par","objective")])
  }
  
  # 2. Initialisation: Valley
  # Beta: Has to be negative, Theta1: High, Theta2: Low, Theta3: High
  # t1: Around first observed inflexion point, t2: Around second observed inflexion point
  indLowerInflexionPoint <- match(min(vecGradients[1:(indMinMiddleMean-1)], na.rm=TRUE), vecGradients[1:(indMinMiddleMean-1)])
  indUpperInflexionPoint <- indMinMiddleMean - 1 + match(max(vecGradients[indMinMiddleMean:(nTimepts-1)], na.rm=TRUE), vecGradients[indMinMiddleMean:(nTimepts-1)])
  vecParamGuess <- c(1,log(vecExpressionMeans[1]+1),log(scaMinMiddleMean+1),log(vecExpressionMeans[nTimepts]+1),
    (vecTimepoints[indLowerInflexionPoint]+vecTimepoints[indLowerInflexionPoint+1])/2,
    (vecTimepoints[indUpperInflexionPoint]+vecTimepoints[indUpperInflexionPoint+1])/2 )
  
  if("optim" %in% optim_method){
    if(strMode=="batch"){
      lsFitValley <- unlist( optim(
        par=vecParamGuess, 
        fn=evalLogLikImpulseBatch_comp, 
        vecX=vecTimepoints,
        matY=matCounts, 
        scaDispEst=scaDispersionEstimate,
        matNormConst=matNormConst,
        matboolObserved=matboolObserved,
        method="BFGS", 
        control=list(maxit=MAXIT,fnscale=-1)
      )[c("par","value","convergence")] )
    }else if(strMode=="timecourses"){
      lsFitValley <- unlist( optim(
        par=vecParamGuess, 
        fn=evalLogLikImpulseByTC_comp, 
        vecX=vecTimepoints,
        matY=matCounts, 
        scaDispEst=scaDispersionEstimate,
        matNormConst=matNormConst,
        matTranslationFactors=matTranslationFactors,
        matboolObserved=matboolObserved,
        method="BFGS", 
        control=list(maxit=MAXIT,fnscale=-1)
      )[c("par","value","convergence")] )
    }else if(strMode=="singlecell"){
      lsFitValley <- unlist( optim(
        par=vecParamGuess, 
        fn=evalLogLikImpulseSC_comp, 
        vecX=vecTimepoints,
        matY=matCounts, 
        scaDispEst=scaDispersionEstimate, 
        scaDropoutEst=scaDropoutEstimate,
        matNormConst=matNormConst,
        matboolObserved=matboolObserved, 
        matboolZero=matboolZero,
        method="BFGS", 
        control=list(maxit=MAXIT,fnscale=-1)
      )[c("par","value","convergence")] )
    } else {
      stop(paste0("ERROR: Unrecognised strMode in fitImpulse(): ",strMode))
    }
  }
  if("nlminb" %in% optim_method){
    stop("switched optimisation to maximisation which is not used here")
    lsFitValley2 <- unlist(nlminb(start=vecParamGuess, objective=evalLogLikImpulseByTC_comp, vecX=vecTimepoints,
      matY=matCounts, scaDispEst=scaDispersionEstimate)[c("par","objective")])
  }
  
  # (IV) Process fits
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
  if(strMode=="batch"){
    lsBestFitSummary <- c(dfFitsByInitialisation[,indBestFit],scaMu,scaLoglikNull)
    names(lsBestFitSummary) <- c("beta","h0","h1","h2","t1","t2",
      "logL_H1","converge_H1","mu","logL_H0")
  } else {
    vecMuTimecoursesWithNA <- array(NA,length(vecboolObservedTimecourse))
    vecMuTimecoursesWithNA[vecboolObservedTimecourse] <- vecMuTimecourses
    lsBestFitSummary <- c(dfFitsByInitialisation[,indBestFit],scaMu,scaLoglikNull,vecMuTimecoursesWithNA)
    nTimecourses <- length(vecboolObservedTimecourse)
    vecColnamesMubyTimecourse <- paste0(rep("muByTimecourse",nTimecourses),c(1:nTimecourses))
    names(lsBestFitSummary) <- c("beta","h0","h1","h2","t1","t2",
      "logL_H1","converge_H1","mu","logL_H0",vecColnamesMubyTimecourse)
  }
  # Remove log scale from count parameters of impulse model
  lsBestFitSummary[c("h0","h1","h2")] <- exp(lsBestFitSummary[c("h0","h1","h2")])
  
  return(lsBestFitSummary)
}

#' Fits impulse models to all genes of a dataset
#' 
#' Fits impulse models to all genes of a dataset. Performs parallelisation.
#' Internal function of \code{fitImpulse}.
#' 
#' @details Maximum number of iterations for optimisation is set within this 
#'    function as \code{MAXIT}. In the case of single condition differential 
#'    expression over time, this function is called once for the case condition.
#'    In the case of case and control condition, this function is called three 
#'    times: data from case, control and combined conditions.
#'    Calls \code{fitImpulse_gene}.
#' 
#' @seealso Called by \code{fitImpulse}. Calls \code{fitImpulse_gene}.
#'   
#' @param ImpulseDE2_arr3DCountData (3D array genes x samples x replicates)
#'    Count data: \code{arr2DCountData} reshaped into a 3D array. For internal use.
#' @param vecDispersions (vector number of genes) Inverse of gene-wise 
#'    negative binomial dispersion coefficients computed by DESeq2.
#' @param matNormConst: (matrix samples x replicates) Normalisation
#'    constants for each replicate. Missing samples are set NA.
#' @param dfAnnotationRedRed (data frame) Reduced version of 
#'    \code{dfAnnotationRedFull}. Lists co-variables of samples: 
#'    Sample, Condition, Time. Time must be numeric. For internal use.
#' @param strCaseName (str) Name of the case condition in \code{dfAnnotationRedFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#'    \code{dfAnnotationRedFull}.
#' @param nProcessesAssigned (scalar) [Default 3] Number of processes for parallelisation. The
#'    specified value is internally changed to \code{min(detectCores() - 1, nProc)} 
#'    using the \code{detectCores} function from the package \code{parallel} to 
#'    avoid overload.
#' @param NPARAM (scalar) [Default 6] Number of parameters of impulse model.
#' 
#' @return lsFitResults_matrix (list length 2) List of two matrices which
#'    contain parameter fits and model values for given time course for the
#'    given condition.
#' \itemize{
#'    \item \code{parameter}(genes x [beta, h0, h1, h2, t1, t2, logL_H1, 
#'        converge_H1, mu, logL_H0, converge_H0]) 
#'        Beta to t2 are parameters of the impulse model, mu is the single 
#'        parameter of the mean model, logL are log likelihoods of full 
#'        (H1) and reduced model (H0) respectively, converge is convergence 
#'        status of numerical optimisation of model fitting by
#'        \code{optim} from \code{stats} of either model.
#'    \item \code{values} (genes x time points) Contains the counts 
#'        predicted by the impulse model at the observed time points.
#' }
#' @export

fitImpulse_matrix <- function(arr3DCountDataCondition, vecDispersions, 
  matNormConst, vecTimepoints,
  strCaseName, strControlName=NULL, strMode="batch", 
  nProcessesAssigned=3, NPARAM=6){
  
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
    # Variables
    assign("arr3DCountDataCondition", arr3DCountDataCondition, envir = my.env)
    assign("vecDispersions", vecDispersions, envir = my.env)
    assign("matNormConst", matNormConst, envir = my.env)
    assign("vecTimepoints", vecTimepoints, envir = my.env)
    assign("strMode", strMode, envir = my.env)
    assign("NPARAM", NPARAM, envir = my.env)
    assign("MAXIT", MAXIT, envir = my.env)
    # Functions
    assign("fitImpulse", fitImpulse, envir = my.env)
    assign("fitImpulse_gene",fitImpulse_gene, envir = my.env)
    assign("fitImpulse_matrix",fitImpulse_matrix, envir = my.env)
    assign("calcImpulse_comp", calcImpulse_comp, envir = my.env)
    assign("evalLogLikImpulseBatch_comp", evalLogLikImpulseBatch_comp, envir = my.env)
    assign("evalLogLikImpulseByTC_comp", evalLogLikImpulseByTC_comp, envir = my.env)
    assign("evalLogLikImpulseSC_comp", evalLogLikImpulseSC_comp, envir = my.env)
    
    clusterExport(cl=cl, varlist=c(
      "arr3DCountDataCondition",
      "vecDispersions", 
      "matNormConst",
      "vecTimepoints",
      "strMode",
      "MAXIT",
      "NPARAM",
      "calcImpulse_comp", 
      "evalLogLikImpulseBatch_comp", 
      "evalLogLikImpulseByTC_comp",
      "evalLogLikImpulseSC_comp", 
      "fitImpulse", 
      "fitImpulse_matrix", 
      "fitImpulse_gene"
    ), envir = my.env)
    
    # Fit impulse model to each gene of matrix and get impulse parameters:
    # clusterApply runs the function impulse_fit_gene_wise
    # The input data are distributed to nodes by lsGeneIndexByCore partitioning
    lsmatFits <- clusterApply(cl, 1:length(lsGeneIndexByCore),
      function(z){ t(sapply( lsGeneIndexByCore[[z]],
        function(x){fitImpulse_gene(
          matCounts=arr3DCountDataCondition[x,,],
          scaDispersionEstimate=vecDispersions[x],
          matNormConst=matNormConst,
          vecTimepoints=vecTimepoints, 
          strMode=strMode,
          NPARAM=NPARAM,
          MAXIT=MAXIT )}
      ))})
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
      fitImpulse_gene(
        matCounts=arr3DCountDataCondition[x,,],
        scaDispersionEstimate=vecDispersions[x],
        matNormConst=matNormConst,
        vecTimepoints=vecTimepoints,
        strMode=strMode,
        NPARAM=NPARAM,
        MAXIT=MAXIT )})
    matFits_temp <- array(NA,c(length(matFits),length(matFits[[1]])))
    for (i in 1:length(matFits)){
      matFits_temp[i,] <- matFits[[i]]
    }
    matFits <- matFits_temp
  }
  
  # Name output columns
  if(strMode=="batch"){
    colnames(matFits) <- c("beta","h0","h1","h2","t1","t2","logL_H1",
      "converge_H1","mu","logL_H0")
  } else {
    nTimecourses <- dim(arr3DCountDataCondition)[3]
    vecColnamesMubyTimecourse <- paste0(rep("muByTimecourse",nTimecourses),c(1:nTimecourses))
    colnames(matFits) <- c("beta","h0","h1","h2","t1","t2","logL_H1",
      "converge_H1","mu","logL_H0",vecColnamesMubyTimecourse)
  }
  rownames(matFits) <- rownames(arr3DCountDataCondition)
  
  ### DAVID: can reduce this if clause?
  # Use obtained impulse parameters to calculate impulse fit values
  matTheta <- matFits[,1:NPARAM]
  matTheta[,2:4] <- log(matTheta[,2:4])
  if(nrow(matFits) == 1){
    # if matrix contains only one gene
    matImpulseValues <- as.data.frame(t(calcImpulse_comp(matTheta,
      unique(sort(vecTimepoints)))),row.names=rownames(matFits))
  } else {                    
    # if matrix contains > 1 genes
    matImpulseValues <- t(apply(matTheta,1,function(x){calcImpulse_comp(x,
      unique(sort(vecTimepoints)))}))
  } 
  colnames(matImpulseValues) = unique(sort(vecTimepoints))
  rownames(matImpulseValues) <- rownames(arr3DCountDataCondition)
  
  # Report results.
  lsFitResults_matrix <- list()
  lsFitResults_matrix[[1]] <- matFits   
  # ["beta","h0","h1","h2","t1","t2","logL_H1",
  # "converge_H1","mu","logL_H0",("muByTimecourse")] (x3 if with control)
  lsFitResults_matrix[[2]] <- matImpulseValues
  # Values of impulse model on input timepoints by gene
  names(lsFitResults_matrix) <- c("parameters","values")
  return(lsFitResults_matrix)
}

#' Fits impulse model to a timecourse dataset
#' 
#' This function coordinates helper functions for fitting:
#'  1. [Helper \code{fitImpulse_gene}] Fit impulse model to a single gene.
#'      Calls cost functions (\code{evalLogLikImpulseBatch_comp},
#'      \code{evalLogLikImpulseByTC_comp} and \code{evalLogLikImpulseSC_comp})
#'      and \code{calcImpulse_comp}.
#'  2. [Helper \code{fitImpulse_matrix}] Fit impulse model to matrix of genes.
#'      Calls \code{fitImpulse_gene}.
#'  3. [This function] Prepare data and fit model. Calls \code{impulse_fit_matrix}.
#'  
#' @seealso Called by \code{runImpulseDE2}. Calls \code{fitImpulse_matrix}.
#'  
#' @param arr3DCountData (3D array genes x samples x replicates)
#'    Count data: \code{arr2DCountData} reshaped into a 3D array.
#'    For internal use.
#' @param vecDispersions (vector number of genes) Inverse of gene-wise 
#'    negative binomial dispersion coefficients computed by DESeq2.
#' @param matNormConst: (matrix samples x replicates) Normalisation
#'    constants for each replicate. Missing samples are set NA.
#' @param dfAnnotationRedRed (data frame) Reduced version of 
#'    \code{dfAnnotationRedFull}. Lists co-variables of samples: 
#'    Sample, Condition, Time. Time must be numeric. For internal use.
#' @param strCaseName (str) Name of the case condition in \code{dfAnnotationRedFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#'    \code{dfAnnotationRedFull}.
#' @param nProcessesAssigned (scalar) [Default 3] Number of processes for parallelisation. The
#'    specified value is internally changed to \code{min(detectCores() - 1, nProc)} 
#'    using the \code{detectCores} function from the package \code{parallel} to 
#'    avoid overload.
#' @param NPARAM (scalar) [Default 6] Number of parameters of impulse model.
#' 
#' @return lsFitResults_all (list length 2 or 6) List of matrices which
#'    contain parameter fits and model values for given time course for the
#'    case condition (and control and combined if control is present).
#'    Each parameter matrix is called parameter_'condition' and has the form
#'    (genes x [beta, h0, h1, h2, t1, t2, logL_H1, converge_H1, mu, logL_H0, 
#'    converge_H0]) where beta to t2 are parameters of the impulse
#'    model, mu is the single parameter of the mean model, logL are
#'    log likelihoods of full (H1) and reduced model (H0) respectively, converge
#'    is convergence status of numerical optimisation of model fitting by
#'    \code{optim} from \code{stats} of either model. Each value matrix is called
#'    value_'condition' and has the form (genes x time points) and contains the
#'    counts predicted by the impulse model at the observed time points.
#' @export

fitImpulse <- function(arr3DCountData, vecDispersions=NULL, 
  matNormConst, dfAnnotationRed,
  strCaseName=NULL, strControlName=NULL, strMode="batch",
  nProcessesAssigned=3, NPARAM=6){
  
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
    lsTimepointsCase <- as.vector( dfAnnotationRed[dfAnnotationRed$Condition %in% strCaseName,"Time"] )
    lsTimepointsCtrl <- as.vector( dfAnnotationRed[dfAnnotationRed$Condition %in% strControlName,"Time"] )
    lsTimepointsAll <- unique(c(lsTimepointsCase,lsTimepointsCtrl))
    # Pad with NA sample if one time point not present in one of the two conditions
    # Define new array dimensions
    nGenes <- dim(arr3DCountData[,,])[1]
    nTimepoints <- length(lsTimepointsAll)
    nRep <- dim(arr3DCountData[,,])[3]
    # Leave non-sampled positions in new array (i.e. entire time points of
    # either case or control) as NA.
    arr3DCountDataAll <- array(NA,c(nGenes,nTimepoints,nRep*2))
    arr3DCountDataAll[,lsTimepointsCase %in% lsTimepointsAll,1:nRep] <- arr3DCountData[,(dfAnnotationRed$Condition %in% strCaseName),]
    arr3DCountDataAll[,lsTimepointsCtrl %in% lsTimepointsAll,(nRep+1):(nRep*2)] <- arr3DCountData[,(dfAnnotationRed$Condition %in% strControlName),]
    
    # (II) Case: Count data of case
    arr3DCountDataCase = arr3DCountData[,(dfAnnotationRed$Condition %in% strCaseName),]
    
    # (III) Control: Count data of control
    arr3DCountDataCtrl = arr3DCountData[,(dfAnnotationRed$Condition %in% strControlName),]
    
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
      colnames(arr3DCountDataCase) <- colnames(arr3DCountData[,(dfAnnotationRed$Condition %in% strCaseName),])
      rownames(arr3DCountDataCtrl) <- rownames(arr3DCountData)
      colnames(arr3DCountDataCtrl) <- colnames(arr3DCountData[,(dfAnnotationRed$Condition %in% strControlName),])
    }
    lsarr3DCountData <- list(arr3DCountDataAll, arr3DCountDataCase, arr3DCountDataCtrl)
    
  } else {
    arr3DCountDataAll = arr3DCountData[,(dfAnnotationRed$Condition %in% strCaseName),]
    # Formatting into 3D breaks if only single gene is selected:
    # If exactly one gene is tested:
    if(dim(arr3DCountDataAll)[1]==1){ 
      arr3DCountDataAll_temp <- array(NA,c(1,dim(arr3DCountData)[2],dim(arr3DCountData)[3]))
      for(z in 1:dim(arr3DCountData)[3]){
        arr3DCountDataAll_temp[1,,z] <- arr3DCountDataAll[,z]
      }
      arr3DCountDataAll <- arr3DCountDataAll_temp
      rownames(arr3DCountDataAll) <- rownames(arr3DCountData)
      colnames(arr3DCountDataAll) <- colnames(arr3DCountData[,(dfAnnotationRed$Condition %in% strCaseName),])
    }
    lsarr3DCountData <- list(arr3DCountDataAll)
  }
  
  # Fitting for different runs
  for (iRuns in 1:nRuns){
    lsFitResults_run <- fitImpulse_matrix(
      arr3DCountDataCondition=lsarr3DCountData[[iRuns]],
      matNormConst=matNormConst,
      vecDispersions=vecDispersions,
      vecTimepoints=as.numeric(as.character(
        dfAnnotationRed[dfAnnotationRed$Sample %in% colnames(lsarr3DCountData[[iRuns]]),"Time"])),
      strCaseName=strCaseName,
      strControlName=strControlName,
      strMode=strMode,
      nProcessesAssigned=nProcessesAssigned,
      NPARAM=NPARAM )
    
    # DAVID: do i still need this?
    #imp_res$impulse_parameters <- (imp_res$impulse_parameters)[g_names,]
    #imp_res$impulse_fits <- (imp_res$impulse_fits)[g_names,]
    
    names(lsFitResults_run) <- paste(c("parameters","values"),lsLabels[iRuns],sep="_")
    lsFitResults_all[[names(lsFitResults_run)[1]]] <- lsFitResults_run[[1]]
    lsFitResults_all[[names(lsFitResults_run)[2]]] <- lsFitResults_run[[2]]
  }
  
  return(lsFitResults_all)
}

