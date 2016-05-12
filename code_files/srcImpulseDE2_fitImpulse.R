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
#' @param vecNormConst: (numeric vector number of replicates) 
#'    Normalisation constants for each replicate.
#' @param vecTimepointAssign (numeric vector number replicates) Timepoints 
#'    assigned to replicates.
#' @param vecTimecourseAssign (numeric vector number replicates) Time courses 
#'    assigned to replicates. NULL if not operating in strMode="timecourses".
#' @param dfAnnotationFull (Table) Lists co-variables of individual replicates:
#'    Replicate, Sample, Condition, Time. Time must be numeric.
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

fitImpulse_gene <- function(vecCounts, scaDispersionEstimate, 
  vecNormConst,
  vecTimepointAssign, vecTimecourseAssign, dfAnnotationFull,
  strMode="batch", NPARAM=6, MAXIT=1000){
  
  optim_method <- "optim"
  #optim_method <- "nlminb"
  #optim_method <- c("optim","nlminb")
  
  # Get boolean observation vectors:
  vecboolObserved <- !is.na(vecCounts)
  vecboolZero <- vecCounts==0
  
  # The vectors vecTimepoints and vecTimecourses are shared
  # between all genes of a condition if no observations
  # are NA. They are computed for every gene to account
  # for time point or time course drop  out, if all replicates
  # corresponding to a time point or time course are NA.
  # Get list of time points observed in this condition and gene
  vecTimepoints <- sort(unique( vecTimepointAssign ))
  nTimepts <- length(vecTimepoints)
  # Get vector of numeric time point assignment indices:
  # One index pointing at an element in vecTimepoints for
  # each element in vecTimepointAssign, i.e. each replicate.
  vecindTimepointAssign <- match(vecTimepointAssign, vecTimepoints)
  # Get list of time courses observed in this condition and gene
  if(strMode=="timecourses"){
    vecTimecourses <- unique( vecTimecourseAssign )
    nTimecourses <- length(vecTimecourses)
  }
  
  # DAVID: todo export this into extra function and handle entire matrix at once
  # (II) Fit mean model
  # Fit mean to normalised counts. Count normalisation 
  # corresponds to model normalisation in impulse model 
  # fitting.
  scaMu <- mean(vecCounts, na.rm=TRUE)
  if(strMode=="batch" | strMode=="singlecell"){
    # Null model (no longitudinal sampling): 
    # Fit one mean for all replicates.
    
    # Fitting is done above: scaMu is used in different parts
    # of code, too. 
    # Evaluate likelihood of null model
    scaLoglikNull <- sum(dnbinom(
      vecCounts[vecboolObserved], 
      mu=scaMu, 
      size=scaDispersionEstimate, 
      log=TRUE))
  }else if(strMode=="timecourses"){
    # Null model (longitudinal sampling): 
    # Fit one mean for each time course (longitudinal series).
    
    # Fit null model.
    vecMuTimecourses <- sapply(vecTimecourses,
      function(tc){mean(vecCounts[vecTimecourseAssign==tc], na.rm=TRUE)})
    names(vecMuTimecourses) <- vecTimecourses
    
    # Evaluate likelihood of null model
    scaLoglikNull <- sum(dnbinom(
      vecCounts[vecboolObserved], 
      mu=(vecMuTimecourses[vecTimecourseAssign])[vecboolObserved], 
      size=scaDispersionEstimate, 
      log=TRUE))
  } else {
    stop(paste0("ERROR: Unrecognised strMode in fitImpulse(): ",strMode))
  }
  
  # Compute translation factors for fitting:
  if(strMode=="timecourses"){
    # Compute translation factors: Normalisation factor
    # to scale impulse model to indivindual time courses.
    # Note: Translation factor is the same for all replicates
    # in a time course.
    # Reference: Overall scaled mean
    scaMuScaled <- mean(vecCounts/vecNormConst, na.rm=TRUE)
    # Scaled means by timecourse
    vecMuTimecoursesScaled <- sapply(vecTimecourses, function(tc){
      mean(vecCounts[vecTimecourseAssign==tc]/vecNormConst[vecTimecourseAssign==tc], na.rm=TRUE)
    })
    names(vecMuTimecoursesScaled) <- vecTimecourses
    # Scaled mean ratio per replicate
    vecTranslationFactors <- vecMuTimecoursesScaled[vecTimecourseAssign]/scaMuScaled
  }
  # Compute statistics for initialisation:
  # Expression means by timepoint
  vecExpressionMeans <- sapply(vecTimepoints,
    function(tp){mean(vecCounts[vecTimepointAssign==tp], na.rm=TRUE)})  
  scaMaxMiddleMean <- max(vecExpressionMeans[2:(nTimepts-1)], na.rm=TRUE)
  scaMinMiddleMean <- min(vecExpressionMeans[2:(nTimepts-1)], na.rm=TRUE)
  # +1 to push indicices up from middle stretch to entire window (first is omitted here)
  indMaxMiddleMean <- match(scaMaxMiddleMean,vecExpressionMeans[2:(nTimepts-1)]) + 1
  indMinMiddleMean <- match(scaMinMiddleMean,vecExpressionMeans[2:(nTimepts-1)]) + 1
  # Gradients between neighbouring points
  vecGradients <- unlist( lapply(c(1:(nTimepts-1)),function(x){
    (vecExpressionMeans[x+1]-vecExpressionMeans[x])/(vecTimepoints[x+1]-vecTimepoints[x])}) )
  
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
        vecY=vecCounts,
        scaDispEst=scaDispersionEstimate,
        vecNormConst=vecNormConst,
        vecindTimepointAssign=vecindTimepointAssign,
        vecboolObserved=vecboolObserved,
        method="BFGS",
        control=list(maxit=MAXIT,fnscale=-1)
      )[c("par","value","convergence")] )
    }else if(strMode=="timecourses"){
      lsFitPeak <- unlist( optim(
        par=vecParamGuess, 
        fn=evalLogLikImpulseByTC_comp, 
        vecX=vecTimepoints,
        vecY=vecCounts, 
        scaDispEst=scaDispersionEstimate,
        vecNormConst=vecNormConst,
        vecTranslationFactors=vecTranslationFactors,
        vecindTimepointAssign=vecindTimepointAssign,
        vecboolObserved=vecboolObserved,
        method="BFGS", 
        control=list(maxit=MAXIT,fnscale=-1)
      )[c("par","value","convergence")] )
    }else if(strMode=="singlecell"){
      lsFitPeak <- unlist( optim(par=vecParamGuess, 
        fn=evalLogLikImpulseSC_comp, 
        vecX=vecTimepoints,
        vecY=vecCounts, 
        scaDispEst=scaDispersionEstimate,
        scaDropoutEst=scaDropoutEstimate,
        vecNormConst=vecNormConst,
        vecindTimepointAssign=vecindTimepointAssign,
        vecboolObserved=vecboolObserved, 
        vecboolZero=vecboolZero,
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
        vecY=vecCounts, 
        scaDispEst=scaDispersionEstimate,
        vecNormConst=vecNormConst,
        vecindTimepointAssign=vecindTimepointAssign,
        vecboolObserved=vecboolObserved,
        method="BFGS", 
        control=list(maxit=MAXIT,fnscale=-1)
      )[c("par","value","convergence")] )
    }else if(strMode=="timecourses"){
      lsFitValley <- unlist( optim(
        par=vecParamGuess, 
        fn=evalLogLikImpulseByTC_comp, 
        vecX=vecTimepoints,
        vecY=vecCounts, 
        scaDispEst=scaDispersionEstimate,
        vecNormConst=vecNormConst,
        vecTranslationFactors=vecTranslationFactors,
        vecindTimepointAssign=vecindTimepointAssign,
        vecboolObserved=vecboolObserved,
        method="BFGS", 
        control=list(maxit=MAXIT,fnscale=-1)
      )[c("par","value","convergence")] )
    }else if(strMode=="singlecell"){
      lsFitValley <- unlist( optim(
        par=vecParamGuess, 
        fn=evalLogLikImpulseSC_comp, 
        vecX=vecTimepoints,
        vecY=vecCounts, 
        scaDispEst=scaDispersionEstimate, 
        scaDropoutEst=scaDropoutEstimate,
        vecNormConst=vecNormConst,
        vecindTimepointAssign=vecindTimepointAssign,
        vecboolObserved=vecboolObserved,
        vecboolZero=vecboolZero,
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
  if(strMode=="batch" | strMode=="singlecell"){
    lsBestFitSummary <- c(dfFitsByInitialisation[,indBestFit],scaMu,scaLoglikNull)
    names(lsBestFitSummary) <- c("beta","h0","h1","h2","t1","t2",
      "logL_H1","converge_H1","mu","logL_H0")
  } else if(strMode=="timecourses") {
    lsBestFitSummary <- c(dfFitsByInitialisation[,indBestFit],scaMu,scaLoglikNull,vecMuTimecourses)
    nTimecourses <- length(vecMuTimecourses)
    vecColnamesMubyTimecourse <- paste0(rep("muByTimecourse",nTimecourses),c(1:nTimecourses))
    names(lsBestFitSummary) <- c("beta","h0","h1","h2","t1","t2",
      "logL_H1","converge_H1","mu","logL_H0",vecColnamesMubyTimecourse)
  }  else {
    stop(paste0("ERROR: Unrecognised strMode in fitImpulse(): ",strMode))
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
#' @param arr2DCountData (2D array genes x replicates) Count data: Reduced 
#'    version of \code{matCountData}. For internal use.
#' @param vecDispersions (vector number of genes) Inverse of gene-wise 
#'    negative binomial dispersion coefficients computed by DESeq2.
#' @param vecNormConst: (numeric vector number of replicates) 
#'    Normalisation constants for each replicate.
#' @param vecTimepointAssign (numeric vector number replicates) Timepoints 
#'    assigned to replicates.
#' @param vecTimecourseAssign (numeric vector number replicates) Time courses 
#'    assigned to replicates. NULL if not operating in strMode="timecourses".
#' @param dfAnnotationFull (Table) Lists co-variables of individual replicates:
#'    Replicate, Sample, Condition, Time. Time must be numeric.
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

fitImpulse_matrix <- function(arr2DCountDataCondition, vecDispersions, 
  vecNormConst, 
  vecTimepointAssign, vecTimecourseAssign, dfAnnotationFull,
  strCaseName, strControlName=NULL, strMode="batch", 
  nProcessesAssigned=3, NPARAM=6){
  
  # Maximum number of iterations for numerical optimisation of
  # likelihood function in MLE fitting of Impulse model:
  MAXIT <- 1000
  
  # Set number of processes to number of cores assigned if available
  nProcesses <- min(detectCores() - 1, nProcessesAssigned)
  if(nrow(arr2DCountDataCondition) > max(2*nProcesses,10)){
    # Use parallelisation if number of genes/centroids to fit is large
    
    # Define partitioning of genes onto nodes: lsGeneIndexByCore
    lsGeneIndexByCore = list()
    bord = floor(nrow(arr2DCountDataCondition)/nProcesses)
    for (i in 1:nProcesses){
      lsGeneIndexByCore[[i]] <- ((i-1)*bord+1):(i*bord)
    }
    if(nProcesses*bord < nrow(arr2DCountDataCondition)){
      # Add remaining genes to last node
      lsGeneIndexByCore[[nProcesses]] <-  c(lsGeneIndexByCore[[nProcesses]],(nProcesses*bord+1):nrow(arr2DCountDataCondition))
    }
    
    cl <- makeCluster(nProcesses, outfile = "ImpulseDE2_ClusterOut.txt")
    #DAVID take fitImpulse out?
    my.env <- new.env()
    # Variables
    assign("arr2DCountDataCondition", arr2DCountDataCondition, envir = my.env)
    assign("vecDispersions", vecDispersions, envir = my.env)
    assign("vecNormConst", vecNormConst, envir = my.env)
    assign("vecTimepointAssign", vecTimepointAssign, envir = my.env)
    assign("vecTimecourseAssign", vecTimecourseAssign, envir = my.env)
    assign("dfAnnotationFull", dfAnnotationFull, envir = my.env)
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
      "arr2DCountDataCondition",
      "vecDispersions", 
      "vecNormConst",
      "vecTimepointAssign",
      "vecTimecourseAssign",
      "dfAnnotationFull",
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
          vecCounts=arr2DCountDataCondition[x,],
          scaDispersionEstimate=vecDispersions[x],
          vecNormConst=vecNormConst,
          vecTimepointAssign=vecTimepointAssign,
          vecTimecourseAssign=vecTimecourseAssign,
          dfAnnotationFull=dfAnnotationFull,
          strMode=strMode,
          NPARAM=NPARAM,
          MAXIT=MAXIT )}
      ))})
    # Give output rownames again, which are lost above
    for(i in 1:length(lsGeneIndexByCore)){
      rownames(lsmatFits[[i]]) <- rownames(arr2DCountDataCondition[lsGeneIndexByCore[[i]],])
    }
    # Concatenate the output objects of each node
    matFits <- do.call(rbind,lsmatFits)
    # Parallelisation ends here
    stopCluster(cl)
    
  } else {
    # Do not use parallelisation if number of genes to fit is small
    
    # Fit impulse model to each gene of matrix and get impulse parameters
    lsGeneIndexByCore_all <- 1:nrow(arr2DCountDataCondition)
    matFits <- lapply(lsGeneIndexByCore_all,function(x){
      fitImpulse_gene(
        vecCounts=arr2DCountDataCondition[x,],
        scaDispersionEstimate=vecDispersions[x],
        vecNormConst=vecNormConst,
        vecTimepointAssign=vecTimepointAssign,
        vecTimecourseAssign=vecTimecourseAssign,
        dfAnnotationFull=dfAnnotationFull,
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
    nTimecourses <- length(unique( vecTimecourseAssign ))
    vecColnamesMubyTimecourse <- paste0(rep("muByTimecourse",nTimecourses),c(1:nTimecourses))
    colnames(matFits) <- c("beta","h0","h1","h2","t1","t2","logL_H1",
      "converge_H1","mu","logL_H0",vecColnamesMubyTimecourse)
  }
  rownames(matFits) <- rownames(arr2DCountDataCondition)
  
  ### DAVID: can reduce this if clause?
  # Use obtained impulse parameters to calculate impulse fit values
  matTheta <- matFits[,1:NPARAM]
  matTheta[,2:4] <- log(matTheta[,2:4])
  if(nrow(matFits) == 1){
    # if matrix contains only one gene
    matImpulseValues <- as.data.frame(t(calcImpulse_comp(matTheta,
      unique(sort(vecTimepointAssign)))),row.names=rownames(matFits))
  } else {                    
    # if matrix contains > 1 genes
    matImpulseValues <- t(apply(matTheta,1,function(x){calcImpulse_comp(x,
      unique(sort(vecTimepointAssign)))}))
  } 
  colnames(matImpulseValues) = unique(sort(vecTimepointAssign))
  rownames(matImpulseValues) <- rownames(arr2DCountDataCondition)
  
  # Report results.
  lsFitResults_matrix <- list()
  lsFitResults_matrix[[1]] <- matFits
  lsFitResults_matrix[[2]] <- matImpulseValues
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
#' @param arr2DCountData (2D array genes x replicates) Count data: Reduced 
#'    version of \code{matCountData}. For internal use.
#' @param vecDispersions (vector number of genes) Inverse of gene-wise 
#'    negative binomial dispersion coefficients computed by DESeq2.
#' @param vecNormConst: (numeric vector number of replicates) 
#'    Normalisation constants for each replicate.
#' @param dfAnnotationFull (Table) Lists co-variables of individual replicates:
#'    Replicate, Sample, Condition, Time. Time must be numeric.
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

fitImpulse <- function(arr2DCountData, vecDispersions=NULL, 
  vecNormConst, dfAnnotationFull,
  strCaseName=NULL, strControlName=NULL, strMode="batch",
  nProcessesAssigned=3, NPARAM=6){
  
  #g_names = rownames(arr3DCountData)
  lsFitResults_all = list()
  
  if(!is.null(strControlName)){
    lsLabels <- c("combined","case","control")
    names(lsLabels) <- c("combined","case","control")
  } else {
    lsLabels <- c("case")
    names(lsLabels) <- c("case")
  }
  
  # Perform 3 runs (combined, case and control) if control timecourse is present
  if(!is.null(strControlName)){    
    lsReplicatesByCond <- list(
      dfAnnotationFull$Replicate,
      dfAnnotationFull[dfAnnotationFull$Condition %in% strCaseName,]$Replicate,
      dfAnnotationFull[dfAnnotationFull$Condition %in% strControlName,]$Replicate )
    names(lsReplicatesByCond) <- c("combined","case","control")   
  } else {
    lsReplicatesByCond <- list(
      dfAnnotationFull$Replicate )
    names(lsReplicatesByCond) <- c("case")
  }
  if(strMode=="timecourses"){
    vecTimecourseAssign <- dfAnnotationFull[match(colnames(arr2DCountData),dfAnnotationFull$Replicate),]$Timecourse
  } else {
    vecTimecourseAssign <- NULL
  }
  
  # Fitting for different runs
  for (label in lsLabels){
    lsFitResults_run <- fitImpulse_matrix(
      arr2DCountDataCondition=arr2DCountData[,lsReplicatesByCond[[label]]],
      vecNormConst=vecNormConst[lsReplicatesByCond[[label]]],
      vecDispersions=vecDispersions,
      vecTimepointAssign=dfAnnotationFull[match(colnames(arr2DCountData[,lsReplicatesByCond[[label]]]),dfAnnotationFull$Replicate),]$Time,
      vecTimecourseAssign=vecTimecourseAssign,
      dfAnnotationFull=dfAnnotationFull,
      strCaseName=strCaseName,
      strControlName=strControlName,
      strMode=strMode,
      nProcessesAssigned=nProcessesAssigned,
      NPARAM=NPARAM )
    
    # DAVID: do i still need this?
    #imp_res$impulse_parameters <- (imp_res$impulse_parameters)[g_names,]
    #imp_res$impulse_fits <- (imp_res$impulse_fits)[g_names,]
    
    names(lsFitResults_run) <- paste(c("parameters","values"),lsLabels[label],sep="_")
    lsFitResults_all[[names(lsFitResults_run)[1]]] <- lsFitResults_run[[1]]
    lsFitResults_all[[names(lsFitResults_run)[2]]] <- lsFitResults_run[[2]]
  }
  
  return(lsFitResults_all)
}

