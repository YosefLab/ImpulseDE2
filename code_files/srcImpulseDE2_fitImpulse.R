#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Impulse model fit    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#findInflextionPoints()
#estimateParamPeak()
#estimateParamValley()
#computeLogLikNull()
#computeTranslationFactors()

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
#' @param vecCounts (count vector number of replicates) Count data.
#' @param scaDispersionEstimate (scalar) Inverse of negative binomial 
#'    dispersion coefficients computed by DESeq2 for given gene.
#' @param vecDropoutRate: (probability vector number of replicates) 
#'    [Default NULL] Dropout rate/mixing probability of zero inflated 
#'    negative binomial mixturemodel for each gene and cell.
#' @param vecProbNB: (probability vector number of replicates) 
#'    Probability of observations to come from negative binomial 
#'    component of mixture model.
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

fitImpulse_gene <- function(vecCounts, 
  scaDispersionEstimate,
  vecDropoutRate=NULL, vecProbNB=NULL,
  vecNormConst,
  vecTimepointAssign, vecTimecourseAssign, 
  dfAnnotationFull,
  strMode="batch", NPARAM=6, MAXIT=1000){

  optim_method <- "optim"
  #optim_method <- "nlminb"
  #optim_method <- c("optim","nlminb")
  
  # Get boolean observation vectors:
  vecboolObserved <- !is.na(vecCounts)
  vecboolZero <- vecCounts==0
  vecboolNotZeroObserved <- !is.na(vecCounts) & (vecCounts!=0)
  
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
  if(strMode=="batch"){
    # Null model (no longitudinal sampling): 
    # Fit one mean for all replicates.
    
    # Fitting null model:
    # Fit mean to normalised counts. Count normalisation 
    # corresponds to model normalisation in impulse model 
    # fitting.
    scaMu <- mean(vecCounts/vecNormConst, na.rm=TRUE)
    
    # Evaluate likelihood of null model:
    # Scale null model according to normalisation factors.
    scaLogLikNull <- sum(dnbinom(
      vecCounts[vecboolObserved], 
      mu=scaMu*vecNormConst[vecboolObserved], 
      size=scaDispersionEstimate, 
      log=TRUE))
  } else if(strMode=="singlecell"){
    # Null model (no longitudinal sampling):
    # Single cells are independent observations
    # of the time series.
    # Fit one mean for all replicates, with mixture model.
    # Note on data scaling: Normalisation factors were
    # computed based on raw data.
    
    # Fit null model:
    # Fit weighted mean to count data: Weight observations
    # by their probability to belong to the negative
    # binomial component of the mixture model.
    scaMu <- sum(vecCounts/vecNormConst*vecProbNB, na.rm=TRUE)/sum(vecProbNB)
    
    # Evaluate likelihood of null model
    # Likelihood of zero counts:
    vecLikNullZeros <- (1-vecDropoutRate[vecboolZero])*
        dnbinom( vecCounts[vecboolZero], 
          mu=scaMu, 
          size=scaDispersionEstimate, 
          log=FALSE) +
      vecDropoutRate[vecboolZero]
    # Replace zero likelihood observation with machine precision
    # for taking log.
    scaLogLikNullZeros <- sum( log(vecLikNullZeros[vecLikNullZeros!=0]) +
        sum(vecLikNullZeros==0)*log(.Machine$double.eps) )
    # Likelihood of non-zero counts:
    vecLikNullNonzeros <- (1-vecDropoutRate[vecboolNotZeroObserved])*
      dnbinom(vecCounts[vecboolNotZeroObserved], 
        mu=scaMu, 
        size=scaDispersionEstimate, 
        log=FALSE)
    # Replace zero likelihood observation with machine precision
    # for taking log.
    scaLogLikNullNonzeros <- sum( log(vecLikNullNonzeros[vecLikNullNonzeros!=0]) +
        sum(vecLikNullNonzeros==0)*log(.Machine$double.eps) )
    
    # Compute likelihood of all data:
    scaLogLikNull <- scaLogLikNullZeros + scaLogLikNullNonzeros
  } else if(strMode=="timecourses"){
    # Null model (longitudinal sampling): 
    # Fit one mean for each time course (longitudinal series).
    
    # Fit null model.
    # Fit mean to normalised counts. Count normalisation 
    # corresponds to model normalisation in impulse model 
    # fitting.
    scaMu <- mean(vecCounts/vecNormConst, na.rm=TRUE)
    vecMuTimecourses <- sapply(vecTimecourses,
      function(tc){mean(vecCounts[vecTimecourseAssign==tc], na.rm=TRUE)})
    names(vecMuTimecourses) <- vecTimecourses
    
    # Evaluate likelihood of null model
    scaLogLikNull <- sum(dnbinom(
      vecCounts[vecboolObserved], 
      mu=(vecMuTimecourses[as.vector(vecTimecourseAssign)])[vecboolObserved]*vecNormConst[vecboolObserved], 
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
    
    # Scaled mean ratio per replicate
    vecTranslationFactors <- vecMuTimecourses[as.vector(vecTimecourseAssign)]/scaMu
    # Note: If all replicates are zero, scaMu is zero and
    # vecTranslationFactors are NA. Genes only with zero counts
    # are removed in processData(). The input data to this function
    # ma still consist entirely of zeros if this function is called
    # on the subset of case or control data (if control condition
    # is given). Those subsets may be only zeros. This exception
    # is caught here and vecTranslationFactors set to 1 from NA.
    # This removes the effect of vecTranslationFactors on fitting.
    # The negative binomial density can still be
    # evaluated as it is initialised with values from an impulse model
    # padded with zero counts.s 
    if(scaMu==0){
      vecTranslationFactors <- array(1,length(vecTranslationFactors))
    }
  }
  # Compute statistics for initialisation:
  # Expression means by timepoint
  if(strMode=="batch" | strMode=="timecourses"){
    vecExpressionMeans <- sapply(vecTimepoints,
      function(tp){mean(vecCounts[vecTimepointAssign==tp], na.rm=TRUE)})
  } else if(strMode=="singlecell"){
    vecExpressionMeans <- unlist(sapply( vecTimepoints,
      function(tp){sum((vecCounts/vecNormConst*vecProbNB)[vecTimepointAssign==tp], na.rm=TRUE)/sum(vecProbNB[vecTimepointAssign==tp])} ))
  } else {
    stop(paste0("ERROR: Unrecognised strMode in fitImpulse(): ",strMode))
  }
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
        vecDropoutRateEst=vecDropoutRate,
        vecNormConst=vecNormConst,
        vecindTimepointAssign=vecindTimepointAssign,
        vecboolNotZeroObserved=vecboolNotZeroObserved, 
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
        vecDropoutRateEst=vecDropoutRate,
        vecNormConst=vecNormConst,
        vecindTimepointAssign=vecindTimepointAssign,
        vecboolNotZeroObserved=vecboolNotZeroObserved,
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
    lsBestFitSummary <- c(dfFitsByInitialisation[,indBestFit],scaMu,scaLogLikNull)
    names(lsBestFitSummary) <- c("beta","h0","h1","h2","t1","t2",
      "logL_H1","converge_H1","mu","logL_H0")
  } else if(strMode=="timecourses") {
    vecTranslationFactorsUnique <- vecTranslationFactors[as.vector(vecTimecourses)]
    lsBestFitSummary <- c(dfFitsByInitialisation[,indBestFit],
      scaMu,
      scaLogLikNull,
      vecMuTimecourses,
      vecTranslationFactorsUnique)
    vecColnamesMubyTimecourse <- paste0(rep("mu_",length(vecMuTimecourses)), names(vecMuTimecourses))
    vecColnameTranslationFactors <- paste0(rep("TranslationFac_",length(vecTranslationFactorsUnique)), names(vecTranslationFactorsUnique))
    names(lsBestFitSummary) <- c("beta","h0","h1","h2","t1","t2",
      "logL_H1","converge_H1",
      "mu",
      "logL_H0",
      vecColnamesMubyTimecourse,
      vecColnameTranslationFactors)
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
#' @param arr2DCountData: (count matrix genes x replicates) Count data.
#' @param vecDispersions: (vector number of genes) Inverse of gene-wise 
#'    negative binomial dispersion coefficients computed by DESeq2.
#' @param matDropoutRate: (probability matrix genes x replicates) Dropout 
#'    rate/mixing probability of zero inflated negative binomial mixture 
#'    model for each gene and cell.
#' @param matProbNB: (probability vector genes x replicates) 
#'    Probability of observations to come from negative binomial 
#'    component of mixture model.
#' @param vecNormConst: (numeric vector number of replicates) 
#'    Normalisation constants for each replicate.
#' @param vecTimepointAssign: (numeric vector number replicates) Timepoints 
#'    assigned to replicates.
#' @param vecTimecourseAssign: (numeric vector number replicates) Time courses 
#'    assigned to replicates. NULL if not operating in strMode="timecourses".
#' @param dfAnnotationFull: (Table) Lists co-variables of individual replicates:
#'    Replicate, Sample, Condition, Time. Time must be numeric.
#' @param strCaseName: (str) Name of the case condition in \code{dfAnnotationRedFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#'    \code{dfAnnotationRedFull}.
#' @param nProcessesAssigned: (scalar) [Default 3] Number of processes for parallelisation. The
#'    specified value is internally changed to \code{min(detectCores() - 1, nProc)} 
#'    using the \code{detectCores} function from the package \code{parallel} to 
#'    avoid overload.
#' @param NPARAM: (scalar) [Default 6] Number of parameters of impulse model.
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

fitImpulse_matrix <- function(arr2DCountDataCondition, 
  vecDispersions, 
  matDropoutRate=NULL, matProbNB=NULL,
  vecNormConst, 
  vecTimepointAssign, vecTimecourseAssign, 
  dfAnnotationFull,
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
    assign("matDropoutRate", matDropoutRate, envir = my.env)
    assign("matProbNB", matProbNB, envir = my.env)
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
      "matDropoutRate",
      "matProbNB",
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
    if(strMode=="singlecell"){
      # Call fitting with single cell parameters
      lsmatFits <- clusterApply(cl, 1:length(lsGeneIndexByCore),
        function(z){ t(sapply( lsGeneIndexByCore[[z]],
          function(x){fitImpulse_gene(
            vecCounts=arr2DCountDataCondition[x,],
            scaDispersionEstimate=vecDispersions[x],
            vecDropoutRate=matDropoutRate[x,],
            vecProbNB=matProbNB[x,],
            vecNormConst=vecNormConst,
            vecTimepointAssign=vecTimepointAssign,
            vecTimecourseAssign=vecTimecourseAssign,
            dfAnnotationFull=dfAnnotationFull,
            strMode=strMode,
            NPARAM=NPARAM,
            MAXIT=MAXIT )}
        ))})
    } else {
      # Call fitting without single cell parameters
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
    }
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
    if(strMode=="singlecell"){
      # Call fitting with single cell parameters
      matFits <- lapply(lsGeneIndexByCore_all,function(x){
        fitImpulse_gene(
          vecCounts=arr2DCountDataCondition[x,],
          scaDispersionEstimate=vecDispersions[x],
          vecDropoutRate=matDropoutRate[x,],
          vecProbNB=matProbNB[x,],
          vecNormConst=vecNormConst,
          vecTimepointAssign=vecTimepointAssign,
          vecTimecourseAssign=vecTimecourseAssign,
          dfAnnotationFull=dfAnnotationFull,
          strMode=strMode,
          NPARAM=NPARAM,
          MAXIT=MAXIT )})
    } else {
      # Call fitting without single cell parameters
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
    }
    # Call fitting without single cell parameters
    matFits <- do.call(rbind,matFits)
    rownames(matFits) <- rownames(arr2DCountDataCondition)
  }
  
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
#' @param arr2DCountData: (count matrix genes x replicates) Count data.
#' @param vecDispersions: (vector number of genes) Inverse of gene-wise 
#'    negative binomial dispersion coefficients computed by DESeq2.
#' @param matDropoutRate: (probability matrix genes x replicates) Dropout 
#'    rate/mixing probability of zero inflated negative binomial mixture 
#'    model for each gene and cell.
#' @param matProbNB: (probability vector genes x replicates) 
#'    Probability of observations to come from negative binomial 
#'    component of mixture model.
#' @param vecNormConst: (numeric vector number of replicates) 
#'    Normalisation constants for each replicate.
#' @param dfAnnotationFull (Table) Lists co-variables of individual replicates:
#'    Replicate, Sample, Condition, Time. Time must be numeric.
#' @param strCaseName (str) Name of the case condition in \code{dfAnnotationRedFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#'    \code{dfAnnotationRedFull}.
#' @param strMode: (str) [Default "batch"] {"batch","timecourses","singlecell"}
#'    Mode of model fitting.
#' @param nProcessesAssigned: (scalar) [Default 1] Number of processes for parallelisation. The
#'    specified value is internally changed to \code{min(detectCores() - 1, nProc)} 
#'    using the \code{detectCores} function from the package \code{parallel} to 
#'    avoid overload.
#' @param NPARAM: (scalar) [Default 6] Number of parameters of impulse model.
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

fitImpulse <- function(arr2DCountData, 
  vecDispersions, 
  matDropoutRate, matProbNB,
  vecNormConst, dfAnnotationFull,
  strCaseName, strControlName=NULL, strMode="batch",
  nProcessesAssigned=3, NPARAM=6){
  
  #g_names = rownames(arr3DCountData)
  lsFitResults_all = list()
  
  # Condition labels to be used in runs:
  # These labels are used to group replicates
  if(!is.null(strControlName)){
    lsLabels <- c("combined","case","control")
    names(lsLabels) <- c("combined","case","control")
  } else {
    lsLabels <- c("case")
    names(lsLabels) <- c("case")
  }
  
  # Create lists of replicates to be used per run
  if(!is.null(strControlName)){    
    lsReplicatesByCond <- list(
      colnames(arr2DCountData),
      (colnames(arr2DCountData))[dfAnnotationFull[match(colnames(arr2DCountData),dfAnnotationFull$Replicate),]$Condition %in% strCaseName],
      (colnames(arr2DCountData))[dfAnnotationFull[match(colnames(arr2DCountData),dfAnnotationFull$Replicate),]$Condition %in% strControlName] )
    names(lsReplicatesByCond) <- c("combined","case","control")   
  } else {
    lsReplicatesByCond <- list(
      colnames(arr2DCountData) )
    names(lsReplicatesByCond) <- c("case")
  }
  # Get time course assignments of replicates
  if(strMode=="timecourses"){
    vecTimecourseAssign <- dfAnnotationFull[match(colnames(arr2DCountData),dfAnnotationFull$Replicate),]$Timecourse
    names(vecTimecourseAssign) <- colnames(arr2DCountData)
  } else {
    vecTimecourseAssign <- NULL
  }
  # Get time point assignments of replicates
  vecTimepointAssign <- dfAnnotationFull[match(colnames(arr2DCountData),dfAnnotationFull$Replicate),]$Time
  names(vecTimepointAssign) <- colnames(arr2DCountData)
  
  # Fitting for different runs
  for (label in lsLabels){
    if(strMode=="singlecell"){
      # Call fitting with single cell parameters
      lsFitResults_run <- fitImpulse_matrix(
        arr2DCountDataCondition=arr2DCountData[,lsReplicatesByCond[[label]]],
        vecNormConst=vecNormConst[lsReplicatesByCond[[label]]],
        vecDispersions=vecDispersions,
        matDropoutRate=matDropoutRate[,lsReplicatesByCond[[label]]],
        matProbNB=matProbNB[,lsReplicatesByCond[[label]]],
        vecTimepointAssign=vecTimepointAssign[lsReplicatesByCond[[label]]],
        vecTimecourseAssign=vecTimecourseAssign[lsReplicatesByCond[[label]]],
        dfAnnotationFull=dfAnnotationFull,
        strCaseName=strCaseName,
        strControlName=strControlName,
        strMode=strMode,
        nProcessesAssigned=nProcessesAssigned,
        NPARAM=NPARAM )
    } else {
      # Call fitting without singlecell parameters
      lsFitResults_run <- fitImpulse_matrix(
        arr2DCountDataCondition=arr2DCountData[,lsReplicatesByCond[[label]]],
        vecNormConst=vecNormConst[lsReplicatesByCond[[label]]],
        vecDispersions=vecDispersions,
        vecTimepointAssign=vecTimepointAssign[lsReplicatesByCond[[label]]],
        vecTimecourseAssign=vecTimecourseAssign[lsReplicatesByCond[[label]]],
        dfAnnotationFull=dfAnnotationFull,
        strCaseName=strCaseName,
        strControlName=strControlName,
        strMode=strMode,
        nProcessesAssigned=nProcessesAssigned,
        NPARAM=NPARAM )
    }
    
    # DAVID: do i still need this?
    #imp_res$impulse_parameters <- (imp_res$impulse_parameters)[g_names,]
    #imp_res$impulse_fits <- (imp_res$impulse_fits)[g_names,]
    
    names(lsFitResults_run) <- paste(c("parameters","values"),lsLabels[label],sep="_")
    lsFitResults_all[[names(lsFitResults_run)[1]]] <- lsFitResults_run[[1]]
    lsFitResults_all[[names(lsFitResults_run)[2]]] <- lsFitResults_run[[2]]
  }
  
  return(lsFitResults_all)
}

