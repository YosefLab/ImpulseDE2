### Impulse and constant model fitting

#' Compute peak and valley impulse model parameter initialisations
#' for data of one gene
#' 
#' [Model fitting function hierarchy: helper to level 3 out of 4]
#' This is a fitting helper function which computes parameter intialisations
#' and does not wrap or execute numerical optimisation.
#' The peak model models a maximum between start and end time,
#' the valley model models a minimum between start and end time.
#' 
#' @seealso Called by \link{fitConstImpulseGene}.
#' 
#' @param vecCounts (numeric vector number of samples)
#' Read count data.
#' @param vecTimepoints (numeric vector length number of samples)
#' Time coordinates of each sample.
#' @param vecSizeFactors (numeric vector number of samples) 
#' Model scaling factors for each sample which take
#' sequencing depth into account (size factors).
#' @param lsvecidxBatch (list length number of confounding variables)
#' List of index vectors. 
#' One vector per confounding variable.
#' Each vector has one entry per sample with the index batch
#' within the given confounding variable of the given sample.
#' Batches are enumerated from 1 to number of batches.
#' 
#' @return (list length 2)
#' \itemize{
#' \item peak (numeric vector length 6)
#' \{beta, h0, h1, h2, t1, t2\}
#' Peak model initialisations of impulse model parameters.
#' \item valley (numeric vector length 6)
#' \{beta, h0, h1, h2, t1, t2\}
#' Valley model initialisations of impulse model parameters.
#' }
#' 
#' @author David Sebastian Fischer
estimateImpulseParam <- function(vecCounts, vecTimepoints, vecSizeFactors, 
                                 lsvecidxBatch) {
    
    # Compute general statistics for initialisation:
    vecTimepointsUnique <- unique(vecTimepoints)
    scaMeanCount <- mean(vecCounts, na.rm = TRUE)
    # Expression means by timepoint
    vecCountsSFcorrected <- vecCounts/vecSizeFactors
    vecCountsSFcorrectedNorm <- vecCountsSFcorrected / scaMeanCount
    if (!is.null(lsvecidxBatch)) {
        # Estimate batch factors
        vecBatchFactors <- array(1, length(vecCounts))
        for (vecidxBatch in lsvecidxBatch) {
            vecBatchFactorsConfounder <- tapply(
                vecCountsSFcorrectedNorm, 
                vecidxBatch, 
                mean, na.rm=TRUE)
            # Catch exception that all observations of a batch are zero or all
            # observations are zero:
            vecBatchFactorsConfounder[is.na(vecBatchFactorsConfounder) | 
                                          vecBatchFactorsConfounder == 0] <- 1
            vecBatchFactors <- vecBatchFactors * 
                vecBatchFactorsConfounder[vecidxBatch]
        }
        vecCountsSFBatchcorrected <- vecCountsSFcorrected/vecBatchFactors
        vecExpressionMeans <- tapply(
            vecCountsSFBatchcorrected, 
            match(vecTimepoints,vecTimepointsUnique), 
            mean, na.rm=TRUE)
    } else {
        vecExpressionMeans <- tapply(
            vecCountsSFcorrected, 
            match(vecTimepoints,vecTimepointsUnique), 
            mean, na.rm=TRUE)
    }
    scaNTimepoints <- length(vecTimepointsUnique)
    vecidxMiddle <- seq(2, scaNTimepoints-1, by=1)
    scaMaxMiddleMean <- max(vecExpressionMeans[vecidxMiddle], 
                            na.rm = TRUE)
    scaMinMiddleMean <- min(vecExpressionMeans[vecidxMiddle], 
                            na.rm = TRUE)
    # +1 to push indicices up from middle stretch to entire window (first is
    # omitted here)
    indMaxMiddleMean <- match(scaMaxMiddleMean, 
                              vecExpressionMeans[vecidxMiddle]) + 1
    indMinMiddleMean <- match(scaMinMiddleMean, 
                              vecExpressionMeans[vecidxMiddle]) + 1
    # Gradients between neighbouring points
    vecDiffExpr <- diff(vecExpressionMeans)
    vecDiffTime <- diff(vecTimepointsUnique)
    vecGradients <- vecDiffExpr / vecDiffTime
    vecGradients[is.na(vecGradients) | !is.finite(vecGradients)] <- 0
    
    # Compute peak initialisation Beta: Has to be negative, Theta1: Low,
    # Theta2: High, Theta3: Low t1: Around first observed inflexion point,
    # t2: Around second observed inflexion point
    vecdidxFirstPart <- seq(1, indMaxMiddleMean-1, by=1)
    vecdidxSecndPart <- seq(indMaxMiddleMean, length(vecGradients), by=1)
    indLowerInflexionPoint <- match(max(
        vecGradients[vecdidxFirstPart], na.rm = TRUE), 
        vecGradients[vecdidxFirstPart])
    indUpperInflexionPoint <- indMaxMiddleMean - 1 + 
        match(min(
            vecGradients[vecdidxSecndPart], na.rm = TRUE), 
            vecGradients[vecdidxSecndPart])
    vecParamGuessPeak <- c(
        1, 
        log(vecExpressionMeans[1] + 1), 
        log(scaMaxMiddleMean + 1), 
        log(vecExpressionMeans[scaNTimepoints] + 1), 
        (vecTimepointsUnique[indLowerInflexionPoint] + 
             vecTimepointsUnique[indLowerInflexionPoint + 1])/2, 
        (vecTimepointsUnique[indUpperInflexionPoint] + 
             vecTimepointsUnique[indUpperInflexionPoint + 1])/2
    )
    
    # Compute valley initialisation Beta: Has to be negative, Theta1: High,
    # Theta2: Low, Theta3: High t1: Around first observed inflexion point,
    # t2: Around second observed inflexion point
    indLowerInflexionPoint <- match(min(
        vecGradients[1:(indMinMiddleMean - 1)], na.rm = TRUE), 
        vecGradients[1:(indMinMiddleMean - 1)])
    indUpperInflexionPoint <- indMinMiddleMean - 1 + match(max(
        vecGradients[indMinMiddleMean:(scaNTimepoints - 1)], na.rm = TRUE), 
        vecGradients[indMinMiddleMean:(scaNTimepoints - 1)])
    vecParamGuessValley <- c(
        1, 
        log(vecExpressionMeans[1] + 1), 
        log(scaMinMiddleMean + 1), 
        log(vecExpressionMeans[scaNTimepoints] + 1), 
        (vecTimepointsUnique[indLowerInflexionPoint] + 
             vecTimepointsUnique[indLowerInflexionPoint + 1])/2, 
        (vecTimepointsUnique[indUpperInflexionPoint] + 
             vecTimepointsUnique[indUpperInflexionPoint + 1])/2
    )
    
    return(list(peak = vecParamGuessPeak, valley = vecParamGuessValley))
}

#' Fit a constant model to data of a gene
#' 
#' [Model fitting function hierarchy: 4 out of 4]
#' This quarterny fitting wrapper performs constant model fitting:
#' This function executes numerical optimisaiton and error-handling
#' thereof.
#' 
#' @seealso Called by \link{fitConstImpulseGene} to fit constant
#' model to samples of one condition and one gene.
#' Calls constant model cost function 
#' \link{evalLogLikMu} within \link{optim}.
#' 
#' @param vecCounts (numeric vector number of samples)
#' Read count data.
#' @param scaDisp (scalar) Gene-wise 
#' negative binomial dispersion hyper-parameter.
#' @param vecSizeFactors (numeric vector number of samples) 
#' Model scaling factors for each sample which take
#' sequencing depth into account (size factors).
#' @param lsvecidxBatch (list length number of confounding variables)
#' List of index vectors. 
#' One vector per confounding variable.
#' Each vector has one entry per sample with the index batch
#' within the given confounding variable of the given sample.
#' Batches are enumerated from 1 to number of batches.
#' @param MAXIT (scalar) [Default 1000] 
#' Maximum number of BFGS iterations for model fitting with \link{optim}.
#' @param RELTOL (scalar) [Default 10^(-8)]
#' Maximum relative change in loglikelihood to reach convergence in
#' numerical optimisation by BFGS in \link{optim}.
#' @param trace (scalar) [Defaul 0]
#' Reporting parameter of \link{optim}.
#' @param REPORT (scalar) [Default 10]
#' Reporting parameter of \link{optim}.
#' 
#' @return (list) List of constant fit parameters and results.
#' \itemize{
#' \item scaMu (scalar) Maximum likelihood estimator of
#' negative binomial mean parameter.
#' \item lsvecBatchFactors (list length number of confounders)
#' List of vectors of scalar batch correction factors for each sample.
#' These are also maximum likelihood estimators.
#' NULL if no confounders given.
#' \item scaDispParam (scalar) Dispersion parameter estimate
#' used in fitting (hyper-parameter).
#' \item scaLL (scalar) Loglikelihood of data under maximum likelihood
#' estimator model.
#' \item scaConvergence (scalar) 
#' Convergence status of optim on constant model.
#' }
#' 
#' @author David Sebastian Fischer
fitConstModel <- function(
    vecCounts, scaDisp, vecSizeFactors, lsvecidxBatch, 
    MAXIT = 1000, RELTOL = 10^(-8), trace = 0, REPORT = 10) {
    
    vecParamGuess <- log(mean(vecCounts, na.rm = TRUE) + 1)
    if (!is.null(lsvecidxBatch)) {
        for (vecidxConfounder in lsvecidxBatch) {
            vecParamGuess <- c(vecParamGuess, 
                               rep(0, length(unique(vecidxConfounder)) - 1))
        }
    }
    
    lsFit <- tryCatch({
        optim(par = vecParamGuess, fn = evalLogLikMu_comp, 
              vecCounts = vecCounts, scaDisp = scaDisp, 
              vecSizeFactors = vecSizeFactors, 
              lsvecidxBatch = lsvecidxBatch, 
              vecboolObserved = !is.na(vecCounts), method = "BFGS", 
              control = list(maxit = MAXIT, reltol = RELTOL, fnscale = -1)
        )[c("par", "value", "convergence")]
    }, error = function(strErrorMsg) {
        print(paste0("ERROR: Fitting null model: fitConstModel().", 
                     " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
        print(paste0("vecParamGuess ", 
                     paste(vecParamGuess, collapse = " ")))
        print(paste0("vecCounts ", 
                     paste(vecCounts, collapse = " ")))
        print(paste0("scaDisp ", 
                     paste(scaDisp, collapse = " ")))
        print(paste0("vecSizeFactors ", 
                     paste(vecSizeFactors, collapse = " ")))
        print(paste0("lsvecidxBatch ", 
                     paste(lsvecidxBatch, collapse = " ")))
        print(paste0("MAXIT ", MAXIT))
        print(strErrorMsg)
        stop(strErrorMsg)
    })
    
    # Extract parameter estimates
    scaMu <- exp(lsFit$par[1])
    # Catch boundary of likelihood domain on mu space:
    if (scaMu < 10^(-10)) {
        scaMu <- 10^(-10)
    }
    scaNParamUsed <- 1
    if (!is.null(lsvecidxBatch)) {
        lsvecBatchFactors <- list()
        for(i in seq(1,length(lsvecidxBatch))) {
            vecidxConfounder <- lsvecidxBatch[[i]]
            scaNBatchFactors <- max(vecidxConfounder) - 1  
            # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining 
            # factors scale based on the first batch.
            vecBatchFactors <- c(1, exp(lsFit$par[
                (scaNParamUsed + 1):
                    (scaNParamUsed + scaNBatchFactors)] ))
            scaNParamUsed <- scaNParamUsed + scaNBatchFactors
            # Catch boundary of likelihood domain on batch factor space:
            vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
            vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
            lsvecBatchFactors[[i]] <- vecBatchFactors
        }
    } else {
        lsvecBatchFactors <- NULL
    }
    
    return(list(scaMu = scaMu, lsvecBatchFactors = lsvecBatchFactors, 
                scaDispParam = scaDisp, 
                scaLL = lsFit$value, scaConvergence = lsFit$convergence))
}

#' Fit an impulse model to data of a gene
#' 
#' [Model fitting function hierarchy: 4 out of 4]
#' This quarterny fitting wrapper performs impulse model fitting:
#' This function executes numerical optimisaiton and error-handling
#' thereof.
#' 
#' @seealso Called by \link{fitConstImpulseGene} to fit impulse
#' model to samples of one condition and one gene.
#' Calls impulse model cost function 
#' \link{evalLogLikImpulse_comp} within \link{optim}.
#' 
#' @param vecImpulseParamGuess (numeric vector length 6)
#' \{beta, h0, h1, h2, t1, t2\}
#' Initialisations of impulse model parameters.
#' @param vecCounts (numeric vector number of samples)
#' Read count data.
#' @param scaDisp (scalar) Gene-wise 
#' negative binomial dispersion hyper-parameter.
#' @param vecSizeFactors (numeric vector number of samples) 
#' Model scaling factors for each sample which take
#' sequencing depth into account (size factors).
#' @param lsvecidxBatch (list length number of confounding variables)
#' List of index vectors. 
#' One vector per confounding variable.
#' Each vector has one entry per sample with the index batch
#' within the given confounding variable of the given sample.
#' Batches are enumerated from 1 to number of batches.
#' @param vecTimepointsUnique
#' (numeric vector length number of unique time points)
#' Unique time points of set of time points of given samples.
#' @param vecidxTimepoint (index vector length number of samples)
#' Index of of time point assigned to each sample in vector
#' vecTimepointsUnique.
#' @param MAXIT (scalar) [Default 1000] 
#' Maximum number of BFGS iterations for model fitting with \link{optim}.
#' @param RELTOL (scalar) [Default 10^(-8)]
#' Maximum relative change in loglikelihood to reach convergence in
#' numerical optimisation by BFGS in \link{optim}.
#' @param trace (scalar) [Defaul 0]
#' Reporting parameter of \link{optim}.
#' @param REPORT (scalar) [Default 10]
#' Reporting parameter of \link{optim}.
#' 
#' @return (list) List of impulse fit parameters and results.
#' \itemize{
#' \item vecImpulseParam (numeric vector length 6)
#' \{beta, h0, h1, h2, t1, t2\}
#' Maximum likelihood estimators of impulse model parameters.
#' \item vecImpulseValue (numeric vector length number of time points)
#' Values of impulse model fit at time points used for fit.
#' \item lsvecBatchFactors (list length number of confounders)
#' List of vectors of scalar batch correction factors for each sample.
#' These are also maximum likelihood estimators.
#' NULL if no confounders given.
#' \item scaDispParam (scalar) Dispersion parameter estimate
#' used in fitting (hyper-parameter).
#' \item scaLL (scalar) Loglikelihood of data under maximum likelihood
#' estimator model.
#' \item scaConvergence (scalar) 
#' Convergence status of optim on impulse model.
#' }
#' 
#' @author David Sebastian Fischer
fitImpulseModel <- function(
    vecImpulseParamGuess, vecCounts, scaDisp, vecSizeFactors, 
    lsvecidxBatch, vecTimepointsUnique, vecidxTimepoint, 
    MAXIT = 1000, RELTOL = 10^(-8), trace = 0, REPORT = 10) {
    
    
    vecParamGuess <- vecImpulseParamGuess
    if (!is.null(lsvecidxBatch)) {
        for (vecidxConfounder in lsvecidxBatch) {
            vecParamGuess <- c(
                vecParamGuess, 
                rep(0, length(unique(vecidxConfounder)) - 1))
        }
    }
    
    lsFit <- tryCatch({
        optim(par = vecParamGuess, fn = evalLogLikImpulse_comp, 
              vecCounts = vecCounts, scaDisp = scaDisp, 
              vecSizeFactors = vecSizeFactors, 
              vecTimepointsUnique = vecTimepointsUnique, 
              vecidxTimepoint = vecidxTimepoint, lsvecidxBatch = lsvecidxBatch, 
              vecboolObserved = !is.na(vecCounts), method = "BFGS", 
              control = list(maxit = MAXIT, reltol = RELTOL, fnscale = -1)
        )[c("par", "value", "convergence")]
    }, error = function(strErrorMsg) {
        print(paste0("ERROR: Fitting impulse model: fitImpulseModel().", 
                     " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
        print(paste0("vecParamGuess ", paste(vecParamGuess, collapse = " ")))
        print(paste0("vecCounts ", paste(vecCounts, collapse = " ")))
        print(paste0("scaDisp ", paste(scaDisp, collapse = " ")))
        print(paste0("vecSizeFactors ", 
                     paste(vecSizeFactors, collapse = " ")))
        print(paste0("vecTimepointsUnique ", 
                     paste(vecTimepointsUnique, collapse = " ")))
        print(paste0("vecidxTimepoint ", 
                     paste(vecidxTimepoint, collapse = " ")))
        print(paste0("lsvecidxBatch ", 
                     paste(lsvecidxBatch, collapse = " ")))
        print(paste0("MAXIT ", MAXIT))
        print(strErrorMsg)
        stop(strErrorMsg)
    })
    
    # Extract parameter estimates
    vecImpulseParam <- lsFit$par[1:6]
    vecImpulseParam[2:4] <- exp(vecImpulseParam[2:4])
    names(vecImpulseParam) <- c("beta", "h0", "h1", "h2", "t1", "t2")
    vecImpulseValue <- evalImpulse_comp(
        vecImpulseParam = vecImpulseParam, 
        vecTimepoints = vecTimepointsUnique)[vecidxTimepoint]
    names(vecImpulseValue) <- names(vecCounts)
    scaNParamUsed <- 6
    if (!is.null(lsvecidxBatch)) {
        lsvecBatchFactors <- list()
        for(i in seq(1,length(lsvecidxBatch))) {
            vecidxConfounder <- lsvecidxBatch[[i]]
            scaNBatchFactors <- max(vecidxConfounder) - 1  
            # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining 
            # factors scale based on the first batch.
            vecBatchFactors <- c(1, exp(lsFit$par[
                (scaNParamUsed + 1):
                    (scaNParamUsed + scaNBatchFactors)] ))
            scaNParamUsed <- scaNParamUsed + scaNBatchFactors
            # Catch boundary of likelihood domain on batch factor space:
            vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
            vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
            lsvecBatchFactors[[i]] <- vecBatchFactors
        }
    } else {
        lsvecBatchFactors <- NULL
    }
    
    return(list(vecImpulseParam = vecImpulseParam, 
                vecImpulseValue = vecImpulseValue, 
                lsvecBatchFactors = lsvecBatchFactors, 
                scaDispParam = scaDisp, scaLL = lsFit$value, 
                scaConvergence = lsFit$convergence))
}

#' Fit an impulse and constant model to a single gene
#' 
#' [Model fitting function hierarchy: 3 out of 4]
#' This tertiary fitting wrapper calls the optimisation wrappers
#' for the individual fitting operations to be performed on the 
#' observations of this gene.
#' Structure of this function:
#' \itemize{
#' \item Fit impulse model
#' \itemize{
#' \item Initialise impulse model parameters (peak and valley)
#' \item Fit impulse model (peak initialisation)
#' \item Fit impulse model (valley initialisation)
#' }
#' \item Select best impulse model fit from initialisations,
#' \item Fit constant model (if constant model is to be fit).
#' }
#' 
#' @seealso Called by \link{fitConstImpulseGene} to fit constant and impulse
#' model to samples of one condition and one gene.
#' Calls impulse parameter initialisation function
#' \link{estimateImpulseParam} and 
#' optimisation wrappers 
#' \link{fitImpulseModel} and \link{fitConstModel}.
#' 
#' @param vecCounts (numeric vector number of samples)
#' Read count data.
#' @param scaDisp (scalar) Gene-wise 
#' negative binomial dispersion hyper-parameter.
#' @param vecSizeFactors (numeric vector number of samples) 
#' Model scaling factors for each sample which take
#' sequencing depth into account (size factors).
#' @param vecTimepointsUnique (numeric vector length number of unique
#' timepoints) Vector of unique time coordinates observed in this condition.
#' @param vecidxTimepoint (numeric vector length number of samples)
#' Index of the time coordinates of each sample (reference is
#' vecTimepointsUnique).
#' @param lsvecidxBatch (idx list length number of confounding variables)
#' List of vectors. 
#' One vector per confounding variable.
#' Each vector has one entry per sample with the index of the batch ID
#' within the given confounding variable of the given sample. Reference
#' is the list of unique batch ids for each confounding variable.
#' @param boolFitConst (bool) Whether to fit a constant model.
#' @param MAXIT (scalar) [Default 1000] 
#' Maximum number of BFGS iterations for model fitting with \link{optim}.
#' 
#' @return (list length 2)
#' Impulse and constant model fit to gene observations.
#' \itemize{
#' \item lsImpulseFit (list) List of impulse fit parameters and results.
#' \itemize{
#' \item vecImpulseParam (numeric vector length 6)
#' \{beta, h0, h1, h2, t1, t2\}
#' Maximum likelihood estimators of impulse model parameters.
#' \item vecImpulseValue (numeric vector length number of time points)
#' Values of impulse model fit at time points used for fit.
#' \item lsvecBatchFactors (list length number of confounders)
#' List of vectors of scalar batch correction factors for each sample.
#' These are also maximum likelihood estimators.
#' NULL if no confounders given.
#' \item scaDispParam (scalar) Dispersion parameter estimate
#' used in fitting (hyper-parameter).
#' \item scaLL (scalar) Loglikelihood of data under maximum likelihood
#' estimator model.
#' \item scaConvergence (scalar) 
#' Convergence status of optim on impulse model.
#' }
#' \item lsConstFit (list) List of constant fit parameters and results.
#' \itemize{
#' \item scaMu (scalar) Maximum likelihood estimator of
#' negative binomial mean parameter.
#' \item lsvecBatchFactors (list length number of confounders)
#' List of vectors of scalar batch correction factors for each sample.
#' These are also maximum likelihood estimators.
#' NULL if no confounders given.
#' \item scaDispParam (scalar) Dispersion parameter estimate
#' used in fitting (hyper-parameter).
#' \item scaLL (scalar) Loglikelihood of data under maximum likelihood
#' estimator model.
#' \item scaConvergence (scalar) 
#' Convergence status of optim on constant model.
#' }
#' }
#' 
#' @author David Sebastian Fischer
fitConstImpulseGene <- function(
    vecCounts, scaDisp, vecSizeFactors, vecTimepointsUnique, 
    vecidxTimepoint, lsvecidxBatch, boolFitConst, MAXIT = 1000) {
    
    # (I) Fit Impulse model 1. Compute initialisations
    lsParamGuesses <- estimateImpulseParam(
        vecCounts = vecCounts, vecTimepoints = vecTimepointsUnique[vecidxTimepoint], 
        lsvecidxBatch = lsvecidxBatch, vecSizeFactors = vecSizeFactors)
    vecParamGuessPeak <- lsParamGuesses$peak
    vecParamGuessValley <- lsParamGuesses$valley
    
    # 2. Initialisation: Peak
    lsFitPeak <- fitImpulseModel(vecImpulseParamGuess = vecParamGuessPeak, 
                                 vecCounts = vecCounts, scaDisp = scaDisp, 
                                 vecSizeFactors = vecSizeFactors, 
                                 vecTimepointsUnique = vecTimepointsUnique, 
                                 vecidxTimepoint = vecidxTimepoint, 
                                 lsvecidxBatch = lsvecidxBatch, MAXIT = MAXIT)
    # 3. Initialisation: Valley
    lsFitValley <- fitImpulseModel(
        vecImpulseParamGuess = vecParamGuessValley, 
        vecCounts = vecCounts, scaDisp = scaDisp, 
        vecSizeFactors = vecSizeFactors, 
        vecTimepointsUnique = vecTimepointsUnique, 
        vecidxTimepoint = vecidxTimepoint, 
        lsvecidxBatch = lsvecidxBatch, MAXIT = MAXIT)
    
    # (II) Select best fit and report fit type
    if (lsFitValley$scaLL > lsFitPeak$scaLL) {
        lsBestImpulseFit <- lsFitValley
    } else {
        lsBestImpulseFit <- lsFitPeak
    }
    
    if (boolFitConst) {
        # (III) Fit null model and compute loglikelihood of null model Only
        # necessary if running case only DE analysis as otherwise two impulse
        # models are compared.
        lsConstFit <- fitConstModel(
            vecCounts = vecCounts, scaDisp = scaDisp, 
            vecSizeFactors = vecSizeFactors, lsvecidxBatch = lsvecidxBatch)
    } else {
        lsConstFit <- NULL
    }
    
    return(list(lsImpulseFit = lsBestImpulseFit, lsConstFit = lsConstFit))
}

#' Fits impulse and constant models to all genes on all samples
#' of a condition
#' 
#' [Model fitting function hierarchy: 2 out of 4]
#' This secondary fitting wrapper performs parralelisation of 
#' model fitting across genes.
#' 
#' @seealso Called by \link{fitModels} to fit constant and impulse
#' model to samples of one condition.
#' Calls \link{fitConstImpulseGene} to perform fitting on each gene.
#' 
#' @param matCountDataProcCondition (matrix genes x samples)
#' Read count data.
#' @param vecSizeFactors (numeric vector number of samples) 
#' Model scaling factors for each sample which take
#' sequencing depth into account (size factors).
#' @param vecDispersions (vector number of genes)  Gene-wise 
#' negative binomial dispersion hyper-parameters.
#' @param vecTimepoints (numeric vector length number of samples)
#' Time coordinates of each sample.
#' @param lsvecBatches (list length number of confounding variables)
#' List of vectors. 
#' One vector per confounding variable.
#' Each vector has one entry per sample with the name of the batch ID
#' within the given confounding variable of the given sample.
#' @param boolFitConst (bool) Whether to fit a constant model.
#' 
#' @return (list length 5)
#' \itemize{
#' \item lsFits (list of lists length number of genes) 
#' List of model fits for each gene.
#' Each gene entry is a list of model fits to the individual models:
#' Impulse model and constant model (if boolFitConst is TRUE).
#' At this level, the sigmoid model fit can be added later.
#' Each model fit per gene is a list of fitting parameters and results.
#' \itemize{
#' \item Gene ID (list length 2)
#' Impulse and constant model fit to gene observations.
#' One entry of this format for all gene IDs.
#' \itemize{
#' \item lsImpulseFit (list) List of impulse fit parameters and results.
#' \itemize{
#' \item vecImpulseParam (numeric vector length 6)
#' \{beta, h0, h1, h2, t1, t2\}
#' Maximum likelihood estimators of impulse model parameters.
#' \item vecImpulseValue (numeric vector length number of time points)
#' Values of impulse model fit at time points used for fit.
#' \item lsvecBatchFactors (list length number of confounders)
#' List of vectors of scalar batch correction factors for each sample.
#' These are also maximum likelihood estimators.
#' NULL if no confounders given.
#' \item scaDispParam (scalar) Dispersion parameter estimate
#' used in fitting (hyper-parameter).
#' \item scaLL (scalar) Loglikelihood of data under maximum likelihood
#' estimator model.
#' \item scaConvergence (scalar) 
#' Convergence status of optim on impulse model.
#' }
#' \item lsConstFit (list) List of constant fit parameters and results.
#' \itemize{
#' \item scaMu (scalar) Maximum likelihood estimator of
#' negative binomial mean parameter.
#' \item lsvecBatchFactors (list length number of confounders)
#' List of vectors of scalar batch correction factors for each sample.
#' These are also maximum likelihood estimators.
#' NULL if no confounders given.
#' \item scaDispParam (scalar) Dispersion parameter estimate
#' used in fitting (hyper-parameter).
#' \item scaLL (scalar) Loglikelihood of data under maximum likelihood
#' estimator model.
#' \item scaConvergence (scalar) 
#' Convergence status of optim on constant model.
#' }
#' }
#' \item vecTimepointsUnique (numeric vector length number of unique
#' timepoints) Vector of unique time coordinates observed in this condition.
#' \item vecidxTimepoint (idx vector length number of samples)
#' Index of the time coordinates of each sample (reference is
#' vecTimepointsUnique).
#' \item lsvecBatchUnique (list number of confounders)
#' List of string vectors. One vector per confounder: vector of unique batches
#' in this confounder.
#' \item lsvecidxBatches (idx list length number of confounding variables)
#'   List of index vectors. 
#'   One vector per confounding variable.
#'   Each vector has one entry per sample with the index of the batch ID
#'   within the given confounding variable of the given sample. Reference
#'   is the list of unique batch ids for each confounding variable.
#' }
#' }
#' 
#' @author David Sebastian Fischer
fitConstImpulse <- function(
    matCountDataProcCondition, vecDispersions, vecSizeFactors, 
    vecTimepoints, lsvecBatches, boolFitConst) {
    
    # Maximum number of iterations for numerical optimisation of likelihood
    # function in MLE fitting of impulse and constant model:
    MAXIT <- 1000
    
    # Get sample group assignment indices: Time and batch
    vecTimepointsUnique <- sort(unique(vecTimepoints))
    vecidxTimepoint <- match(vecTimepoints, vecTimepointsUnique)
    # Get batch assignments
    if (length(lsvecBatches) > 0) {
        lsvecBatchUnique <- list()
        lsvecidxBatch <- list()
        for (batchfactor in lsvecBatches) {
            lsvecBatchUnique[[length(lsvecBatchUnique) + 1]] <- 
                unique(batchfactor)
            lsvecidxBatch[[length(lsvecidxBatch) + 1]] <- 
                match(batchfactor, unique(batchfactor))
        }
        names(lsvecBatchUnique) <- names(lsvecBatches)
        names(lsvecidxBatch) <- names(lsvecBatches)
    } else {
        lsvecBatchUnique <- NULL
        lsvecidxBatch <- NULL
    }
    
    lsFits <- bplapply(rownames(matCountDataProcCondition), function(x) {
        fitConstImpulseGene(
            vecCounts = matCountDataProcCondition[x, ], 
            scaDisp = vecDispersions[x], vecSizeFactors = vecSizeFactors, 
            vecTimepointsUnique = vecTimepointsUnique, 
            vecidxTimepoint = vecidxTimepoint, 
            lsvecidxBatch = lsvecidxBatch, boolFitConst = boolFitConst, 
            MAXIT = MAXIT)
    })
    names(lsFits) <- rownames(matCountDataProcCondition)
    
    return(list(lsFits = lsFits, vecTimepointsUnique = vecTimepointsUnique, 
                vecidxTimepoint = vecidxTimepoint, lsvecBatchUnique = lsvecBatchUnique, 
                lsvecidxBatch = lsvecidxBatch))
}

#' Fits impulse and constant models to a timecourse dataset
#' 
#' [Model fitting function hierarchy: 1 out of 4]
#' This primary wrapper coordinates fitting of impulse and constant model
#' to separate conditions according to the differential expression
#' mode (case-only or case-control).
#' 
#' @seealso Calls \link{fitConstImpulse}
#' once for each condition with the appropriate parameters and samples.
#' 
#' @param objectImpulseDE2 (object class ImpulseDE2Object)
#' Object to be fit.
#' @param vecConfounders (vector of strings number of confounding variables)
#' Factors to correct for during batch correction.
#' Names refer to columns in dfAnnotation.
#' @param boolCaseCtrl (bool) 
#' Whether to perform case-control analysis. Does case-only
#' analysis if FALSE.
#' 
#' @return objectImpulseDE2 (object class ImpulseDE2Object)
#' Object with sigmoidal fit added: objectImpulseDE2@lsModelFits
#' is updated to:
#' lsModelFits (list length number of conditions fit (1 or 3) +1)
#' \{'case'\} or \{'case', 'control', 'combined'\}
#' One model fitting object for each condition:
#' In case-only DE analysis, only the condition \{'case'\} is fit.
#' In case-control DE analysis, the conditions 
#' \{'case', 'control','combined\} are fit.
#' Each condition entry is a list of model fits for each gene.
#' Each gene entry is a list of model fits to the individual models:
#' Impulse model and constant model (if boolFitConst is TRUE).
#' At this level, the sigmoid model fit can be added later.
#' Each model fit per gene is a list of fitting parameters and results.
#' \itemize{
#' \item IdxGroups (list length number of conditions)
#' Samples grouped by time points and by batches and time point vectors. 
#' Sample groups are stored in the form of index vectors in which
#' samples of the same time point or batch have the same index.
#' \itemize{
#' \item Condition ID (list length 5)
#' List of index vectors and time points.
#' One entry of this format for each condition.
#' \itemize{
#' \item vecTimepointsUnique (numeric vector length number of unique
#' timepoints) Vector of unique time coordinates observed in this condition.
#' \item vecidxTimepoint (idx vector length number of samples)
#' Index of the time coordinates of each sample (reference is
#' vecTimepointsUnique).
#' \item lsvecBatchUnique (list number of confounders)
#' List of string vectors. One vector per confounder: vector of unique batches
#' in this confounder.
#' \item lsvecidxBatches (idx list length number of confounding variables)
#' List of index vectors. 
#' One vector per confounding variable.
#' Each vector has one entry per sample with the index of the batch ID
#' within the given confounding variable of the given sample. Reference
#' is the list of unique batch ids for each confounding variable.
#'   \item vecSamples (vector number of samples) Names of samples fit
#' for this condition in same order as index vectors above.
#' }
#' }
#' \item Condition ID (list length number of genes)
#' List of fits for each gene to the samples of this condition.
#' One entry of this format for all conditions fit.
#' \itemize{
#' \item Gene ID (list length 2)
#' Impulse and constant model fit to gene observations.
#' One entry of this format for all gene IDs.
#' \itemize{
#' \item lsImpulseFit (list) List of impulse fit parameters and results.
#' \itemize{
#' \item vecImpulseParam (numeric vector length 6)
#' \{beta, h0, h1, h2, t1, t2\}
#' Maximum likelihood estimators of impulse model parameters.
#' \item vecImpulseValue (numeric vector length number of time points)
#' Values of impulse model fit at time points used for fit.
#' \item lsvecBatchFactors (list length number of confounders)
#' List of vectors of scalar batch correction factors for each sample.
#' These are also maximum likelihood estimators.
#' NULL if no confounders given.
#' \item scaDispParam (scalar) Dispersion parameter estimate
#' used in fitting (hyper-parameter).
#' \item scaLL (scalar) Loglikelihood of data under maximum likelihood
#' estimator model.
#' \item scaConvergence (scalar) 
#' Convergence status of optim on impulse model.
#' }
#' \item lsConstFit (list) List of constant fit parameters and results.
#' \itemize{
#' \item scaMu (scalar) Maximum likelihood estimator of
#' negative binomial mean parameter.
#' \item lsvecBatchFactors (list length number of confounders)
#' List of vectors of scalar batch correction factors for each sample.
#' These are also maximum likelihood estimators.
#' NULL if no confounders given.
#' \item scaDispParam (scalar) Dispersion parameter estimate
#' used in fitting (hyper-parameter).
#' \item scaLL (scalar) Loglikelihood of data under maximum likelihood
#' estimator model.
#' \item scaConvergence (scalar) 
#' Convergence status of optim on constant model.
#' }
#' }
#' }
#' }
#' 
#' @author David Sebastian Fischer
fitModels <- function(objectImpulseDE2, vecConfounders, boolCaseCtrl) {
    
    lsFitResults_all = list()
    dfAnnot <- get_dfAnnotationProc(obj=objectImpulseDE2)
    
    # Conditions to be fitted separately
    if (boolCaseCtrl) {
        vecLabels <- c("combined", "case", "control")
    } else {
        vecLabels <- c("case")
    }
    
    # Create lists of samples per condition
    if (boolCaseCtrl) {
        lsSamplesByCond <- list(
            combined = dfAnnot$Sample, 
            case = dfAnnot[dfAnnot$Condition == "case", ]$Sample, 
            control = dfAnnot[dfAnnot$Condition == "control", ]$Sample)
    } else {
        lsSamplesByCond <- list(
            case = dfAnnot[dfAnnot$Condition == "case", ]$Sample)
    }
    
    # Get batch assignments of samples
    if (!is.null(vecConfounders)) {
        lsvecBatches <- lapply(vecConfounders, function(confounder) {
            vecBatches <- dfAnnot[, confounder]
            names(vecBatches) <- dfAnnot$Sample
            return(vecBatches)
        })
        names(lsvecBatches) <- vecConfounders
    } else {
        lsvecBatches <- NULL
    }
    
    # Get time point assignments of samples
    vecTimepoints <- dfAnnot$Time
    names(vecTimepoints) <- colnames(get_matCountDataProc(obj=objectImpulseDE2))
    
    # Fitting for different runs Fit constant model only if doing case-only
    # analysis for which the constant model is the null model or in the
    # combined case of case-ctrl to derive a inferred mean parameter for
    # reference.
    lsFitResultsByCond <- lapply(vecLabels, function(label) {
        lsFitResults <- fitConstImpulse(
            matCountDataProcCondition = get_matCountDataProc(
                obj=objectImpulseDE2)[,lsSamplesByCond[[label]]], 
            vecSizeFactors = get_vecSizeFactors(
                obj=objectImpulseDE2)[lsSamplesByCond[[label]]], 
            vecDispersions = get_vecDispersions(obj=objectImpulseDE2), 
            vecTimepoints = vecTimepoints[lsSamplesByCond[[label]]], 
            lsvecBatches = lapply(lsvecBatches, function(confounder) 
                confounder[lsSamplesByCond[[label]]] ), 
            boolFitConst = TRUE)
        # Note: Don't need constant fit with case-control other than for results
        # table and transient identification. It is however fast and the
        # inferred means may be of interest - do for all conditions.
        return(lsFitResults)
    })
    lsModelFitsByCondFormat <- lapply(
        lsFitResultsByCond, function(res) res$lsFits)
    names(lsModelFitsByCondFormat) <- vecLabels
    lsModelFitsByCondFormat$IdxGroups <- lapply(
        lsFitResultsByCond, function(res) {
            list(vecTimepointsUnique = res$vecTimepointsUnique, 
                 vecidxTimepoint = res$vecidxTimepoint, 
                 lsvecBatchUnique = res$lsvecBatchUnique, 
                 lsvecidxBatch = res$lsvecidxBatch)
        })
    names(lsModelFitsByCondFormat$IdxGroups) <- vecLabels
    for (label in vecLabels){ 
        lsModelFitsByCondFormat$IdxGroups[[label]]$Samples <- 
            lsSamplesByCond[[label]] 
    }
    
    objectImpulseDE2 <- set_lsModelFits(obj=objectImpulseDE2, 
                                        element=lsModelFitsByCondFormat)
    return(objectImpulseDE2)
}
