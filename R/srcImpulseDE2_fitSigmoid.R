### Sigmoidal  model fit

#' Compute up and down sigmoid model parameter initialisations
#' for data of one gene
#' 
#' [Model fitting function hierarchy: helper to level 2 out of 3]
#' This is a fitting helper function which computes parameter intialisations
#' and does not wrap or execute numerical optimisation.
#' The up model models a sigmoidal expression increase over time,
#' the down model a sigmoidal decrease over time.
#' 
#' @seealso Called by \link{fitSigmoidGene}.
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
#' \{beta, h0, h1, t\}
#' Up model initialisations of sigmoidal model parameters.
#' \item valley (numeric vector length 6)
#' \{beta, h0, h1, t\}
#' Down model initialisations of sigmoidal model parameters.
#' }
#' 
#' @author David Sebastian Fischer
estimateSigmoidParam <- function(
    vecCounts, vecTimepoints, vecSizeFactors, 
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
            vecBatchFactors <- vecBatchFactors * vecBatchFactorsConfounder[vecidxBatch]
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
    idxMiddleTP <- round(scaNTimepoints/2)
    vecdidxFirstPart <- seq(1, idxMiddleTP-1, by=1)
    vecdidxSecndPart <- seq(idxMiddleTP, scaNTimepoints, by=1)
    scaMaxEarlyMean <- max(vecExpressionMeans[vecdidxFirstPart], na.rm = TRUE)
    scaMinEarlyMean <- min(vecExpressionMeans[vecdidxFirstPart], na.rm = TRUE)
    scaMaxLateMean <- max(vecExpressionMeans[vecdidxSecndPart], 
                          na.rm = TRUE)
    scaMinLateMean <- min(vecExpressionMeans[vecdidxSecndPart], 
                          na.rm = TRUE)
    
    # Compute up initialisation
    vecParamGuessUp <- c(
        1, 
        log(scaMinEarlyMean + 1), 
        log(scaMaxLateMean + 1), 
        vecTimepointsUnique[idxMiddleTP] )
    
    # Compute down initialisation
    vecParamGuessDown <- c(
        1, 
        log(scaMaxEarlyMean + 1), 
        log(scaMinLateMean + 1), 
        vecTimepointsUnique[idxMiddleTP] )
    
    return(list(up = vecParamGuessUp, down = vecParamGuessDown))
}


#' Fit a sigmoidal model to data of a gene
#' 
#' [Model fitting function hierarchy: 3 out of 3]
#' This tertiary fitting wrapper performs sigmoidal model fitting:
#' This function executes numerical optimisaiton and error-handling
#' thereof.
#' 
#' @seealso Called by \link{fitSigmoidGene} to fit sigmoidal
#' model to samples of one condition and one gene.
#' Calls sigmoidal model cost function 
#' \link{evalLogLikSigmoid_comp} within \link{optim}.
#' 
#' @param vecSigmoidParamGuess (numeric vector length 4)
#' \{beta, h0, h1, t\}
#' Up model initialisations of sigmoidal model parameters.
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
#' @return (list) List of sigmoid fit parameters and results.
#' \itemize{
#' \item vecSigmoidParam (numeric vector length 4)
#' \{beta, h0, h1, t\}
#' Maximum likelihood estimators of sigmoidal model parameters.
#' \item vecSigmoidValue (numeric vector length number of time points)
#' Values of sigmoid model fit at time points used for fit.
#' \item lsvecBatchFactors (list length number of confounders)
#' List of vectors of scalar batch correction factors for each sample.
#' These are also maximum likelihood estimators.
#' NULL if no confounders given.
#' \item scaDispParam (scalar) Dispersion parameter estimate
#' used in fitting (hyper-parameter).
#' \item scaLL (scalar) Loglikelihood of data under maximum likelihood
#' estimator model.
#' \item scaConvergence (scalar) 
#' Convergence status of optim for sigmoid model.
#' }
#' 
#' @author David Sebastian Fischer
fitSigmoidModel <- function(
    vecSigmoidParamGuess, vecCounts, scaDisp, vecSizeFactors, 
    lsvecidxBatch, vecTimepointsUnique, vecidxTimepoint, 
    MAXIT = 1000, RELTOL = 10^(-8), trace = 0, REPORT = 10) {
    
    
    vecParamGuess <- vecSigmoidParamGuess
    if (!is.null(lsvecidxBatch)) {
        for (vecidxConfounder in lsvecidxBatch) {
            vecParamGuess <- c( vecParamGuess, 
                                rep(0, length(unique(vecidxConfounder)) - 1)) 
        }
    }
    
    lsFit <- tryCatch({
        optim(par = vecParamGuess, fn = evalLogLikSigmoid_comp, 
              vecCounts = vecCounts, scaDisp = scaDisp, 
              vecSizeFactors = vecSizeFactors, 
              vecTimepointsUnique = vecTimepointsUnique, 
              vecidxTimepoint = vecidxTimepoint, lsvecidxBatch = lsvecidxBatch, 
              vecboolObserved = !is.na(vecCounts), method = "BFGS", 
              control = list(maxit = MAXIT, reltol = RELTOL, fnscale = -1)
        )[c("par", "value", "convergence")]
    }, error = function(strErrorMsg) {
        print(paste0("ERROR: Fitting sigmoid model: fitSigmoidModel().", 
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
    vecSigmoidParam <- lsFit$par[1:4]
    vecSigmoidParam[2:3] <- exp(vecSigmoidParam[2:3])
    vecSigmoidParam[2:3][vecSigmoidParam[2:3] < 10^(-10)] <- 10^(-10)
    vecSigmoidParam[2:3][vecSigmoidParam[2:3] > 10^(10)] <- 10^(10)
    names(vecSigmoidParam) <- c("beta", "h0", "h1", "t")
    vecSigmoidValue <- evalSigmoid_comp(
        vecSigmoidParam = vecSigmoidParam, 
        vecTimepoints = vecTimepointsUnique)[vecidxTimepoint]
    names(vecSigmoidValue) <- names(vecCounts)
    scaNParamUsed <- 4
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
    
    return(list(
        vecSigmoidParam = vecSigmoidParam, vecSigmoidValue = vecSigmoidValue, 
        lsvecBatchFactors = lsvecBatchFactors, scaDispParam = scaDisp, 
        scaLL = lsFit$value, scaConvergence = lsFit$convergence))
}

#' Fit a sigmoidal model to a single gene
#' 
#' [Model fitting function hierarchy: 2 out of 3]
#' This secondary fitting wrapper calls the optimisation wrappers
#' for the individual fitting operations to be performed on the 
#' observations of this gene.
#' Structure of this function:
#' \itemize{
#' \item Fit sigmoidal model
#' \itemize{
#' \item Initialise sigmoidal model parameters (up and down)
#' \item Fit sigmoidal model (up initialisation)
#' \item Fit sigmoidal model (down initialisation)
#' }
#' \item Select best sigmoidal model fit from initialisations,
#' }
#' 
#' @seealso Called by \link{fitSigmoidModels} to fit
#' sigmoidal model to samples of one condition and one gene.
#' Calls sigmoidal parameter initialisation function
#' \link{estimateSigmoidParam} and 
#' optimisation wrapper \link{fitSigmoidModel}.
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
#' @param vecidxTimepoint (idx vector length number of samples)
#' Index of the time coordinates of each sample (reference is
#' vecTimepointsUnique).
#' @param lsvecidxBatch (idx list length number of confounding variables)
#' List of vectors. 
#' One vector per confounding variable.
#' Each vector has one entry per sample with the index of the batch ID
#' within the given confounding variable of the given sample. Reference
#' is the list of unique batch ids for each confounding variable.
#' @param MAXIT (scalar) [Default 1000] 
#' Maximum number of BFGS iterations for model fitting with \link{optim}.
#' 
#' @return (list) 
#' List of sigmoidal fit parameters and results.
#' \itemize{
#' \item vecSigmoidParam (numeric vector length 4)
#' \{beta, h0, h1, t\}
#' Maximum likelihood estimators of sigmoidal model parameters.
#' \item vecSigmoidValue (numeric vector length number of time points)
#' Values of sigmoid model fit at time points used for fit.
#' \item lsvecBatchFactors (list length number of confounders)
#' List of vectors of scalar batch correction factors for each sample.
#' These are also maximum likelihood estimators.
#' NULL if no confounders given.
#' \item scaDispParam (scalar) Dispersion parameter estimate
#' used in fitting (hyper-parameter).
#' \item scaLL (scalar) Loglikelihood of data under maximum likelihood
#' estimator model.
#' \item scaConvergence (scalar) 
#' Convergence status of optim on sigmoidal model.
#' }
#' 
#' @author David Sebastian Fischer
fitSigmoidGene <- function(
    vecCounts, scaDisp, vecSizeFactors, vecTimepointsUnique, 
    vecidxTimepoint, lsvecidxBatch, MAXIT = 1000) {
    
    # (I) Fit sigmoidal model 1. Compute initialisations
    lsParamGuesses <- estimateSigmoidParam(
        vecCounts = vecCounts, 
        vecTimepoints = vecTimepointsUnique[vecidxTimepoint], 
        lsvecidxBatch = lsvecidxBatch, 
        vecSizeFactors = vecSizeFactors)
    vecParamGuessUp <- lsParamGuesses$up
    vecParamGuessDown <- lsParamGuesses$down
    
    # 2. Initialisation: Up
    lsFitUp <- fitSigmoidModel(
        vecSigmoidParamGuess = vecParamGuessUp, vecCounts = vecCounts, 
        scaDisp = scaDisp, vecSizeFactors = vecSizeFactors, 
        vecTimepointsUnique = vecTimepointsUnique, 
        vecidxTimepoint = vecidxTimepoint, lsvecidxBatch = lsvecidxBatch, 
        MAXIT = MAXIT)
    # 3. Initialisation: Down
    lsFitDown <- fitSigmoidModel(
        vecSigmoidParamGuess = vecParamGuessDown, vecCounts = vecCounts, 
        scaDisp = scaDisp, vecSizeFactors = vecSizeFactors, 
        vecTimepointsUnique = vecTimepointsUnique, 
        vecidxTimepoint = vecidxTimepoint, lsvecidxBatch = lsvecidxBatch, 
        MAXIT = MAXIT)
    
    # (II) Select best fit and report fit type
    if (lsFitDown$scaLL > lsFitUp$scaLL) {
        lsbestSigmoidFit <- lsFitDown
    } else {
        lsbestSigmoidFit <- lsFitUp
    }
    
    return(lsbestSigmoidFit)
}

#' Fits sigmoidal models to all genes on all all samples
#' of a condition
#' 
#' [Model fitting function hierarchy: 1 out of 3]
#' This primary fitting wrapper performs parralelisation of 
#' model fitting across genes.
#' 
#' @seealso Calls \link{fitSigmoidGene} to perform fitting on each gene.
#' 
#' @param objectImpulseDE2 (object class ImpulseDE2Object)
#' Object to be fit with sigmoidal model. Needs to be fitted with impulse 
#' model before.
#' @param vecConfounders (vector of strings number of confounding variables)
#' Factors to correct for during batch correction.
#' Names refer to columns in dfAnnotation.
#' @param strCondition (str)
#' Name of condition entry in lsModelFits for which sigmoidal
#' models are to be fit to each gene.
#' 
#' @return objectImpulseDE2 (object class ImpulseDE2Object)
#' Object with sigmoidal fit added: objectImpulseDE2@@lsModelFits
#' is updated to:
#' lsModelFits (list length number of conditions fit (1 or 3) +1)
#' \{'case'\} or \{'case', 'control', 'combined'\}
#' This is the lsModelFits object handed to this function with additional
#' sigmoid fit entries for every gene for the given condition.
#' One model fitting object for each condition:
#' In case-only DE analysis, only the condition \{'case'\} is fit.
#' In case-control DE analysis, the conditions 
#' \{'case', 'control','combined\} are fit.
#' Each condition entry is a list of model fits for each gene.
#' Each gene entry is a list of model fits to the individual models:
#' Impulse model, constant model and sigmoidal fit.
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
#' Impulse, constant and sigmoidal model fit to gene observations.
#' One entry of this format for all gene IDs.
#' \itemize{
#' \item lsImpulseFit (list) List of impulse fit parameters and results.
#' For details, read the annotation of \link{fitModels}.
#' \item lsConstFit (list) List of constant fit parameters and results.
#' For details, read the annotation of \link{fitModels}.
#' \item ls SigmoidFit (list) List of sigmoidal fit parameters and results.
#' \itemize{
#' \item vecSigmoidParam (numeric vector length 4)
#' \{beta, h0, h1, t\}
#' Maximum likelihood estimators of sigmoidal model parameters.
#' \item vecSigmoidValue (numeric vector length number of time points)
#' Values of sigmoid model fit at time points used for fit.
#' \item lsvecBatchFactors (list length number of confounders)
#' List of vectors of scalar batch correction factors for each sample.
#' These are also maximum likelihood estimators.
#' NULL if no confounders given.
#' \item scaDispParam (scalar) Dispersion parameter estimate
#' used in fitting (hyper-parameter).
#' \item scaLL (scalar) Loglikelihood of data under maximum likelihood
#' estimator model.
#' \item scaConvergence (scalar) 
#' Convergence status of optim on sigmoidal model.
#' }
#' }
#' }
#' }
#'
#' @examples
#' lsSimulatedData <- simulateDataSetImpulseDE2(
#' vecTimePointsA   = rep(seq(1,8),3),
#' vecTimePointsB   = NULL,
#' vecBatchesA      = NULL,
#' vecBatchesB      = NULL,
#' scaNConst        = 0,
#' scaNImp          = 20,
#' scaNLin          = 10,
#' scaNSig          = 20)
#' objectImpulseDE2 <- runImpulseDE2(
#' matCountData    = lsSimulatedData$matObservedCounts, 
#' dfAnnotation    = lsSimulatedData$dfAnnotation,
#' boolCaseCtrl    = FALSE,
#' vecConfounders  = NULL,
#' boolIdentifyTransients = FALSE,
#' scaNProc        = 1 )
#' # You could have used boolIdentifyTransients=TRUE
#' # to avoid the following post wrapper fitting.
#' objectImpulseDE2 <- fitSigmoidModels(
#' objectImpulseDE2 = objectImpulseDE2,
#' vecConfounders   = NULL,
#' strCondition     = 'case')
#' objectImpulseDE2 <- updateDEAnalysis(
#' objectImpulseDE2=objectImpulseDE2,
#' scaQThresTransients=0.001)
#' head(objectImpulseDE2$dfImpulseDE2Results)
#' # dfImpulseDE2Results now contain 'transients-analysis'.
#' 
#' @author David Sebastian Fischer
#' 
#' @export
fitSigmoidModels <- function(objectImpulseDE2, vecConfounders, strCondition) {
    
    # Load objects from output class
    matCountDataProc <- objectImpulseDE2@matCountDataProc
    dfAnnotationProc <- objectImpulseDE2@dfAnnotationProc
    lsModelFits <- objectImpulseDE2@lsModelFits
    vecSizeFactors <- objectImpulseDE2@vecSizeFactors
    vecDispersions <- objectImpulseDE2@vecDispersions
    
    vecSamplesCond <- dfAnnotationProc[
        dfAnnotationProc$Condition == strCondition, ]$Sample
    
    # Get batch assignments of samples
    lsvecidxBatchCond <- lsModelFits$IdxGroups[[strCondition]]$lsvecidxBatch
    # Get time point assignments of samples
    vecTimepointsUniqueCond <- lsModelFits$IdxGroups[[strCondition]]$
        vecTimepointsUnique
    vecidxTimepointCond <- lsModelFits$IdxGroups[[strCondition]]$
        vecidxTimepoint
    
    # Developmental note: Compared to impulse/constant fitting, this
    # function does not iterate over conditions as this is likely only used
    # for one condition (case). Therefore merge two wrappers used for
    # impulse/const fit into one here.
    
    # Maximum number of iterations for numerical optimisation of likelihood
    # function in MLE fitting of sigmoidal model:
    MAXIT <- 1000
    
    lsSigmoidFits <- bplapply(rownames(matCountDataProc), function(x) {
        fitSigmoidGene(
            vecCounts = matCountDataProc[x, vecSamplesCond], 
            scaDisp = vecDispersions[x], 
            vecSizeFactors = vecSizeFactors[vecSamplesCond], 
            vecTimepointsUnique = vecTimepointsUniqueCond, 
            vecidxTimepoint = vecidxTimepointCond, 
            lsvecidxBatch = lsvecidxBatchCond, MAXIT = MAXIT)
    })
    names(lsSigmoidFits) <- rownames(matCountDataProc)
    
    # Add sigmoid fits into model fit data structure to preexisting impulse
    # (and constant) fits.
    for (x in rownames(matCountDataProc)) {
        lsModelFits[[strCondition]][[x]]$lsSigmoidFit <- lsSigmoidFits[[x]]
    }
    
    objectImpulseDE2@lsModelFits <- lsModelFits
    return(objectImpulseDE2)
}
