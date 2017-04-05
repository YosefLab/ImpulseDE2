###  Simulate a data set

#' Simulate a data set for ImpulseDE2
#' 
#' Simulates a data set with genes with constant and impulse
#' expression traces. Expression strength and variation in impulse
#' like traces are parameterised and random. All temporary files
#' are saved into dirOutSimulation and only the objects necessary
#' for running ImpulseDE2 (the count matrix and the annotation table
#' are returned). The remaining objects representing hidden
#' parameters can be used to evaluate parameter estimates.
#' 
#' @seealso Called by separately by user.
#' 
#' @param vecTimePointsA (numeric vector number of time points)
#' Number of time points in batch A.
#' @param vecTimePointsB (numeric vector number of time points)
#' Number of time points in batch B.
#' @param vecBatchesA (str vector number of samples in vecTimePointsA) 
#' [Default NULL]
#' Batch IDs of each sample in condition A. Set to NULL if 
#' simulating without batch effects.
#' @param vecBatchesB (str vector number of samples in vecTimePointsB) 
#' [Default NULL]
#' Batch IDs of each sample in condition B. Set to NULL if 
#' simulating without batch effects.
#' @param scaNConst (scalar) Number of constant genes in data set.
#' @param scaNImp (scalar) Number of impulse distributed genes in data set.
#' @param scaNLin (scalar) Number of linear distributed genes in data set.
#' @param scaNSig (scalar) Number of sigmoid distributed genes in data set.
#' @param scaNRand (scalar) [Default NULL] Number of random 
#' distributed genes in data set.
#' @param scaSeedInit (scalar) [Default 1] Scalar based on which seeds are chosen.
#' One vlaue correspond sto a unique set of seeds for all random number generations.
#' @param scaMumax (scalar) [Default 1000]
#' Maximum expression mean parameter to be used.
#' @param boolOneConstMu (bool) [Default False]
#' Don't sample constant trajectories from uniform [0,scaMumax]
#' but set all to scaMumax
#' @param scaSDExpressionChange (scalar) [Default 1]
#' Standard deviation of normal distribution from which the 
#' amplitude change within an impulse trace is drawn.
#' @param scaSDRand (scalar) [Default 0]
#' Standard deviation of normal distribution from which the 
#' random deviations are drawn.
#' @param scaMuSizeEffect (numeric vector number of genes) [Default NULL]
#' Mean of normal distribution of which scaNLing factor for 
#' size effects per sample are drawn.
#' @param scaSDSizeEffect (numeric vector number of genes) [Default NULL]
#' Standard deviation of normal distribution of which scaling factor for 
#' size effects per sample are drawn.
#' @param scaMuBatchEffect (numeric vector number of genes) [Default NULL]
#' Mean of normal distribution of which scaling factor for 
#' batch effects per gene are drawn (reference is batch A).
#' @param scaSDBatchEffect (numeric vector number of genes) [Default NULL]
#' Standard deviation of normal distribution of which scaling factor for 
#' batch effects per gene are drawn (reference is batch A).
#' @param dirOutSimulation (directory) [Default NULL]
#' Directory to which simulated parameter objects are 
#' saved to.
#' 
#' @return list (length 2)
#' \itemize{
#' \item dfAnnotation (data frame samples x covariates) 
#' {Sample, Condition, Time (numeric), TimeCateg (str)
#' (and confounding variables if given).}
#' Annotation table with covariates for each sample.
#' \item matSampledCountsObserved (matrix genes x cells)
#' Sampled count data of all cells after drop-out.
#' }
#' 
#' @examples
#' lsSimulatedData <- simulateDataSetImpulseDE2(
#' vecTimePointsA   = rep(seq(1,8),3),
#' vecTimePointsB   = NULL,
#' vecBatchesA      = NULL,
#' vecBatchesB      = NULL,
#' scaNConst        = 30,
#' scaNImp          = 10,
#' scaNLin          = 10,
#' scaNSig          = 10)
#' 
#' @author David Sebastian Fischer
#' 
#' @export
simulateDataSetImpulseDE2 <- function(
    vecTimePointsA, vecTimePointsB, vecBatchesA, 
    vecBatchesB, scaNConst, scaNImp, scaNLin, scaNSig, scaNRand = 0, 
    scaSeedInit = 1, scaMumax = 1000, boolOneConstMu = FALSE, 
    scaSDExpressionChange = 1, scaSDRand = NULL, scaMuSizeEffect = 1, 
    scaSDSizeEffect = 0.1, scaMuBatchEffect = NULL, 
    scaSDBatchEffect = NULL, dirOutSimulation = NULL) {
    
    #### Internal functions Evalute impulse model at time points
    evalImpulse <- function(t, beta, t1, t2, h0, h1, h2) {
        return(1/h1 * (h0 + (h1 - h0) * 1/(1 + exp(-beta * (t - t1)))) * 
                   (h2 + (h1 - h2) * 1/(1 + exp(beta * (t - t2)))))
    }
    # Evalute sigmoid model at time points
    evalSigmoid <- function(t, beta, t1, h0, h1) {
        return(h0 + (h1 - h0)/(1 + exp(-beta * (t - t1))))
    }
    
    #### Simulate data (case only) Create annotation and sample naming
    vecSamplesA <- vecTimePointsA
    names(vecSamplesA) <- sapply(seq(1, length(vecSamplesA)), function(i) {
        paste0("A_", vecSamplesA[i], "_Rep", 
               match(i, which(vecSamplesA == vecSamplesA[i])))
    })
    vecSamplesB <- vecTimePointsB
    if (!is.null(vecSamplesB)) {
        names(vecSamplesB) <- sapply(seq(1, length(vecSamplesB)), function(i) {
            paste0("B_", vecSamplesB[i], "_Rep", 
                   match(i, which(vecSamplesB == vecSamplesB[i])))
        })
        boolCaseCtrl <- TRUE
    } else {
        boolCaseCtrl <- FALSE
    }
    vecSamples <- c(vecSamplesA, vecSamplesB)
    scaNSamples <- length(vecSamples)
    vecTimePointsUnique <- unique(vecSamples)
    vecindTimePointAssign <- match(vecSamples, vecTimePointsUnique)
    vecTimePointsUniqueB <- unique(vecSamplesB)
    vecindTimePointAssignB <- match(vecSamplesB, vecTimePointsUniqueB)
    
    if (is.null(vecBatchesA)) {
        print("Setting no batch structure.")
        vecBatchesA <- rep("B_NULL", length(vecSamplesA))
        vecBatchesB <- rep("B_NULL", length(vecSamplesB))
    }
    dfAnnotation <- data.frame(Sample = names(vecSamples), 
                               Condition = rep("case", length(vecSamples)), 
                               Time = vecSamples, 
                               Batch = c(vecBatchesA, vecBatchesB), 
                               stringsAsFactors = FALSE)
    rownames(dfAnnotation) <- dfAnnotation$Sample
    if (boolCaseCtrl) {
        dfAnnotation[dfAnnotation$Sample %in% names(vecSamplesB), ]$Condition <-
            rep("control", sum(dfAnnotation$Sample %in% names(vecSamplesB)))
    }
    
    # 1. Create hidden data set
    if (scaNConst > 0) {
        vecConstIDs <- paste0(rep("gene", scaNConst), c(1:scaNConst))
    } else {
        vecConstIDs <- NULL
    }
    if (scaNImp > 0) {
        vecImpulseIDs <- paste0(rep("gene", scaNImp), 
                                c((scaNConst + 1):(scaNConst + scaNImp)))
    } else {
        vecImpulseIDs <- NULL
    }
    if (scaNLin > 0) {
        vecLinIDs <- paste0(rep("gene", scaNLin), 
                            c((scaNConst + scaNImp + 1):
                                  (scaNConst + scaNImp + scaNLin)))
    } else {
        vecLinIDs <- NULL
    }
    if (scaNSig > 0) {
        vecSigIDs <- paste0(rep("gene", scaNSig), 
                            c((scaNConst + scaNImp + scaNLin + 1):
                                  (scaNConst + scaNImp + scaNLin + scaNSig)))
    } else {
        vecSigIDs <- NULL
    }
    if (scaNRand > 0) {
        vecRandIDs <- paste0(rep("gene", scaNRand), 
                             c((scaNConst + scaNImp + scaNLin + scaNSig + 1):
                                   (scaNConst + scaNImp + scaNLin + 
                                        scaNSig + scaNRand)))
    } else {
        vecRandIDs <- NULL
    }
    
    scaNGenes <- scaNConst + scaNImp + scaNLin + scaNSig + scaNRand
    scaEps <- 1e-05
    scaSeedsUsed <- 0
    
    # a. Draw means from uniform (first half of genes): one mean per gene
    if (scaNConst > 0) {
        set.seed(scaSeedInit + scaSeedsUsed)
        scaSeedsUsed <- scaSeedsUsed + 1
        if (boolOneConstMu) {
            vecMuConstHidden <- rep(scaMumax, scaNConst) 
        } else {
            vecMuConstHidden <- runif(scaNConst) * scaMumax
        }
        
        vecMuConstHidden[vecMuConstHidden < scaEps] <- scaEps
        matMuConstHidden <- matrix(vecMuConstHidden, nrow = scaNConst, 
                                   ncol = scaNSamples, byrow = FALSE)
        rownames(matMuConstHidden) <- vecConstIDs
        colnames(matMuConstHidden) <- names(vecSamples)
    } else {
        matMuConstHidden <- NULL
    }
    
    # b. Draw means from impulse model
    if (scaNImp > 0) {
        set.seed(scaSeedInit + scaSeedsUsed)
        scaSeedsUsed <- scaSeedsUsed + 1
        beta <- runif(scaNImp) * 2 + 0.5
        set.seed(scaSeedInit + scaSeedsUsed)
        scaSeedsUsed <- scaSeedsUsed + 1
        ta <- runif(scaNImp) * max(vecTimePointsUnique)
        set.seed(scaSeedInit + scaSeedsUsed)
        scaSeedsUsed <- scaSeedsUsed + 1
        tb <- runif(scaNImp) * max(vecTimePointsUnique)
        t1 <- ta
        t1[tb < ta] <- tb[tb < ta]
        t2 <- tb
        t2[tb < ta] <- ta[tb < ta]
        set.seed(scaSeedInit + scaSeedsUsed)
        scaSeedsUsed <- scaSeedsUsed + 1
        h0 <- runif(scaNImp) * scaMumax
        set.seed(scaSeedInit + scaSeedsUsed)
        scaSeedsUsed <- scaSeedsUsed + 1
        h1 <- h0 * abs(rnorm(n = scaNImp, mean = 1, 
                             sd = scaSDExpressionChange))
        set.seed(scaSeedInit + scaSeedsUsed)
        scaSeedsUsed <- scaSeedsUsed + 1
        h2 <- h0 * abs(rnorm(n = scaNImp, mean = 1, 
                             sd = scaSDExpressionChange))
        h0[h0 < scaEps] <- scaEps
        h1[h1 < scaEps] <- scaEps
        h2[h2 < scaEps] <- scaEps
        lsMuImpulseHidden <- lapply(seq(1, scaNImp), function(gene) {
            evalImpulse(t = vecTimePointsUnique, 
                        beta = beta[gene], t1 = t1[gene], 
                        t2 = t2[gene], h0 = h0[gene], h1 = h1[gene], 
                        h2 = h2[gene])[vecindTimePointAssign]
        })
        matImpulseModelHidden <- cbind(beta, h0, h1, h2, t1, t2)
        matMuImpulseHidden <- do.call(rbind, lsMuImpulseHidden)
        rownames(matImpulseModelHidden) <- vecImpulseIDs
        rownames(matMuImpulseHidden) <- vecImpulseIDs
        colnames(matMuImpulseHidden) <- names(vecSamples)
        if (boolCaseCtrl & scaNImp > 1) {
            # Keep start point the same, this doesn't work well for impulse model
            scaNImpDE <- round(scaNImp/2)
            set.seed(scaSeedInit + scaSeedsUsed)
            scaSeedsUsed <- scaSeedsUsed + 1
            betaDE <- runif(scaNImpDE) * 2 + 0.5
            set.seed(scaSeedInit + scaSeedsUsed)
            scaSeedsUsed <- scaSeedsUsed + 1
            ta <- runif(scaNImpDE) * max(vecTimePointsUniqueB)
            set.seed(scaSeedInit + scaSeedsUsed)
            scaSeedsUsed <- scaSeedsUsed + 1
            tb <- runif(scaNImpDE) * max(vecTimePointsUniqueB)
            t1DE <- ta
            t1DE[tb < ta] <- tb[tb < ta]
            t2DE <- tb
            t2DE[tb < ta] <- ta[tb < ta]
            set.seed(scaSeedInit + scaSeedsUsed)
            scaSeedsUsed <- scaSeedsUsed + 1
            h1DE <- h0[1:scaNImpDE] * abs(rnorm(n = scaNImpDE, mean = 1, 
                                                sd = scaSDExpressionChange))
            set.seed(scaSeedInit + scaSeedsUsed)
            scaSeedsUsed <- scaSeedsUsed + 1
            h2DE <- h0[1:scaNImpDE] * abs(rnorm(n = scaNImpDE, mean = 1, 
                                                sd = scaSDExpressionChange))
            h1DE[h1 < scaEps] <- scaEps
            h2DE[h2 < scaEps] <- scaEps
            lsMuImpulseHiddenDE <- lapply(seq(1, scaNImpDE), function(gene) {
                evalImpulse(t = vecTimePointsUniqueB, beta = betaDE[gene], 
                            t1 = t1DE[gene], t2 = t2DE[gene], h0 = h0[gene], 
                            h1 = h1DE[gene], 
                            h2 = h2DE[gene])[vecindTimePointAssignB]
            })
            matMuImpulseHiddenDE <- do.call(rbind, lsMuImpulseHiddenDE)
            rownames(matMuImpulseHiddenDE) <- vecImpulseIDs[1:scaNImpDE]
            colnames(matMuImpulseHiddenDE) <- names(vecSamplesB)
            matMuImpulseHidden[seq(1, scaNImpDE), names(vecSamplesB)] <- 
                matMuImpulseHiddenDE
        }
    } else {
        matMuImpulseHidden <- NULL
        matImpulseModelHidden <- NULL
    }
    
    
    # c. Linear functions Draw linear model parameters
    if (scaNLin > 0) {
        set.seed(scaSeedInit + scaSeedsUsed)
        scaSeedsUsed <- scaSeedsUsed + 1
        vecInitialLevel <- runif(scaNLin) * scaMumax
        set.seed(scaSeedInit + scaSeedsUsed)
        scaSeedsUsed <- scaSeedsUsed + 1
        vecFinalLevel <- vecInitialLevel * 
            abs(rnorm(n = scaNLin, mean = 1, sd = scaSDExpressionChange))
        vecInitialLevel[vecInitialLevel < scaEps] <- scaEps
        vecFinalLevel[vecFinalLevel < scaEps] <- scaEps
        # Evaluate linear functions
        scaDeltaTtot <- max(vecTimePointsUnique) - min(vecTimePointsUnique)
        matMuLinHidden <- do.call(rbind, lapply(seq(1, scaNLin), function(i) {
            (vecInitialLevel[i] + (vecFinalLevel[i] - vecInitialLevel[i])/
                 scaDeltaTtot * (vecTimePointsUnique - min(vecTimePointsUnique))
            )[vecindTimePointAssign]
        }))
        rownames(matMuLinHidden) <- vecLinIDs
        colnames(matMuLinHidden) <- names(vecSamples)
        if (boolCaseCtrl & scaNLin > 1) {
            # Keep start point the same, this doesn't work well for impulse model
            scaNLinDE <- round(scaNLin/2)
            set.seed(scaSeedInit + scaSeedsUsed)
            scaSeedsUsed <- scaSeedsUsed + 1
            vecFinalLevelDE <- vecInitialLevel[1:scaNLinDE] * 
                abs(rnorm(n = scaNLinDE, mean = 1, sd = scaSDExpressionChange))
            vecFinalLevelDE[vecFinalLevelDE < scaEps] <- scaEps
            # Evaluate linear functions
            scaDeltaTtot <- max(vecTimePointsUniqueB) - min(vecTimePointsUniqueB)
            matMuLinHiddenDE <- do.call(rbind, lapply(
                seq(1, scaNLinDE), 
                function(i) { 
                    (vecInitialLevel[i] + (vecFinalLevelDE[i] - 
                                               vecInitialLevel[i])/
                         scaDeltaTtot * (vecTimePointsUniqueB - 
                                             min(vecTimePointsUniqueB))
                    )[vecindTimePointAssignB]
                }))
            rownames(matMuLinHiddenDE) <- vecLinIDs[1:scaNLinDE]
            colnames(matMuLinHiddenDE) <- names(vecSamplesB)
            matMuLinHidden[seq(1, scaNImpDE), names(vecSamplesB)] <- 
                matMuLinHiddenDE
        }
    } else {
        matMuLinHidden <- NULL
    }
    
    # d. Sigmoid functions
    if (scaNSig > 0) {
        # Draw sigmoid model parameters
        set.seed(scaSeedInit + scaSeedsUsed)
        scaSeedsUsed <- scaSeedsUsed + 1
        vech0 <- runif(scaNSig) * scaMumax
        set.seed(scaSeedInit + scaSeedsUsed)
        scaSeedsUsed <- scaSeedsUsed + 1
        vech1 <- vech0 * abs(rnorm(n = scaNSig, mean = 1, 
                                   sd = scaSDExpressionChange))
        set.seed(scaSeedInit + scaSeedsUsed)
        scaSeedsUsed <- scaSeedsUsed + 1
        vecBeta <- runif(scaNSig) * 4 + 0.5
        set.seed(scaSeedInit + scaSeedsUsed)
        scaSeedsUsed <- scaSeedsUsed + 1
        vecT1 <- runif(scaNSig) * (max(vecTimePointsUnique) - 
                                       min(vecTimePointsUnique)) + 
            min(vecTimePointsUnique)
        vech0[vech0 < scaEps] <- scaEps
        vech1[vech1 < scaEps] <- scaEps
        # Evaluate sigmoid functions
        matMuSigHidden <- do.call(rbind, lapply(
            seq(1, scaNSig), function(i) {
                evalSigmoid(t = vecTimePointsUnique, 
                            beta = vecBeta[i], t1 = vecT1[i], 
                            h0 = vech0[i], h1 = vech1[i])[vecindTimePointAssign]
            }))
        rownames(matMuSigHidden) <- vecSigIDs
        colnames(matMuSigHidden) <- names(vecSamples)
        if (boolCaseCtrl & scaNSig > 1) {
            # Keep start point the same, this doesn't work well for impulse model
            scaNSigDE <- round(scaNSig/2)
            set.seed(scaSeedInit + scaSeedsUsed)
            scaSeedsUsed <- scaSeedsUsed + 1
            vech1DE <- vech0[1:scaNSigDE] * 
                abs(rnorm(n = scaNSigDE, mean = 1, 
                          sd = scaSDExpressionChange))
            set.seed(scaSeedInit + scaSeedsUsed)
            scaSeedsUsed <- scaSeedsUsed + 1
            vecBetaDE <- runif(scaNSigDE) * 4 + 0.5
            set.seed(scaSeedInit + scaSeedsUsed)
            scaSeedsUsed <- scaSeedsUsed + 1
            vecT1DE <- runif(scaNSigDE) * (max(
                vecTimePointsUniqueB) - min(vecTimePointsUniqueB)) + 
                min(vecTimePointsUniqueB)
            vech1DE[vech1DE < scaEps] <- scaEps
            # Evaluate sigmoid functions
            matMuSigHiddenDE <- do.call(rbind, lapply(
                seq(1, scaNSigDE), function(i)
                    evalSigmoid(t = vecTimePointsUniqueB, beta = vecBetaDE[i], 
                                t1 = vecT1DE[i], h0 = vech0[i], 
                                h1 = vech1DE[i])[vecindTimePointAssignB]
            ))
            rownames(matMuSigHiddenDE) <- vecSigIDs[1:scaNSigDE]
            colnames(matMuSigHiddenDE) <- names(vecSamplesB)
            matMuSigHidden[seq(1, scaNSigDE), names(vecSamplesB)] <- 
                matMuSigHiddenDE
        }
    } else {
        matMuSigHidden <- NULL
    }
    
    # e. Random data: random fluctuations around mean
    if (scaNRand > 0) {
        set.seed(scaSeedInit + scaSeedsUsed)
        scaSeedsUsed <- scaSeedsUsed + 1
        vecMu <- runif(scaNRand) * scaMumax
        vecMu[vecMu < scaEps] <- scaEps
        if (is.null(scaSDRand)) {
            scaSDRand <- scaSDExpressionChange
        }
        # Evaluate sigmoid functions
        matMuRandHidden <- do.call(rbind, lapply(seq(1, scaNRand), function(i) {
            set.seed(scaSeedInit + scaSeedsUsed)
            scaSeedsUsed <- scaSeedsUsed + 1
            vecRands <- vecMu[i] * abs(rnorm(n = length(vecSamples), mean = 1, 
                                             sd = scaSDRand))
            vecRands[vecRands < scaEps] <- scaEps
            return(vecRands)
        }))
        rownames(matMuRandHidden) <- vecRandIDs
        colnames(matMuRandHidden) <- names(vecSamples)
    } else {
        matMuRandHidden <- NULL
    }
    
    # f. Merge data
    matMuHidden <- do.call(rbind, list(
        matMuConstHidden, matMuImpulseHidden, 
        matMuLinHidden, matMuSigHidden, matMuRandHidden))
    if (boolCaseCtrl) {
        vecCaseCtrlDEIDs <- c()
        if (scaNImp > 1) {
            vecCaseCtrlDEIDs <- c(vecCaseCtrlDEIDs, vecImpulseIDs[1:scaNImpDE])
        }
        if (scaNLin > 1) {
            vecCaseCtrlDEIDs <- c(vecCaseCtrlDEIDs, vecLinIDs[1:scaNLinDE])
        }
        if (scaNSig > 1) {
            vecCaseCtrlDEIDs <- c(vecCaseCtrlDEIDs, vecSigIDs[1:scaNSigDE])
        }
    }
    
    # Add scaling factors a) Sample size factors
    set.seed(scaSeedInit + scaSeedsUsed)
    scaSeedsUsed <- scaSeedsUsed + 1
    vecSizeFactorsHidden <- rnorm(n = scaNSamples, mean = scaMuSizeEffect, 
                                  sd = scaSDSizeEffect)
    vecSizeFactorsHidden[vecSizeFactorsHidden < 0.1] <- 0.1
    vecSizeFactorsHidden[vecSizeFactorsHidden > 10] <- 10
    names(vecSizeFactorsHidden) <- names(vecSamples)
    # Scale by size factors
    matMuHiddenScaled <- matMuHidden * matrix(
        vecSizeFactorsHidden, nrow = dim(matMuHidden)[1], 
        ncol = dim(matMuHidden)[2], byrow = TRUE)
    # b) Batch factors
    vecBatchesUnique <- unique(dfAnnotation$Batch)
    vecindBatches <- match(dfAnnotation$Batch, vecBatchesUnique)
    # Check that batches dont coincide with case-ctrl split of samples
    if (boolCaseCtrl & length(vecBatchesUnique) == 2) {
        if (setequal(dfAnnotation[dfAnnotation$Batch == 
                                  vecBatchesUnique[1], ]$Sample, 
                     dfAnnotation[dfAnnotation$Condition == 
                                  "case", ]$Sample) | 
            setequal(dfAnnotation[dfAnnotation$Batch == 
                                  vecBatchesUnique[1], ]$Sample, 
                     dfAnnotation[dfAnnotation$Condition == 
                                  "ctrl", ]$Sample)) {
            stop("Batch structure coincides with case-control structure.")
        }
    }
    if (length(vecBatchesUnique) > 1) {
        matBatchFactorsUnqiueHidden <- do.call(rbind, lapply(
            seq(1, scaNGenes), function(i) {
                set.seed(scaSeedInit + scaSeedsUsed)
                scaSeedsUsed <- scaSeedsUsed + 1
                return(c(1, rnorm(n = length(vecBatchesUnique) - 1, 
                                  mean = scaMuBatchEffect, 
                                  sd = scaSDBatchEffect)))
            }))
        
        matBatchFactorsUnqiueHidden[matBatchFactorsUnqiueHidden < 0.1] <- 0.1
        matBatchFactorsUnqiueHidden[matBatchFactorsUnqiueHidden > 10] <- 10
        matBatchFactorsHidden <- matBatchFactorsUnqiueHidden[, vecindBatches]
        rownames(matBatchFactorsHidden) <- rownames(matMuHidden)
        # Scale
        matMuHiddenScaled <- matMuHiddenScaled * matBatchFactorsHidden
    }
    rownames(matMuHiddenScaled) <- rownames(matMuHidden)
    colnames(matMuHiddenScaled) <- colnames(matMuHidden)
    
    # e. draw dispersions by gene - one per condition in case-control
    set.seed(scaSeedInit + scaSeedsUsed)
    scaSeedsUsed <- scaSeedsUsed + 1
    vecDispHidden <- rnorm(dim(matMuHiddenScaled)[1], mean = 1, sd = 0.1) * 
        10
    names(vecDispHidden) <- rownames(matMuHiddenScaled)
    vecDispHidden[vecDispHidden < 1] <- 1
    
    # f. add noise - draw from negative binomial
    matObservedData <- do.call(rbind, lapply(
        seq(1, dim(matMuHiddenScaled)[1]), 
        function(gene) {
            sapply(names(vecSamples), function(j) {
                rnbinom(n = 1, mu = matMuHiddenScaled[gene, j], 
                        size = vecDispHidden[gene])
            })
        }))
    rownames(matObservedData) <- rownames(matMuHiddenScaled)
    colnames(matObservedData) <- colnames(matMuHiddenScaled)
    
    # Counts
    matObservedCounts <- round(matObservedData)
    
    # Save simulation
    if (!is.null(dirOutSimulation)) {
        save(vecSamples, file = file.path(dirOutSimulation, 
                                          "Simulation_vecPT.RData"))
        save(vecConstIDs, file = file.path(dirOutSimulation, 
                                           "Simulation_vecConstIDs.RData"))
        save(vecImpulseIDs, file = file.path(dirOutSimulation, 
                                             "Simulation_vecImpulseIDs.RData"))
        save(vecLinIDs, file = file.path(dirOutSimulation, 
                                         "Simulation_vecLinIDs.RData"))
        save(vecSigIDs, file = file.path(dirOutSimulation, 
                                         "Simulation_vecSigIDs.RData"))
        save(vecSigIDs, file = file.path(dirOutSimulation, 
                                         "Simulation_vecRandIDs.RData"))
        
        save(vecSizeFactorsHidden, 
             file = file.path(dirOutSimulation, 
                              "Simulation_vecSizeFactorsHidden.RData"))
        if (length(vecBatchesUnique) > 1) {
            save(matBatchFactorsHidden, 
                 file = file.path(dirOutSimulation, 
                                  "Simulation_matBatchFactorsHidden.RData"))
        }
        if (boolCaseCtrl) {
            save(vecCaseCtrlDEIDs, 
                 file = file.path(dirOutSimulation, 
                                  "Simulation_vecCaseCtrlDEIDs.RData"))
        }
        
        save(vecDispHidden, 
             file = file.path(dirOutSimulation, 
                              "Simulation_vecDispHidden.RData"))
        save(matImpulseModelHidden, 
             file = file.path(dirOutSimulation, 
                              "Simulation_matImpulseModelHidden.RData"))
        save(matMuHidden, 
             file = file.path(dirOutSimulation, 
                              "Simulation_matMuHidden.RData"))
        save(matMuHiddenScaled, 
             file = file.path(dirOutSimulation, 
                              "Simulation_matMuHiddenScaled.RData"))
        save(matObservedData, 
             file = file.path(dirOutSimulation, 
                              "Simulation_matObservedData.RData"))
        save(matObservedCounts, 
             file = file.path(dirOutSimulation, 
                              "Simulation_matObservedCounts.RData"))
    }
    
    return(list(dfAnnotation = dfAnnotation, 
                matObservedCounts = matObservedCounts))
}
