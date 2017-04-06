#### Perform differential expression analysis and identify transient genes

#' Perform differential expression analysis and identification
#' of transiently activated or deactivated genes.
#' 
#' Performs model selction based on loglikelihood ratio tests.
#' The primary model selection is the differential expression analysis.
#' The secondary model selection is the selection between a sigmoidal
#' and an impulse fit for differentially expressed genes which is used
#' to define transiently activated or deactivated genes.
#' 
#' @seealso Called by \link{runImpulseDE2}.
#' 
#' @param objectImpulseDE2 (object class ImpulseDE2Object)
#' Object containing fits to be evaluated.
#' @param boolCaseCtrl (bool) 
#' Whether to perform case-control analysis. Does case-only
#' analysis if FALSE.
#' @param boolIdentifyTransients (bool) [Defaul FALSE]
#' Whether to identify transiently activated or deactivated 
#' genes. This involves an additional fitting of sigmoidal models
#' and hypothesis testing between constant, sigmoidal and impulse model.
#' @param scaQThresTransients (scalar) [Default 0.001]
#' FDR-corrected p-value threshold for hypothesis tests between
#' impulse, sigmoidal and constant model used to identify transiently
#' regulated genes.
#' 
#' @return (ImpulseDE2Object)
#' Input object with dfDEAnalysis updated to:
#' dfDEAnalysis (data frame samples x reported characteristics) 
#' Summary of fitting procedure and 
#' differential expression results for each gene.
#' \itemize{
#' \item Gene: Gene ID.
#' \item p: P-value for differential expression.
#' \item padj: Benjamini-Hochberg false-discovery rate corrected p-value
#' for differential expression analysis.
#' \item loglik_full: Loglikelihood of full model.
#' \item loglik_red: Loglikelihood of reduced model.
#' \item df_full: Degrees of freedom of full model.
#' \item df_red: Degrees of freedom of reduced model
#' \item mean: Inferred mean parameter of constant model of first batch.
#' From combined samples in case-ctrl. 
#' \item allZero (bool) Whether there were no observed non-zero 
#' observations of this gene.
#' If TRUE, fitting and DE analsysis were skipped and entry is NA.
#' }
#' Entries only present in case-only DE analysis:
#' \itemize{
#' \item converge_impulse: Convergence status of optim for 
#' impulse model fit (full model).
#' \item converge_const: Convergence status of optim for 
#' constant model fit (reduced model).
#' }
#' Entries only present in case-control DE analysis:
#' \itemize{
#' \item converge_combined: Convergence status of optim for 
#' impulse model fit to case and control samples combined (reduced model).
#' \item converge_case: Convergence status of optim for 
#' impulse model fit to samples of case condition (full model 1/2).
#' \item converge_control: Convergence status of optim for 
#' impulse model fit to samples of control condition (full model 2/2).
#' }
#' Entries only present if boolIdentifyTransients is TRUE:
#' \itemize{
#' \item converge_sigmoid: Convergence status of optim for 
#' sigmoid model fit to samples of case condition.
#' \item impulseTOsigmoid_p: P-value of loglikelihood ratio test
#' impulse model fit versus sigmoidal model on samples of case condition.
#' \item impulseTOsigmoid_padj: Benjamini-Hochberg 
#' false-discovery rate corrected p-value of loglikelihood ratio test
#' impulse model fit versus sigmoid model on samples of case condition.
#' \item sigmoidTOconst_p: P-value of loglikelihood ratio test
#' sigmoidal model fit versus constant model on samples of case condition.
#' \item sigmoidTOconst_padj: Benjamini-Hochberg 
#' false-discovery rate corrected p-value of loglikelihood ratio test
#' sigmoidal model fit versus constant model on samples of case condition.
#' \item isTransient (bool) Whether gene is transiently
#' activated or deactivated and differentially expressed.
#' \item isMonotonous (bool) Whether gene is not transiently
#' activated or deactivated and differentially expressed. This scenario
#' corresponds to a montonous expression level increase or decrease.
#' }
#' 
#' @author David Sebastian Fischer
runDEAnalysis <- function(
    objectImpulseDE2, boolCaseCtrl, 
    boolIdentifyTransients, scaQThresTransients = 0.001) {
    
    dfAnnot <- get_dfAnnotationProc(obj = objectImpulseDE2)
    
    if (!boolCaseCtrl) {
        # Case-only: Take Values from the only fitting performed: case The full
        # model is the impulse fit.
        vecLogLikFull <- sapply(get_lsModelFits(objectImpulseDE2)$case, 
                                function(fit) fit$lsImpulseFit$scaLL)
        # The reduced model is the mean fit.
        vecLogLikRed <- sapply(get_lsModelFits(objectImpulseDE2)$case, 
                               function(fit) fit$lsConstFit$scaLL)
        # Mean inferred expression:
        vecMuCase <- sapply(get_lsModelFits(objectImpulseDE2)$case, 
                            function(fit) fit$lsConstFit$scaMu)
        if (!is.null(get_vecConfounders(objectImpulseDE2))) {
            scaNBatchFactors <- sum(sapply(
                get_vecConfounders(objectImpulseDE2), function(confounder)
                length(unique(dfAnnot[
                    , confounder])) - 1
            ))
        } else {
            scaNBatchFactors <- 0
        }
        # 6 impulse model parameters, 1 dispersion estimate and 1 batch factor
        # for each batch (except for the first one) for each confounder.
        scaDegFreedomFull <- 6 + 1 + scaNBatchFactors
        # 1 constant model parameter1, 1 dispersion estimate and 1 batch factor
        # for each batch (except for the first one) for each confounder.
        scaDegFreedomRed <- 1 + 1 + scaNBatchFactors
        
        vecConvergenceImpulse <- sapply(
            get_lsModelFits(objectImpulseDE2)$case, 
            function(fit) fit$lsImpulseFit$scaConvergence)
        vecConvergenceConst <- sapply(
            get_lsModelFits(objectImpulseDE2)$case, 
            function(fit) fit$lsConstFit$scaConvergence)
    } else {
        # Case-control: Full model: Case and control model separate: The log
        # likelihood of the full model is the sum of the log likelihoods of
        # case and control fits.
        vecLogLikFull <- sapply(get_lsModelFits(objectImpulseDE2)$case, 
                                function(fit) fit$lsImpulseFit$scaLL) + 
            sapply(get_lsModelFits(objectImpulseDE2)$control, 
                   function(fit) fit$lsImpulseFit$scaLL)
        # The reduced model is the combined data fit.
        vecLogLikRed <- sapply(get_lsModelFits(obj=objectImpulseDE2)$combined, 
                               function(fit) fit$lsImpulseFit$scaLL)
        # Mean inferred expression: On combined data
        vecMuCombined <- sapply(get_lsModelFits(obj=objectImpulseDE2)$combined, 
                                function(fit) fit$lsConstFit$scaMu)
        if (!is.null(get_vecConfounders(obj=objectImpulseDE2))) {
            scaNBatchFactorsFull <- sum(sapply(
                get_vecConfounders(obj=objectImpulseDE2), function(confounder) {
                    length(unique(dfAnnot[dfAnnot$Condition == "case", 
                                          confounder])) - 1 + 
                        length(unique(dfAnnot[dfAnnot$Condition == "control", 
                                              confounder])) - 1
                }))
            scaNBatchFactorsRed <- sum(sapply(
                get_vecConfounders(obj=objectImpulseDE2), function(confounder) {
                length(unique(dfAnnot[, confounder])) - 1
            }))
        } else {
            scaNBatchFactorsFull <- 0
            scaNBatchFactorsRed <- 0
        }
        # 6 impulse model parameters for each case and control, 1 dispersion
        # estimate and 1 batch factor for each batch (except for the first one)
        # for each confounder in each condition.
        scaDegFreedomFull <- 6 * 2 + 1 + scaNBatchFactorsFull
        # 6 impulse model parameters, 1 dispersion estimate and 1 batch factor
        # for each batch (except for the first one) for each confounder.
        scaDegFreedomRed <- 6 + 1 + scaNBatchFactorsRed
        vecConvergenceImpulseCombined <- sapply(
            get_lsModelFits(obj=objectImpulseDE2)$combined, 
            function(fit) fit$lsImpulseFit$scaConvergence)
        vecConvergenceImpulseCase <- sapply(
            get_lsModelFits(obj=objectImpulseDE2)$case, 
            function(fit) fit$lsImpulseFit$scaConvergence)
        vecConvergenceImpulseControl <- sapply(
            get_lsModelFits(obj=objectImpulseDE2)$control, 
            function(fit) fit$lsImpulseFit$scaConvergence)
    }
    
    # Compute difference in degrees of freedom between null model and
    # alternative model.
    scaDeltaDegFreedom <- scaDegFreedomFull - scaDegFreedomRed
    # Compute test statistic: Deviance
    vecDeviance <- 2 * (vecLogLikFull - vecLogLikRed)
    # Get p-values from Chi-square distribution (assumption about null
    # model)
    vecPvalue <- pchisq(vecDeviance, 
                        df = scaDeltaDegFreedom, lower.tail = FALSE)
    # Multiple testing correction (Benjamini-Hochberg)
    vecPvalueBH <- p.adjust(vecPvalue, method = "BH")
    
    scaNGenes <- dim(get_matCountDataProc(obj=objectImpulseDE2))[1]
    if (!boolCaseCtrl) {
        # Case-only
        dfDEAnalysis = data.frame(
            Gene = row.names(get_matCountDataProc(obj=objectImpulseDE2)), 
            p = vecPvalue, padj = vecPvalueBH, loglik_full = vecLogLikFull, 
            loglik_red = vecLogLikRed, 
            df_full = rep(scaDegFreedomFull, scaNGenes), 
            df_red = rep(scaDegFreedomRed, scaNGenes), 
            mean = vecMuCase, 
            converge_impulse = vecConvergenceImpulse, 
            converge_const = vecConvergenceConst, 
            stringsAsFactors = FALSE)
    } else {
        # Case-control
        dfDEAnalysis = data.frame(
            Gene = row.names(get_matCountDataProc(obj=objectImpulseDE2)), 
            p = as.numeric(vecPvalue), padj = as.numeric(vecPvalueBH), 
            loglik_full = vecLogLikFull, loglik_red = vecLogLikRed, 
            df_full = rep(scaDegFreedomFull, scaNGenes), 
            df_red = rep(scaDegFreedomRed, scaNGenes), 
            mean = vecMuCombined, 
            converge_combined = vecConvergenceImpulseCombined, 
            converge_case = vecConvergenceImpulseCase, 
            converge_control = vecConvergenceImpulseControl, 
            stringsAsFactors = FALSE)
    }
    rownames(dfDEAnalysis) <- row.names(get_matCountDataProc(
        obj=objectImpulseDE2))
    
    # Add entries if sigmoidal fit was performed and transient tranjectories
    # are to be differentiated from monotonous ones. Do this in case
    # condition.
    if (boolIdentifyTransients) {
        # Take the sigmoid model as reference if transient are tested versus
        # montonous model
        vecLogLikImpulse <- sapply(get_lsModelFits(obj=objectImpulseDE2)$case, 
                                   function(fit) fit$lsImpulseFit$scaLL)
        vecLogLikSigmoid <- sapply(get_lsModelFits(obj=objectImpulseDE2)$case, 
                                   function(fit) fit$lsSigmoidFit$scaLL)
        vecLogLikConst <- sapply(get_lsModelFits(obj=objectImpulseDE2)$case, 
                                 function(fit) fit$lsConstFit$scaLL)
        if (boolCaseCtrl) {
            # Get number of batch parameters only in case condition
            if(!is.null(get_vecConfounders(obj=objectImpulseDE2))) {
                scaNBatchFactors <- sum(sapply(
                    get_vecConfounders(obj=objectImpulseDE2), function(confounder) {
                        length(unique(dfAnnot[dfAnnot$Condition == "case", 
                                              confounder])) - 1
                    }))
            } else {
                scaNBatchFactors <- 0
            }
        }
        scaDegFreedomImpulse <- 6 + 1 + scaNBatchFactors
        scaDegFreedomSigmoid <- 4 + 1 + scaNBatchFactors
        scaDegFreedomConst <- 1 + 1 + scaNBatchFactors
        
        # Compare impulse to sigmoid
        vecDevianceImpulseSigmoid <- 2 * (vecLogLikImpulse - vecLogLikSigmoid)
        vecPvalueImpulseSigmoid <- pchisq(
            vecDevianceImpulseSigmoid, df = scaDegFreedomImpulse - 
                scaDegFreedomSigmoid, lower.tail = FALSE)
        vecPvalueImpulseSigmoidBH <- p.adjust(vecPvalueImpulseSigmoid, 
                                              method = "BH")
        
        # Compare impulse to constant This is the same as the differential
        # expression p-value if DE analysis is case-only
        vecDevianceImpulseConst <- 2 * (vecLogLikImpulse - vecLogLikConst)
        vecPvalueImpulseConst <- pchisq(
            vecDevianceImpulseConst, 
            df = scaDegFreedomImpulse - scaDegFreedomConst, lower.tail = FALSE)
        vecPvalueImpulseConstBH <- p.adjust(vecPvalueImpulseConst, 
                                            method = "BH")
        
        # Compare sigmoid to constant Trajectory is montonous if this is
        # significant but vecPvalueImpulseSigmoidBH is not
        vecDevianceSigmoidConst <- 2 * (vecLogLikSigmoid - vecLogLikConst)
        vecPvalueSigmoidConst <- pchisq(
            vecDevianceSigmoidConst, 
            df = scaDegFreedomSigmoid - scaDegFreedomConst, lower.tail = FALSE)
        vecPvalueSigmoidConstBH = p.adjust(vecPvalueSigmoidConst, 
                                           method = "BH")
        
        # Add entries into DE table
        dfDEAnalysis$converge_sigmoid <- sapply(
            get_lsModelFits(obj=objectImpulseDE2)$case, 
            function(fit) fit$lsSigmoidFit$scaConvergence)
        dfDEAnalysis$impulseTOsigmoid_p <- 
            as.numeric(vecPvalueImpulseSigmoid)
        dfDEAnalysis$impulseTOsigmoid_padj <- 
            as.numeric(vecPvalueImpulseSigmoidBH)
        dfDEAnalysis$sigmoidTOconst_p <- 
            as.numeric(vecPvalueSigmoidConst)
        dfDEAnalysis$sigmoidTOconst_padj <- 
            as.numeric(vecPvalueSigmoidConstBH)
        # Classify trajectories as tranient change or monotonous change
        # (transition).  Note that significant impulse vs sigmoid hits include
        # monotonous fits which are better fit by impulse than by sigmoid, this
        # is corrected for here.
        vecTimePointsCase <- sort(
            get_lsModelFits(obj=objectImpulseDE2)$IdxGroups$case$vecTimepointsUnique, 
            decreasing = FALSE)
        vecboolMonotonousImpulseTraject <- sapply(
            get_lsModelFits(obj=objectImpulseDE2)$case, 
            function(fit) {
                # Do not use stored per sample fits from 
                # get_lsModelFits(obj=objectImpulseDE2) here: 
                # These fits contain batch factors and we are only interested in
                # the raw structure of the fit.
                vecImpulseValues <- evalImpulse_comp(
                    vecImpulseParam = fit$lsImpulseFit$vecImpulseParam, 
                    vecTimepoints = vecTimePointsCase)
                boolMonotonous <- 
                    (max(vecImpulseValues[2:length((vecImpulseValues) - 1)]) <=
                         max(vecImpulseValues[c(1, length(vecImpulseValues))])) & 
                    (min(vecImpulseValues[2:length((vecImpulseValues) - 1)]) >= 
                         min(vecImpulseValues[c(1, length(vecImpulseValues))]))
                return(boolMonotonous)
            })
        # Is transient if significantly better fit by impulse than sigmoidal
        # model and not monotonous
        dfDEAnalysis$isTransient <- vecPvalueImpulseSigmoidBH <= 
            scaQThresTransients & 
            !vecboolMonotonousImpulseTraject
        # Is transition (monotonous) if (1) not significantly better fit by
        # impulse than sigmoidal model or is monotonous trajectory (which is
        # better fit by impulse than sigmoid) and (2) is differentially
        # expressed, defined as significantly better fit by sigmoid than
        # constant.
        dfDEAnalysis$isMonotonous <- (
            vecPvalueImpulseSigmoidBH > scaQThresTransients | 
                vecboolMonotonousImpulseTraject) & 
            vecPvalueSigmoidConstBH <= 
            scaQThresTransients
    }
    
    dfDEAnalysis <- dfDEAnalysis[match(get_vecAllIDs(obj=objectImpulseDE2), 
                                       dfDEAnalysis$Gene), ]
    rownames(dfDEAnalysis) <- get_vecAllIDs(obj=objectImpulseDE2)
    vecboolAllZero <- !(get_vecAllIDs(obj=objectImpulseDE2) %in% 
                            rownames(get_matCountDataProc(
                                obj=objectImpulseDE2)))
    dfDEAnalysis$allZero <- vecboolAllZero
    
    objectImpulseDE2 <- set_dfImpulseDE2Results(obj=objectImpulseDE2, 
                                                element=dfDEAnalysis)
    return(objectImpulseDE2)
}

#' Update dfImpulseDE2Results after sigmoids have been fit 
#' through external call
#' 
#' This is a userfriendly wrapper of runDEAnalysis for this update
#' scenario.
#' 
#' @seealso Called by separately by user.
#' 
#' @param objectImpulseDE2 (object class ImpulseDE2Object)
#' Object containing fits to be evaluated.
#' @param scaQThresTransients (scalar) [Default 0.001]
#' FDR-corrected p-value threshold for hypothesis tests between
#' impulse, sigmoidal and constant model used to identify transiently
#' regulated genes.
#' 
#' @return objectImpulseDE2 (ImpulseDE2Object)
#' Input object with dfDEAnalysis updated to:
#' dfDEAnalysis (data frame samples x reported characteristics) 
#' Summary of fitting procedure and 
#' differential expression results for each gene.
#' \itemize{
#' \item Gene: Gene ID.
#' \item p: P-value for differential expression.
#' \item padj: Benjamini-Hochberg false-discovery rate corrected p-value
#' for differential expression analysis.
#' \item loglik_full: Loglikelihood of full model.
#' \item loglik_red: Loglikelihood of reduced model.
#' \item df_full: Degrees of freedom of full model.
#' \item df_red: Degrees of freedom of reduced model
#' \item mean: Inferred mean parameter of constant model of first batch.
#' From combined samples in case-ctrl. 
#' \item allZero (bool) Whether there were no observed non-zero 
#' observations of this gene.
#' If TRUE, fitting and DE analsysis were skipped and entry is NA.
#' }
#' Entries only present in case-only DE analysis:
#' \itemize{
#' \item converge_impulse: Convergence status of optim for 
#' impulse model fit (full model).
#' \item converge_const: Convergence status of optim for 
#' constant model fit (reduced model).
#' }
#' Entries only present in case-control DE analysis:
#' \itemize{
#' \item converge_combined: Convergence status of optim for 
#' impulse model fit to case and control samples combined (reduced model).
#' \item converge_case: Convergence status of optim for 
#' impulse model fit to samples of case condition (full model 1/2).
#' \item converge_control: Convergence status of optim for 
#' impulse model fit to samples of control condition (full model 2/2).
#' }
#' Entries only present if boolIdentifyTransients is TRUE:
#' \itemize{
#' \item converge_sigmoid: Convergence status of optim for 
#' sigmoid model fit to samples of case condition.
#' \item impulseTOsigmoid_p: P-value of loglikelihood ratio test
#' impulse model fit versus sigmoidal model on samples of case condition.
#' \item impulseTOsigmoid_padj: Benjamini-Hochberg 
#' false-discovery rate corrected p-value of loglikelihood ratio test
#' impulse model fit versus sigmoid model on samples of case condition.
#' \item sigmoidTOconst_p: P-value of loglikelihood ratio test
#' sigmoidal model fit versus constant model on samples of case condition.
#' \item sigmoidTOconst_padj: Benjamini-Hochberg 
#' false-discovery rate corrected p-value of loglikelihood ratio test
#' sigmoidal model fit versus constant model on samples of case condition.
#' \item isTransient (bool) Whether gene is transiently
#' activated or deactivated and differentially expressed.
#' \item isMonotonous (bool) Whether gene is not transiently
#' activated or deactivated and differentially expressed. This scenario
#' corresponds to a montonous expression level increase or decrease.
#' }
#'  
#' @examples
#' lsSimulatedData <- simulateDataSetImpulseDE2(
#' vecTimePointsA   = rep(seq(1,8),3),
#' vecTimePointsB   = NULL,
#' vecBatchesA      = NULL,
#' vecBatchesB      = NULL,
#' scaNConst        = 0,
#' scaNImp          = 50,
#' scaNLin          = 0,
#' scaNSig          = 50)
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
updateDEAnalysis <- function(objectImpulseDE2, scaQThresTransients = 0.001) {
    
    objectImpulseDE2 <- runDEAnalysis(
        objectImpulseDE2 = objectImpulseDE2, 
        boolCaseCtrl = get_boolCaseCtrl(obj=objectImpulseDE2), 
        boolIdentifyTransients = TRUE, scaQThresTransients = scaQThresTransients)
    
    return(objectImpulseDE2)
}
