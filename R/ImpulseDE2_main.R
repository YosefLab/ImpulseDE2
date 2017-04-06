#' @import Biobase
#' @import BiocParallel
#' @import circlize
#' @importFrom compiler cmpfun
#' @import ComplexHeatmap
#' @importFrom cowplot plot_grid
#' @import DESeq2
#' @import ggplot2
#' @importFrom grDevices dev.off graphics.off pdf
#' @import knitr
#' @import Matrix
#' @import methods
#' @importFrom stats dnbinom median optim p.adjust pchisq rnbinom rnorm runif sd time
#' @import SummarizedExperiment
#' @importFrom utils packageDescription
NULL

### Main function / Wrapper function

#' ImpulseDE2 wrapper
#'
#' Wrapper to run ImpulseDE2 on bulk omics count data.
#' This wrapper can perform the entire analysis pipeline of 
#' ImpulseDE2 on its own if the right parameters are supplied.
#' To run ImpulseDE2 on bulk omics count data, use the minimal
#' parameter set:
#' \itemize{
#' \item matCountData
#' \item dfAnnotation
#' \item boolCaseCtrl
#' \item vecConfounders
#' }
#' Additionally, you can provide:
#' \itemize{
#' \item scaNProc to set the number of processes for parallelisation.
#' \item scaQThres to set the cut off for your DE gene list. 
#' \item vecDispersionsExternal to supply external dispersion parameters
#' which may be necessary depending on your confounding factors (runImpulseDE2
#' will tell you if it is necessary).
#' \item vecSizeFactorsExternal to supply external size factors.
#' \item boolVerbose to control stdout output.
#' }
#' 
#' @details ImpulseDE2 is based on the impulse model proposed by
#' Chechik and Koller (Chechik and Koller, 2009).
#' The computational complexity of ImpulseDE2 is linear in the
#' number of genes and linear in the number of samples.
#' 
#' @aliases ImpulseDE2 wrapper
#' 
#' @seealso Calls the following functions:
#' \link{processData}, 
#' \link{runDESeq2},
#' \link{computeNormConst},
#' \link{fitModels},
#' \link{fitSigmoidModels},
#' \link{runDEAnalysis}. 
#' The following functions are additionally available to the user:
#' \link{fitSigmoidModels},
#' \link{plotGenes},  
#' \link{plotHeatmap},
#' \link{runDEAnalysis},   
#' \link{simulateDataSetImpulseDE2}.
#' 
#' @param matCountData (matrix genes x samples) [Default NULL] 
#' Read count data, unobserved entries are NA. 
#' Can be SummarizedExperiment object.
#' @param dfAnnotation (data frame samples x covariates) 
#' {Sample, Condition, Time (numeric), TimeCateg (str)
#' (and confounding variables if given).}
#' Annotation table with covariates for each sample.
#' @param boolCaseCtrl (bool) [Default FALSE]
#' Whether to perform case-control analysis. Does case-only
#' analysis if FALSE.
#' @param vecConfounders (vector of strings number of confounding variables)
#' Factors to correct for during batch correction. Have to 
#' supply dispersion factors if more than one is supplied.
#' Names refer to columns in dfAnnotation.
#' @param scaNProc (scalar) [Default 1] Number of processes for 
#' parallelisation.
#' @param scaQThres (scalar) [Default NULL] 
#' FDR-corrected p-value cutoff for significance.
#' @param vecDispersionsExternal (vector length number of
#' genes in matCountData) [Default NULL]
#' Externally generated list of gene-wise dispersion factors
#' which overides DESeq2 generated dispersion factors.
#' @param vecSizeFactorsExternal (vector length number of
#' cells in matCountData) [Default NULL]
#' Externally generated list of size factors which override
#' size factor computation in ImpulseDE2.
#' @param boolIdentifyTransients (bool) [Defaul FALSE]
#' Whether to identify transiently activated or deactivated 
#' genes. This involves an additional fitting of sigmoidal models
#' and hypothesis testing between constant, sigmoidal and impulse model.
#' @param boolVerbose (bool) [Default TRUE] Whether to print
#' progress to stdout.
#' 
#' @return (object of class ImpulseDE2Object)
#' This object can be treated as a list with 2 elements:
#' (list length 2)
#' \itemize{
#' \item vecDEGenes (list number of genes) Genes IDs identified
#' as differentially expressed by ImpulseDE2 at threshold \code{scaQThres}.
#' \item dfDEAnalysis (data frame samples x reported characteristics) 
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
#' objectImpulseDE2 <- runImpulseDE2(
#' matCountData    = lsSimulatedData$matObservedCounts, 
#' dfAnnotation    = lsSimulatedData$dfAnnotation,
#' boolCaseCtrl    = FALSE,
#' vecConfounders  = NULL,
#' scaNProc        = 1 )
#' head(objectImpulseDE2$dfImpulseDE2Results)
#' 
#' @author David Sebastian Fischer
#' 
#' @export
runImpulseDE2 <- function(
    matCountData = NULL, dfAnnotation = NULL, boolCaseCtrl = FALSE, 
    vecConfounders = NULL, scaNProc = 1, scaQThres = NULL, 
    vecDispersionsExternal = NULL, 
    vecSizeFactorsExternal = NULL, boolIdentifyTransients = FALSE, 
    boolVerbose = TRUE) {
    
    strMessage <- paste0("ImpulseDE2 for count data, v", 
                         packageDescription("ImpulseDE2", fields = "Version"))
    if (boolVerbose) { message(strMessage) }
    strReport <- strMessage
    
    tm_runImpulseDE2 <- system.time({
        # 1. Process input data
        strMessage <- "# Process input"
        if (boolVerbose) { message(strMessage) }
        strReport <- paste0(strReport, "\n", strMessage)
        # Extract count matrix if handed SummarizedExperiment
        if (class(matCountData) == "SummarizedExperiment"){ 
            matCountData <- assay(matCountData)
        }
        lsProcessedData <- processData(
            dfAnnotation = dfAnnotation, matCountData = matCountData, 
            boolCaseCtrl = boolCaseCtrl, vecConfounders = vecConfounders, 
            vecDispersionsExternal = vecDispersionsExternal, 
            vecSizeFactorsExternal = vecSizeFactorsExternal)
        
        matCountDataProc <- lsProcessedData$matCountDataProc
        dfAnnotationProc <- lsProcessedData$dfAnnotationProc
        vecSizeFactorsExternalProc <- 
            lsProcessedData$vecSizeFactorsExternalProc
        vecDispersionsExternalProc <- 
            lsProcessedData$vecDispersionsExternalProc
        if (boolVerbose) { 
            write(lsProcessedData$strReportProcessing, file = "", ncolumns = 1)
        }
        strReport <- paste0(strReport, lsProcessedData$strReportProcessing)
        
        # Set parallelisation
        if (scaNProc > 1) {
            register(MulticoreParam(workers = scaNProc))
        } else {
            register(SerialParam())
        }
        
        # 2. Run DESeq2 Use dispersion factors from DESeq2 if
        if (is.null(vecDispersionsExternal)) {
            strMessage <- paste0("# Run DESeq2: Using dispersion factors",
                                 "computed by DESeq2.")
            if (boolVerbose) { message(strMessage) }
            strReport <- paste0(strReport, "\n", strMessage)
            tm_runDESeq2 <- system.time({
                vecDispersions <- runDESeq2(
                    dfAnnotationProc = dfAnnotationProc, 
                    matCountDataProc = matCountDataProc, 
                    boolCaseCtrl = boolCaseCtrl, 
                    vecConfounders = vecConfounders)
            })
            strMessage <- paste0("Consumed time: ", 
                                 round(tm_runDESeq2["elapsed"]/60,2), " min.")
            if (boolVerbose) { message(strMessage) }
            strReport <- paste0(strReport, "\n", strMessage)
        } else {
            # Use externally provided dispersions and reorder
            strMessage <- "# Using externally supplied dispersion factors."
            if (boolVerbose) { message(strMessage) }
            strReport <- paste0(strReport, "\n", strMessage)
            vecDispersions <- vecDispersionsExternalProc
        }
        
        # 3. Compute size factors
        strMessage <- "# Compute size factors"
        if (boolVerbose) { message(strMessage) }
        strReport <- paste0(strReport, "\n", strMessage)
        vecSizeFactors <- computeNormConst(
            matCountDataProc = matCountDataProc, 
            vecSizeFactorsExternal = vecSizeFactorsExternalProc)
        
        # 4. Create instance of ImpulseDE2Object Create ImpulseDE2 object
        objectImpulseDE2 <- new(
            "ImpulseDE2Object", dfImpulseDE2Results = NULL, 
            vecDEGenes = NULL, lsModelFits = NULL, 
            matCountDataProc = matCountDataProc, 
            vecAllIDs = rownames(matCountData), 
            dfAnnotationProc = dfAnnotationProc, 
            vecSizeFactors = vecSizeFactors, vecDispersions = vecDispersions, 
            boolCaseCtrl = boolCaseCtrl, vecConfounders = vecConfounders, 
            scaNProc = scaNProc, scaQThres = scaQThres, strReport = strReport)
        
        # 5. Fit null and alternative model to each gene
        strMessage <- "# Fitting null and alternative model to the genes"
        if (boolVerbose) { message(strMessage) }
        objectImpulseDE2 <- append_strReport(obj = objectImpulseDE2,
                                             s = strMessage)
        tm_fitImpulse <- system.time({
            objectImpulseDE2 <- fitModels(
                objectImpulseDE2 = objectImpulseDE2, 
                vecConfounders = vecConfounders, boolCaseCtrl = boolCaseCtrl)
        })
        strMessage <- paste0("Consumed time: ", 
                             round(tm_fitImpulse["elapsed"]/60, 2), " min.")
        if (boolVerbose) { message(strMessage) }
        objectImpulseDE2 <- append_strReport(obj = objectImpulseDE2,
                                             s = strMessage)
        
        # 6. Fit sigmoid model to case condition if desired
        if (boolIdentifyTransients) {
            strMessage <- "# Fitting sigmoid model to case condition"
            if (boolVerbose) { message(strMessage) }
            objectImpulseDE2 <- append_strReport(obj = objectImpulseDE2,
                                                 s = strMessage)
            tm_fitSigmoid <- system.time({
                objectImpulseDE2 <- fitSigmoidModels(
                    objectImpulseDE2 = objectImpulseDE2, 
                    vecConfounders = vecConfounders, strCondition = "case")
            })
            strMessage <- paste0("Consumed time: ", 
                                 round(tm_fitSigmoid["elapsed"]/60, 2), " min.")
            if (boolVerbose) { message(strMessage) }
            objectImpulseDE2 <- append_strReport(obj = objectImpulseDE2,
                                                 s = strMessage)
        }
        
        # 7. Differentially expression analysis based on model fits
        strMessage <- "# Differentially expression analysis based on model fits"
        if (boolVerbose) { message(strMessage) }
        objectImpulseDE2 <- append_strReport(obj = objectImpulseDE2,
                                             s = strMessage)
        objectImpulseDE2 <- runDEAnalysis(
            objectImpulseDE2 = objectImpulseDE2, 
            boolCaseCtrl = get_boolCaseCtrl(obj=objectImpulseDE2),
            boolIdentifyTransients = boolIdentifyTransients)
        
        if (!is.null(scaQThres)) {
            vecDEGenes <- as.vector(objectImpulseDE2$dfImpulseDE2Results[
                as.numeric(objectImpulseDE2$dfImpulseDE2Results$padj) <= 
                    scaQThres, "Gene"])
            strMessage <- paste0("Found ", length(vecDEGenes), " DE genes", 
                                 " at a FDR corrected p-value cut off of ", 
                                 scaQThres, ".")
            if (boolVerbose) { message(strMessage) }
            objectImpulseDE2 <- append_strReport(obj = objectImpulseDE2,
                                                 s = strMessage)
        } else {
            vecDEGenes <- NULL
        }
        objectImpulseDE2 <- set_vecDEGenes(obj=objectImpulseDE2,
                                           element=vecDEGenes)
    })
    strMessage <- "Finished running ImpulseDE2."
    if (boolVerbose) { message(strMessage) }
    objectImpulseDE2 <- append_strReport(obj = objectImpulseDE2,
                                         s = strMessage)
    strMessage <- paste0("TOTAL consumed time: ", 
                         round(tm_runImpulseDE2["elapsed"]/60, 2), " min.")
    if (boolVerbose) { message(strMessage) }
    objectImpulseDE2 <- append_strReport(obj = objectImpulseDE2,
                                         s = strMessage)
    
    return(objectImpulseDE2)
}
