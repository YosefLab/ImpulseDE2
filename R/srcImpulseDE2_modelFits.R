### Return model fits

#' Return a gene x sample and gene x time point (for each condition) matrix with model fits.
#' 
#' @seealso Called seperately by user.
#' 
#' @param objectImpulseDE2 (instance of class ImpulseDE2Object)
#' ImpulseDE2 output object to create heatmap from.
#' 
#' @return 
#' \itemize{
#' \item sample (data.frame)
#' Data frame with model fit by gene and sample.
#' \item case (data.frame)
#' Data frame with model fit by gene and time point in case condition
#' \item control (data.frame)
#' Data frame with model fit by gene and time point in control condition
#' }
#' 
#' @examples
#' lsSimulatedData <- simulateDataSetImpulseDE2(
#' vecTimePointsA   = rep(seq(1,8),3),
#' vecTimePointsB   = NULL,
#' vecBatchesA      = do.call(c, lapply(c("A", "B", "C"), function(x) rep(x, 8))),
#' vecBatchesB      = NULL,
#' scaNConst        = 0,
#' scaNImp          = 50,
#' scaNLin          = 0,
#' scaNSig          = 50,
#' scaMuBatchEffect = 1, 
#' scaSDBatchEffect = 1)
#' objectImpulseDE2 <- runImpulseDE2(
#' matCountData    = lsSimulatedData$matObservedCounts, 
#' dfAnnotation    = lsSimulatedData$dfAnnotation,
#' boolCaseCtrl    = FALSE,
#' vecConfounders  = c("Batch"),
#' boolIdentifyTransients = TRUE,
#' scaNProc        = 1 )
#' modelFits <- computeModelFits(objectImpulseDE2=objectImpulseDE2)
#' head(modelFits$case)
#' 
#' @author David Sebastian Fischer
#' 
#' @export
computeModelFits <- function(
    objectImpulseDE2
    ) {
    
    dfAnnot <- get_dfAnnotationProc(obj=objectImpulseDE2)
    
    # Compute by time point and by sample:
    vecTimePointsToEval <- dfAnnot$Time
    matImpulseValueSample <- do.call(rbind, lapply(objectImpulseDE2@vecAllIDs, function(x) {
        vecSfValues <- get_vecSizeFactors(obj=objectImpulseDE2)
        vecImpulseValuesCase <- ImpulseDE2:::evalImpulse_comp(
            vecImpulseParam = 
                get_lsModelFits(obj=objectImpulseDE2)$case[[x]]$
                lsImpulseFit$vecImpulseParam, 
            vecTimepoints = dfAnnot$Time
        )
        vecModelValues <- vecImpulseValuesCase * vecSfValues
        if (!is.null(objectImpulseDE2@lsModelFits$IdxGroups$case$lsvecidxBatch)) {
            vecBatchFactorsCase <- get_lsModelFits(
                obj=objectImpulseDE2)$case[[x]]$lsImpulseFit$lsvecBatchFactors
            for (i in seq(length(vecBatchFactorsCase))) {
                vecModelValues <- vecModelValues * vecBatchFactorsCase[[i]][
                    objectImpulseDE2@lsModelFits$IdxGroups$case$lsvecidxBatch[[i]]
                ]
            }
        }
        
        if (!is.null(objectImpulseDE2@lsModelFits$IdxGroups$control)) {
            vecImpulseValuesCtrl <- ImpulseDE2:::evalImpulse_comp(
                vecImpulseParam = 
                    get_lsModelFits(obj=objectImpulseDE2)$control[[x]]$
                    lsImpulseFit$vecImpulseParam, 
                vecTimepoints = dfAnnot$Time
            )
            vecModelValuesCtrl <- vecImpulseValuesCtrl * vecSfValues
            if (!is.null(objectImpulseDE2@lsModelFits$IdxGroups$case$lsvecidxBatch)) {
                vecBatchFactorsCtrl <- get_lsModelFits(
                    obj=objectImpulseDE2)$control[[x]]$lsImpulseFit$lsvecBatchFactors
                for (i in seq(length(vecBatchFactorsCase))) {
                    vecModelValuesCtrl <- vecModelValuesCtrl * vecBatchFactorsCtrl[[i]][
                        objectImpulseDE2@lsModelFits$IdxGroups$control$lsvecidxBatch[[i]]
                    ]
                }
            }
            vecModelValues[dfAnnot$Condition == "control"] <- 
                vecModelValuesCtrl[dfAnnot$Condition == "control"]
        }
        return(vecModelValues)
    }))
    rownames(matImpulseValueSample) <- dfAnnot$vecAllIDs
    colnames(matImpulseValueSample) <- dfAnnot$Sample

    # Compute by time point in case condition:
    matImpulseValueCase <- do.call(rbind, lapply(objectImpulseDE2@vecAllIDs, function(x) {
        evalImpulse_comp(
            vecImpulseParam = 
                get_lsModelFits(obj=objectImpulseDE2)$case[[x]]$
                lsImpulseFit$vecImpulseParam, 
            vecTimepoints = sort(unique(dfAnnot$Time), decreasing = FALSE)
        )
    }))
    rownames(matImpulseValueCase) <- objectImpulseDE2@vecAllIDs
    colnames(matImpulseValueCase) <- sort(unique(dfAnnot$Time), decreasing = FALSE)
    
    # Compute by time point in control condition:
    if (!is.null(objectImpulseDE2@lsModelFits$IdxGroups$control)) {
        matImpulseValueCtrl <- do.call(rbind, lapply(objectImpulseDE2@vecAllIDs, function(x) {
            evalImpulse_comp(
                vecImpulseParam = 
                    get_lsModelFits(obj=objectImpulseDE2)$control[[x]]$
                    lsImpulseFit$vecImpulseParam, 
                vecTimepoints = sort(unique(dfAnnot$Time), decreasing = FALSE)
            )
        }))
        rownames(matImpulseValueCtrl) <- objectImpulseDE2@vecAllIDs
        colnames(matImpulseValueCtrl) <- sort(unique(dfAnnot$Time), decreasing = FALSE)
    } else {
        matImpulseValueCtrl <- NULL
    }
    
    return(list(sample = matImpulseValueSample, 
                case = matImpulseValueCase, 
                control = matImpulseValueCtrl))
}
