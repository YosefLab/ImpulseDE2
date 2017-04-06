### Plot impulse fits

#' @include srcImpulseDE2_evalImpulse.R
#' @include srcImpulseDE2_classImpulseDE2Object.R
NULL

#' Plots the impulse fits and data
#' 
#' Plots the impulse fits and data to pdf and return a list of gplots. Points
#' are size factor normalised data. Consider using boolSimplePlot=TRUE
#' if the plot seems to crowded.
#' 
#' @seealso Called by separately by user.
#' 
#' @param vecGeneIDs (string vector) [Default NULL]
#' Gene names to be plotted. Must be in rownames of 
#' objectImpulseDE2@matCountDataProc.
#' Supply either vecGeneIDs or scaNTopIDs.
#' @param scaNTopIDs (int) [Default NULL]
#' Number of top differentially expressed (by q-value) genes to 
#' be plotted
#' Supply either vecGeneIDs or scaNTopIDs.
#' @param objectImpulseDE2 (ImpulseDE2 object)
#' Object previously fitted to be used for plotting.
#' @param boolCaseCtrl (bool) Whether to create case-ctrl plot.
#' @param dirOut (dir) [Default NULL]
#' Directory into which pdf is printed.
#' @param strFileName (str) [Default 'ImpulseDE2_Trajectories.pdf']
#' File name of pdf with plots.
#' @param boolMultiplePlotsPerPage (bool) [Default TRUE]
#' Whether to create grid with multiple plots on each page of pdf.
#' @param boolSimplePlot (bool) [Default TRUE]
#' Whether to omit batch structure in plotting of model fits
#' and only plot fit to first batch/all data (if no confounders were given).
#' This strongly simplifies plots and is recommended e.g. for case-ctrl data.
#' @param vecRefPval (vector length vecGeneIDs) [Default NULL]
#' P/Q-values to be displayed alongside ImpulseDE2 q-value
#' for differential expression in plot titles.
#' @param strNameRefMethod (str) [Default NULL]
#' Name of reference method used to generate vecRefPval.
#' Mentioned in plot titles.
#' 
#' @return lsgplotsID (gplot list length vecGeneIDs)
#' List of gplots for IDs in vecGeneIDs. This is secondary 
#' output next to the .pdf and can be used to extract 
#' single plots or assemble plots differently.
#' 
#' @examples
#' lsSimulatedData <- simulateDataSetImpulseDE2(
#' vecTimePointsA   = rep(seq(1,8),3),
#' vecTimePointsB   = NULL,
#' vecBatchesA      = NULL,
#' vecBatchesB      = NULL,
#' scaNConst        = 0,
#' scaNImp          = 40,
#' scaNLin          = 20,
#' scaNSig          = 40)
#' objectImpulseDE2 <- runImpulseDE2(
#' matCountData    = lsSimulatedData$matObservedCounts, 
#' dfAnnotation    = lsSimulatedData$dfAnnotation,
#' boolCaseCtrl    = FALSE,
#' vecConfounders  = NULL,
#' boolIdentifyTransients = FALSE,
#' scaNProc        = 1 )
#' lsgplotsID <- plotGenes(
#' scaNTopIDs=5,
#' objectImpulseDE2=objectImpulseDE2,
#' boolCaseCtrl=FALSE,
#' boolMultiplePlotsPerPage=TRUE,
#' boolSimplePlot=FALSE)
#' lsgplotsID[[1]]
#' 
#' @author David Sebastian Fischer
#' 
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' 
#' @export
plotGenes <- function(
    vecGeneIDs = NULL, scaNTopIDs = NULL, objectImpulseDE2, 
    boolCaseCtrl, dirOut = NULL, strFileName = "ImpulseDE2_Trajectories.pdf", 
    boolMultiplePlotsPerPage = TRUE, boolSimplePlot = FALSE, vecRefPval = NULL, 
    strNameRefMethod = NULL) {
    
    dfAnnot <- get_dfAnnotationProc(obj=objectImpulseDE2)
    # Set graphical parameters
    scaNPlotsPerPage <- 4
    # Colour-blind palette
    cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                   "#0072B2", "#D55E00", "#CC79A7")
    
    # Check input
    if (is.null(get_lsModelFits(objectImpulseDE2))) {
        error(paste0("objectImpulseDE2 does not contain model fits. ",
                     "Run ImpulseDE2_main first."))
    }
    if (is.null(vecGeneIDs) & is.null(scaNTopIDs)) {
        stop("Supply either vecGeneIDs or scaNTopIDs.")
    }
    if (!is.null(vecGeneIDs) & !is.null(scaNTopIDs)) { 
        stop("Only one of the two: vecGeneIDs or scaNTopIDs.")
    }
    if (is.null(get_lsModelFits(objectImpulseDE2)$
                IdxGroups$case$lsvecBatchUnique) & 
        !is.null(boolSimplePlot)) {
        print("Setting boolSimplePlot=TRUE as no batch structure was found.")
        boolSimplePlot <- TRUE
    }
    if (!is.null(dirOut) && !file.exists(dirOut)) {
        stop("Output directory dirOut not available.")
    }
    if (!is.null(vecRefPval) && (names(vecRefPval) != vecGeneIDs)) {
        stop("Names of vecRefPval have to be IDs from vecGeneIDs.")
    }
    if (!boolSimplePlot & boolCaseCtrl & boolMultiplePlotsPerPage) {
        warning(paste0("Plots are likely overloaded. ",
                       "Consider switching to boolSimplePlot=TRUE ",
                       "or boolMultiplePlotsPerPage=FALSE."))
    }
    # Use top scaNTopIDs IDs if no specific IDs were supplied.
    if (is.null(vecGeneIDs)) 
        vecGeneIDs <- objectImpulseDE2$dfImpulseDE2Results[
            with(objectImpulseDE2$dfImpulseDE2Results, 
                 order(padj)), ]$Gene[1:scaNTopIDs]
    scaNIDs <- length(vecGeneIDs)
    lsgplotsID <- list()
    for (id in vecGeneIDs) {
        vecTimePointsFit <- seq(
            min(get_dfAnnotationProc(obj=objectImpulseDE2)$Time), 
                                max(get_dfAnnotationProc(obj=objectImpulseDE2)$Time), 
                                length.out = 100)
        vecCaseImpulseParam <- get_lsModelFits(
            obj=objectImpulseDE2)$case[[id]]$lsImpulseFit$vecImpulseParam
        vecCaseImpulseValue <- evalImpulse_comp(
            vecImpulseParam = vecCaseImpulseParam, 
            vecTimepoints = vecTimePointsFit)
        
        # Can only resolve one level of batch effects as this is plotted as
        # different function fits to the indididual batches
        if (!is.null(get_vecConfounders(obj=objectImpulseDE2))) {
            vecCaseBatchFactors <- get_lsModelFits(
                obj=objectImpulseDE2)$case[[id]]$lsImpulseFit$lsvecBatchFactors[[1]]
            vecBatchLabelsCase <- get_lsModelFits(
                obj=objectImpulseDE2)$IdxGroups$case$lsvecBatchUnique[[1]]
        } else {
            vecCaseBatchFactors <- 1
            vecBatchLabelsCase <- " "
        }
        
        vecValueToPlotCase <- do.call(c, lapply(
            vecCaseBatchFactors, function(f) vecCaseImpulseValue * f))
        
        if (boolCaseCtrl) {
            vecSamples <- dfAnnot$Sample
            dfRaw <- data.frame(
                normCounts = get_matCountDataProc(
                    obj=objectImpulseDE2)[id,vecSamples]/
                    get_vecSizeFactors(objectImpulseDE2)[vecSamples], 
                time = dfAnnot[vecSamples, ]$Time, 
                Batch = dfAnnot[
                    vecSamples, get_vecConfounders(obj=objectImpulseDE2)[1]], 
                Condition = dfAnnot[vecSamples, ]$Condition)
            
            vecControlImpulseParam <- get_lsModelFits(
                obj=objectImpulseDE2)$control[[id]]$lsImpulseFit$vecImpulseParam
            vecControlImpulseValue <- evalImpulse_comp(
                vecImpulseParam = vecControlImpulseParam, 
                vecTimepoints = vecTimePointsFit)
            vecCombinedImpulseParam <-get_lsModelFits(
                obj=objectImpulseDE2)$combined[[id]]$lsImpulseFit$vecImpulseParam
            vecCombinedImpulseValue <- evalImpulse_comp(
                vecImpulseParam = vecCombinedImpulseParam, 
                vecTimepoints = vecTimePointsFit)
            
            vecControlBatchFactors <- get_lsModelFits(
                obj=objectImpulseDE2)$control[[id]]$lsImpulseFit$lsvecBatchFactors[[1]]
            vecCombinedBatchFactors <- get_lsModelFits(
                obj=objectImpulseDE2)$combined[[id]]$lsImpulseFit$lsvecBatchFactors[[1]]
            
            vecBatchLabelsCtrl <- get_lsModelFits(
                obj=objectImpulseDE2)$IdxGroups$control$lsvecBatchUnique[[1]]
            vecBatchLabelsComb <- get_lsModelFits(
                obj=objectImpulseDE2)$IdxGroups$combined$lsvecBatchUnique[[1]]
            
            vecValueToPlotCtrl <- do.call(c, lapply(
                vecControlBatchFactors, function(f) vecControlImpulseValue * f))
            vecValueToPlotComb <- do.call(c, lapply(
                vecCombinedBatchFactors, function(f) vecCombinedImpulseValue * f))
            
            if (boolSimplePlot) {
                dfFit <- data.frame(
                    time = rep(vecTimePointsFit, 3), 
                    value = c(vecCaseImpulseValue, 
                              vecControlImpulseValue, 
                              vecCombinedImpulseValue), 
                    Condition = c(rep("case", length(vecTimePointsFit)), 
                                  rep("control", length(vecTimePointsFit)), 
                                  rep("combined", length(vecTimePointsFit))))
                if (!is.null(vecBatchLabelsComb)) {
                    gplotGene <- ggplot() + 
                        geom_point(data = dfRaw, 
                                   aes(x = time, y = normCounts, 
                                                     colour = Condition, 
                                                     shape = Batch)) + 
                        geom_line(data = dfFit, 
                                  aes(x = time, y = value, 
                                                    colour = Condition))
                } else {
                    gplotGene <- ggplot() + 
                        geom_point(data = dfRaw, 
                                   aes(x = time, y = normCounts, 
                                                     colour = Condition)) + 
                        geom_line(data = dfFit, 
                                  aes(x = time, y = value, 
                                                    colour = Condition))
                }
            } else {
                dfFit <- data.frame(
                    time=c(rep(vecTimePointsFit, 
                               length(vecCaseBatchFactors)),
                           rep(vecTimePointsFit, 
                               length(vecControlBatchFactors)),
                           rep(vecTimePointsFit, 
                               length(vecCombinedBatchFactors))),
                    value=c(vecValueToPlotCase,
                            vecValueToPlotCtrl,
                            vecValueToPlotComb),
                    Condition=c(rep(rep("case", length(vecTimePointsFit)), 
                                    length(vecCaseBatchFactors)),
                                rep(rep("control", length(vecTimePointsFit)), 
                                    length(vecControlBatchFactors)),
                                rep(rep("combined", length(vecTimePointsFit)), 
                                    length(vecCombinedBatchFactors))),
                    BatchFit=c(sapply(vecBatchLabelsCase, function(label) 
                        rep(label, length(vecTimePointsFit)) ),
                        sapply(vecBatchLabelsCtrl, function(label) 
                            rep(label, length(vecTimePointsFit)) ),
                        sapply(vecBatchLabelsComb, function(label) 
                            rep(label, length(vecTimePointsFit)) ) )
                )
                gplotGene <- ggplot() + 
                    geom_point(data = dfRaw, aes(x = time, y = normCounts, 
                                                 colour = Condition, shape = Batch)) + 
                    geom_line(data = dfFit, aes(x = time, y = value, 
                                                colour = Condition, 
                                                linetype = BatchFit))
            }
        } else {
            vecSamples <- dfAnnot[
                dfAnnot$Condition == "case", ]$Sample
            if (!is.null(get_vecConfounders(obj=objectImpulseDE2))) {
                vecBatches <- dfAnnot[
                    vecSamples, get_vecConfounders(obj=objectImpulseDE2)[1]]
            } else {
                vecBatches <- rep(" ", length(vecSamples))
            }
            
            dfRaw <- data.frame(
                normCounts = get_matCountDataProc(
                    obj=objectImpulseDE2)[id, vecSamples]/
                    get_vecSizeFactors(obj=objectImpulseDE2)[vecSamples], 
                time = dfAnnot[vecSamples, ]$Time, 
                Batch = vecBatches, 
                condition = dfAnnot[vecSamples, ]$Condition)
            
            if (boolSimplePlot) {
                dfFit <- data.frame(
                    time = vecTimePointsFit, value = vecCaseImpulseValue, 
                    Batch = rep(vecBatchLabelsCase[1], length(vecTimePointsFit)))
                if (!is.null(vecBatchLabelsCase)) {
                    gplotGene <- ggplot() + 
                        geom_point(data = dfRaw, aes(x = time,  y = normCounts, 
                                                     shape = Batch)) + 
                        geom_line(data = dfFit, 
                                  aes(x = time, y = value, linetype = Batch))
                } else {
                    gplotGene <- ggplot() + 
                        geom_point(data = dfRaw, aes(x = time, y = normCounts)) + 
                        geom_line(data = dfFit, aes(x = time, 
                                                    y = value))
                }
            } else {
                dfFit <- data.frame(
                    time = rep(vecTimePointsFit, length(vecCaseBatchFactors)), 
                    value = vecValueToPlotCase, 
                    condition = rep(rep("case", length(vecTimePointsFit)), 
                                    length(vecCaseBatchFactors)), 
                    BatchFit = as.vector(sapply(
                        vecBatchLabelsCase, function(label) 
                            rep(label, length(vecTimePointsFit)))))
                gplotGene <- ggplot() + 
                    geom_point(data = dfRaw, 
                               aes(x = time, y = normCounts, 
                                                 colour = Batch, 
                                                 shape = Batch)) + 
                    geom_line(data = dfFit, 
                              aes(x = time, y = value, 
                                  colour = BatchFit, 
                                  linetype = BatchFit))
            }
        }
        gplotID <- gplotGene + scale_colour_manual(values = cbPalette) + 
            xlab("time [hours]") + ylab("Read counts")
        if (!is.null(strNameRefMethod)) {
            gplotID <- gplotID + 
                labs(title = paste0(
                    id, ": ImpulseDE2 [", 
                    round(log(objectImpulseDE2$dfImpulseDE2Results[id, ]$padj)/
                              log(10), 2), 
                    "] ", strNameRefMethod, " [", 
                    round(log(vecRefPval[id])/log(10), 2), 
                    "]"))
        } else {
            gplotID <- gplotID + labs(title = paste0(
                id, ": log10 q = ", 
                round(log(objectImpulseDE2$dfImpulseDE2Results[id, ]$padj)/
                          log(10), 2)))
        }
        lsgplotsID[[length(lsgplotsID) + 1]] <- gplotID
    }
    
    # Print to file
    if (!is.null(dirOut)) {
        dirFileOut <- paste0(dirOut, strFileName)
        print(paste0("Creating ", dirFileOut))
        graphics.off()
        if (boolMultiplePlotsPerPage) {
            pdf(dirFileOut)
            scaNPages <- scaNIDs%/%scaNPlotsPerPage
            if (scaNIDs%%scaNPlotsPerPage == 0) 
                scaNPages <- scaNPages - 1
            for (p in seq(0, scaNPages)) {
                if (p < scaNIDs%/%scaNPlotsPerPage) {
                    vecidxPlots <- seq((p * scaNPlotsPerPage + 1), 
                                       ((p + 1) * (scaNPlotsPerPage)))
                } else {
                    vecidxPlots <- seq((p * scaNPlotsPerPage + 1), scaNIDs)
                }
                print(plot_grid(plotlist = lsgplotsID[vecidxPlots], align = "h", 
                                nrow = scaNPlotsPerPage/2, ncol = 2, 
                                rel_widths = c(1, 1), rel_heights = c(1, 1, 1)))
            }
        } else {
            pdf(dirFileOut)
            for (p in seq(1, scaNIDs)) {
                print(lsgplotsID[[p]])
            }
        }
        dev.off()
        graphics.off()
    }
    
    # Return raw list of gplots
    return(lsgplotsID)
}
