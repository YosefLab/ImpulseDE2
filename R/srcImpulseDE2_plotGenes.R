#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Plot impulse fits    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Plots the impulse fits and data
#' 
#' Plots the impulse fits and data to pdf and return a list of gplots.
#' 
#' @seealso Called by separately by user.
#' 
#' @param vecGeneIDs: (string vector) Gene names to be plotted,
#'     Must be in rownames of matCountDataProc.
#' @param objectImpulseDE2: (ImpulseDE2 object)
#'    Object previously fitted to be used for plotting.
#' @param boolCaseCtrl: (bool) Whether to create case-ctrl plot.
#' @param dirOut: (dir) Directory into which pdf is printed.
#' @param strFileName: (str) [Default "ImpulseDE2_Trajectories.pdf"]
#'    File name of pdf with plots.
#' @param vecRefPval: (vector length vecGeneIDs) [Default NULL]
#'    P/Q-values to be displayed alongside ImpulseDE2 q-value
#'    for differential expression in plot titles.
#' @param strNameRefMethod: (str) [Default NULL]
#'    Name of reference method used to generate vecRefPval.
#'    Mentioned in plot titles.
#' 
#' @return lsgplotsID: (gplot list length vecGeneIDs)
#'    List of gplots for IDs in vecGeneIDs. This is secondary 
#'    output next to the .pdf and can be used to extract 
#'    single plots or assemble plots differently.
#' 
#' @author David Sebastian Fischer
#' 
#' @export
plotGenes <- function(vecGeneIDs, 
                      objectImpulseDE2,
                      boolCaseCtrl,
                      dirOut,
                      strFileName="ImpulseDE2_Trajectories.pdf",
                      vecRefPval=NULL, 
                      strNameRefMethod=NULL){
  
  # Load objects from output class
  matCountDataProc <- objectImpulseDE2@matCountDataProc
  dfAnnotationProc <- objectImpulseDE2@dfAnnotationProc
  lsModelFits <- objectImpulseDE2@lsModelFits
  vecSizeFactors <- objectImpulseDE2@vecSizeFactors
  vecDispersions <- objectImpulseDE2@vecDispersions
  vecConfounders <- objectImpulseDE2@vecConfounders
  
  # Set graphical parameters
  scaNPlotsPerPage <- 6
  scaRefTextSize <- 14
  scaResScale <- 2
  scaLabelsScale <- 2
  
  # Check input
  if(is.null(dirOut) | !file.exists(dirOut)) stop("Output directory dirOut not available.")
  if(!is.null(vecRefPval)) if(names(vecRefPval) != vecGeneIDs) stop("Names of vecRefPval have to be IDs from vecGeneIDs.")
  
  scaNIDs <- length(vecGeneIDs)
  lsgplotsID <- list()
  for(id in vecGeneIDs){
    
    vecTimePointsFit <- seq(min(dfAnnotationProc$Time), max(dfAnnotationProc$Time), length.out=100)
    vecCaseImpulseParam <- objectImpulseDE2@lsModelFits$case[[id]]$lsImpulseFit$vecImpulseParam
    vecCaseImpulseValue <- evalImpulse_comp(vecImpulseParam=vecCaseImpulseParam,
                                            vecTimepoints=vecTimePointsFit)
    
    # Can only resolve one level of batch effects as this is plotted
    # as different function fits to the indididual batches
    if(!is.null(vecConfounders)){
      vecCaseBatchFactors <- objectImpulseDE2@lsModelFits$case[[id]]$lsvecBatchFactors[[1]]
      vecBatchLabelsCase <- objectImpulseDE2@lsModelFits$IdxGroups$case$lsvecBatchUnique[[1]]
    } else {
      vecCaseBatchFactors <- 1
      vecBatchLabelsCase <- " "
    }
    
    vecValueToPlotCase <- do.call(c, lapply(vecCaseBatchFactors, function(f) vecCaseImpulseValue*f) )
    
    if(boolCaseCtrl){
      vecSamples <- dfAnnotationProc$Sample
      dfRaw <- data.frame(
        normCounts=matCountDataProc[id,vecSamples]/vecSizeFactors[vecSamples],
        time=dfAnnotationProc[vecSamples,]$Time,
        Batch=dfAnnotationProc[vecSamples,vecConfounders[1]],
        Condition=dfAnnotationProc[vecSamples,]$Condition
      )
      
      vecControlImpulseParam <- objectImpulseDE2@lsModelFits$control[[id]]$lsImpulseFit$vecImpulseParam
      vecControlImpulseValue <- evalImpulse_comp(vecImpulseParam=vecControlImpulseParam,
                                                 vecTimepoints=vecTimePointsFit)
      vecCombinedImpulseParam <- objectImpulseDE2@lsModelFits$combined[[id]]$lsImpulseFit$vecImpulseParam
      vecCombinedImpulseValue <- evalImpulse_comp(vecImpulseParam=vecCombinedImpulseParam,
                                                  vecTimepoints=vecTimePointsFit)
      
      vecControlBatchFactors <- objectImpulseDE2@lsModelFits$ctrl[[id]]$lsvecBatchFactors[[1]]
      vecCombinedBatchFactors <- objectImpulseDE2@lsModelFits$comb[[id]]$lsvecBatchFactors[[1]]
      
      vecBatchLabelsCtrl <- objectImpulseDE2@lsModelFits$IdxGroups$ctrl$lsvecBatchUnique[[1]]
      vecBatchLabelsComb <- objectImpulseDE2@lsModelFits$IdxGroups$comb$lsvecBatchUnique[[1]]
      
      vecValueToPlotCtrl <- do.call(c, lapply(vecCtrlBatchFactors, function(f) vecControlImpulseValue*f) )
      vecValueToPlotComb <- do.call(c, lapply(vecCombBatchFactors, function(f) vecCombinedImpulseValue*f) )
      
      dfFit <- data.frame(
        time=c(rep(vecTimePointsFit, length(vecCaseBatchFactors)),
               rep(vecTimePointsFit, length(vecCtrlBatchFactors)),
               rep(vecTimePointsFit, length(vecCombBatchFactors))),
        value=c(vecValueToPlotCase,
                vecValueToPlotCtrl,
                vecValueToPlotComb),
        Condition=c(rep(rep("case", length(vecTimePointsFit)), length(vecCaseBatchFactors)),
                    rep(rep("ctrl", length(vecTimePointsFit)), length(vecCtrlBatchFactors)),
                    rep(rep("comb", length(vecTimePointsFit)), length(vecCombBatchFactors))),
        Fit=c(sapply(vecBatchLabelsCase, function(label) rep(label, length(vecTimePointsFit)) ),
                sapply(vecBatchLabelsCtrl, function(label) rep(label, length(vecTimePointsFit)) ),
                sapply(vecBatchLabelsComb, function(label) rep(label, length(vecTimePointsFit)) ) )
      )
      gplotGene <- ggplot() +
        geom_point(data=dfRaw, aes(x=time, y=normCounts, colour=Condition, shape=Batch), size=2) +
        geom_line(data=dfFit, aes(x=time, y=value, colour=Condition, linetype=Fit), size=1)
    } else {
      vecSamples <- dfAnnotationProc[dfAnnotationProc$Condition=="case",]$Sample
      if(!is.null(vecConfounders)){ vecBatches <- dfAnnotationProc[vecSamples,vecConfounders[1]]
      } else { vecBatches <- rep(" ", length(vecSamples)) }
      
      dfRaw <- data.frame(
        normCounts=matCountDataProc[id,vecSamples]/vecSizeFactors[vecSamples],
        time=dfAnnotationProc[vecSamples,]$Time,
        Batch=vecBatches,
        condition=dfAnnotationProc[vecSamples,]$Condition
      )
      
      dfFit <- data.frame(
        time=rep(vecTimePointsFit, length(vecCaseBatchFactors)),
        value=vecValueToPlotCase,
        condition=rep(rep("case", length(vecTimePointsFit)), length(vecCaseBatchFactors)),
        Batch=as.vector(sapply(vecBatchLabelsCase, function(label) rep(label, length(vecTimePointsFit)) ))
      )
      gplotGene <- ggplot() +
        geom_point(data=dfRaw, aes(x=time, y=normCounts, colour=Batch, shape=Batch), size=2) +
        geom_line(data=dfFit, aes(x=time, y=value, colour=Batch, linetype=Batch), size=1)
    }
    gplotID <- gplotGene + 
      xlab("time [hours]") +
      ylab("Read counts")
    if(!is.null(strNameRefMethod)){
      gplotID <- gplotID +
        labs(title=paste0(id, ": ImpulseDE2 [",
                          round(log(objectImpulseDE2$dfImpulseDE2Results[id,]$padj)/log(10), 2), 
                          "] ", strNameRefMethod," [", 
                          round(log(vecRefPval[id])/log(10), 2), "]" ),
             colour="Batch", linetype="Batch", shape="Batch")
    } else {
      gplotID <- gplotID +
        labs(title=paste0(id, ": ImpulseDE2 [",
                          round(log(objectImpulseDE2$dfImpulseDE2Results[id,]$padj)/log(10), 2), 
                          "] "))
    }
    lsgplotsID[[length(lsgplotsID)+1]] <- gplotID
  }
  
  # Printing to file
  dirFileOut <- paste0(dirOut, strFileName)
  print(paste0("Creating ", dirFileOut))
  graphics.off()
  pdf(dirFileOut)
  for(p in seq(0,scaNIDs %/% scaNPlotsPerPage)){
    if(p < scaNIDs %/% scaNPlotsPerPage){ vecidxPlots <- seq((p*scaNPlotsPerPage+1),((p+1)*(scaNPlotsPerPage)))
    } else { vecidxPlots <- seq((p*scaNPlotsPerPage+1),scaNIDs) }
    print(plot_grid(plotlist=lsgplotsID[vecidxPlots],
                    align="h",
                    nrow=3, ncol=2,
                    rel_widths=c(1,1),
                    rel_heights=c(1,1,1)))
  }
  dev.off()
  graphics.off()
  
  # Return raw list of gplots
  return(lsgplotsID)
}