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
  if(!is.null(vecRefPval) & names(vecRefPval) != vecGeneIDs) stop("Names of vecRefPval have to be IDs from vecGeneIDs.")
  
  scaNIDs <- length(vecGeneIDs)
  lsgplotsID <- list()
  for(id in vecGeneIDs){
    
    vecTimePointsFit <- seq(min(dfAnnotation$Time), max(dfAnnotation$Time), length.out=100)
    vecCaseImpulseParam <- lsResImpulseDE2$lsModelFits$case[[id]]$lsImpulseFit$vecImpulseParam
    vecCaseImpulseValue <- evalImpulse_comp(vecImpulseParam=vecCaseImpulseParam,
                                            vecTimepoints=vecTimePointsFit)
    
    # Can only resolve one level of batch effects as this is plotted
    # as different function fits to the indididual batches
    vecCaseBatchFactors <- lsResImpulseDE2$lsModelFits$case[[id]]$lsvecBatchFactors[[1]]
    vecBatchLabelsCase <- lsResImpulseDE2$lsModelFits$IdxGroups$case$lsvecBatchUnique[[1]]
    
    vecValueToPlotCase <- do.call(c, lapply(vecCaseBatchFactors, function(f) vecCaseImpulseValue*f) )
    
    if(boolCaseCtrl){
      vecSamples <- dfAnnotationProc$Gene
      dfRaw <- data.frame(
        normCounts=matCountDataProc[id,vecSamples]/vecSizeFactors[vecSamples],
        time=dfAnnotationProc[vecSamples,]$Time,
        batch=dfAnnotation[vecSamples,vecConfounders[1]],
        condition=dfAnnotation[vecSamples,]$Condition
      )
      
      vecControlImpulseParam <- lsResImpulseDE2$lsModelFits$control[[id]]$lsImpulseFit$vecImpulseParam
      vecControlImpulseValue <- evalImpulse_comp(vecImpulseParam=vecControlImpulseParam,
                                                 vecTimepoints=vecTimePointsFit)
      vecCombinedImpulseParam <- lsResImpulseDE2$lsModelFits$combined[[id]]$lsImpulseFit$vecImpulseParam
      vecCombinedImpulseValue <- evalImpulse_comp(vecImpulseParam=vecCombinedImpulseParam,
                                                  vecTimepoints=vecTimePointsFit)
      
      vecControlBatchFactors <- lsResImpulseDE2$lsModelFits$ctrl[[id]]$lsvecBatchFactors[[1]]
      vecCombinedBatchFactors <- lsResImpulseDE2$lsModelFits$comb[[id]]$lsvecBatchFactors[[1]]
      
      vecBatchLabelsCtrl <- lsResImpulseDE2$lsModelFits$IdxGroups$ctrl$lsvecBatchUnique[[1]]
      vecBatchLabelsComb <- lsResImpulseDE2$lsModelFits$IdxGroups$comb$lsvecBatchUnique[[1]]
      
      vecValueToPlotCtrl <- do.call(c, lapply(vecCtrlBatchFactors, function(f) vecControlImpulseValue*f) )
      vecValueToPlotComb <- do.call(c, lapply(vecCombBatchFactors, function(f) vecCombinedImpulseValue*f) )
      
      dfFit <- data.frame(
        time=c(rep(vecTimePointsFit, length(vecCaseBatchFactors)),
               rep(vecTimePointsFit, length(vecCtrlBatchFactors)),
               rep(vecTimePointsFit, length(vecCombBatchFactors))),
        value=c(vecValueToPlotCase,
                vecValueToPlotCtrl,
                vecValueToPlotComb),
        condition=c(rep(rep("case", length(vecTimePointsFit)), length(vecCaseBatchFactors)),
                    rep(rep("ctrl", length(vecTimePointsFit)), length(vecCtrlBatchFactors)),
                    rep(rep("comb", length(vecTimePointsFit)), length(vecCombBatchFactors))),
        batchfit=c(sapply(vecBatchLabelsCase, function(label) rep(label, length(vecTimePointsFit)) ),
                sapply(vecBatchLabelsCtrl, function(label) rep(label, length(vecTimePointsFit)) ),
                sapply(vecBatchLabelsComb, function(label) rep(label, length(vecTimePointsFit)) ) )
      )
      gplotGene <- ggplot() +
        geom_point(data=dfRaw, aes(x=time, y=normCounts, colour=condition, shape=batch), size=6) +
        geom_line(data=dfFit, aes(x=time, y=value, colour=condition, linetype=batchfit), size=2)
    } else {
      vecSamples <- dfAnnotationProc[dfAnnotationProc$Condition=="case",]$Gene
      dfRaw <- data.frame(
        normCounts=matCountDataProc[id,vecSamples]/vecSizeFactors[vecSamples],
        time=dfAnnotationProc[vecSamples,]$Time,
        batch=dfAnnotation[vecSamples,vecConfounders[1]],
        condition=dfAnnotation[vecSamples,]$Condition
      )
      
      dfFit <- data.frame(
        time=rep(vecTimePointsFit, length(vecCaseBatchFactors)),
        value=vecValueToPlotCase,
        condition=rep(rep("case", length(vecTimePointsFit)), length(vecCaseBatchFactors)),
        batchfit=sapply(vecBatchLabelsCase, function(label) rep(label, length(vecTimePointsFit)) )
      )
      gplotGene <- ggplot() +
        geom_point(data=dfRaw, aes(x=time, y=normCounts, colour=batch, shape=batch), size=6) +
        geom_line(data=dfFit, aes(x=time, y=value, colour=batch, linetype=batchfit), size=2)
    }
    lsgplotsID[[length(lsgplotsID)+1]] <- gplotGene + 
      labs(title=paste0(id, ": ImpulseDE2 [",
                        round(log(matPval[id,]$ImpulseDE2)/log(10), 2), 
                        "] ", strRefMethod," [", 
                        round(log(matPval[id,strRefMethod])/log(10), 2), "]" )) +
      xlab("time [hours]") +
      ylab("Read counts")
  }
  
  ### Put pages together
  # Add margins
  marginObject <- theme(plot.margin = unit(c(1,1,1,2), "cm")) #t,r,b,l
  lsgplotsID <- lapply(lsgplotsID, "+", marginObject)
  
  # Add text and legend formats
  themeObjectStandard <- theme(axis.text=element_text(size=scaRefTextSize*scaResScale),
                               axis.title=element_text(size=scaRefTextSize*scaResScale,face="bold"),
                               plot.title=element_text(size=rel(scaResScale),face="bold"),
                               legend.text=element_text(size=scaRefTextSize*scaResScale),
                               legend.title=element_text(size=scaRefTextSize*scaResScale),
                               legend.key.size = unit(scaResScale, 'lines'))
  # themeObjectNoLegend is currently not used. Add in lsthemeObjects if displaying many plots
  themeObjectNoLegend <- theme(axis.text=element_text(size=scaRefTextSize*scaResScale),
                               axis.title=element_text(size=scaRefTextSize*scaResScale,face="bold"),
                               plot.title=element_text(size=rel(scaResScale),face="bold"),
                               legend.position="none")
  lsthemeObjects <- list(themeObjectStandard,themeObjectStandard,
                         themeObjectStandard,themeObjectStandard,
                         themeObjectStandard,themeObjectStandard )
  lsgplotsID <- lapply(seq(1, length(lsgplotsID)), function(g){
    lsgplotsID[[g]] + lsthemeObjects[[g]]
  })
  
  # Printing to file
  dirFileOut <- paste0(dirOut, strFileName)
  print(paste0("Creating ", dirFileOut))
  graphics.off()
  pdf(dirFileOut, width=480*2*2, height=480*2*3)
  for(p in seq(1,scaNIDs %/% scaNPlotsPerPage)){
    if(p < scaNIDs %/% scaNPlotsPerPage){ vecidxPlots <- seq((p*scaNPlotsPerPage+1),(p*(scaNPlotsPerPage+1)))
    } else { vecidxPlots <- seq((p*scaNPlotsPerPage+1),scaNIDs) }
    print(plot_grid(plotlist=lsgplotsID[[vecidxPlots]],
                    align="h",
                    nrow=3, ncol=2,
                    rel_widths=c(1,1),
                    rel_heights=c(1,1,1),
                    labels="AUTO",
                    label_size=scaRefTextSize*scaResScale*scaLabelsScale))
  }
  dev.off()
  graphics.off()
  
  # Return raw list of gplots
  return(lsgplotsID)
}