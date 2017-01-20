#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Plot z-value heatmaps    ++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Plot structured z-value heatmaps of differentially expressed genes
#' 
#' Creates a complexHeatmap heatmap structured into subsets of genes according
#' to their behaviour and sorted by peak time for raw counts and for the
#' fitted signal.
#' 
#' @seealso Called seperately by used.
#' 
#' @param objectImpulseDE2 (instance of class ImpulseDE2Object)
#'    ImpulseDE2 output object to create heatmap from.
#' @param strCondition (str) {"case","control","combined}
#'    Heatmap is created from samples of this condition.
#' @param boolIdentifyTransients (bool) 
#'    Whether to structure heatmap into transient and transition
#'    trajectories, only possible if sigmoids were fit to the
#'    indicated condition.
#' @param scaQThres (scalar) FDR-corrected p-value threshold
#'    for calling differentially expressed genes: Only genes
#'    below this threshold are included in the heatmap.
#'    
#' @return (list length 3)
#'    \itemize{
#'      \item complexHeatmapRaw (complexHeatmap plot)
#'      Heatmap of raw data by time point: Average of the
#'      size factor (and batch factor) normalised counts 
#'      per time point and gene.
#'      Plot with draw(complexHeatmapRaw).
#'      \item complexHeatmapFit (complexHeatmap plot)
#'      Heatmap of impulse-fitted data by time point.
#'      Plot with draw(complexHeatmapFit).
#'      \item lsvecGeneGroups (list)
#'      List of gene ID vectors: One per heatmap group 
#'      with all gene IDs of the the profiles displayed
#'      in the heatmap.
#'    }
#'    
#' @author David Sebastian Fischer
#' 
#' @import ComplexHeatmap
#' @import circlize
#' 
#' @export
plotHeatmap <- function(objectImpulseDE2,
                        strCondition,
                        boolIdentifyTransients,
                        scaQThres=0.01){
  
  # Load objects from output class
  matCountDataProc <- objectImpulseDE2@matCountDataProc
  dfAnnotationProc <- objectImpulseDE2@dfAnnotationProc
  lsModelFits <- objectImpulseDE2@lsModelFits
  dfImpulseDE2Results  <- objectImpulseDE2@dfImpulseDE2Results
  vecSizeFactors <- objectImpulseDE2@vecSizeFactors
  
  scaNGenes <- dim(matCountDataProc)[1]
  # Order genes by time of extremum (peak/valley)
  vecSignificantIDs <- rownames(dfImpulseDE2Results[!is.na(dfImpulseDE2Results$padj) & 
                                                      dfImpulseDE2Results$padj<scaQThres,])
  vecTimePointsToEval <- sort(unique(dfAnnotationProc$Time), decreasing=FALSE)
  scaNTPtoEvaluate <- length(vecTimePointsToEval)
  matImpulseValue <- do.call(rbind, lapply(vecSignificantIDs, function(x){
    evalImpulse_comp(vecImpulseParam=lsModelFits[[strCondition]][[x]]$lsImpulseFit$vecImpulseParam,
                     vecTimepoints=vecTimePointsToEval)
  }))
  rownames(matImpulseValue) <- vecSignificantIDs
  matidxMaxTimeSort <-t(apply(matImpulseValue, 1, function(genevalues){
    sort(genevalues, decreasing=TRUE, index.return=TRUE)$ix
  }))
  vecMaxTime <- vecTimePointsToEval[matidxMaxTimeSort[,1]]
  matidxMinTimeSort <- t(apply(matImpulseValue, 1, function(genevalues){
    sort(genevalues, decreasing=FALSE, index.return=TRUE)$ix
  }))
  
  if(boolIdentifyTransients){
    # Group into transients and monotonous
    vecidxTransient <- which(dfImpulseDE2Results[vecSignificantIDs,]$isTransient)
    #vecidxTransientSort <- vecidxTransient[do.call(order, as.data.frame(matidxMaxTimeSort[vecidxTransient,1]))]
    
    vecidxMonotonous <- which(dfImpulseDE2Results[vecSignificantIDs,]$isMonotonous)
    #vecidxMonotonousSort <- vecidxMonotonous[do.call(order, as.data.frame(matidxMaxTimeSort[vecidxMonotonous,1]))]
    
    # Fine sort montonous transition signals into up/down
    if(length(vecidxMonotonous) == 1) matImpulseMonot <- t(matImpulseValue[vecidxMonotonous,])
    else matImpulseMonot <- matImpulseValue[vecidxMonotonous,]
    vecboolMonotonousUp <- apply(matImpulseMonot, 1, function(gene){
      gene[1] < gene[scaNTPtoEvaluate]
    })
    vecboolMonotonousDown <- !vecboolMonotonousUp
    
    vecidxMonotonousUp <- vecidxMonotonous[vecboolMonotonousUp]
    vecidxMonotonousUpSort <- vecidxMonotonousUp[do.call(order, as.data.frame(matidxMinTimeSort[vecidxMonotonousUp,1]))]
    vecidxMonotonousDown <- vecidxMonotonous[vecboolMonotonousDown]
    vecidxMonotonousDownSort <- vecidxMonotonousDown[do.call(order, as.data.frame(matidxMinTimeSort[vecidxMonotonousDown,1]))]
    
    # Fine sort transitive signals into off/on
    if(length(vecidxTransient) == 1) matImpulseTransient <- t(matImpulseValue[vecidxTransient,])
    else matImpulseTransient <- matImpulseValue[vecidxTransient,]
    vecboolTransientValley <- apply(matImpulseTransient, 1, function(genevalues){
      boolValley <- any(genevalues[2:(scaNTPtoEvaluate-1)] < genevalues[1] &
                          genevalues[2:(scaNTPtoEvaluate-1)] < genevalues[scaNTPtoEvaluate])
      return(boolValley)
    })
    vecboolTransientPeak <- !vecboolTransientValley
    
    vecidxTransientPeak <- vecidxTransient[vecboolTransientPeak]
    vecidxTransientPeakSort <- vecidxTransientPeak[do.call(order, as.data.frame(matidxMaxTimeSort[vecidxTransientPeak,1]))]
    vecidxTransientValley <- vecidxTransient[vecboolTransientValley]
    vecidxTransientValleySort <- vecidxTransientValley[do.call(order, as.data.frame(matidxMinTimeSort[vecidxTransientValley,1]))]
    
    vecidxAllSort <- c(vecidxMonotonousUpSort,
                       vecidxMonotonousDownSort,
                       vecidxTransientPeakSort,
                       vecidxTransientValleySort)
    
    vecTrajectoryType <- c(rep("up", length(vecidxMonotonousUpSort)), 
                           rep("down", length(vecidxMonotonousDownSort)), 
                           rep("*up", length(vecidxTransientPeakSort)),
                           rep("*down", length(vecidxTransientValleySort)) )
    lsvecGeneGroups <- list( transition_up   = vecSignificantIDs[vecidxMonotonousUpSort],
                             transition_down = vecSignificantIDs[vecidxMonotonousDownSort],
                             transient_up    = vecSignificantIDs[vecidxTransientPeakSort],
                             transient_down  = vecSignificantIDs[vecidxTransientValleySort] )
  } else {
    vecidxAllSort <- sort(vecMaxTime, decreasing=FALSE, index.return=TRUE)$ix
    vecTrajectoryType <- rep(" ", length(vecidxAllSort))
    lsvecGeneGroups <- list(all = vecSignificantIDs[vecidxAllSort])
  }
  
  #vecidxAllSort <- vecidxMaxSort
  #vecTrajectoryType <- rep("up", length(vecidxAllSort))
  
  vecUniqueTP <- unique(dfAnnotationProc$Time)
  vecSizeFactors <- computeNormConst(
    matCountDataProc=matCountDataProc,
    vecSizeFactorsExternal=vecSizeFactors )
  matSizefactors <- matrix(vecSizeFactors,
                           nrow=length(vecSignificantIDs),
                           ncol=dim(matCountDataProc)[2],
                           byrow=TRUE)
  
  # 1. Plot raw data
  matDataNorm <- matCountDataProc[vecSignificantIDs,]/matSizefactors
  #if(length(objectImpulseDE2@lsModelFits[[strCondition]][[1]]$lsImpulseFit$lsvecBatchFactors)>0){
  #  for(vecBatchFactors in )
  #}
  matDataHeat <- do.call(cbind, lapply(vecUniqueTP, function(tp){
    vecidxCols <- which(dfAnnotationProc$Time %in% tp)
    if(length(vecidxCols)>1){ return(rowMeans(matDataNorm[,vecidxCols], na.rm=TRUE))
    } else { return(matDataNorm[,vecidxCols]) }
  }))
  colnames(matDataHeat) <- vecUniqueTP
  rownames(matDataHeat) <- NULL
  # Row normalisation: z-scores
  matDataHeatZ <- do.call(rbind, lapply(seq(1,dim(matDataHeat)[1]), function(i){
    (matDataHeat[i,]-mean(matDataHeat[i,], na.rm=TRUE))/sd(matDataHeat[i,], na.rm=TRUE)
  }))
  
  # Create heatmap of raw data, ordered by fits
  complexHeatmapRaw <- Heatmap(matDataHeatZ[vecidxAllSort,],
                               gap                  =unit(5, "mm"),
                               split                = vecTrajectoryType,
                               cluster_rows         = FALSE,
                               cluster_columns      = FALSE,
                               heatmap_legend_param = list(title = "z-score") )
  
  # 2. Plot fitted data
  matDataHeat <- matImpulseValue
  colnames(matDataHeat) <- vecUniqueTP
  rownames(matDataHeat) <- NULL
  # Row normalisation: z-scores
  matDataHeatZ <- do.call(rbind, lapply(seq(1,dim(matDataHeat)[1]), function(i){
    (matDataHeat[i,]-mean(matDataHeat[i,], na.rm=TRUE))/sd(matDataHeat[i,], na.rm=TRUE)
  }))
  
  # Create heatmap of raw data, ordered by fits
  complexHeatmapFit <- Heatmap(matDataHeatZ[vecidxAllSort,],
                               gap                  = unit(5, "mm"),
                               split                = vecTrajectoryType,
                               cluster_rows         = FALSE,
                               cluster_columns      = FALSE,
                               heatmap_legend_param = list(title = "z-score") )
  
  return(list(
    complexHeatmapRaw = complexHeatmapRaw,
    complexHeatmapFit = complexHeatmapFit,
    lsvecGeneGroups   = lsvecGeneGroups
  ))
}