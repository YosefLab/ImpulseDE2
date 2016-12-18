
plotHeatmap <- function(matCountDataProc,
                        dfAnnotationProc,
                        lsModelFits,
                        dfImpulseDE2Results,
                        boolIdentifyTransients,
                        strCondition,
                        vecSizeFactors=NULL,
                        scaQThres=0.001){
  
  scaNTPtoEvaluate <- 100 # Number of time points to evalute for each model
  scaNGenes <- dim(matCountDataProc)[1]
  # Order genes by time of extremum (peak/valley)
  vecSignificantIDs <- rownames(dfImpulseDE2Results[dfImpulseDE2Results$padj<scaQThres,])
  vecTimePointsToEval <- seq(min(dfAnnotationProc$Time), max(dfAnnotationProc$Time), length.out=scaNTPtoEvaluate)
  matImpulseValue <- do.call(rbind, lapply(vecSignificantIDs, function(x){
    evalImpulse_comp(vecImpulseParam=lsModelFits[[strCondition]][[x]]$lsImpulseFit$vecImpulseParam,
                                      vecTimepoints=vecTimePointsToEval)
  }))
  rownames(matImpulseValue) <- vecSignificantIDs
  vecMaxTime <- apply(matImpulseValue, 1, function(genevalues){
    idxPeakTime <- sort(genevalues, decreasing=TRUE, index.return=TRUE)$ix[1]
    return(vecTimePointsToEval[idxPeakTime])
  })
  vecMinTime <- apply(matImpulseValue, 1, function(genevalues){
    idxPeakTime <- sort(genevalues, decreasing=FALSE, index.return=TRUE)$ix[1]
    return(vecTimePointsToEval[idxPeakTime])
  })
  
  if(boolIdentifyTransients){
    vecboolValleyUnSig <- apply(matImpulseValue, 1, function(genevalues){
      # Assume no complete constant fit initially
      boolValley <- any(genevalues[2:(scaNTPtoEvaluate-1)] < genevalues[1] &
                          genevalues[2:(scaNTPtoEvaluate-1)] < genevalues[scaNTPtoEvaluate])
      return(boolValley)
    })
    vecboolValley <- vecboolValleyUnSig & dfImpulseDE2Results[vecSignificantIDs,]$isTransient
    vecboolPeak <- !vecboolValley
    
    vecidxPeak <- which(vecboolPeak)
    vecidxPeakSort <- vecidxPeak[sort(vecMaxTime[vecidxPeak], decreasing=FALSE, index.return=TRUE)$ix]
    vecidxValley <- which(vecboolValley)
    vecidxValleySort <- vecidxValley[sort(vecMinTime[vecidxValley], decreasing=FALSE, index.return=TRUE)$ix]
    
    vecidxAllSort <- c(vecidxPeakSort, vecidxValleySort)
    
    vecTrajectoryType <- c(rep("peak", length(vecidxPeakSort)), 
                           rep("valley", length(vecidxValleySort)))
  } else {
    vecidxAllSort <- sort(vecMaxTime, decreasing=FALSE, index.return=TRUE)$ix
    vecTrajectoryType <- rep(" ", length(vecidxAllSort))
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
  matDataHeat <- do.call(cbind, lapply(vecUniqueTP, function(tp){
    vecidxCols <- which(dfAnnotationProc$Time %in% tp)
    if(length(vecidxCols)>1){ return(rowMeans(matDataNorm[,vecidxCols], na.rm=TRUE))
    } else { return(matDataNorm[,vecidxCols]) }
  }))
  colnames(matDataHeat) <- vecUniqueTP
  rownames(matDataHeat) <- NULL
  matDataHeat <- matDataHeat[apply(matDataHeat, 1, function(gene) any(gene>0) ),]
  # Row normalisation: z-scores
  matDataHeatZ <- do.call(rbind, lapply(seq(1,dim(matDataHeat)[1]), function(i){
    (matDataHeat[i,]-mean(matDataHeat[i,], na.rm=TRUE))/sd(matDataHeat[i,], na.rm=TRUE)
  }))
  
  # Create heatmap of raw data, ordered by fits
  complexHeatmapRaw <- Heatmap(matDataHeatZ[vecidxAllSort,],
                               #row_order                = ,
                               split                    = vecTrajectoryType,
                               cluster_rows             = FALSE,
                               #show_row_dend            = FALSE,
                               cluster_columns          = FALSE,
                               heatmap_legend_param = list(title = "z-score") )
  draw(complexHeatmapRaw)
  
  # 2. Plot fitted data
  vecidxTPMatch <- sapply(vecUniqueTP, function(tp) which.min(abs(vecTimePointsToEval-tp)) )
  matDataHeat <- matImpulseValue[,vecidxTPMatch]
  colnames(matDataHeat) <- vecUniqueTP
  rownames(matDataHeat) <- NULL
  matDataHeat <- matDataHeat[apply(matDataHeat, 1, function(gene) any(gene>0) ),]
  # Row normalisation: z-scores
  matDataHeatZ <- do.call(rbind, lapply(seq(1,dim(matDataHeat)[1]), function(i){
    (matDataHeat[i,]-mean(matDataHeat[i,], na.rm=TRUE))/sd(matDataHeat[i,], na.rm=TRUE)
  }))
  
  # Create heatmap of raw data, ordered by fits
  complexHeatmapFit <- Heatmap(matDataHeatZ,
                               row_order                = vecidxAllSort,
                               split                    = vecTrajectoryType,
                               #km                       = 6, # Number of k-means clusters
                               cluster_rows             = FALSE, #hclustObject,
                               #show_row_dend            = FALSE,
                               cluster_columns          = FALSE,
                               heatmap_legend_param = list(title = "z-score") )
  draw(complexHeatmapFit)
  
  return(list(
    complexHeatmapRaw=complexHeatmapRaw,
    complexHeatmapFit=complexHeatmapFit
  ))
}