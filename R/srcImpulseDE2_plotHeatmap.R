
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
  vecMinTime <- vecTimePointsToEval[matidxMinTimeSort[,1]]
  
  if(boolIdentifyTransients){
    # Group into transients and monotonous
    vecidxTransient <- which(dfImpulseDE2Results[vecSignificantIDs,]$isTransient)
    #vecidxTransientSort <- vecidxTransient[do.call(order, as.data.frame(matidxMaxTimeSort[vecidxTransient,1]))]
    
    vecidxMonotonous <- which(dfImpulseDE2Results[vecSignificantIDs,]$isMonotonous)
    #vecidxMonotonousSort <- vecidxMonotonous[do.call(order, as.data.frame(matidxMaxTimeSort[vecidxMonotonous,1]))]
    
    # Finbe sort montonous transition signals into up/down
    vecboolMonotonousUp <- apply(matImpulseValue[vecidxMonotonous,], 1, function(gene){
      gene[1] < gene[scaNTPtoEvaluate]
    })
    vecboolMonotonousDown <- !vecboolMonotonousUp
    
    vecidxMonotonousUp <- vecidxMonotonous[vecboolMonotonousUp]
    vecidxMonotonousUpSort <- vecidxMonotonousUp[do.call(order, as.data.frame(matidxMinTimeSort[vecidxMonotonousUp,1]))]
    vecidxMonotonousDown <- vecidxMonotonous[vecboolMonotonousDown]
    vecidxMonotonousDownSort <- vecidxMonotonousDown[do.call(order, as.data.frame(matidxMinTimeSort[vecidxMonotonousDown,1]))]
    
    # Fine sort transitive signals into off/on
    vecboolTransientValley <- apply(matImpulseValue[vecidxTransient,], 1, function(genevalues){
      # Assume no complete constant fit initially
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
  matDataHeat <- matImpulseValue
  colnames(matDataHeat) <- vecUniqueTP
  rownames(matDataHeat) <- NULL
  matDataHeat <- matDataHeat[apply(matDataHeat, 1, function(gene) any(gene>0) ),]
  # Row normalisation: z-scores
  matDataHeatZ <- do.call(rbind, lapply(seq(1,dim(matDataHeat)[1]), function(i){
    (matDataHeat[i,]-mean(matDataHeat[i,], na.rm=TRUE))/sd(matDataHeat[i,], na.rm=TRUE)
  }))
  
  # Create heatmap of raw data, ordered by fits
  complexHeatmapFit <- Heatmap(matDataHeatZ[vecidxAllSort,],
                               #row_order                = vecidxAllSort,
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