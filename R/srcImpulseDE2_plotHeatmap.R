
plotHeatmap <- function(matCountDataProc,
                        lsModelFits,
                        dfAnnotationProc,
                        strCondition){
  
  scaEvaluationRes <- 100 # Number of time points to evalute for each model
  # Order genes by time of extremum (peak/valley)
  # Uses model fits for this
  vecTimePoints <- seq(min(dfAnnotationProc$Time), max(dfAnnotationProc$Time), length.out=scaEvaluationRes)
  vecPeakTime <- sapply(lsModelFits[[strCondition]], function(genefit){
    
  })
  
  # Create heatmap of raw data, ordered by fits
  complexHeatmapRaw <- Heatmap(matDataHeatZ,
                               row_order                = vecidxSortByCluster,
                               split                    = letters[vecidxClustersSorted],
                               #km                       = 6, # Number of k-means clusters
                               cluster_rows             = FALSE, #hclustObject,
                               #show_row_dend            = FALSE,
                               cluster_columns          = FALSE,
                               heatmap_legend_param = list(title = "z-score", color_bar = "discrete") )
  
  
  # Create heatmap of fitted models, ordered by fits
  return(list(
    complexHeatmapRaw=complexHeatmapRaw,
    complexHeatmapFit=complexHeatmapFit
  ))
}