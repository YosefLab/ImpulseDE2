#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++     Plot ZINB fits    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Plots the zero inflated negative binomial fits and data to pdf.
#' 
#' Plots the zero inflated negative binomial fits and data to pdf.
#' 
#' @seealso Called by \code{runPseudoDE}.
#' 
#' @param vecGeneIDs: (string vector) Gene names to be plotted. 
#'    Elements must be rownames of matCounts.
#' @param matCounts: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param matMuCluster:
#' @param vecDispersions: (vector number of genes) Gene-wise 
#'    negative binomial dispersion coefficients.
#' @param matProbNB: (probability matrix genes x samples) 
#'    Probability of observations to come from negative binomial 
#'    component of mixture model.
#' @param vecClusterAssignments:
#' @param lsResultsClustering (list {"Assignments","Centroids","K"})
#'    \itemize{
#'      \item   Assignments: (integer vector length number of
#'        cells) Index of cluster assigned to each cell.
#'      \item   Centroids: 1D Coordinates of cluster centroids,
#'        one scalar per centroid.
#'      \item   K: (scalar) Number of clusters selected.
#'      }
#' @param dfAnnotationImpulseDE2: (Table)
#'    Annotation table. Lists covariates of samples: 
#'    Sample, Condition, Time. Time must be numeric. No longitudinal
#'    series given as scRNA-seq data are not longitudinal.
#' @param strPDFname: (str) Name of .pdf with plots.
#' 
#' @return NULL
#' @export

plotZINBfits <- function(vecGeneIDs, 
  matCounts,
  matMuCluster, 
  vecDispersions,
  matProbNB,
  vecClusterAssignments,
  lsResultsClustering,
  dfAnnotation, 
  strPDFname="PseudoDE_ZINBfits.pdf"){
  
  # Only for batch/singlecell plotting:
  # Width of negative binomial pdf (ylim) in time units for plotting
  # E.g. PDF_WIDTH <- 1 means that the inferred negative binomial
  # density of a sample will be scaled into the interval [0,1]
  # and plotted at sample x time units vertically in the interval
  # [x,x+1] time units, with the y axis (the counts) being the
  # horizontal axis of the pdf.
  #PDF_WIDTH <- 1
  
  # Scale count data by size factors for plotting:
  # The impulse models are fit based on normalised means.
  # Therefore, the model curves follow the normalised
  # count data and not the raw count data. However, 
  # fitting was still performed based on raw count data.
  #matCounts <- matCounts / matrix(vecNormConst, 
  #  nrow=dim(matCounts)[1], 
  #  ncol=dim(matCounts)[2], 
  #  byrow=TRUE)
  
  vecClusters <- sort(unique(vecClusterAssignments))
  vecindClusterAssignments <- match(vecClusterAssignments,vecClusters)
  
  # Define grid for printing plots
  par(mfrow=c(3,3), xpd=TRUE)
  scaLegendInset <- -0.65
  
  # Create list of plots
  lsPlots <- list()
  for (geneID in vecGeneIDs){
    # Plot observed points in blue - all time courses in same colour
    if(FALSE){
      scaYlim_lower <- min(matCounts[geneID,],na.rm=TRUE)
      scaYlim_upper <- max(matCounts[geneID,],na.rm=TRUE)
      plot(vecindClusterAssignments,
        matCounts[geneID,],
        col="blue",pch=3,
        xlim=c(0,max(vecindClusterAssignments,na.rm=TRUE)+PDF_WIDTH), 
        ylim=c(scaYlim_lower,scaYlim_upper),
        xlab="Cluster", ylab="Counts",
        main=geneID )
      
      # Plot inferred negative binomial pdf at each time point in black (vertical)
      vecXCoordPDF <- seq(round(scaYlim_lower),round(scaYlim_upper), by=1 )
      for(cl in vecClusters){
        vecYCoordPDF <- dnbinom(vecXCoordPDF,mu=matMuCluster[geneID,cl],
          size=vecDispersions[geneID])
        # Scale Y_coord to uniform peak heights of 1
        # This translates into width of one time unit in plot
        vecYCoordPDF <- vecYCoordPDF * PDF_WIDTH/max(vecYCoordPDF)
        # Plot pdf vertically at time point
        lines(x=match(cl,vecClusters)+vecYCoordPDF,y=vecXCoordPDF,col="black")
      }
      
      # Plot mean of each cluster
      points(seq(1,length(vecClusters)),
        matMuCluster[geneID,],
        col="black",pch=1)
    } else{
      for(cl in vecClusters){
        vecXCoordPDF <- seq(0,max(matCounts[geneID,vecindClusterAssignments==match(cl,vecClusters)],na.rm=TRUE))
        vecYCoordPDF <- dnbinom(vecXCoordPDF,mu=matMuCluster[geneID,cl],
          size=vecDispersions[geneID])
        # Plot both
        vecMixtures <- c("Negative binomial","Dropout")
        dfData <- data.frame( counts=as.vector(matCounts[geneID,vecindClusterAssignments==match(cl,vecClusters)] ),
          mixture=vecMixtures[as.numeric(matProbNB[geneID,vecindClusterAssignments==match(cl,vecClusters)] < 0.5)+1]  )
        dfDensity <- data.frame( counts=vecXCoordPDF, density=vecYCoordPDF*(dim(dfData)[1]))
        lsPlots[[length(lsPlots)+1]] <- suppressMessages( ggplot( ) +
            geom_histogram( data=dfData, aes(x=counts, y=..count.., fill=mixture), alpha=0.5 ) +
            geom_line( data=dfDensity, aes(x=counts, y=density), col = "blue") +
            labs(title=paste0(geneID," Cluster ", cl, " Observed cells: ",length(dfData$counts),
              "\n Non-zero counts: ", sum(dfData$counts!=0), " Inferred mean: ", round(matMuCluster[geneID,cl]))) +
            xlab("Counts") +
            ylab("Frequency") )
      }
    }
  }
  
  # Write plots into .pdf
  pdf(strPDFname,height=6.0,width=9.0)
  for(plotTF in lsPlots){
    suppressMessages( print(plotTF) )
  }
  dev.off()
}