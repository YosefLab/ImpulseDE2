#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++     Plot ZINB fits    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Plots the zero inflated negative binomial fits and data to pdf
#' 
#' Plots the zero inflated negative binomial fits and data to pdf
#' 
#' @seealso Called by \code{runPseudoDE}.
#' 
#' @param lsGeneIDs (character vector) of gene names to be plotted; must be
#     trownames of arr3DCountData.
#' @param arr2DCountData (2D array genes x replicates) Count data: Reduced 
#'    version of \code{matCountData}. For internal use.
#' @param vecNormConst: (numeric vector number of replicates) 
#'    Normalisation constants for each replicate.
#' @param dfAnnotation (Table) Lists co-variables of individual replicates:
#'    Replicate, Sample, Condition, Time. Time must be numeric.
#' @param lsImpulseFits (list length 2 or 6) List of matrices which
#'    contain parameter fits and model values for given time course for the
#'    case condition (and control and combined if control is present).
#'    Each parameter matrix is called parameter_'condition' and has the form
#'    (genes x [beta, h0, h1, h2, t1, t2, logL_H1, converge_H1, mu, logL_H0, 
#'    ) where beta to t2 are parameters of the impulse
#'    model, mu is the single parameter of the mean model, logL are
#'    log likelihoods of full (H1) and reduced model (H0) respectively, converge
#'    is convergence status of numerical optimisation of model fitting by
#'    \code{optim} from \code{stats}. Each value matrix is called
#'    value_'condition' and has the form (genes x time points) and contains the
#'    counts predicted by the impulse model at the observed time points.
#' @param dfDEAnalysis (data frame genes x fitting characteristics) 
#'    Summary of fitting procedure for each gene.
#' @param vecRefPval (vec length genes) Method 2 (DESeq2) adjusted p-values
#' @param strCaseName (str) Name of the case condition in \code{dfAnnotationFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#'    \code{dfAnnotationFull}.
#' @param strMode: (str) [Default "batch"] {"batch","timecourses","singlecell"}
#'    Mode of model fitting.
#' @param NPARAM (scalar) [Default 6] Number of parameters of impulse model.
#' @param strFileNameSuffix (character string) [Default ""] File extention.
#' @param title_string (character string) [Default ""] Title for each plot.
#' @param strPlotSubtitle (character string) [Default ""] Subtitle for each plot.
#' @param NPARAM (scalar) [Default 6] Number of parameters of impulse model.
#' 
#' @return NULL
#' @export

plotZINBfits <- function(lsGeneIDs, arr2DCountData,
  matClusterMeans, vecDispersions, matProbNB,
  vecClusterAssignments, lsResultsClustering,
  dfAnnotation, strPDFname="PseudoDE_ZINBfits.pdf"){
  
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
  #arr2DCountData <- arr2DCountData / matrix(vecNormConst, 
  #  nrow=dim(arr2DCountData)[1], 
  #  ncol=dim(arr2DCountData)[2], 
  #  byrow=TRUE)
  
  # Name genes if not previously named (already done if this function is
  # called within the wrapper function).
  if(length(grep("[a-zA-Z]",rownames(arr2DCountData))) == 0 ){
    rownames(arr2DCountData) <- paste(rownames(arr2DCountData),"G", sep = "_")
    lsGeneIDs <- paste(lsGeneIDs, "G", sep = "_")
  }
  
  vecClusters <- sort(unique(vecClusterAssignments))
  vecindClusterAssignments <- match(vecClusterAssignments,vecClusters)
  
  # Define grid for printing plots
  par(mfrow=c(3,3), xpd=TRUE)
  scaLegendInset <- -0.65
  
  lsPlots <- list()
  for (geneID in lsGeneIDs){
    # Plot observed points in blue - all time courses in same colour
    if(FALSE){
      scaYlim_lower <- min(arr2DCountData[geneID,],na.rm=TRUE)
      scaYlim_upper <- max(arr2DCountData[geneID,],na.rm=TRUE)
      plot(vecindClusterAssignments,
        arr2DCountData[geneID,],
        col="blue",pch=3,
        xlim=c(0,max(vecindClusterAssignments,na.rm=TRUE)+PDF_WIDTH), 
        ylim=c(scaYlim_lower,scaYlim_upper),
        xlab="Cluster", ylab="Counts",
        main=geneID )
      
      # Plot inferred negative binomial pdf at each time point in black (vertical)
      vecXCoordPDF <- seq(round(scaYlim_lower),round(scaYlim_upper), by=1 )
      for(cl in vecClusters){
        vecYCoordPDF <- dnbinom(vecXCoordPDF,mu=matClusterMeans[geneID,match(cl,vecClusters)],
          size=vecDispersions[geneID])
        # Scale Y_coord to uniform peak heights of 1
        # This translates into width of one time unit in plot
        vecYCoordPDF <- vecYCoordPDF * PDF_WIDTH/max(vecYCoordPDF)
        # Plot pdf vertically at time point
        lines(x=match(cl,vecClusters)+vecYCoordPDF,y=vecXCoordPDF,col="black")
      }
      
      # Plot mean of each cluster
      points(seq(1,length(vecClusters)),
        matClusterMeans[geneID,],
        col="black",pch=1)
    } else{
      for(cl in vecClusters){
        vecXCoordPDF <- seq(0,max(arr2DCountData[geneID,vecindClusterAssignments==match(cl,vecClusters)],na.rm=TRUE))
        vecYCoordPDF <- dnbinom(vecXCoordPDF,mu=matClusterMeans[geneID,match(cl,vecClusters)],
          size=vecDispersions[geneID])
        # Plot both
        vecMixtures <- c("Negative binomial","Dropout")
        dfData <- data.frame( counts=as.vector(arr2DCountData[geneID,vecindClusterAssignments==match(cl,vecClusters)] ),
          mixture=vecMixtures[as.numeric(matProbNB[geneID,vecindClusterAssignments==match(cl,vecClusters)] < 0.5)+1]  )
        dfDensity <- data.frame( counts=vecXCoordPDF, density=vecYCoordPDF*(dim(dfData)[1]))
        lsPlots[[length(lsPlots)+1]] <- suppressMessages( ggplot( ) +
          geom_histogram( data=dfData, aes(x=counts, y=..count.., fill=mixture), alpha=0.5 ) +
          geom_line( data=dfDensity, aes(x=counts, y=density), col = "blue") +
          labs(title=paste0(geneID," Cluster ", cl, " Observed cells: ",length(dfData$counts),
            "\n Non-zero counts: ", sum(dfData$counts!=0), " Inferred mean: ", round(matClusterMeans[geneID,match(cl,vecClusters)]))) +
          xlab("Counts") +
          ylab("Density") )
      }
    }
  }
  
  # Open .pdf
  pdf(strPDFname,height=6.0,width=9.0)
  for(plotTF in lsPlots){
    suppressMessages( print(plotTF) )
  }
  # Close .pdf
  dev.off()
}