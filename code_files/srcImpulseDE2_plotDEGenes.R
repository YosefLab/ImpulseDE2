#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Plot impulse fits    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Plots the impulse fits and data to pdf
#' 
#' Plots the impulse fits and data to pdf.
#' 
#' @seealso Called by \code{runImpulseDE2}.
#' 
#' @param vecGeneIDs (string vector) Gene names to be plotted,
#      must be rownames of matCountDataProc
#' @param matCountDataProc: (matrix genes x samples)
#'    Count data: Reduced version of \code{matCountData}. 
#'    For internal use.
#' @param dfAnnotationProc: (Table) Processed annotation table. 
#'    Lists co-variables of samples: 
#'    Sample, Condition, Time (numeric), TimeCateg (str)
#'    (and LongitudinalSeries). For internal use.
#' @param matNormConst: (numeric matrix genes x samples) 
#'    Model scaling factors for each observation: Take
#'    sequencing depth and longitudinal time series mean
#'    within a gene into account (size and translation
#'    factors).
#' @param matSizeFactors: (numeric matrix genes x samples) 
#'    Model scaling factors for each observation only accounting
#'    for sequencing depth into account (size factors).
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
#' @param strCaseName (str) Name of the case condition in \code{dfAnnotationProcFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#'    \code{dfAnnotationProcFull}.
#' @param strMode: (str) [Default "batch"] {"batch","longitudinal","singlecell"}
#'    Mode of model fitting.
#' @param NPARAM (scalar) [Default 6] Number of parameters of impulse model.
#' @param strFileNameSuffix (character string) [Default ""] File extention.
#' @param title_string (character string) [Default ""] Title for each plot.
#' @param strPlotSubtitle (character string) [Default ""] Subtitle for each plot.
#' @param NPARAM (scalar) [Default 6] Number of parameters of impulse model.
#' 
#' @return NULL
#' @export

plotDEGenes <- function(vecGeneIDs, 
  matCountDataProc, 
  matNormConst, matSizeFactors,
  dfAnnotationProc, lsImpulseFits, matClusterMeansFitted=NULL,
  dfImpulseResults, vecRefPval=NULL, 
  strCaseName, strControlName=NULL, strMode="batch",
  strFileNameSuffix = "", strPlotTitleSuffix = "", strPlotSubtitle = "",
  strNameMethod1="ImpulseDE2", strNameMethod2=NULL,
  boolSimplePlot=FALSE, boolLogPlot=FALSE,
  NPARAM=6){
  
  # Reconstruct translation factors
  matTranslationFactors <- matNormConst / matSizeFactors
  rownames(matTranslationFactors) <- rownames(matNormConst)
  # Scale count data by size factors for plotting:
  # The impulse models are fit based on normalised means.
  # Therefore, the model curves follow the normalised
  # count data and not the raw count data. However, 
  # fitting was still performed based on raw count data.
  matCountDataProc <- matCountDataProc / matSizeFactors
  if(boolLogPlot){
    # Add Pseudocount for plotting convenience
    matCountDataProc <- log(matCountDataProc+1)/log(10)
    matClusterMeansFitted <- log(matClusterMeansFitted)/log(10)
    strLogPlot <- "log_10"
  } else {
    strLogPlot <- NULL
  }
  
  # Name genes if not previously named (already done if this function is
  # called within the wrapper function).
  if(length(grep("[a-zA-Z]",rownames(matCountDataProc))) == 0 ){
    rownames(matCountDataProc) <- paste(rownames(matCountDataProc),"G", sep = "_")
    vecGeneIDs <- paste(vecGeneIDs, "G", sep = "_")
  }
  
  vecConditions <- as.vector( dfAnnotationProc[match(
    colnames(matCountDataProc),
    dfAnnotationProc$Sample),
    ]$Condition )
  vecTimepointAssign <- as.vector( dfAnnotationProc[match(
    colnames(matCountDataProc),
    dfAnnotationProc$Sample),
    ]$Time )
  vecTimepoints <- sort(unique(vecTimepointAssign))
  if(!is.null(strControlName)){
    vecidxSamplesCase <- vecConditions %in% strCaseName
    vecidxSamplesCtrl <- vecConditions %in% strControlName
  }
  if(strMode=="longitudinal"){
    vecLongitudinalSeriesAssign <- dfAnnotationProc[match(
      colnames(matCountDataProc),
      dfAnnotationProc$Sample),
      ]$LongitudinalSeries
    vecLongitudinalSeries <- unique(vecLongitudinalSeriesAssign)
    vecindLongitudinalSeriesAssign <- match(vecLongitudinalSeriesAssign, vecLongitudinalSeries)
    vecindLongitudinalSeriesAssignUnique <- sapply(vecLongitudinalSeries,
      function(longser){min(which(vecLongitudinalSeriesAssign %in% longser))})
  }
  
  # Only for batch/singlecell plotting:
  # Width of negative binomial pdf (ylim) in time units for plotting
  # E.g. PDF_WIDTH <- 1 means that the inferred negative binomial
  # density of a sample will be scaled into the interval [0,1]
  # and plotted at sample x time units vertically in the interval
  # [x,x+1] time units, with the y axis (the counts) being the
  # horizontal axis of the pdf.
  PDF_WIDTH <- max(vecTimepoints)/length(vecTimepoints)
  
  # Open .pdf
  pdf(paste("ImpulseDE2_",strFileNameSuffix,".pdf",sep=""),height=6.0,width=9.0)
  
  # Define grid for printing plots
  if (length(vecGeneIDs) == 1){
    par(mfrow=c(1,1), xpd=TRUE)
    scaLegendInset <- -0.15
  } else if (length(vecGeneIDs) <= 4){
    par(mfrow=c(2,2), xpd=TRUE)
    scaLegendInset <- -0.35
  } else if (length(vecGeneIDs) <= 6){
    par(mfrow=c(2,3), xpd=TRUE)
    scaLegendInset <- -0.3
  } else {
    par(mfrow=c(3,3), xpd=TRUE)
    scaLegendInset <- -0.65
  }
  # Time points for plotting of impulse model
  vecX <- seq(min(vecTimepoints),max(vecTimepoints),min(0.1,(max(vecTimepoints)-min(vecTimepoints))/100))
  # Find elements in vecX corresponding (closest) to observed time point
  indVecXObs <- unlist(lapply( 
    vecTimepoints, function(t){match( min(abs(vecX-t)), abs(vecX-t) )} 
  ))
  for (geneID in vecGeneIDs){
    if(is.null(strControlName)){
      # Without control  data 
      
      # Convert h0,h1,h2 to log space again
      lsImpulseParamCaseLog <- lsImpulseFits$parameters_case[geneID,1:NPARAM]
      lsImpulseParamCaseLog[c("h0","h1","h2")] <- log( lsImpulseParamCaseLog[c("h0","h1","h2")] )
      vecCaseValues <- calcImpulse_comp(lsImpulseParamCaseLog,vecX)
      if(boolLogPlot){ vecCaseValues <- log(vecCaseValues+1)/log(10) }
      
      pval_Impulse <- round( log(dfImpulseResults[geneID,]$adj.p)/log(10), 2 )
      if(!is.null(vecRefPval)){
        pval_Method2 <- round( log(vecRefPval[geneID])/log(10), 2 )
      }
      
      if(boolSimplePlot){
        # Only plot points and impulse trace
        # Plot observed points in blue - all time courses in same colour
        scaYlim_lower <- min( min(matCountDataProc[geneID,],na.rm=TRUE), min(vecCaseValues[indVecXObs]) )
        scaYlim_upper <- max( max(matCountDataProc[geneID,],na.rm=TRUE), max(vecCaseValues[indVecXObs]) )
        strPvalImpulse <- paste0(strNameMethod1," ",pval_Impulse)
        if(!is.null(vecRefPval)){
          strPvalMethod2 <- paste0(" ",strNameMethod2," ",pval_Method2)
        } else {
          strPvalMethod2 <- NULL
        }
        plot(vecTimepointAssign,
          matCountDataProc[geneID,],
          col="blue",pch=3,
          xlim=c(0,max(vecTimepoints,na.rm=TRUE)+PDF_WIDTH), 
          ylim=c(scaYlim_lower,scaYlim_upper),
          xlab="Time", ylab=paste0(strLogPlot," Impulse fit and count data"),
          main=paste0(geneID," ",strPlotTitleSuffix," log(Pval):\n ",strPvalImpulse,
            strPvalMethod2 ),sub=strPlotSubtitle)
        # Plot impulse model within boundaries of observed points
        vecCaseValuesToPlot <- vecCaseValues
        indImpulseValToPlot <- vecCaseValuesToPlot >= scaYlim_lower & vecCaseValuesToPlot <= scaYlim_upper
        vecCaseValuesToPlot[!indImpulseValToPlot] <- NA
        points(vecX, vecCaseValuesToPlot,col="black", type="l")
        
        # Plot mean of each time point
        points(vecTimepoints,
          sapply(vecTimepoints,function(tp){mean(matCountDataProc[geneID,vecTimepointAssign==tp],na.rm=TRUE)}),
          col="black",pch=1)
        if(strMode=="singlecell"){
          # Plot inferred mean of each time point
          points(vecTimepoints,
            matClusterMeansFitted[geneID,],
            col="red",pch=1)
        }
      } else {
        # Plot points, impulse trace and
        # batch: inferred negative binomials
        # longitudinal: longitudinal series scaled impulse traces
        
        if(strMode=="batch" | strMode=="singlecell"){
          # Plot observed points in blue - all time courses in same colour
          scaYlim_lower <- min( min(matCountDataProc[geneID,],na.rm=TRUE), min(vecCaseValues[indVecXObs]) )
          scaYlim_upper <- max( max(matCountDataProc[geneID,],na.rm=TRUE), max(vecCaseValues[indVecXObs]) )
          strPvalImpulse <- paste0(strNameMethod1," ",pval_Impulse)
          if(!is.null(vecRefPval)){
            strPvalMethod2 <- paste0(" ",strNameMethod2," ",pval_Method2)
          } else {
            strPvalMethod2 <- NULL
          }
          plot(vecTimepointAssign,
            matCountDataProc[geneID,],
            col="blue",pch=3,
            xlim=c(0,max(vecTimepoints,na.rm=TRUE)+PDF_WIDTH), 
            ylim=c(scaYlim_lower,scaYlim_upper),
            xlab="Time", ylab=paste0(strLogPlot," Impulse fit and count data"),
            main=paste0(geneID," ",strPlotTitleSuffix," log(Pval):\n ",strPvalImpulse,
              strPvalMethod2 ),sub=strPlotSubtitle)
          # Plot impulse model within boundaries of observed points
          vecCaseValuesToPlot <- vecCaseValues
          indImpulseValToPlot <- vecCaseValuesToPlot >= scaYlim_lower & vecCaseValuesToPlot <= scaYlim_upper
          vecCaseValuesToPlot[!indImpulseValToPlot] <- NA
          points(vecX, vecCaseValuesToPlot,col="blue", type="l")
          
          # Plot inferred negative binomial pdf at each time point in black (vertical)
          vecXCoordPDF <- seq(round(scaYlim_lower),round(scaYlim_upper), by=1 )
          # Get means of negative binomial at each time point
          vecCaseValueAtTP <- calcImpulse_comp(lsImpulseParamCaseLog,vecTimepoints)
          for(tp in vecTimepoints){
            vecYCoordPDF <- dnbinom(vecXCoordPDF,mu=vecCaseValueAtTP[match(tp,vecTimepoints)],
              size=as.numeric(as.vector(dfImpulseResults[geneID,]$size)) )
            # Scale Y_coord to uniform peak heights of 1
            # This translates into width of one time unit in plot
            vecYCoordPDF <- vecYCoordPDF * PDF_WIDTH/max(vecYCoordPDF)
            # Plot pdf vertically at time point
            lines(x=tp+vecYCoordPDF,y=vecXCoordPDF,col="black")
          }
          
          # Plot mean of each time point
          points(vecTimepoints,
            sapply(vecTimepoints,function(tp){mean(matCountDataProc[geneID,vecTimepointAssign==tp],na.rm=TRUE)}),
            col="black",pch=1)
          
          legend(x="bottomright",
            legend=c(strCaseName),
            fill=c("blue"), 
            cex=0.6, 
            inset=c(0,scaLegendInset))
        } else if(strMode=="longitudinal"){
          vecTranslationFactors <- matTranslationFactors[geneID,vecindLongitudinalSeriesAssignUnique]
          
          # Create colour vector
          vecCol <- rainbow(n=length(vecLongitudinalSeries))
          
          # Create plot and plot data of first time course
          scaYlim_lower <- min( min(matCountDataProc[geneID,],na.rm=TRUE), 
            min(vecCaseValues[indVecXObs]*min(vecTranslationFactors)) )
          scaYlim_upper <- max( max(matCountDataProc[geneID,],na.rm=TRUE), 
            max(vecCaseValues[indVecXObs]*max(vecTranslationFactors)) )
          strPvalImpulse <- paste0(strNameMethod1," ",pval_Impulse)
          if(!is.null(vecRefPval)){
            strPvalMethod2 <- paste0(" ",strNameMethod2," ",pval_Method2)
          } else {
            strPvalMethod2 <- NULL
          }
          plot(vecTimepointAssign[vecindLongitudinalSeriesAssign==1],
            matCountDataProc[geneID, vecindLongitudinalSeriesAssign==1],
            col=vecCol[1],pch=3,
            xlim=c(0,max(vecTimepoints,na.rm=TRUE)), 
            ylim=c(scaYlim_lower,scaYlim_upper),
            xlab="Time", ylab=paste0(strLogPlot," Impulse fit and count data"),
            main=paste0(geneID," ",strPlotTitleSuffix," log(Pval):\n ",strPvalImpulse,
              strPvalMethod2 ),sub=strPlotSubtitle)
          # Plot impulse fit to time course
          vecCaseValuesToPlotTC <- vecCaseValues*vecTranslationFactors[1]
          indImpulseValToPlot <- vecCaseValuesToPlotTC >= scaYlim_lower & vecCaseValuesToPlotTC <= scaYlim_upper
          vecCaseValuesToPlotTC[!indImpulseValToPlot] <- NA
          points(vecX, vecCaseValuesToPlotTC,col=vecCol[1], type="l")
          
          # Plot remaining time courses
          if(length(vecLongitudinalSeries)>1){
            for(longser in 2:length(vecLongitudinalSeries)){
              # Plot data of time course
              points(x=vecTimepointAssign[vecindLongitudinalSeriesAssign==longser],
                y=matCountDataProc[geneID, vecindLongitudinalSeriesAssign==longser],
                col=vecCol[longser],pch=3)
              # Plot impulse within boundaries of observed points
              vecCaseValuesToPlotTC <- vecCaseValues*vecTranslationFactors[longser]
              indImpulseValToPlot <- vecCaseValuesToPlotTC >= scaYlim_lower & vecCaseValuesToPlotTC <= scaYlim_upper
              vecCaseValuesToPlotTC[!indImpulseValToPlot] <- NA
              points(vecX, 
                vecCaseValuesToPlotTC,
                col=vecCol[longser], type="l")
            }
          }
          
          # Plot mean of each time point
          points(vecTimepoints,
            sapply(vecTimepoints,function(tp){mean(matCountDataProc[geneID,vecTimepointAssign==tp],na.rm=TRUE)}),
            col="black",pch=1)
          
          legend(x="bottomright",
            legend=vecLongitudinalSeries,
            fill=vecCol, 
            cex=0.6, 
            inset=c(0,scaLegendInset-0.04*length(vecLongitudinalSeries)))
        } else {
          stop(paste0("ERROR: Unrecognised strMode in plotDEGenes(): ",strMode))
        }
      }
      
    } else {
      # With control data:
      # Do not distinguish batch/singlecell/longitudinals plotting:
      # Plots become to crowded if more than 3 traces (case, control
      # combined) and two point clouds (case, control) are plotted.
      
      # Convert h0,h1,h2 to log space again
      lsImpulseParamCaseLog <- lsImpulseFits$parameters_case[geneID,1:NPARAM]
      lsImpulseParamCaseLog[c("h0","h1","h2")] <- log( lsImpulseParamCaseLog[c("h0","h1","h2")] )
      vecCaseValues <- calcImpulse_comp(lsImpulseParamCaseLog,vecX)
      
      lsImpulseParamCtrlLog <- lsImpulseFits$parameters_control[geneID,1:NPARAM]
      lsImpulseParamCtrlLog[c("h0","h1","h2")] <- log( lsImpulseParamCtrlLog[c("h0","h1","h2")] )
      lsCtrlValues <- calcImpulse_comp(lsImpulseParamCtrlLog,vecX)
      
      lsImpulseParamCombLog <- lsImpulseFits$parameters_combined[geneID,1:NPARAM]
      lsImpulseParamCombLog[c("h0","h1","h2")] <- log( lsImpulseParamCombLog[c("h0","h1","h2")] )
      lsCombValues <- calcImpulse_comp(lsImpulseParamCombLog,vecX)
      
      pval_Impulse <- round( log(dfImpulseResults[geneID,]$adj.p)/log(10), 2 )
      if(!is.null(vecRefPval)){
        pval_Method2 <- round( log(vecRefPval[geneID])/log(10), 2 )
      }
      
      # Plot observed points: Case green, control red
      scaYlim_lower <- min( min(matCountDataProc[geneID,],na.rm=TRUE),
        min(vecCaseValues[indVecXObs]),
        min(lsCtrlValues[indVecXObs]),
        min(lsCombValues[indVecXObs]) )
      scaYlim_upper <- max( max(matCountDataProc[geneID,],na.rm=TRUE),
        max(vecCaseValues[indVecXObs]),
        max(lsCtrlValues[indVecXObs]),
        max(lsCombValues[indVecXObs]) )
      strPvalImpulse <- paste0(strNameMethod1," ",pval_Impulse)
      if(!is.null(vecRefPval)){
        strPvalMethod2 <- paste0(" ",strNameMethod2," ",pval_Method2)
      } else {
        strPvalMethod2 <- NULL
      }
      # Plot points case
      plot(vecTimepointAssign[vecidxSamplesCase],
        (matCountDataProc[geneID,])[vecidxSamplesCase],
        col="green",pch=3,
        xlim=c(0,max(vecTimepoints,na.rm=TRUE)+PDF_WIDTH),
        ylim=c(scaYlim_lower,scaYlim_upper),
        xlab="Time", ylab=paste0(strLogPlot," Impulse fit and count data"),
        main=paste0(geneID," ",strPlotTitleSuffix," log(Pval):\n ",
          strPvalImpulse," ",strPvalMethod2),sub=strPlotSubtitle)
      # Plot points control
      points(vecTimepointAssign[vecidxSamplesCtrl],
        (matCountDataProc[geneID,])[vecidxSamplesCtrl],
        col="red",pch=3, type="p")
      # Plot impulse models within boundaries of observed points for 
      # case, control and combined.
      # Plot impulse model case:
      vecCaseValuesToPlot <- vecCaseValues
      indImpulseValToPlot <- vecCaseValuesToPlot >= scaYlim_lower & vecCaseValuesToPlot <= scaYlim_upper
      vecCaseValuesToPlot[!indImpulseValToPlot] <- NA
      points(vecX, vecCaseValuesToPlot,col="green", type="l")
      # Plot impulse model control:
      lsCtrlValuesToPlot <- lsCtrlValues
      indImpulseValToPlot <- lsCtrlValuesToPlot >= scaYlim_lower & lsCtrlValuesToPlot <= scaYlim_upper
      lsCtrlValuesToPlot[!indImpulseValToPlot] <- NA
      points(vecX, lsCtrlValuesToPlot,col="red", type="l")
      # Plot impulse model combined:
      lsCombValuesToPlot <- lsCombValues
      indImpulseValToPlot <- lsCombValuesToPlot >= scaYlim_lower & lsCombValuesToPlot <= scaYlim_upper
      lsCombValuesToPlot[!indImpulseValToPlot] <- NA
      points(vecX, lsCombValuesToPlot,col="black", type="l")
      
      # Plot mean of each time point
      points(vecTimepoints,
        sapply(vecTimepoints,function(tp){mean(matCountDataProc[geneID,vecTimepointAssign==tp],na.rm=TRUE)}),
        col="black",pch=1)
      
      legend(x="bottomright",
        legend=c(strCaseName,strControlName),
        fill=c("green","red"), 
        cex=0.6, 
        inset=c(0,scaLegendInset))
      
    }
  }
  # Close .pdf
  dev.off()
}