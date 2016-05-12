#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Plot impulse fits    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Plots the impulse fits and data to pdf
#' 
#' Plots the impulse fits and data to pdf.
#' 
#' @seealso Called by \code{runImpulseDE2}.
#' 
#' @param lsGeneIDs (character vector) of gene names to be plotted; must be
#     trownames of arr3DCountData.
#' @param arr3DCountData (3D array genes x samples x replicates)
#'    Count data: \code{arr2DCountData} reshaped into a 3D array. For internal use.
#' @param matNormConst: (matrix samples x replicates) Normalisation
#'    constants for each replicate. Missing samples are set NA.
#' @param dfAnnotationRed (data frame) Reduced version of 
#'    \code{dfAnnotationFull}. Lists co-variables of samples: 
#'    Sample, Condition, Time. Time must be numeric. For internal use.
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
#' @param vecMethod2Results (vec length genes) Method 2 (DESeq2) adjusted p-values
#' @param strCaseName (str) Name of the case condition in \code{dfAnnotationRedFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#'    \code{dfAnnotationRedFull}.
#' @param NPARAM (scalar) [Default 6] Number of parameters of impulse model.
#' @param strFileNameSuffix (character string) [Default ""] File extention.
#' @param title_string (character string) [Default ""] Title for each plot.
#' @param strPlotSubtitle (character string) [Default ""] Subtitle for each plot.
#' @param NPARAM (scalar) [Default 6] Number of parameters of impulse model.
#' 
#' @return NULL
#' @export

plotDEGenes <- function(lsGeneIDs, arr3DCountData, matNormConst,
  dfAnnotationRed, lsImpulseFits, 
  dfImpulseResults, vecMethod2Results, 
  strCaseName, strControlName=NULL, strMode="batch",
  strFileNameSuffix = "", strPlotTitleSuffix = "", strPlotSubtitle = "",
  strNameMethod1="ImpulseDE2", strNameMethod2="DESeq2",
  NPARAM=6){
  
  # Only for batch/singlecell plotting:
  # Width of negative binomial pdf (ylim) in time units for plotting
  # E.g. PDF_WIDTH <- 1 means that the inferred negative binomial
  # density of a sample will be scaled into the interval [0,1]
  # and plotted at sample x time units vertically in the interval
  # [x,x+1] time units, with the y axis (the counts) being the
  # horizontal axis of the pdf.
  PDF_WIDTH <- 1
  
  # Scale count data by size factors for plotting:
  # The impulse models are fit based on normalised means.
  # Therefore, the model curves follow the normalised
  # count data and not the raw count data. However, 
  # fitting was still performed based on raw count data.
  for(gene in 1:dim(arr3DCountData)[1]){ arr3DCountData[gene,,] <- arr3DCountData[gene,,]/matNormConst }
  
  # Name genes if not previously named (already done if this function is
  # called within the wrapper function).
  if(length(grep("[a-zA-Z]",rownames(arr3DCountData))) == 0 ){
    rownames(arr3DCountData) <- paste(rownames(arr3DCountData),"G", sep = "_")
    lsGeneIDs <- paste(lsGeneIDs, "G", sep = "_")
  }
  
  if(!is.null(strControlName)){
    arr3DCountData_Case <- arr3DCountData[,as.vector(dfAnnotationRed[dfAnnotationRed$Condition %in% strCaseName,]$Sample),]
    arr3DCountData_Ctrl <- arr3DCountData[,as.vector(dfAnnotationRed[dfAnnotationRed$Condition %in% strControlName,]$Sample),]
    vecTimepoints_Case <- as.numeric(as.character(dfAnnotationRed[dfAnnotationRed$Condition %in% strCaseName,]$Time))
    vecTimepoints_Ctrl <- as.numeric(as.character(dfAnnotationRed[dfAnnotationRed$Condition %in% strControlName,]$Time))
  }
  vecTimepoints_Comb <- as.numeric(as.character(dfAnnotationRed[dfAnnotationRed$Sample %in% colnames(arr3DCountData),"Time"]))
  
  # Expand vecTimepoints_Comb vector to size of expression array
  if(!is.null(strControlName)){
    arrTimepoints_Case <- t(matrix(
      rep(vecTimepoints_Case,dim(arr3DCountData)[3]),
      length(vecTimepoints_Case),
      dim(arr3DCountData)[3] ))
    arrTimepoints_Ctrl <- t(matrix(
      rep(vecTimepoints_Ctrl,dim(arr3DCountData)[3]),
      length(vecTimepoints_Ctrl),
      dim(arr3DCountData)[3]))
  }
  arrTimepoints_Comb <- t(matrix(rep(vecTimepoints_Comb,dim(arr3DCountData)[3]),length(vecTimepoints_Comb),dim(arr3DCountData)[3]))
  
  # Open .pdf
  pdf(paste("ImpulseDE2_",strFileNameSuffix,".pdf",sep=""),height=6.0,width=9.0)
  
  # Define grid for printing plots
  if (length(lsGeneIDs) == 1){
    par(mfrow=c(1,1), xpd=TRUE)
    scaLegendInset <- -0.15
  } else if (length(lsGeneIDs) <= 4){
    par(mfrow=c(2,2), xpd=TRUE)
    scaLegendInset <- -0.35
  } else if (length(lsGeneIDs) <= 6){
    par(mfrow=c(2,3), xpd=TRUE)
    scaLegendInset <- -0.3
  } else {
    par(mfrow=c(3,3), xpd=TRUE)
    scaLegendInset <- -0.65
  }
  # Time points for plotting of impulse model
  vecX <- seq(0,max(vecTimepoints_Comb),0.1)
  # Find elements in vecX corresponding (closest) to observed time point
  indVecXObs <- unlist(lapply( 
    vecTimepoints_Comb, function(t){match( min(abs(vecX-t)), abs(vecX-t) )} 
  ))
  # Identify columns containing time course mean
  vecindMuByTimecourse <- grep("muByTimecourse",colnames(lsImpulseFits$parameters_case))
  for (geneID in lsGeneIDs){
    if(is.null(strControlName)){
      # Without control  data 
      
      # Convert h0,h1,h2 to log space again
      lsImpulseParamCaseLog <- lsImpulseFits$parameters_case[geneID,1:NPARAM]
      lsImpulseParamCaseLog[c("h0","h1","h2")] <- log( lsImpulseParamCaseLog[c("h0","h1","h2")] )
      vecCaseValues <- calcImpulse_comp(lsImpulseParamCaseLog,vecX)
      
      pval_Impulse <- round( log(dfImpulseResults[geneID,]$adj.p)/log(10), 2 )
      pval_Method2 <- round( log(vecMethod2Results[geneID])/log(10), 2 )
      
      if(strMode=="batch" | strMode=="singlecell"){
        # Plot observed points in blue - all time courses in same colour
        scaYlim_lower <- min( min(arr3DCountData[geneID,,],na.rm=TRUE), min(vecCaseValues[indVecXObs]) )
        scaYlim_upper <- max( max(arr3DCountData[geneID,,],na.rm=TRUE), max(vecCaseValues[indVecXObs]) )
        plot(arrTimepoints_Comb,t(arr3DCountData[geneID,,]),col="blue",pch=3,
          xlim=c(0,max(vecTimepoints_Comb,na.rm=TRUE)+PDF_WIDTH), ylim=c(scaYlim_lower,scaYlim_upper),
          xlab="Time", ylab="Impulse fit and expression values",
          main=paste0(geneID," ",strPlotTitleSuffix," log(Pval):\n ",strNameMethod1," ",pval_Impulse,
            " ",strNameMethod2," ",pval_Method2),sub=strPlotSubtitle)
        # Plot impulse within boundaries of observed points
        vecCaseValuesToPlot <- vecCaseValues
        indImpulseValToPlot <- vecCaseValuesToPlot >= scaYlim_lower & vecCaseValuesToPlot <= scaYlim_upper
        vecCaseValuesToPlot[!indImpulseValToPlot] <- NA
        points(vecX, vecCaseValuesToPlot,col="blue", type="l")
        
        # Plot inferred negative binomial pdf at each time point in black (vertical)
        vecXCoordPDF <- seq(round(scaYlim_lower),round(scaYlim_upper), by=1 )
        # Get means of negative binomial at each time point
        vecCaseValueAtTP <- calcImpulse_comp(lsImpulseParamCaseLog,arrTimepoints_Comb[1,])
        for(tp in arrTimepoints_Comb[1,]){
          vecYCoordPDF <- dnbinom(vecXCoordPDF,mu=vecCaseValueAtTP[match(tp,arrTimepoints_Comb[1,])],
            size=as.numeric(as.vector(dfImpulseResults[geneID,]$size)) )
          # Scale Y_coord to uniform peak heights of 1
          # This translates into width of one time unit in plot
          vecYCoordPDF <- vecYCoordPDF * PDF_WIDTH/max(vecYCoordPDF)
          # Plot pdf vertically at time point
          lines(x=tp+vecYCoordPDF,y=vecXCoordPDF,col="black")
        }
      } else if(strMode=="timecourses"){
        vecTranslationFactors <- lsImpulseFits$parameters_case[geneID,vecindMuByTimecourse]/
          lsImpulseFits$parameters_case[geneID,"mu"]
        # Create colour vector
        vecCol <- rainbow(n=dim(arr3DCountData)[3])
        
        # Create plot and plot data of first time course
        scaYlim_lower <- min( min(arr3DCountData[geneID,,],na.rm=TRUE), min(vecCaseValues[indVecXObs]*min(vecTranslationFactors)) )
        scaYlim_upper <- max( max(arr3DCountData[geneID,,],na.rm=TRUE), max(vecCaseValues[indVecXObs]*max(vecTranslationFactors)) )
        plot(arrTimepoints_Comb[1,],t(arr3DCountData[geneID,,1]),col=vecCol[1],pch=3,
          xlim=c(0,max(vecTimepoints_Comb,na.rm=TRUE)+PDF_WIDTH), ylim=c(scaYlim_lower,scaYlim_upper),
          xlab="Time", ylab="Impulse fit and expression values",
          main=paste0(geneID," ",strPlotTitleSuffix," log(Pval):\n ",strNameMethod1," ",pval_Impulse,
            " ",strNameMethod2," ",pval_Method2),sub=strPlotSubtitle)
        # Plot impulse fit to time course
        vecCaseValuesToPlotTC <- vecCaseValues*vecTranslationFactors[1]
        indImpulseValToPlot <- vecCaseValuesToPlotTC >= scaYlim_lower & vecCaseValuesToPlotTC <= scaYlim_upper
        vecCaseValuesToPlotTC[!indImpulseValToPlot] <- NA
        points(vecX, vecCaseValuesToPlotTC,col=vecCol[1], type="l")
        
        # Plot remaining time courses
        if(dim(arr3DCountData)[3]>1){
          for(timecourse in 2:dim(arr3DCountData)[3]){
            # Plot data of time course
            points(x=arrTimepoints_Comb[1,],y=t(arr3DCountData[geneID,,timecourse]),col=vecCol[timecourse],pch=3)
            # Plot impulse within boundaries of observed points
            vecCaseValuesToPlotTC <- vecCaseValues*vecTranslationFactors[timecourse]
            indImpulseValToPlot <- vecCaseValuesToPlotTC >= scaYlim_lower & vecCaseValuesToPlotTC <= scaYlim_upper
            vecCaseValuesToPlotTC[!indImpulseValToPlot] <- NA
            points(vecX, vecCaseValuesToPlotTC,col=vecCol[timecourse], type="l")
          }
        }
      } else {
        stop(paste0("ERROR: Unrecognised strMode in plotDEGenes(): ",strMode))
      }
      
      legend(x="bottomright",c(strCaseName),fill=c("blue"), cex=0.6, inset=c(0,scaLegendInset))
      # Plot mean of each time point
      points(arrTimepoints_Comb[1,],(apply(arr3DCountData[geneID,,],1,function(x){mean(x,na.rm=TRUE)})),col="black",pch=1)
      
    } else {
      # With control data
      
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
      pval_Method2 <- round( log(vecMethod2Results[geneID])/log(10), 2 )
      
      if(strMode=="batch" | strMode=="singlecell"){
        # Plot observed points: Case green, control red
        scaYlim_lower <- min( min(arr3DCountData[geneID,,],na.rm=TRUE),
          min(vecCaseValues[indVecXObs]),
          min(lsCtrlValues[indVecXObs]),
          min(lsCombValues[indVecXObs]) )
        scaYlim_upper <- max( max(arr3DCountData[geneID,,],na.rm=TRUE),
          max(vecCaseValues[indVecXObs]),
          max(lsCtrlValues[indVecXObs]),
          max(lsCombValues[indVecXObs]) )
        # case
        plot(arrTimepoints_Case,t(arr3DCountData_Case[geneID,,]),
          col="green",pch=3,
          xlim=c(0,max(vecTimepoints_Comb,na.rm=TRUE)+PDF_WIDTH),
          ylim=c(scaYlim_lower,scaYlim_upper),
          xlab="Time", ylab="Impulse fit and expression values",
          main=paste0(geneID," ",strPlotTitleSuffix," log(Pval):\n ",
            strNameMethod1," ",pval_Impulse," ",
            strNameMethod2," ",pval_Method2),sub=strPlotSubtitle)
        points(arrTimepoints_Ctrl,t(arr3DCountData_Ctrl[geneID,,]),
          col="red",pch=3, type="p")
        # Plot impulses within boundaries of observed points: case, control, combined
        # case
        vecCaseValuesToPlot <- vecCaseValues
        indImpulseValToPlot <- vecCaseValuesToPlot >= scaYlim_lower & vecCaseValuesToPlot <= scaYlim_upper
        vecCaseValuesToPlot[!indImpulseValToPlot] <- NA
        points(vecX, vecCaseValuesToPlot,col="green", type="l")
        # control
        lsCtrlValuesToPlot <- lsCtrlValues
        indImpulseValToPlot <- lsCtrlValuesToPlot >= scaYlim_lower & lsCtrlValuesToPlot <= scaYlim_upper
        lsCtrlValuesToPlot[!indImpulseValToPlot] <- NA
        points(vecX, lsCtrlValuesToPlot,col="red", type="l")
        # combined
        lsCombValuesToPlot <- lsCombValues
        indImpulseValToPlot <- lsCombValuesToPlot >= scaYlim_lower & lsCombValuesToPlot <= scaYlim_upper
        lsCombValuesToPlot[!indImpulseValToPlot] <- NA
        points(vecX, lsCombValuesToPlot,col="black", type="l")
        
      } else if(strMode=="timecourses"){
        vecTranslationFactors <- lsImpulseFits$parameters_case[geneID,vecindMuByTimecourse]/lsImpulseFits$parameters_case[geneID,"mu"]
        # Create colour vector
        vecCol <- rainbow(n=dim(arr3DCountData)[3])
        
        # Create plot and plot data of first time course
        scaYlim_lower <- min( min(arr3DCountData[geneID,,],na.rm=TRUE), min(vecCaseValues[indVecXObs]*min(vecTranslationFactors)) )
        scaYlim_upper <- max( max(arr3DCountData[geneID,,],na.rm=TRUE), max(vecCaseValues[indVecXObs]*max(vecTranslationFactors)) )
        plot(arrTimepoints_Comb[1,],t(arr3DCountData[geneID,,1]),col=vecCol[1],pch=3,
          xlim=c(0,max(vecTimepoints_Comb,na.rm=TRUE)+PDF_WIDTH), ylim=c(scaYlim_lower,scaYlim_upper),
          xlab="Time", ylab="Impulse fit and expression values",
          main=paste0(geneID," ",strPlotTitleSuffix," log(Pval):\n ",strNameMethod1," ",pval_Impulse,
            " ",strNameMethod2," ",pval_Method2),sub=strPlotSubtitle)
        # Plot impulse fit to time course
        vecCaseValuesToPlotTC <- vecCaseValues*vecTranslationFactors[1]
        indImpulseValToPlot <- vecCaseValuesToPlotTC >= scaYlim_lower & vecCaseValuesToPlotTC <= scaYlim_upper
        vecCaseValuesToPlotTC[!indImpulseValToPlot] <- NA
        points(vecX, vecCaseValuesToPlotTC,col=vecCol[1], type="l")
        
      } else {
        stop(paste0("ERROR: Unrecognised strMode in plotDEGenes(): ",strMode))
      }
      
      legend(x="bottomright",c(strCaseName,strControlName),fill=c("green","red"), cex=0.6, inset=c(0,scaLegendInset))
      # Plot mean of each time point
      
      points(arrTimepoints_Comb[1,],sapply(arrTimepoints_Comb[1,],function(tp){
        mean(arr3DCountData[geneID, as.vector(dfAnnotationRed[dfAnnotationRed$Time %in% tp,]$Sample),],na.rm=TRUE)}),
        col="black",pch=1)
      
    }
  }
  # Close .pdf
  dev.off()
}