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
#' @param dfAnnotationRed (data frame) Reduced version of 
#'    \code{dfAnnotationFull}. Lists co-variables of samples: 
#'    Sample, Condition, Time. Time must be numeric. For internal use.
#' @param lsImpulseFits (list length 2 or 6) List of matrices which
#'    contain parameter fits and model values for given time course for the
#'    case condition (and control and combined if control is present).
#'    Each parameter matrix is called parameter_'condition' and has the form
#'    (genes x [beta, h0, h1, h2, t1, t2, logL_H1, converge_H1, mu, logL_H0, 
#'    converge_H0]) where beta to t2 are parameters of the impulse
#'    model, mu is the single parameter of the mean model, logL are
#'    log likelihoods of full (H1) and reduced model (H0) respectively, converge
#'    is convergence status of numerical optimisation of model fitting by
#'    \code{optim} from \code{stats} of either model. Each value matrix is called
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

plotDEGenes <- function(lsGeneIDs, arr3DCountData, dfAnnotationRed,
  lsImpulseFits, dfImpulseResults, vecMethod2Results, 
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
  
  # Name genes if not previously named (already done if this function is
  # called within the wrapper function).
  if(length(grep("[a-zA-Z]",rownames(arr3DCountData))) == 0 ){
    rownames(arr3DCountData) <- paste(rownames(arr3DCountData),"G", sep = "_")
    lsGeneIDs <- paste(lsGeneIDs, "G", sep = "_")
  }
  
  if(!is.null(strControlName)){
    arr3DCountData_Case <- arr3DCountData[,dfAnnotationRed$Condition %in% strCaseName,]
    arr3DCountData_Ctrl <- arr3DCountData[,dfAnnotationRed$Condition %in% strControlName,]
    lsTimepoints_Case <-  as.numeric(as.character(dfAnnotationRed[colnames(arr3DCountData_Case),"Time"]))
    lsTimepoints_Ctrl <- as.numeric(as.character(dfAnnotationRed[colnames(arr3DCountData_Ctrl),"Time"]))
  }
  lsTimepoints_All <- as.numeric(as.character(dfAnnotationRed[dfAnnotationRed$Sample %in% colnames(arr3DCountData),"Time"]))
  
  # Expand lsTimepoints_All vector to size of expression array
  if(!is.null(strControlName)){
    arrTimepoints_Case <- t(matrix(rep(lsTimepoints_Case,dim(arr3DCountData)[3]),length(lsTimepoints_Case),dim(arr3DCountData)[3]))
    arrTimepoints_Ctrl <- t(matrix(rep(lsTimepoints_Ctrl,dim(arr3DCountData)[3]),length(lsTimepoints_Ctrl),dim(arr3DCountData)[3]))
  }
  arrTimepoints_All <- t(matrix(rep(lsTimepoints_All,dim(arr3DCountData)[3]),length(lsTimepoints_All),dim(arr3DCountData)[3]))
  
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
    scaLegendInset <- -0.6
  }
  # Time points for plotting of impulse model
  vecX <- seq(0,max(lsTimepoints_All),0.1)
  # Find elements in vecX corresponding (closest) to observed time point
  indVecXObs <- unlist(lapply( 
    lsTimepoints_All, function(t){match( min(abs(vecX-t)), abs(vecX-t) )} 
  ))
  # Identify columns containing time course mean
  vecindMuByTimecourse <- grep("muByTimecourse",colnames(lsImpulseFits$parameters_case))
  for (geneID in lsGeneIDs){
    # Without control  data 
    if(is.null(strControlName)){
      # Convert h0,h1,h2 to log space again
      lsImpulseParamCaseLog <- lsImpulseFits$parameters_case[geneID,1:NPARAM]
      lsImpulseParamCaseLog[c("h0","h1","h2")] <- log( lsImpulseParamCaseLog[c("h0","h1","h2")] )
      lsCaseValues <- calcImpulse_comp(lsImpulseParamCaseLog,vecX)
      
      pval_Impulse <- round( log(dfImpulseResults[geneID,]$adj.p)/log(10), 2 )
      pval_Method2 <- round( log(vecMethod2Results[geneID])/log(10), 2 )
      
      if(strMode=="batch" | strMode=="singlecell"){
        # Plot observed points in blue - all time courses in same colour
        scaYlim_lower <- min( min(arr3DCountData[geneID,,],na.rm=TRUE), min(lsCaseValues[indVecXObs]) )
        scaYlim_upper <- max( max(arr3DCountData[geneID,,],na.rm=TRUE), max(lsCaseValues[indVecXObs]) )
        plot(arrTimepoints_All,t(arr3DCountData[geneID,,]),col="blue",pch=3,
          xlim=c(0,max(lsTimepoints_All,na.rm=TRUE)+PDF_WIDTH), ylim=c(scaYlim_lower,scaYlim_upper),
          xlab="Time", ylab="Impulse fit and expression values",
          main=paste0(geneID," ",strPlotTitleSuffix," log(Pval):\n ",strNameMethod1," ",pval_Impulse,
            " ",strNameMethod2," ",pval_Method2),sub=strPlotSubtitle)
        # Plot impulse within boundaries of observed points
        lsCaseValuesToPlot <- lsCaseValues
        indImpulseValToPlot <- lsCaseValuesToPlot >= scaYlim_lower & lsCaseValuesToPlot <= scaYlim_upper
        lsCaseValuesToPlot[!indImpulseValToPlot] <- NA
        points(vecX, lsCaseValuesToPlot,col="blue", type="l")
        
        # Plot inferred negative binomial pdf at each time point in black (vertical)
        vecXCoordPDF <- seq(round(scaYlim_lower),round(scaYlim_upper), by=1 )
        # Get means of negative binomial at each time point
        vecCaseValueAtTP <- calcImpulse_comp(lsImpulseParamCaseLog,arrTimepoints_All[1,])
        for(tp in arrTimepoints_All[1,]){
          vecYCoordPDF <- dnbinom(vecXCoordPDF,mu=vecCaseValueAtTP[match(tp,arrTimepoints_All[1,])],
            size=as.numeric(as.vector(dfImpulseResults[geneID,]$size)) )
          # Scale Y_coord to uniform peak heights of 1
          # This translates into width of one time unit in plot
          vecYCoordPDF <- vecYCoordPDF * PDF_WIDTH/max(vecYCoordPDF)
          # Plot pdf vertically at time point
          lines(x=tp+vecYCoordPDF,y=vecXCoordPDF,col="black")
        }
      } else if(strMode=="timecourses"){
        vecTranslationFactors <- lsImpulseFits$parameters_case[geneID,vecindMuByTimecourse]/lsImpulseFits$parameters_case[geneID,"mu"]
        # Create colour vector
        vecCol <- rainbow(n=dim(arr3DCountData)[3])
        
        # Create plot and plot data of first time course
        scaYlim_lower <- min( min(arr3DCountData[geneID,,],na.rm=TRUE), min(lsCaseValues[indVecXObs]*min(vecTranslationFactors)) )
        scaYlim_upper <- max( max(arr3DCountData[geneID,,],na.rm=TRUE), max(lsCaseValues[indVecXObs]*max(vecTranslationFactors)) )
        plot(arrTimepoints_All[1,],t(arr3DCountData[geneID,,1]),col=vecCol[1],pch=3,
          xlim=c(0,max(lsTimepoints_All,na.rm=TRUE)+PDF_WIDTH), ylim=c(scaYlim_lower,scaYlim_upper),
          xlab="Time", ylab="Impulse fit and expression values",
          main=paste0(geneID," ",strPlotTitleSuffix," log(Pval):\n ",strNameMethod1," ",pval_Impulse,
            " ",strNameMethod2," ",pval_Method2),sub=strPlotSubtitle)
        legend(x="bottomright",as.character(dfAnnotationRed[1,"Condition"]),fill=c("blue"), cex=0.6, inset=c(0,scaLegendInset))
        # Plot impulse fit to time course
        lsCaseValuesToPlotTC <- lsCaseValues*vecTranslationFactors[1]
        indImpulseValToPlot <- lsCaseValuesToPlotTC >= scaYlim_lower & lsCaseValuesToPlotTC <= scaYlim_upper
        lsCaseValuesToPlotTC[!indImpulseValToPlot] <- NA
        points(vecX, lsCaseValuesToPlotTC,col=vecCol[1], type="l")
        
        # Plot remaining time courses
        if(dim(arr3DCountData)[3]>1){
          for(timecourse in 2:dim(arr3DCountData)[3]){
            # Plot data of time course
            points(x=arrTimepoints_All[1,],y=t(arr3DCountData[geneID,,timecourse]),col=vecCol[timecourse],pch=3)
            # Plot impulse within boundaries of observed points
            lsCaseValuesToPlotTC <- lsCaseValues*vecTranslationFactors[timecourse]
            indImpulseValToPlot <- lsCaseValuesToPlotTC >= scaYlim_lower & lsCaseValuesToPlotTC <= scaYlim_upper
            lsCaseValuesToPlotTC[!indImpulseValToPlot] <- NA
            points(vecX, lsCaseValuesToPlotTC,col=vecCol[timecourse], type="l")
          }
        }
      } else {
        stop(paste0("ERROR: Unrecognised strMode in plotDEGenes(): ",strMode))
      }
      
      # Plot mean of each time point
      points(arrTimepoints_All[1,],(apply(arr3DCountData[geneID,,],1,function(x){mean(x,na.rm=TRUE)})),col="black",pch=1)
      
      # With control data
    } else {
      
      if(TRUE %in% is.na(lsImpulseFits$parameters_case[geneID,])){
        lsCaseValues = lsImpulseFits$values_case[geneID,1]
        boolStatusCase = FALSE
      } else {
        lsImpulseParamCaseLog <- lsImpulseFits$parameters_case[geneID,1:NPARAM]
        lsImpulseParamCaseLog[c("h0","h1","h2")] <- log( lsImpulseParamCaseLog[c("h0","h1","h2")] )
        lsCaseValues = calcImpulse_comp(lsImpulseParamCaseLog,vecX)
        boolStatusCase = TRUE
      }
      if(TRUE %in% is.na(lsImpulseFits$parameters_control[geneID,])){
        lsCtrlValues = lsImpulseFits$values_control[geneID,1]
        boolStatusCtrl = FALSE
      } else {
        lsImpulseParamCtrlLog <- lsImpulseFits$parameters_case[geneID,1:NPARAM]
        lsImpulseParamCtrlLog[c("h0","h1","h2")] <- log( lsImpulseParamCtrlLog[c("h0","h1","h2")] )
        lsCtrlValues = calcImpulse_comp(lsImpulseParamCtrlLog,vecX)
        boolStatusCtrl = TRUE
      }
      if(TRUE %in% is.na(lsImpulseFits$parameters_combined[geneID,])){
        lsCombValues = lsImpulseFits$values_combined[geneID,1]
        boolStatusComb = FALSE
      } else {
        lsImpulseParamCombLog <- lsImpulseFits$parameters_case[geneID,1:NPARAM]
        lsImpulseParamCombLog[c("h0","h1","h2")] <- log( lsImpulseParamCombLog[c("h0","h1","h2")] )
        lsCombValues = calcImpulse_comp(lsImpulseParamCombLog,vecX)
        boolStatusComb = TRUE
      }
      
      plot(arrTimepoints_Case,(t(arr3DCountData_Case[geneID,,])),col="blue",pch=3,xlim=c(0,max(lsTimepoints_All)),
        ylim=c(min(c(as.numeric(arr3DCountData[geneID,,]),as.numeric(lsCaseValues), as.numeric(lsCtrlValues), as.numeric(lsCombValues)),na.rm=TRUE),
          max(c(as.numeric(arr3DCountData[geneID,,]),as.numeric(lsCaseValues), as.numeric(lsCtrlValues), as.numeric(lsCombValues)),na.rm=TRUE)),
        xlab="Time", ylab="Impulse fit and expression values",
        main=paste0(geneID," ",strPlotTitleSuffix," log(Pval):\n ",strNameMethod1," ",pval_Impulse,
          " ",strNameMethod2," ",pval_Method2),sub=strPlotSubtitle)
      
      points(lsTimepoints_Ctrl,arr3DCountData_Ctrl[geneID,],col="red",pch=4)
      
      if(boolStatusCase == FALSE){
        abline(h = lsCaseValues , col = "blue")
      } else {
        points(vecX,lsCaseValues, col = "blue", type="l")
      }
      if(boolStatusCtrl == FALSE){
        abline(h = lsCtrlValues , col = "red")
      } else {
        points(vecX,lsCtrlValues, col = "red", type="l")
      }
      if(boolStatusComb == FALSE){
        abline(h = lsCombValues , col = "grey")
      } else {
        points(vecX,lsCombValues, col = "grey", type="l")
      }
      legend(x="bottomright",c(as.character(dfAnnotationRed$Condition[dfAnnotationRed$Condition != control_name][1]),control_name,"combined"),fill=c("blue","red","grey"), cex=0.6)
    }
  }
  # Close .pdf
  dev.off()
}