#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Plot impulse fits    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# Plots the impulse fits to timecourse data and also the control data if present

# INPUT:
#   lsGeneIDs: (character vector) of gene names to be plotted; must be part
#       of the rownames of arr3DCountData.
#   arr3DCountData: (Numeric 3D array genes x samples x replicates)
#       Contains expression values or similar locus-specific read-outs.
#   dfAnnotationRed: (Table samples x 2[time and condition]) 
#       Co-variables for the samples including condition and time points.
#       Time points must be numeric.
#   lsImpulseFits: (list runs*2 ["parameters","impulse_fits"])
#       parameters: (matrix genes x (NPARAM+1)*runs) +1 is value of 
#         objective of optimisation (i.e. sum of squares or weighted SS)
#       impulse_fits: (matrix genes x timepoints) model values for gene data
#   strControlName: (str) Name of the control condition in annotation_table.
#   strCaseName (str) Name of the control condition in annotation_table.
#   strFileNameSuffix: (character string) File extention.
#   title_string: (character string) Title for each plot.
#   strPlotSubtitle: (character string) Subtitle for each plot.
# OUTPUT:
#   -

plotDEGenes <- function(lsGeneIDs, arr3DCountData, dfAnnotationRed,
  lsImpulseFits, dfImpulseResults = NULL, dfDESeq2Results=NULL, 
  strCaseName=NULL, strControlName=NULL,
  strFileNameSuffix = "", strPlotTitleSuffix = "", strPlotSubtitle = "",
  NPARAM=6){

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
  
  # expand lsTimepoints_All vector to size of expression array
  if(!is.null(strControlName)){
    arrTimepoints_Case <- t(matrix(rep(lsTimepoints_Case,dim(arr3DCountData)[3]),length(lsTimepoints_Case),dim(arr3DCountData)[3]))
    arrTimepoints_Ctrl <- t(matrix(rep(lsTimepoints_Ctrl,dim(arr3DCountData)[3]),length(lsTimepoints_Ctrl),dim(arr3DCountData)[3]))
  }
  arrTimepoints_All <- t(matrix(rep(lsTimepoints_All,dim(arr3DCountData)[3]),length(lsTimepoints_All),dim(arr3DCountData)[3]))
  
  # Open .pdf
  pdf(paste("ImpulseDE_",strFileNameSuffix,".pdf",sep=""),height=6.0,width=9.0)
  
  # Define grid for printing plots
  if (length(lsGeneIDs) == 1){
    par(mfrow=c(1,1))
  } else if (length(lsGeneIDs) <= 4){
    par(mfrow=c(2,2))
  } else if (length(lsGeneIDs) <= 6){
    par(mfrow=c(2,3))
  } else {
    par(mfrow=c(3,3))
  }
  # Time points for plotting of impulse model
  vecX <- seq(0,max(lsTimepoints_All),0.1)
  for (geneID in lsGeneIDs){
    # Without control  data 
    if(is.null(strControlName)){
      # Chose impulse fit if parameters of fitted model are not NAN,
      # only plot first timepoint otherwise
      if(TRUE %in% is.na(lsImpulseFits$parameters_case[geneID,])){
        lsCaseValues <- lsImpulseFits$values_case[geneID,1]
      } else {
        lsCaseValues <- calcImpulse_comp(lsImpulseFits$parameters_case[geneID,1:NPARAM],vecX)
      }
      pval_DEseq <- round( log(dfDESeq2Results[geneID,]$padj)/log(10), 2 )
      pval_Impulse <- round( log(dfImpulseResults[geneID,]$adj.p)/log(10), 2 )

      
      plot(arrTimepoints_All,(t(arr3DCountData[geneID,,])),col="blue",pch=3,xlim=c(0,max(lsTimepoints_All,na.rm=TRUE)),
        ylim=c(min(c(as.numeric(arr3DCountData[geneID,,]),as.numeric(lsCaseValues)),na.rm=TRUE),
          max(c(as.numeric(arr3DCountData[geneID,,]),as.numeric(lsCaseValues)),na.rm=TRUE)),
        xlab="Time", ylab="Impulse fit and expression values",
        main=paste0(geneID," ",strPlotTitleSuffix," log(Pval):\n DESeq2 ",pval_DEseq,
          " ImpulseDE2 ",pval_Impulse),sub=strPlotSubtitle)
      
      points(arrTimepoints_All[1,],(apply(arr3DCountData[geneID,,],1,function(x){mean(x,na.rm=TRUE)})),col="red",pch=1)
      
      if(TRUE %in% is.na(lsImpulseFits$parameters_case[geneID,])){
        abline(h = lsCaseValues , col = "blue")
      } else {
        points(vecX, lsCaseValues, col = "blue", type="l")
      } 
      legend(x="bottomright",as.character(dfAnnotationRed[1,"Condition"]),fill=c("blue"), cex=0.6)
      
    # With control data
    } else {
      
      if(TRUE %in% is.na(lsImpulseFits$parameters_case[geneID,])){
        lsCaseValues = lsImpulseFits$values_case[geneID,1]
        boolStatusCase = FALSE
      } else {
        lsCaseValues = calcImpulse_comp(lsImpulseFits$parameters_case[geneID,1:NPARAM],vecX)
        boolStatusCase = TRUE
      }
      if(TRUE %in% is.na(lsImpulseFits$parameters_control[geneID,])){
        lsCtrlValues = lsImpulseFits$values_control[geneID,1]
        boolStatusCtrl = FALSE
      } else {
        lsCtrlValues = calcImpulse_comp(lsImpulseFits$parameters_control[geneID,1:NPARAM],vecX)
        boolStatusCtrl = TRUE
      }
      if(TRUE %in% is.na(lsImpulseFits$parameters_combined[geneID,])){
        lsCombValues = lsImpulseFits$values_combined[geneID,1]
        boolStatusComb = FALSE
      } else {
        lsCombValues = calcImpulse_comp(lsImpulseFits$parameters_combined[geneID,1:NPARAM],vecX)
        boolStatusComb = TRUE
      }
      
      plot(arrTimepoints_Case,(t(arr3DCountData_Case[geneID,,])),col="blue",pch=3,xlim=c(0,max(lsTimepoints_All)),
        ylim=c(min(c(as.numeric(arr3DCountData[geneID,,]),as.numeric(lsCaseValues), as.numeric(lsCtrlValues), as.numeric(lsCombValues)),na.rm=TRUE),
          max(c(as.numeric(arr3DCountData[geneID,,]),as.numeric(lsCaseValues), as.numeric(lsCtrlValues), as.numeric(lsCombValues)),na.rm=TRUE)),
        xlab="Time", ylab="Impulse fit and expression values",
        main=paste0(geneID," ",strPlotTitleSuffix," log(Pval):\n DESeq2 ",pval_DEseq,
          " ImpulseDE2 ",pval_Impulse),sub=strPlotSubtitle)
      
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