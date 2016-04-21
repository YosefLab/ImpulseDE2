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
  lsImpulseFits, ImpulseDE_res = NULL, DESeq2_res=NULL, 
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
    arr3DCountData_Case <- arr3DCountData[,data_annotation$Condition %in% strCaseName,]
    arr3DCountData_Ctrl <- arr3DCountData[,data_annotation$Condition %in% strControlName,]
    lsTimepoints_Case <-  as.numeric(as.character(data_annotation[colnames(arr3DCountData_Case),"Time"]))
    lsTimepoints_Ctrl <- as.numeric(as.character(data_annotation[colnames(arr3DCountData_Ctrl),"Time"]))
  }
  lsTimepoints_All <- as.numeric(as.character(data_annotation[colnames(arr3DCountData),"Time"]))
  
  # expand lsTimepoints_All vector to size of expression array
  if(!is.null(strControlName)){
    arrTimepoints_Case <- t(matrix(rep(lsTimepoints_Case,dim(arr3DCountData)[3]),length(lsTimepoints_Case),dim(arr3DCountData)[3]))
    arrTimepoints_Ctrl <- t(matrix(rep(lsTimepoints_Ctrl,dim(arr3DCountData)[3]),length(lsTimepoints_Ctrl),dim(arr3DCountData)[3]))
  }
  arrTimepoints_All <- t(matrix(rep(lsTimepoints_All,dim(arr3DCountData)[3]),length(lsTimepoints_All),dim(arr3DCountData)[3]))
  
  # Open .pdf
  pdf(paste("impulse_fit_genes_",strFileNameSuffix,".pdf",sep=""),height=6.0,width=9.0)
  
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
        calc_case <- lsImpulseFits$values_case[geneID,1]
      } else {
        calc_case <- calc_impulse_comp(lsImpulseFits$parameters_case[geneID,1:NPARAM],vecX)
      }
      pval_DEseq <- round( log(DESeq2_res[geneID,]$padj)/log(10), 2 )
      pval_Impulse <- round( log(ImpulseDE_res[geneID,]$adj.p)/log(10), 2 )

      plot(arrTimepoints_All,(t(arr3DCountData[geneID,,])),col="blue",pch=3,xlim=c(0,max(lsTimepoints_All)),
        ylim=c((min(c(as.numeric(arr3DCountData[geneID,,]),as.numeric(calc_case)))-0.5),
          (max(c(as.numeric(arr3DCountData[geneID,,]),as.numeric(calc_case)))+0.5)),
        xlab="Time", ylab="Impulse fit and expression values",
        main=paste0(geneID," ",strPlotTitleSuffix," log(Pval):\n DESeq2 ",pval_DEseq,
          " ImpulseDE2 ",pval_Impulse),sub=strPlotSubtitle)
      
      points(arrTimepoints_All[1,],(apply(arr3DCountData[geneID,,],1,mean)),col="red",pch=1)
      
      if(TRUE %in% is.na(lsImpulseFits$parameters_case[geneID,])){
        abline(h = calc_case , col = "blue")
      } else {
        points(vecX, calc_case, col = "blue", type="l")
      } 
      legend(x="bottomright",as.character(data_annotation[1,"Condition"]),fill=c("blue"), cex=0.6)
      
    # With control data
    } else {
      
      if(TRUE %in% is.na(lsImpulseFits$parameters_case[geneID,])){
        calc_case = lsImpulseFits$values_case[geneID,1]
        status_case = FALSE
      } else {
        calc_case = calc_impulse_comp(lsImpulseFits$parameters_case[geneID,1:NPARAM],vecX)
        status_case = TRUE
      }
      if(TRUE %in% is.na(lsImpulseFits$parameters_control[geneID,])){
        calc_ctrl = lsImpulseFits$values_control[geneID,1]
        status_ctrl = FALSE
      } else {
        calc_ctrl = calc_impulse_comp(lsImpulseFits$parameters_control[geneID,1:NPARAM],vecX)
        status_ctrl = TRUE
      }
      if(TRUE %in% is.na(lsImpulseFits$parameters_combined[geneID,])){
        calc_comb = lsImpulseFits$values_combined[geneID,1]
        status_comb = FALSE
      } else {
        calc_comb = calc_impulse_comp(lsImpulseFits$parameters_combined[geneID,1:NPARAM],vecX)
        status_comb = TRUE
      }
      
      plot(arrTimepoints_Case,(t(arr3DCountData_Case[geneID,,])),col="blue",pch=3,xlim=c(0,max(lsTimepoints_All)),
        ylim=c((min(c(as.numeric(arr3DCountData[geneID,,]),as.numeric(calc_case), as.numeric(calc_ctrl), as.numeric(calc_comb)))-0.5),
          (max(c(as.numeric(arr3DCountData[geneID,,]),as.numeric(calc_case), as.numeric(calc_ctrl), as.numeric(calc_comb)))+0.5)),
        xlab="Time", ylab="Impulse fit and expression values",
        main=paste0(geneID," ",strPlotTitleSuffix," log(Pval):\n DESeq2 ",pval_DEseq,
          " ImpulseDE2 ",pval_Impulse),sub=strPlotSubtitle)
      
      points(lsTimepoints_Ctrl,arr3DCountData_Ctrl[geneID,],col="red",pch=4)
      
      if(status_case == FALSE){
        abline(h = calc_case , col = "blue")
      } else {
        points(vecX,calc_case, col = "blue", type="l")
      }
      if(status_ctrl == FALSE){
        abline(h = calc_ctrl , col = "red")
      } else {
        points(vecX,calc_ctrl, col = "red", type="l")
      }
      if(status_comb == FALSE){
        abline(h = calc_comb , col = "grey")
      } else {
        points(vecX,calc_comb, col = "grey", type="l")
      }
      legend(x="bottomright",c(as.character(data_annotation$Condition[data_annotation$Condition != control_name][1]),control_name,"combined"),fill=c("blue","red","grey"), cex=0.6)
    }
  }
  # Close .pdf
  dev.off()
}