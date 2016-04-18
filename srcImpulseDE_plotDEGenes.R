#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Plot impulse fits    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# Plots the impulse fits to timecourse data and also the control data if present

# INPUT:
#   gene_IDs: (character vector) of gene names to be plotted; must be part
#       of the rownames of data_array.
#   data_array: (Numeric 3D array genes x samples x replicates)
#       Contains expression values or similar locus-specific read-outs.
#   data_annotation: (Table samples x 2[time and condition]) 
#       Co-variables for the samples including condition and time points.
#       Time points must be numeric.
#   imp_fit_genes: (list runs*2 ["impulse_parameters","impulse_fits"])
#       impulse_parameters: (matrix genes x (NPARAM+1)*runs) +1 is value of 
#         objective of optimisation (i.e. sum of squares or weighted SS)
#       impulse_fits: (matrix genes x timepoints) model values for gene data
#   background: (vector number of F-scores simulated under H0)  Empirical 
#       F-values simulated under H0.
#   control_timecourse: (bool) [Default FALSE] Control time timecourse is 
#       part of the data set (TRUE) or not (FALSE).
#   control_name: (str) Name of the control condition in annotation_table.
#   case_name (str) Name of the control condition in annotation_table.
#   file_name_part: (character string) File extention.
#   title_string: (character string) Title for each plot.
#   sub_line: (character string) Subtitle for each plot.
# OUTPUT:
#   -

plotDEGenes <- function(gene_IDs, data_array, data_annotation,imp_fit_genes,
  control_timecourse = FALSE, control_name = NULL, case_name = NULL,
  file_name_part = "", title_line = "", sub_line = "",
  ImpulseDE_res = NULL, DESeq2_res=NULL){
  #pvals_impulse_deseq=NULL){
  # DAVID DEVELOPMENTAL NOTE : remove pval_deseq later again?
  
  NPARAM=6
  
  if(length(grep("[a-zA-Z]",rownames(data_array))) == 0 ){
    rownames(data_array) <- paste(rownames(data_array),"G", sep = "_")
    gene_IDs <- paste(gene_IDs, "G", sep = "_")
  }
  
  ### if control timecourse is present split data into case and control data
  if(control_timecourse == TRUE){
    farr_case <- data_array[,!(data_annotation$Condition %in% control_name),]
    farr_ctrl <- data_array[,data_annotation$Condition %in% control_name,]
    timep_case <-  as.numeric(as.character(data_annotation[colnames(farr_case),"Time"]))
    timep_ctrl <- as.numeric(as.character(data_annotation[colnames(farr_ctrl),"Time"]))
  }
  timep <- as.numeric(as.character(data_annotation[colnames(data_array),"Time"]))
  # expand timep vector to size of expression array
  timep_arr <- t(matrix(rep(timep,dim(data_array)[3]),length(timep),dim(data_array)[3]))
  if(control_timecourse == TRUE){
    timep_case_arr <- t(matrix(rep(timep_case,dim(data_array)[3]),length(timep_case),dim(data_array)[3]))
    timep_ctrl_arr <- t(matrix(rep(timep_ctrl,dim(data_array)[3]),length(timep_ctrl),dim(data_array)[3]))
  }
  
  pdf(paste("impulse_fit_genes_",file_name_part,".pdf",sep=""),height=6.0,width=9.0)
  if (length(gene_IDs) == 1){
    par(mfrow=c(1,1))
  } else if (length(gene_IDs) <= 4){
    par(mfrow=c(2,2))
  } else if (length(gene_IDs) <= 6){
    par(mfrow=c(2,3))
  } else {
    par(mfrow=c(3,3))
  }
  x_vec <- seq(0,max(timep),0.1)
  for (gene_ID in gene_IDs){
    # If there is no control data plot only case time course data and fit
    if(control_timecourse == FALSE){
      # Chose impulse fit if parameters of fitted model are not NAN,
      # only plot first timepoint otherwise
      if(TRUE %in% is.na(imp_fit_genes$impulse_parameters_case[gene_ID,])){
        calc_case <- imp_fit_genes$impulse_fits_case[gene_ID,1]
      } else {
        calc_case <- calc_impulse_comp(imp_fit_genes$impulse_parameters_case[gene_ID,1:NPARAM],x_vec)
      }
      pval_DEseq <- round( log(DESeq2_res[gene_ID,]$padj)/log(10), 2 )
      pval_Impulse <- round( log(ImpulseDE_res[gene_ID,]$adj.p)/log(10), 2 )
      #pval_DEseq <- round( log(pvals_impulse_deseq[gene_ID,"DESeq"]))/log(10), 2 )
      #pval_Impulse <- round( log(pvals_impulse_deseq[gene_ID,"Impulse"])/log(10), 2 )
      plot(timep_arr,(t(data_array[gene_ID,,])),col="blue",pch=3,xlim=c(0,max(timep)),
        ylim=c((min(c(as.numeric(data_array[gene_ID,,]),as.numeric(calc_case)))-0.5),
          (max(c(as.numeric(data_array[gene_ID,,]),as.numeric(calc_case)))+0.5)),
        xlab="Time", ylab="Impulse fit and expression values",
        main=paste0(gene_ID," ",title_line," log(Pval):\n DESeq2 ",pval_DEseq,
          " ImpulseDE2 ",pval_Impulse),sub=sub_line)
      points(timep_arr[1,],(apply(data_array[gene_ID,,],1,mean)),col="red",pch=1)
      if(TRUE %in% is.na(imp_fit_genes$impulse_parameters_case[gene_ID,])){
        abline(h = calc_case , col = "blue")
      } else {
        points(x_vec, calc_case, col = "blue", type="l")
      } 
      legend(x="bottomright",as.character(data_annotation[1,"Condition"]),fill=c("blue"), cex=0.6)
      
      ### if there is control data
    } else if(control_timecourse == TRUE){
      
      if(TRUE %in% is.na(imp_fit_genes$impulse_parameters_case[gene_ID,])){
        calc_case = imp_fit_genes$impulse_fits_case[gene_ID,1]
        status_case = FALSE
      } else {
        calc_case = calc_impulse_comp(imp_fit_genes$impulse_parameters_case[gene_ID,1:NPARAM],x_vec)
        status_case = TRUE
      }
      if(TRUE %in% is.na(imp_fit_genes$impulse_parameters_control[gene_ID,])){
        calc_ctrl = imp_fit_genes$impulse_fits_control[gene_ID,1]
        status_ctrl = FALSE
      } else {
        calc_ctrl = calc_impulse_comp(imp_fit_genes$impulse_parameters_control[gene_ID,1:NPARAM],x_vec)
        status_ctrl = TRUE
      }
      if(TRUE %in% is.na(imp_fit_genes$impulse_parameters_combined[gene_ID,])){
        calc_comb = imp_fit_genes$impulse_fits_combined[gene_ID,1]
        status_comb = FALSE
      } else {
        calc_comb = calc_impulse_comp(imp_fit_genes$impulse_parameters_combined[gene_ID,1:NPARAM],x_vec)
        status_comb = TRUE
      }
      
      plot(timep_case_arr,(t(farr_case[gene_ID,,])),col="blue",pch=3,xlim=c(0,max(timep)),
        ylim=c((min(c(as.numeric(data_array[gene_ID,,]),as.numeric(calc_case), as.numeric(calc_ctrl), as.numeric(calc_comb)))-0.5),
          (max(c(as.numeric(data_array[gene_ID,,]),as.numeric(calc_case), as.numeric(calc_ctrl), as.numeric(calc_comb)))+0.5)),
        xlab="Time", ylab="Impulse fit or log2 expression value",
        main=paste(gene_ID," ",title_line,sep=""),sub=sub_line)
      
      points(timep_ctrl,farr_ctrl[gene_ID,],col="red",pch=4)
      
      if(status_case == FALSE){
        abline(h = calc_case , col = "blue")
      } else {
        points(x_vec,calc_case, col = "blue", type="l")
      }
      if(status_ctrl == FALSE){
        abline(h = calc_ctrl , col = "red")
      } else {
        points(x_vec,calc_ctrl, col = "red", type="l")
      }
      if(status_comb == FALSE){
        abline(h = calc_comb , col = "grey")
      } else {
        points(x_vec,calc_comb, col = "grey", type="l")
      }
      legend(x="bottomright",c(as.character(data_annotation$Condition[data_annotation$Condition != control_name][1]),control_name,"combined"),fill=c("blue","red","grey"), cex=0.6)
    }
  }
  dev.off()
}