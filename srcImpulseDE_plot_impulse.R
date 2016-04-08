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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Plot impulse model fits
#'
#' Plots impulse model fits for the specified gene IDs. In the case of two
#' time courses, the fits for the combined, case and control data are plotted.
#' @aliases plot_impulse
#' @param gene_IDs character vector of gene names to be plotted; must be part
#' of the \code{rownames} of \code{data_array}
#' @param data_array numeric matrix of expression values; genes should be in
#' rows, samples in columns. Data should be properly normalized and
#' log2-transformed as well as filtered for present or variable genes.
#' @param data_annotation table providing co-variables for the samples including
#' condition and time points. Time points must be numeric numbers.
#' @param imp_fit_genes list of  fitted impulse model values and parameters as
#' produced by \code{impulse_DE} as list element \code{impulse_fit_results};
#' another possibility is load the saved fitting object
#' \code{load("impulse_fit_genes.RData")}.
#' @param control_timecourse logical indicating whether a control time
#' timecourse is part of the data set (\code{TRUE}) or not (\code{FALSE}).
#' Default is \code{FALSE}.
#' @param control_name character string specifying the name of the control
#' condition in \code{annotation_table}.
#' @param case_name character string specifying the name of the case
#' condition in \code{annotation_table}. Should be set if more than two
#' conditions are present in \code{data_annotation}.
#' @param file_name_part character string to be used as file extention.
#' @param title_line character string to be used as title for each plot.
#' @param sub_line character string to be used as subtitle for each plot.
#' @export
plot_impulse <- function(gene_IDs, data_array, data_annotation,imp_fit_genes,
                         control_timecourse = FALSE, control_name = NULL, case_name = NULL,
                         file_name_part = "", title_line = "", sub_line = ""){
  
  #print("---Plotting genes")
  
  if(length(grep("[a-zA-Z]",rownames(data_array))) == 0 ){
    rownames(data_array) <- paste(rownames(data_array),"G", sep = "_")
    gene_IDs <- paste(gene_IDs, "G", sep = "_")
  }
  #print(gene_IDs)
  
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
  for (i in 1:length(gene_IDs)){
    # If there is no control data plot only case time course data and fit
    if(control_timecourse == FALSE){
      # Chose impulse fit if parameters of fitted model are not NAN,
      # only plot first timepoint otherwise
      if(TRUE %in% is.na(imp_fit_genes$impulse_parameters_case[gene_IDs[i],])){
        calc_case <- imp_fit_genes$impulse_fits_case[gene_IDs[i],1]
      } else {
        calc_case <- calc_impulse_comp(imp_fit_genes$impulse_parameters_case[gene_IDs[i],1:6],x_vec)
      }
      
      # Plot points
      if(FALSE){
      plot(timep_arr,2^(t(data_array[gene_IDs[i],,])),col="blue",pch=3,xlim=c(0,max(timep)),
           ylim=c(2^(min(c(as.numeric(data_array[gene_IDs[i],,]),as.numeric(calc_case)))-0.5),
                  2^(max(c(as.numeric(data_array[gene_IDs[i],,]),as.numeric(calc_case)))+0.5)),
           xlab="Time", ylab="Impulse fit and expression values",
           main=paste(gene_IDs[i]," ",title_line,sep=""),sub=sub_line)
      points(timep_arr[1,],2^(apply(data_array[gene_IDs[i],,],1,mean)),col="red",pch=1)
      if(TRUE %in% is.na(imp_fit_genes$impulse_parameters_case[gene_IDs[i],])){
        abline(h = 2^calc_case , col = "blue")
      } else {
        points(x_vec, 2^calc_case, col = "blue", type="l")
      }
      }
      if(TRUE){
      # no log
      plot(timep_arr,(t(data_array[gene_IDs[i],,])),col="blue",pch=3,xlim=c(0,max(timep)),
        ylim=c((min(c(as.numeric(data_array[gene_IDs[i],,]),as.numeric(calc_case)))-0.5),
          (max(c(as.numeric(data_array[gene_IDs[i],,]),as.numeric(calc_case)))+0.5)),
        xlab="Time", ylab="Impulse fit and expression values",
        main=paste(gene_IDs[i]," ",title_line,sep=""),sub=sub_line)
      points(timep_arr[1,],(apply(data_array[gene_IDs[i],,],1,mean)),col="red",pch=1)
      if(TRUE %in% is.na(imp_fit_genes$impulse_parameters_case[gene_IDs[i],])){
        abline(h = calc_case , col = "blue")
      } else {
        points(x_vec, calc_case, col = "blue", type="l")
      }
      }
      legend(x="bottomright",as.character(data_annotation[1,"Condition"]),fill=c("blue"), cex=0.6)
      
      ### if there is control data
    } else if(control_timecourse == TRUE){
      
      if(TRUE %in% is.na(imp_fit_genes$impulse_parameters_case[gene_IDs[i],])){
        calc_case = imp_fit_genes$impulse_fits_case[gene_IDs[i],1]
        status_case = FALSE
      } else {
        calc_case = calc_impulse_comp(imp_fit_genes$impulse_parameters_case[gene_IDs[i],1:6],x_vec)
        status_case = TRUE
      }
      if(TRUE %in% is.na(imp_fit_genes$impulse_parameters_control[gene_IDs[i],])){
        calc_ctrl = imp_fit_genes$impulse_fits_control[gene_IDs[i],1]
        status_ctrl = FALSE
      } else {
        calc_ctrl = calc_impulse_comp(imp_fit_genes$impulse_parameters_control[gene_IDs[i],1:6],x_vec)
        status_ctrl = TRUE
      }
      if(TRUE %in% is.na(imp_fit_genes$impulse_parameters_combined[gene_IDs[i],])){
        calc_comb = imp_fit_genes$impulse_fits_combined[gene_IDs[i],1]
        status_comb = FALSE
      } else {
        calc_comb = calc_impulse_comp(imp_fit_genes$impulse_parameters_combined[gene_IDs[i],1:6],x_vec)
        status_comb = TRUE
      }
      
      plot(timep_case_arr,2^(t(farr_case[gene_IDs[i],,])),col="blue",pch=3,xlim=c(0,max(timep)),
           ylim=c(2^(min(c(as.numeric(data_array[gene_IDs[i],,]),as.numeric(calc_case), as.numeric(calc_ctrl), as.numeric(calc_comb)))-0.5),
                  2^(max(c(as.numeric(data_array[gene_IDs[i],,]),as.numeric(calc_case), as.numeric(calc_ctrl), as.numeric(calc_comb)))+0.5)),
           xlab="Time", ylab="Impulse fit or log2 expression value",
           main=paste(gene_IDs[i]," ",title_line,sep=""),sub=sub_line)
      
      points(timep_ctrl,farr_ctrl[gene_IDs[i],],col="red",pch=4)
      
      if(status_case == FALSE){
        abline(h = 2^calc_case , col = "blue")
      } else {
        points(x_vec,2^calc_case, col = "blue", type="l")
      }
      if(status_ctrl == FALSE){
        abline(h = 2^calc_ctrl , col = "red")
      } else {
        points(x_vec,2^calc_ctrl, col = "red", type="l")
      }
      if(status_comb == FALSE){
        abline(h = 2^calc_comb , col = "grey")
      } else {
        points(x_vec,2^calc_comb, col = "grey", type="l")
      }
      legend(x="bottomright",c(as.character(data_annotation$Condition[data_annotation$Condition != control_name][1]),control_name,"combined"),fill=c("blue","red","grey"), cex=0.6)
    }
  }
  dev.off()
}