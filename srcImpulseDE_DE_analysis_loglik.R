#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++     DE analysis   ++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# Detect differentially expressed genes over time

# INPUT:
#   data_array: (Numeric 3D array genes x samples x replicates)
#       Contains expression values or similar locus-specific read-outs.
#   data_annotation: (Table samples x 2[time and condition]) 
#       Co-variables for the samples including condition and time points.
#       Time points must be numeric.
#   weight_mat:
#   imp_fit_results: (list runs*2 ["impulse_parameters","impulse_fits"])
#       impulse_parameters: (matrix genes x (NPARAM+1)*runs) +1 is value of 
#         objective of optimisation (i.e. sum of squares or weighted SS)
#       impulse_fits: (matrix genes x timepoints) model values for gene data
#   background: (vector number of F-scores simulated under H0)  Empirical 
#       F-values simulated under H0.
#   control_timecourse: (bool) [Default FALSE] Control time timecourse is 
#       part of the data set (TRUE) or not (FALSE).
#   control_name: (str) name of the control condition in annotation_table.
#   Q: (scalar) [Defaul 0.01]
# OUTPUT:
#   res2:

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

DE_analysis_loglik <- function(data_array,data_annotation,weight_mat,
  impulse_fit_results,background, control_timecourse = FALSE,
  control_name = NULL,Q = 0.01, dispersion_vector=NULL,NPARAM=6){
  
  ### if control timecourse is present split data into case and control data
  if(control_timecourse == TRUE){
    farr_case <- as.matrix(data_array[,!(data_annotation$Condition %in%
        control_name),])
    farr_ctrl <- as.matrix(data_array[,data_annotation$Condition %in%
        control_name,])
    timep_case <-  as.numeric(as.character(data_annotation[colnames(farr_case),
      "Time"]))
    timep_ctrl <- as.numeric(as.character(data_annotation[colnames(farr_ctrl),
      "Time"]))
  }
  timep <- as.numeric(as.character(data_annotation[colnames(data_array),
    "Time"]))
  
  if(control_timecourse == FALSE){
    loglik_fullmodel <- impulse_fit_results$impulse_parameters_case[,"logL_H1"]
    loglik_redmodel <- impulse_fit_results$impulse_parameters_case[,"logL_H0"]
    mu <- impulse_fit_results$impulse_parameters_case[,"mu"]
    df_fullmodel <- NPARAM
    df_redmodel <- 1
  }
  # Compute difference in degrees of freedom between null model and alternative model.
  df <- df_fullmodel - df_redmodel
  # Compute test statistic: Deviance
  deviance <- 2*(loglik_fullmodel - loglik_redmodel)
  # Get p-values from Chi-square distribution (assumption about null model)
  p <- pchisq(deviance,df=df,lower.tail=FALSE)
  # Multiple testing correction (Benjamini-Hochberg)
  p_BH = p.adjust(p, method = "BH")
  
  # Summarise results
  if(control_timecourse == FALSE){
    result =   as.data.frame(cbind(
      "Gene" = row.names(data_array),
      "p"=as.numeric(p),
      "adj.p"=as.numeric(p_BH),
      "loglik_full"=round(loglik_fullmodel),
      "loglik_red"=round(loglik_redmodel),
      "deviance"=round(deviance),
      "mu"=round(mu),
      "size"=round(dispersion_vector),
      "converge_H1"=impulse_fit_results$impulse_parameters_case[,"converge_H1"],
      "converge_H0"=impulse_fit_results$impulse_parameters_case[,"converge_H0"],
      stringsAsFactors = FALSE))
  }
  result$adj.p <- as.numeric(as.character(result$adj.p))
  result = result[order(result$adj.p),]
  print(paste("Found ",nrow(result[result$adj.p <= Q,])," DE genes",sep=""))
  
  return(result)
}