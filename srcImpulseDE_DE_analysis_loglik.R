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
#   e_type:
#   Q: (scalar) [Defaul 0.01]
# OUTPUT:
#   res2:

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

DE_analysis_loglik <- function(data_array,data_annotation,weight_mat,
  impulse_fit_results,background, control_timecourse = FALSE,
  control_name = NULL,e_type = "Array", Q = 0.01, dispersion_vector=NULL){
  
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
    loglik_fullmodel <- -impulse_fit_results$impulse_parameters_case[,"objective"]
    loglik_redmodel <- -impulse_fit_results$impulse_parameters_case[,"nullfit"]
    df_fullmodel <- 6
    df_redmodel <- 1
  }
  df <- df_fullmodel - df_redmodel
  deviance <- 2*(loglik_fullmodel - loglik_redmodel)
  p <- pchisq(deviance,df=df,lower.tail=FALSE)
  
  p_BH = p.adjust(p, method = "BH")
  
  # SS1 Flag
  ### DAVID: SS flag still valid under my model?
  ### Could leave like this or find new threshold for weighted SS
  #SS_1_s = apply(RESD_1_s, 1, function(x){sum(x^2)})
  #SS1_flag = SS_1_s < 20
  SS1_flag <- array(1,length(p_BH))
  error_index = unlist(lapply(SS1_flag,function(x){if(x == FALSE){
    "Fit stability measure not build yet, srcImpulseDE_DE_analysis"} else {""}}))
  
  # Create warning if no non-DE genes were detected
  if(!(TRUE %in% (p_BH <= Q))){
    warning("No DE genes were detected. Maybe amount of background genes is
              too low.")
  }  
  
  # Summarise results
  result =   as.data.frame(cbind(
    "Gene" = row.names(data_array),
    "adj.p"=as.numeric(p_BH),
    "loglik_full"=round(loglik_fullmodel),
    "loglik_red"=round(loglik_redmodel),
    "deviance"=round(deviance),
    "mu"=round(apply(data_array,1,mean)),
    "size"=round(dispersion_vector),
    stringsAsFactors = FALSE))
  result$adj.p <- as.numeric(as.character(result$adj.p))
  result = result[order(result$adj.p),]
  print(paste("Found ",nrow(result[result$adj.p <= Q,])," DE genes",sep=""))
  
  
  if(control_timecourse == TRUE){
    write.table(as.data.frame(cbind("Gene" = row.names(data_array),
      "adj.p"=p_BH,
      "prediction error" = error_index)),"pvals_and_flags.txt",
      quote = FALSE, sep ="\t", row.names = TRUE, col.names = NA)
  } else {
    if(e_type == "Array"){
      write.table(as.data.frame(cbind("Gene" = row.names(data_array),
        "adj.p"=p_BH,"prediction error" = error_index)),
        "pvals_and_flags.txt", quote = FALSE, sep ="\t", row.names = TRUE,
        col.names = NA)
    } else {
      write.table(as.data.frame(cbind("Gene" = row.names(data_array),
        "adj.p"=p_BH,"prediction error" = error_index)),
        "pvals_and_flags.txt", quote = FALSE, sep ="\t",
        row.names = TRUE, col.names = NA)
    }
  }
  return(result)
}