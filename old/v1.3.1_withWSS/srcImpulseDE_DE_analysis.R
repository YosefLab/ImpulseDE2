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

DE_analysis <- function(data_array,data_annotation,weight_mat,
  impulse_fit_results,background, control_timecourse = FALSE,
  control_name = NULL,e_type = "Array", Q = 0.01){
  
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
  
  ### correct background values
  ### may occur as model fitting not convex
  background[background < 0] = 0
  #background = background[background < 100]
  
  if(control_timecourse == TRUE){
    y_i_null = impulse_fit_results[[2]][rownames(farr_case),as.character(timep)]
    
    RESD_1_s = array(NA,c(dim(farr_case)[1],dim(data_array)[2],dim(data_array)[3]))
    RESD_0_s = array(NA,c(dim(farr_case)[1],dim(data_array)[2],dim(data_array)[3]))
    for (z in 1:dim(data_array)[3]){
      RESD_1_s[,,z] = cbind(
        (farr_case[,,z] -
            impulse_fit_results[[4]][rownames(farr_case),as.character(timep_case)]),
        (farr_ctrl[rownames(farr_case),,z] -
            impulse_fit_results[[6]][rownames(farr_case),as.character(timep_ctrl)]))
      
      RESD_0_s[,,z] = data_array[,,z] - y_i_null
    }
    RESD_1_s = RESD_1_s[,colnames(data_array)]
  } else if(control_timecourse == FALSE){
    y_i_null = rowMeans(data_array)
    
    RESD_1_s = array(NA,dim(data_array))
    RESD_0_s = array(NA,dim(data_array))
    for (z in 1:dim(data_array)[3]){
      RESD_1_s[,,z] = data_array[,,z] - impulse_fit_results[[2]][rownames(data_array),
        as.character(timep)]
      RESD_0_s[,,z] = data_array[,,z] - y_i_null
    }
  }
  
  # Compute F-like-scores. The F-like score is a weighted
  # sum of squares not corrected for by degrees of freedom:
  # It does not assume a parameteric form of the residual
  # distribution (c.f. F-score assumes normal errors).
  #F_stats_rep = array(0,c(dim(RESD_1_s)[1],dim(RESD_1_s)[3]))
  #for (z in 1:dim(F_stats_rep)[2]){
  #  F_stats_rep[,z] = rowSums( ((RESD_0_s[,,z])^2 - (RESD_1_s[,,z])^2) *
  #      1/(weight_mat^2) )
  #}
  WSS0_s_rep <- array(0,c(dim(RESD_1_s)[1],dim(RESD_1_s)[3]))
  WSS1_s_rep <- array(0,c(dim(RESD_1_s)[1],dim(RESD_1_s)[3]))
  for (z in 1:dim(RESD_1_s)[3]){
    WSS0_s_rep[,z] <- rowSums( (RESD_0_s[,,z])^2 * 1/(weight_mat^2) )
    WSS1_s_rep[,z] <- rowSums( (RESD_1_s[,,z])^2 * 1/(weight_mat^2) )
  }
  WSS0_s <- rowSums(WSS0_s_rep)
  WSS1_s <- rowSums(WSS1_s_rep)
  #F_background = rowSums(F_background_rep)
  F_stats <- (WSS0_s - WSS1_s)/WSS1_s
  #F_stats = rowSums(F_stats_rep)
  
  # Bootstrap the p-values
  p  <-  unlist(lapply(F_stats,function(x){length(which(background >=
      x))/length(background)}))
  p_scaled = p.adjust(p, method = "BH")
  p_scaled_orig = p_scaled
  
  ### DAVID: deprecate this?
  # FCs timepoint vs earlier
  if(control_timecourse == FALSE){
    # Calculate fold changes based on means over replicates
    data_array_mean <- apply(data_array,c(1,2),mean)
    if(e_type == "Array"){
      means <- t(apply(data_array_mean, 1,function(x){
        tmp = NULL
        for(i in 1:length(unique(timep))){
          ### DAVID: log space
          tmp = c(tmp,mean(x[which(data_annotation[colnames(data_array),
            "Time"] == sort(unique(timep))[i])]))          
          #tmp = c(tmp,mean(2^x[which(data_annotation[colnames(data_array),
          #  "Time"] == sort(unique(timep))[i])]))
        }
        return(tmp)
      }))
      colnames(means) <- sort(unique(timep))
      # Fold change matrix: genes x number of time transitions (timepoints -1)
      ### DAVID code readability: this is equivalent to matrix(0,nrow(data_array),length(unique(timep))-1))
      Ratios_TPs <- matrix(rep(0,nrow(data_array_mean)*(length(unique(timep))-1)),
        nrow(data_array_mean),(length(unique(timep))-1))
      for(tt in 1:(length(unique(timep))-1)){
        # Get ratios between neighbouring time points
        Ratios_TPs[,tt] <- means[,tt+1]/means[,tt]
      }
      
      FC_TP_DEindex <- apply(Ratios_TPs,1,function(x){
        # Label genes as DE if have at least 2 transitions of factor 2
        if(length(which(x > 2 || x < 0.5)) >= 2 ){ "DE" } else {"not_DE"}
      })
    }
  }
  
  if(control_timecourse == TRUE){
    # Calculate fold changes based on means over replicates
    data_array_mean <- apply(data_array,c(1,2),mean)
    farr_case_mean <- apply(farr_case,c(1,2),mean)
    farr_ctrl_mean <- apply(farr_ctrl,c(1,2),mean)
    
    # FCs case vs control
    means_case <- t(apply(farr_case_mean, 1,function(x){
      tmp = NULL
      for(i in 1:length(unique(timep_case))){
        tmp = c(tmp, mean(2^x[which(data_annotation[colnames(farr_case),
          "Time"] == sort(unique(timep_case))[i])]))
      }
      return(tmp)
    }))
    colnames(means_case) <- sort(unique(timep_case))
    
    means_control <- t(apply(farr_ctrl_mean, 1,function(x){
      tmp = NULL
      for(i in 1:length(unique(timep_ctrl))){
        tmp = c(tmp, mean(2^x[which(data_annotation[colnames(farr_ctrl),
          "Time"] == sort(unique(timep_ctrl))[i])]))
      }
      return(tmp)
    }))
    colnames(means_control) <- sort(unique(timep_ctrl))
    
    common_timepoints = sort(unique(timep_case))[sort(unique(timep_case)) %in%
        unique(timep_ctrl)]
    if(length(common_timepoints) == 0){
      FC_DEindex = rep("not_detectable",nrow(means_case))
    } else {
      meanRatios <- t(apply(cbind(means_case[,as.character(common_timepoints)],
        means_control[,as.character(common_timepoints)]),1,function(x){
          tmp <- x[1:length(unique(common_timepoints))]/
            (x[(length(unique(common_timepoints))+1):(2*
                length(unique(common_timepoints)))])
          return(tmp)
        }))
      FCs <- t(apply(meanRatios,1,function(x){
        # translate FC <1 into FC >1 by taking -inverse
        y = x
        y[x<1] = -1/x[x<1]
        return(y)
      }))
      FC_DEindex <- apply(FCs,1,function(x){
        if(length(which(abs(x)> 2)) >= 2 ){ "DE" } else {"not_DE"}
      })
    }
    
    # ANOVA
    anovas <- t(apply(data_array_mean,1,function(x){summary(aov(as.numeric(x) ~
        data_annotation$Time + data_annotation$Condition +
        data_annotation$Time * data_annotation$Condition))[[1]][2:3,
          "Pr(>F)"]}))
    anovas_FDR =  cbind(p.adjust(anovas[,1], method = "BH"),
      p.adjust(anovas[,2], method = "BH"))
    anovas_pre = anovas_FDR[,1] < Q | anovas_FDR[,2] < Q
    anovas_index <- unlist(lapply(anovas_pre, function(x){if(x == TRUE){"DE"
    } else {"not_DE"}}))
  }
  
  # SS1 Flag
  ### DAVID: SS flag still valid under my model?
  ### Could leave like this or find new threshold for weighted SS
  #SS_1_s = apply(RESD_1_s, 1, function(x){sum(x^2)})
  #SS1_flag = SS_1_s < 20
  SS1_flag <- array(1,length(p_scaled))
  error_index = unlist(lapply(SS1_flag,function(x){if(x == FALSE){
    "Fit stability measure not build yet, srcImpulseDE_DE_analysis"} else {""}}))
  
  ### exit the function without error if no DE genes are detected
  if(!(TRUE %in% (p_scaled <= Q))){
    warning("No DE genes were detected. Maybe amount of background genes is
              too low.")
    return(list(NULL,F_stats))
    
    ### if DE genes are detected, finish FDR correction by using the cutoff
  } else {
    
    ### if control data is present but not as a timecourse, add t-test
    ### p-values to the results
    result =   as.data.frame(cbind("Gene" = row.names(data_array),
      "adj.p"=as.numeric(p_scaled_orig), stringsAsFactors = FALSE))
    result$adj.p <- as.numeric(as.character(result$adj.p))
    result = result[order(result$adj.p),]
    print(paste("Found ",nrow(result[result$adj.p <= Q,])," DE genes",sep=""))
    
    
    if(control_timecourse == TRUE){
      write.table(as.data.frame(cbind("Gene" = row.names(data_array),
        "adj.p"=p_scaled_orig, "at least 2 TPs for case vs. control" =
          FC_DEindex, "ANOVA condition or time*condition" = anovas_index,
        "prediction error" = error_index)),"pvals_and_flags.txt",
        quote = FALSE, sep ="\t", row.names = TRUE, col.names = NA)
    } else {
      if(e_type == "Array"){
        write.table(as.data.frame(cbind("Gene" = row.names(data_array),
          "adj.p"=p_scaled_orig, "at least 2 consecutive TP Ratios" =
            FC_TP_DEindex, "prediction error" = error_index)),
          "pvals_and_flags.txt", quote = FALSE, sep ="\t", row.names = TRUE,
          col.names = NA)
      } else {
        write.table(as.data.frame(cbind("Gene" = row.names(data_array),
          "adj.p"=p_scaled_orig,"prediction error" = error_index)),
          "pvals_and_flags.txt", quote = FALSE, sep ="\t",
          row.names = TRUE, col.names = NA)
      }
    }
    return(list(result[as.numeric(result$adj.p) <= Q,],F_stats))
  }
}