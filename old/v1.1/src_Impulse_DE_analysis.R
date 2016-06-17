#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++     DE analysis   ++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# Detect differentially expressed genes over time

# INPUT:
#   data_table: (Numeric matrix genes x samples) Expression values 
#   data_annotation: (Table samples x 2[time and condition]) 
#       Co-variables for the samples including condition and time points.
#       Time points must be numeric.
#   impulse_fit_results: (list runs*2 ["impulse_parameters","impulse_fits"])
#       impulse_parameters: (matrix genes x (NPARAM+1)*runs) +1 is value of 
#       objective of optimisation (i.e. sum of squares or weighted SS)
#       impulse_fits: (matrix genes x runs) model values for gene data
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

DE_analysis <- function(data_table,data_annotation,impulse_fit_results,
                        background, control_timecourse = FALSE, control_name = NULL,
                        e_type = "Array", Q = 0.01){
  
  ### if control timecourse is present split data into case and control data
  if(control_timecourse == TRUE){
    fmat_case <- as.matrix(data_table[,!(data_annotation$Condition %in%
                                           control_name)])
    fmat_ctrl <- as.matrix(data_table[,data_annotation$Condition %in%
                                        control_name])
    timep_case <-  as.numeric(as.character(data_annotation[colnames(fmat_case),
                                                           "Time"]))
    timep_ctrl <- as.numeric(as.character(data_annotation[colnames(fmat_ctrl),
                                                          "Time"]))
  }
  timep <- as.numeric(as.character(data_annotation[colnames(data_table),
                                                   "Time"]))
  
  ### correct background values
  background[background < 0] = 0
  background = background[background < 100]
  
  if(control_timecourse == TRUE){
    y_i_null = impulse_fit_results[[2]][rownames(fmat_case),as.character(timep)]
    RESD_1_s = cbind(fmat_case - impulse_fit_results[[4]][rownames(fmat_case),
                                                          as.character(timep_case)], (fmat_ctrl[rownames(fmat_case),] -
                                                                                        impulse_fit_results[[6]][rownames(fmat_case),as.character(timep_ctrl)]))
    RESD_1_s = RESD_1_s[,colnames(data_table)]
    
    ### if there is no control data
  } else if(control_timecourse == FALSE){
    y_i_null = rowMeans(data_table)
    RESD_1_s = data_table - impulse_fit_results[[2]][rownames(data_table),
                                                     as.character(timep)]
  }
  
  ### DAVID: F-scores calculated here: Change 3 - different F-score
  RESD_0_s = data_table - y_i_null
  SS_0_s = apply(RESD_0_s, 1, function(x){sum(x^2)})
  SS_1_s = apply(RESD_1_s, 1, function(x){sum(x^2)})
  
  F_stats =  (SS_0_s - SS_1_s) / SS_1_s
  ### avoid cases where SS_1 = 0  & SS_0 = 0--> Quotient will be NaN
  F_stats[ which(SS_0_s == 0 & SS_1_s == 0)] = 0
  
  ### bootstrap the p-values
  ### DAVID: Bootstrap - changed way background is collected previously (change 4) - leave this section
  p  =  unlist(lapply(F_stats,function(x){length(which(background >=
                                                         x))/length(background)}))
  p_scaled = p.adjust(p, method = "BH")
  p_scaled_orig = p_scaled
  
  
  # calculate FCs
  ### DAVID: fold change, not affected by me, arbitrary cut offs for selecting DE
  
  # FCs timepoint vs earlier
  if(control_timecourse == FALSE){
    if(e_type == "Array"){
      means <- t(apply(data_table, 1,function(x){
        tmp = NULL
        for(i in 1:length(unique(timep))){
          tmp = c(tmp,mean(2^x[which(data_annotation[colnames(data_table),
                                                     "Time"] == sort(unique(timep))[i])]))
        }
        return(tmp)
      }))
      colnames(means) <- sort(unique(timep))
      Ratios_TPs <- matrix(rep(0,nrow(data_table)*(length(unique(timep))-1)),
                           nrow(data_table),(length(unique(timep))-1))
      for(tt in 1:(length(unique(timep))-1)){
        Ratios_TPs[,tt] <- means[,tt+1]/means[,tt]
      }
      
      FC_TP_DEindex <- apply(Ratios_TPs,1,function(x){
        if(length(which(x > 2 || x < 0.5)) >= 2 ){ "DE" } else {"not_DE"}
      })
    }
    
    
    if(control_timecourse == TRUE){
      
      # FCs case vs control
      means_case <- t(apply(fmat_case, 1,function(x){
        tmp = NULL
        for(i in 1:length(unique(timep_case))){
          tmp = c(tmp, mean(2^x[which(data_annotation[colnames(fmat_case),
                                                      "Time"] == sort(unique(timep_case))[i])]))
        }
        return(tmp)
      }))
      colnames(means_case) <- sort(unique(timep_case))
      
      means_control <- t(apply(fmat_ctrl, 1,function(x){
        tmp = NULL
        for(i in 1:length(unique(timep_ctrl))){
          tmp = c(tmp, mean(2^x[which(data_annotation[colnames(fmat_ctrl),
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
          y = x
          y[x<1] = -1/x[x<1]
          return(y)
        }))
        FC_DEindex <- apply(FCs,1,function(x){
          if(length(which(abs(x)> 2)) >= 2 ){ "DE" } else {"not_DE"}
        })
      }
      
      # ANOVA
      anovas <- t(apply(data_table,1,function(x){summary(aov(as.numeric(x) ~
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
    SS1_flag = SS_1_s < 20
    error_index = unlist(lapply(SS1_flag,function(x){if(x == FALSE){"fit might be
                                                                    unstable due to large variation"} else {""}}))
    
    ### exit the function without error if no DE genes are detected
    if(!(TRUE %in% (p_scaled <= Q))){
      warning("No DE genes were detected. Maybe amount of background genes is
              too low.")
      return(NULL)
      
      ### if DE genes are detected, finish FDR correction by using the cutoff
    } else {
      
      ### if control data is present but not as a timecourse, add t-test
      ### p-values to the results
      result =   as.data.frame(cbind("Gene" = row.names(data_table),
                                     "adj.p"=as.numeric(p_scaled_orig), stringsAsFactors = FALSE))
      result$adj.p <- as.numeric(as.character(result$adj.p))
      result = result[order(result$adj.p),]
      print(paste("Found ",nrow(result[result$adj.p <= Q,])," DE genes",sep=""))
      
      
      if(control_timecourse == TRUE){
        write.table(as.data.frame(cbind("Gene" = row.names(data_table),
                                        "adj.p"=p_scaled_orig, "at least 2 TPs for case vs. control" =
                                          FC_DEindex, "ANOVA condition or time*condition" = anovas_index,
                                        "prediction error" = error_index)),"pvals_and_flags.txt",
                    quote = FALSE, sep ="\t", row.names = TRUE, col.names = NA)
      } else {
        if(e_type == "Array"){
          write.table(as.data.frame(cbind("Gene" = row.names(data_table),
                                          "adj.p"=p_scaled_orig, "at least 2 consecutive TP Ratios" =
                                            FC_TP_DEindex, "prediction error" = error_index)),
                      "pvals_and_flags.txt", quote = FALSE, sep ="\t", row.names = TRUE,
                      col.names = NA)
        } else {
          write.table(as.data.frame(cbind("Gene" = row.names(data_table),
                                          "adj.p"=p_scaled_orig,"prediction error" = error_index)),
                      "pvals_and_flags.txt", quote = FALSE, sep ="\t",
                      row.names = TRUE, col.names = NA)
        }
      }
      return(result[as.numeric(result$adj.p) <= Q,])
    }
    }