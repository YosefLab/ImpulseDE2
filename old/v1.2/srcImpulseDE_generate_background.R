#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Background generation   +++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Generates the background for the DE analysis. 
### Algorithm outlined in description of impulse_bg.
### Can be parallelized onto 'n_proc' nodes.

###   0. ["Pre-script"] Reduce expression array to array of variable genes
###   1. [Subfunction impulse_bg] Generate list of F-scores under H0.
###   2. ["Script body"] Parallelise generation of emprirical F-score 
###         distribution.

# INPUT:
#   data_array: (Numeric 3D array genes x samples x replicates).
#       Contains expression values or similar locus-specific read-outs.
#   data_annotation: (Table samples x 2[time and condition]) 
#       Co-variables for the samples including condition and time points.
#       Time points must be numeric.
#   sigma_mat:
#   n_iter: (scalar) [Default 100] Number of iterations, which are performed 
#       to fit the impulse model to the clusters.
#   imp_fit_genes: (list runs*2 ["impulse_parameters","impulse_fits"])
#       impulse_parameters: (matrix genes x (NPARAM+1)*runs) +1 is value of 
#         objective of optimisation (i.e. sum of squares or weighted SS)
#       impulse_fits: (matrix genes x timepoints) model values for gene data
#   control_timecourse: (bool) [Default FALSE] Control time timecourse is 
#       part of the data set (TRUE) or not (FALSE).
#   control_name: (str) name of the control condition in annotation_table.
#   no_of_clus:
#   n_rands:
#   n_proc: (scalar) [Default 4] number of processes, which can be used on 
#       the machine to run the background calculation in parallel
# OUTPUT:
#   res2: (vector n_rands) Empirical F-values simulated under H0.

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

generate_background <- function(data_array, data_annotation, sigma_mat, 
    n_iter = 100, imp_fit_genes, control_timecourse = FALSE, control_name = NULL, 
    no_of_clus = NULL, n_rands = 50000, n_proc = 4){
  
  # Get number of pre- and fine clusters as input from the genes
  # Exclude genes for which the model was not fitted due to almost no variation
  if(control_timecourse == TRUE){
    no_of_clus = no_of_clus[c(3,4,7,8,11,12)]
    # Select genes which were fitted, i.e. not low variation in any case
    data_array <- data_array[!(is.na(imp_fit_genes[[1]][,"beta"]) | is.na(imp_fit_genes[[3]][,"beta"]) | is.na(imp_fit_genes[[5]][,"beta"])),,]
    sigma_mat <- sigma_mat[!(is.na(imp_fit_genes[[1]][,"beta"]) | is.na(imp_fit_genes[[3]][,"beta"]) | is.na(imp_fit_genes[[5]][,"beta"])),]
  } else if(control_timecourse == FALSE & is.null(control_name)){
    no_of_clus = no_of_clus[c(3,4)]
    # Select genes which were fitted, i.e. not low variation
    data_array <- data_array[!(is.na(imp_fit_genes[[1]][,"beta"])),,]
    sigma_mat <- sigma_mat[!(is.na(imp_fit_genes[[1]][,"beta"])),]
  }
  # Get fitting results for all fitted genes (not low variation)
  imp_fit_genes <- lapply(imp_fit_genes,function(x){x <- x[rownames(data_array),]})
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~ 1. Fit background ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # Generates background for subset of total samples to be processed:
  # This function operates on  n_rand <- round(n_randoms_tot/no_proc) samples
  # 1.  n_rand samples are drawn with replacement from all genes.
  # 2.  Calculate residuals for random samples.
  # 3.  Draw random residuals and add them to the null model:
  #     Sample n (timepoints) residuals with replacement from 
  #     n residuals of gene i to add to null model of gene i.
  # 4.  Cluster data simulated under H0. Number of clusters as
  #     determined on observed data.
  # 5.  Fit impulse model to the clusters of data simulated under H0.
  # 6.  Fit impulse model to the genes of data simulated under H0.
  # 7.  Compute sum of squares for H0 and H1.
  # 8.  Compute F-scores.
  
  # INPUT:
  #   data_arr:(Numeric 3D array genes x samples x replicates).
  #       Contains expression values or similar locus-specific read-outs.
  #   data_annot: (Table samples x 2[time and condition]) 
  #       Co-variables for the samples including condition and time points.
  #       Time points must be numeric.
  #   sigma_mat:
  #   n_it: (scalar) [Defaul 100] Number of iterations, which are performed 
  #       to fit the impulse model to the clusters.
  #   ctrl_tk: (bool) [Default FALSE] Control time timecourse is 
  #       part of the data set (TRUE) or not (FALSE).
  #   ctrl_name: (str) name of the control condition in annotation_table.
  #   n_clust:
  #   n_randoms_tot:
  #   n_proc: (scalar) [Default 4] number of processes, which can be used on 
  #       the machine to run the background calculation in parallel
  # OUTPUT:
  #   F_background: (vector round(n_randoms_tot/no_proc)) Empirical F-values
  #       under H0.
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  impulse_bg <- function(data_arr, data_annot, sigma_mat,  n_it = 100, 
    imp_fit, ctrl_tk = FALSE, ctrl_name = NULL, n_clust = NULL, 
    n_randoms_tot = 50000, no_proc = 4){
    
    n_rand <- round(n_randoms_tot/no_proc)  ### divide randomizations by no. of nodes
    
    timep <- as.numeric(as.character(data_annot[colnames(data_arr),"Time"]))
    NPARAM = 6
    
    # 1.  n_rand samples are drawn with replacement from all genes.
    i_s <- sample(1:nrow(data_arr),n_rand, replace = TRUE)
    data_arr_rand <- data_arr[i_s,,]
    colnames(data_arr_rand) <- colnames(data_arr)
    rownames(data_arr_rand) <- paste("bg_",c(1:length(i_s)))
    sigma_mat_rand <- sigma_mat[i_s,]
    colnames(sigma_mat_rand) <- colnames(data_arr_rand)
    rownames(sigma_mat_rand) <- rownames(data_arr_rand)
    
    #############################
    # I.  With control data
    if(ctrl_tk == TRUE){
      data_list <- list(data_arr_rand, 
                        data_arr_rand[,!(data_annot$Condition %in% ctrl_name),],
                        data_arr_rand[,data_annot$Condition %in% ctrl_name,])
      sigma_list <- list(sigma_mat_rand,
                        sigma_mat_rand[,!(data_annot$Condition %in% ctrl_name)],
                        sigma_mat_rand[,data_annot$Condition %in% ctrl_name])
            
      # 2.  Calculate residuals for random samples.
      # Get combined data (H0)
      y_i_null_s = (imp_fit[[2]][rownames(data_arr),])[i_s, as.character(data_annot[colnames(data_list[[1]]),"Time"])]
      # Get case time courses (H1)
      calc_impulse_comp_data_arr_rand_case <- 
          (imp_fit[[4]][rownames(data_arr),])[i_s,as.character(data_annot[colnames(data_list[[2]]),"Time"])]
      # Get control time courses (H1)
      calc_impulse_comp_data_arr_rand_ctrl <- 
          (imp_fit[[6]][rownames(data_arr),])[i_s, as.character(data_annot[colnames(data_list[[3]]),"Time"])]
      
      # Compute residuals for all replicates and normalise residuals to be 
      # i.i.d., i.e. correct for non-uniform variance.
      nRESD_1_s = array(NA,dim(data_arr_rand))
      for(z in 1:dim(data_arr_rand)[3]){
        nRESD_1_s[,,z] = cbind( (data_list[[2]][,,z] - calc_impulse_comp_data_arr_rand_case) / sigma_list[[2]],
                                (data_list[[3]][,,z] - calc_impulse_comp_data_arr_rand_ctrl) / sigma_list[[3]] )
      }
      nRESD_1_s = nRESD_1_s[,colnames(data_list[[1]]),]
      
      # 3.  Draw random residuals and add them to the null model:
      #     Sample n (timepoints) residuals with replacement from 
      #     n residuals of gene i to add to null model of gene i.
      
      # Sampled and scaled up residuals
      newRESD_1_s = array(NA,dim(nRESD_1_s))
      # Bootstrapped samples under H0: Mean + bootstrapped residual
      Y_star_s = array(NA,dim(nRESD_1_s))
      for(z in 1:dim(nRESD_1_s)[3]){
        # Normalised residuals are drawn across time points for each gene x and replicate
        # z. The residuals are then scaled by the standard deviation of the timrpoint
        # of the gene they have been sampled into. The returned array has the form
        # (genes x timepoints)
        newRESD_1_s = t(apply(nRESD_1_s[,,z],1,function(x){x[sample(1:ncol(nRESD_1_s),
            ncol(nRESD_1_s),replace=TRUE)]})) * sigma_mat_rand
        # Generate bootstrapped samples: (genes x timepoints) combined +
        # (genes x timepoints) bootstrapped residuals.
        Y_star_s[,,z] = y_i_null_s + newRESD_1_s
      }
      colnames(Y_star_s) <- colnames(data_arr_rand)
      rownames(Y_star_s) <- paste("bg_",c(1:length(i_s)))
      
      # Compute standard deviations over bootstrapped replicates
      sigma_mat_bg <- compute_stdvs(observed_data=Y_star_s)
      
      # 4.  Cluster data simulated under H0. Number of clusters as
      #     determined on observed data. Cluster on means over 
      #     bootstrapped replicates.
      tm_c <- system.time({
        clustering_results_background <- cluster_genes_for_impulse(
            apply(Y_star_s,c(1,2),mean), data_annot,ctrl_tk, ctrl_name, 
            plot_clusters = TRUE, no_of_clusters = n_clust, 
            n_genes = nrow(data_arr))
      })
      print(paste("Consumed time bg clustering: ",round(tm_c["elapsed"]/60,2)," min",sep=""))
      
      # 5. Fit impulse model to the clusters of data simulated under H0.
      tm_f <- system.time({
        impulse_fit_clusters_bg <- impulse_fit(clustering_results_background,
            data_annot,sigma_mat_bg,n_it, control_timecourse = ctrl_tk,
            ctrl_name, fit_backg = TRUE)
      })
      print(paste("Consumed time bg clus fit: ",round(tm_f["elapsed"]/60,2)," min",sep=""))
      
      # 6.  Fit impulse model to the genes of data simulated under H0.
      tm_fg <- system.time({
        impulse_fit_genes_bg <- impulse_fit(Y_star_s, data_annot, 
            sigma_mat_bg, 1, ctrl_tk, ctrl_name,clustering_results_background, 
            impulse_fit_clusters_bg, fit_backg = TRUE )
      })
      #      save(impulse_fit_genes_bg, file=file.path(getwd(),paste("impulse_fit_genes_bg_", round(runif(1)*10000),".RData",sep="")))
      print(paste("Consumed time bg gene fit: ",round(tm_fg["elapsed"]/60,2)," min",sep=""))

      # 7.  Compute residuals for H0 and H1.
      RESD_1_s_b = array(NA,dim(Y_star_s))
      RESD_0_s_b = array(NA,dim(Y_star_s))
      for (z in 1:dim(Y_star_s)[3]){
        RESD_1_s_b[,,z] = cbind((Y_star_s[rownames(data_list[[1]]),colnames(data_list[[2]]),z] -
          impulse_fit_genes_bg[[4]][rownames(data_list[[1]]),
          as.character(data_annot[colnames(data_list[[2]]),"Time"])]),
          (Y_star_s[rownames(data_list[[1]]),colnames(data_list[[3]]),z] -
          impulse_fit_genes_bg[[6]][rownames(data_list[[1]]),
          as.character(data_annot[colnames(data_list[[3]]),"Time"])]))
        
        RESD_0_s_b[,,z] = Y_star_s[rownames(data_list[[1]]),colnames(data_list[[1]]),z] -
          impulse_fit_genes_bg[[2]][rownames(data_list[[1]]),
          as.character(data_annot[colnames(data_list[[1]]),"Time"])]
      }
      colnames(RESD_1_s_b) <- c(colnames(data_list[[2]]), colnames(data_list[[3]]))
      RESD_1_s_b = RESD_1_s_b[,colnames(data_list[[1]]),]
      
    ###############################
    # II. No control data
    } else {
      # 2.  Calculate residuals for random samples.
      # Under the assumption of equal replicate means, the residuals
      # are computed relative to the same overall row mean.
      
      # DEVELOPMENT NOTE: In v1.3, the model will be conditioned on patient-gene
      # means by setting them to 0. Therefore the overall mean will 
      # be 0 as well and the following still valid.
      
      y_i_null_s = rowMeans(data_arr_rand)
      calc_impulse_comp_data_arr_rand <- (imp_fit[[2]][rownames(data_arr),])[i_s,as.character(data_annot[,"Time"])]

      # Compute residuals for all replicates and normalise residuals to be 
      # i.i.d., i.e. correct for non-uniform variance.
      nRESD_1_s = array(NA,dim(data_arr_rand))
      for(z in 1:dim(data_arr_rand)[3]){
        nRESD_1_s[,,z] = (data_arr_rand[,,z] - calc_impulse_comp_data_arr_rand) / sigma_mat_rand
      }
      
      # 3.  Draw random residuals and add them to the null model:
      #     Sample n (timepoints) residuals with replacement from 
      #     n residuals of gene i to add to null model of gene i.
      #     The null model data + residuals (data simulated under
      #     H0) are in Y_star_s. Residuals are drawn for each
      #     replicate independently to create simulated data set
      #     with the same number of replicates.
      
      # Sampled and scaled up residuals
      newRESD_1_s = array(NA,dim(nRESD_1_s))
      # Bootstrapped samples under H0: Mean + bootstrapped residual
      Y_star_s = array(NA,dim(nRESD_1_s))
      for(z in 1:dim(nRESD_1_s)[3]){
        # Normalised residuals are drawn across time points for each gene x and replicate
        # z. The residuals are then scaled by the standard deviation of the timrpoint
        # of the gene they have been sampled into. The returned array has the form
        # (timepoints x genes)
        newRESD_1_s = apply(nRESD_1_s[,,z],1,function(x){x[sample(1:ncol(nRESD_1_s),
            ncol(nRESD_1_s),replace=TRUE)]}) * t(sigma_mat_rand)
        # Generate bootstrapped samples: (1 x genes) means + (timepoints x genes)
        # bootstrapped residuals. Take transpose to transform into 
        # (genes x timepoints)
        Y_star_s[,,z] = t(y_i_null_s + newRESD_1_s)
      }
      colnames(Y_star_s) <- colnames(data_arr_rand)
      rownames(Y_star_s) <- paste("bg_",c(1:length(i_s)))
      
      # Compute standard deviations over bootstrapped replicates
      sigma_mat_bg <- compute_stdvs(observed_data=Y_star_s)
      
      # 4.  Cluster data simulated under H0. Number of clusters as
      #     determined on observed data. Cluster on means over 
      #     bootstrapped replicates.
      tm_c <- system.time({
          clustering_results_background <- cluster_genes_for_impulse(
              apply(Y_star_s,c(1,2),mean), data_annot,ctrl_tk, ctrl_name, 
              plot_clusters = FALSE, no_of_clusters = n_clust, 
              n_genes = nrow(data_arr))
      })
      print(paste("Consumed time bg clustering: ",round(tm_c["elapsed"]/60,2)," min",sep=""))

      # 5. Fit impulse model to the clusters of data simulated under H0.
      tm_f <- system.time({
        impulse_fit_clusters_bg <- impulse_fit(clustering_results_background, data_annot, 
            sigma_mat_bg, n_it, control_timecourse = ctrl_tk, ctrl_name, fit_backg = TRUE)
      })
      print(paste("Consumed time bg clust fit: ",round(tm_f["elapsed"]/60,2)," min",sep=""))
      
      # 6.  Fit impulse model to the genes of data simulated under H0.
      ### DAVID changed n_iter from 1 to n_it
      tm_fg <- system.time({
        impulse_fit_genes_bg <- impulse_fit(Y_star_s, data_annot, sigma_mat_bg, 
            n_it, ctrl_tk, ctrl_name, clustering_results_background, impulse_fit_clusters_bg, 
            fit_backg = TRUE )
      })
      #      save(impulse_fit_genes_bg, file=file.path(getwd(),paste("impulse_fit_genes_bg_", round(runif(1)*10000),".RData",sep="")))
      print(paste("Consumed time bg gene fit: ",round(tm_fg["elapsed"]/60,2)," min",sep=""))
      
      
      # 7.  Compute residuals for H0 and H1.
      RESD_1_s_b = array(NA,dim(Y_star_s))
      RESD_0_s_b = array(NA,dim(Y_star_s))
      for (z in 1:dim(Y_star_s)[3]){
        RESD_1_s_b[,,z] = Y_star_s[,,z] - impulse_fit_genes_bg[[2]][rownames(Y_star_s),
          as.character(data_annot[colnames(Y_star_s),"Time"])]
        RESD_0_s_b[,,z] = Y_star_s[,,z] - rowMeans(Y_star_s)
      }
    }

    # 8.  Compute F-like-scores. The F-like score is a weighted
    #     sum of squares not corrected for by degrees of freedom:
    #     It does not assume a parameteric form of the residual
    #     distribution (c.f. F-score assumes normal errors).
    F_background_rep = array(0,c(dim(RESD_1_s_b)[1],dim(RESD_1_s_b)[3]))
    for (z in 1:dim(RESD_1_s_b)[3]){
      F_background_rep[,z] = rowSums( (RESD_0_s_b[,,z] - RESD_1_s_b[,,z])^2 *
          1/(sigma_mat_bg^2) )
    }

    F_background = rowSums(F_background_rep)
    
    print(dim(F_background_rep))
    print(length(F_background))
    
    return(F_background)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~ 2. Parallelise background fitting ~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # set number of nodes to maximum - 1
  mc <- min(detectCores() - 1, n_proc)

  # Define a new environment, which is used for parallelization
  my.env <- new.env()
  assign("n_rands", n_rands, envir = my.env)
  assign("mc", mc, envir = my.env)
  assign("data_array", data_array, envir = my.env)
  assign("data_annotation", data_annotation, envir = my.env)
  assign("sigma_mat", sigma_mat, envir = my.env)
  assign("n_iter", n_iter, envir = my.env)
  assign("control_timecourse", control_timecourse, envir = my.env)
  assign("control_name", control_name, envir = my.env)
  assign("imp_fit_genes", imp_fit_genes, envir = my.env)
  assign("no_of_clus", no_of_clus, envir = my.env)
  assign("cluster_genes_for_impulse", cluster_genes_for_impulse, envir = my.env)
  assign("impulse_fit", impulse_fit, envir = my.env)
  assign("calc_impulse_comp", calc_impulse_comp, envir = my.env)
  assign("two_impulses_WLS_comp", two_impulses_WLS_comp, envir = my.env)
  assign("two_impulses_OLS_comp", two_impulses_OLS_comp, envir = my.env)
  assign("compute_stdvs", compute_stdvs, envir = my.env)
  environment(impulse_bg) <- my.env
  
  # Create Cluster
  cl <- makeCluster(mc, outfile="cluster_out_random.txt")
  # Export variables and functions to each node
  clusterExport(cl=cl, varlist=c("impulse_bg","n_rands","mc","data_array","data_annotation",
                                 "sigma_mat","n_iter","control_timecourse","control_name",
                                 "imp_fit_genes","mc", "no_of_clus","calc_impulse_comp",
                                 "cluster_genes_for_impulse","impulse_fit","two_impulses_WLS_comp",
                                  "two_impulses_OLS_comp","compute_stdvs"), envir = my.env)
  # Load libraries onto each node:
  junk <- clusterEvalQ(cl,library(boot))
  ### DAVID Bug fix need amap,parallel on each node
  # Load amap for Kmeans
  junk <- clusterEvalQ(cl,library(amap))
  junk <- clusterEvalQ(cl,library(parallel))
  clusterSetRNGStream(cl,123)
  # Parallelized generation of background F-scores: is independent
  # and parallelization can be carried out on arbitrary many nodes
  # operating independent of each other. All nodes sample a number
  # genes at random, simulate data under H0 and compute the F-scores
  # which are merged from all nodes to give the background.
  res <- clusterEvalQ(cl,impulse_bg(data_arr=data_array, data_annot=data_annotation, 
      sigma_mat=sigma_mat, n_it=n_iter,imp_fit=imp_fit_genes, 
      ctrl_tk=control_timecourse, ctrl_name=control_name, n_clust=no_of_clus, 
      n_randoms_tot=n_rands, no_proc=mc))
  environment(res) <- my.env
  res2 <- unlist(res)
  # Close cluster
  stopCluster(cl)
  return(res2)

  ## if not paralellized
  #print("Unparallelized:")
  #res <- impulse_bg(data_array, data_annotation, sigma_mat, n_iter, imp_fit_genes,
  #       control_timecourse, control_name, no_of_clus, n_rands, 1)
  #  return(res)
}
