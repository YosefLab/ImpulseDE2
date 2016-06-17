#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Background generation   +++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Generates the background for the DE analysis. 
### Algorithm outlined in description of impulse_bg.
### Can be parallelized onto 'n_proc' nodes.

###   1. [Subfunction impulse_bg] Generate list of F-scores under H0.
###   2. ["Script body"] Parallelise generation of emprirical F-score 
###         distribution.

# INPUT:
#   data_table: (Numeric matrix genes x samples) Expression values 
#   data_annotation: (Table samples x 2[time and condition]) 
#       Co-variables for the samples including condition and time points.
#       Time points must be numeric.
#   sigma_vec:
#   n_iter: (scalar) [Default 100] Number of iterations, which are performed 
#       to fit the impulse model to the clusters.
#   imp_fit_genes: (list runs*2 ["impulse_parameters","impulse_fits"])
#       impulse_parameters: (matrix genes x (NPARAM+1)*runs) +1 is value of 
#       objective of optimisation (i.e. sum of squares or weighted SS)
#       impulse_fits: (matrix genes x runs) model values for gene data
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

generate_background <- function(data_table, data_annotation, sigma_vec, n_iter = 100, imp_fit_genes,
                                control_timecourse = FALSE, control_name = NULL, no_of_clus = NULL,
                                n_rands = 50000, n_proc = 4){
  
  ### get number of pre- and fine clusters as input from the genes
  ### also exclude genes for which the model was not fitted due to almost no variation
  if(control_timecourse == TRUE){
    no_of_clus = no_of_clus[c(3,4,7,8,11,12)]
    data_table <- data_table[!(is.na(imp_fit_genes[[1]][,"beta"]) | is.na(imp_fit_genes[[3]][,"beta"]) | is.na(imp_fit_genes[[5]][,"beta"])),]
  } else if(control_timecourse == FALSE & is.null(control_name)){
    no_of_clus = no_of_clus[c(3,4)]
    data_table <- data_table[!(is.na(imp_fit_genes[[1]][,"beta"])),]
  }
  imp_fit_genes <- lapply(imp_fit_genes,function(x){x <- x[rownames(data_table),]})
  
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
  #   data_table: (Numeric matrix genes x samples) Expression values 
  #   data_annot: (Table samples x 2[time and condition]) 
  #       Co-variables for the samples including condition and time points.
  #       Time points must be numeric.
  #   sigma_vec:
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
  #   background: (vector round(n_randoms_tot/no_proc)) Empirical F-values
  #       under H0.
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  impulse_bg <- function(data_tab, data_annot, sigma_vec,  n_it = 100, imp_fit, ctrl_tk = FALSE, ctrl_name = NULL,
                         n_clust = NULL, n_randoms_tot = 50000, no_proc = 4){
    
    n_rand <- round(n_randoms_tot/no_proc)  ### divide randomizations by no. of nodes
    
    background = rep(0,n_rand)
    timep <- as.numeric(as.character(data_annot[colnames(data_tab),"Time"]))
    NPARAM = 6
    
    # 1.  n_rand samples are drawn with replacement from all genes.
    i_s <- sample(1:nrow(data_tab),n_rand, replace = TRUE)
    data_tab_rand <- data_tab[i_s,]
    rownames(data_tab_rand) <- paste("bg_",c(1:length(i_s)))
    
    ### ---> if control data is present
    ### DAVID: not control data for our case
    if(ctrl_tk == TRUE){
      data_list <- list(data_tab_rand, data_tab_rand[,!(data_annot$Condition %in% ctrl_name)],
                        data_tab_rand[,data_annot$Condition %in% ctrl_name])
      
      # 2.  Calculate residuals for random samples.
      y_i_null_s = (imp_fit[[2]][rownames(data_tab),])[i_s, as.character(data_annot[colnames(data_list[[1]]),"Time"])]
      calc_impulse_comp_data_tab_rand_case <- 
          (imp_fit[[4]][rownames(data_tab),])[i_s,as.character(data_annot[colnames(data_list[[2]]),"Time"])]
      
      ### if real control timecourse is present
      calc_impulse_comp_data_tab_rand_ctrl <- 
          (imp_fit[[6]][rownames(data_tab),])[i_s, as.character(data_annot[colnames(data_list[[3]]),"Time"])]
      RESD_1_s = cbind((data_list[[2]] - calc_impulse_comp_data_tab_rand_case),
                       (data_list[[3]] - calc_impulse_comp_data_tab_rand_ctrl))
      RESD_1_s = RESD_1_s[,colnames(data_list[[1]])]
      
      # 3.  Draw random residuals and add them to the null model:
      #     Sample n (timepoints) residuals with replacement from 
      #     n residuals of gene i to add to null model of gene i.
      Y_star_s = y_i_null_s + t(apply(RESD_1_s,
          1,function(x){x[sample(1:ncol(RESD_1_s), ncol(data_tab_rand),replace=TRUE)]}))
      colnames(Y_star_s) <- colnames(data_tab_rand)
      rownames(Y_star_s) <- paste("bg_",c(1:length(i_s)))
      
      #      write.table(data_tab_rand,"data_tab_rand.txt", sep="\t",quote=FALSE)
      #      write.table(Y_star_s,paste("Y_star_s_",round(runif(1)*10000),".txt", sep = ""), sep="\t",quote=FALSE, col.names = FALSE, row.names = FALSE)
      
      # 4.  Cluster data simulated under H0. Number of clusters as
      #     determined on observed data.
      tm_c <- system.time({
        clustering_results_background <- cluster_genes_for_impulse(Y_star_s, data_annot,
            ctrl_tk, ctrl_name, plot_clusters = TRUE, no_of_clusters = n_clust, n_genes = nrow(data_tab))
      })
      print(paste("Consumed time bg clustering: ",round(tm_c["elapsed"]/60,2)," min",sep=""))
      
      # 5. Fit impulse model to the clusters of data simulated under H0.
      tm_f <- system.time({
        impulse_fit_clusters_bg <- impulse_fit(clustering_results_background,data_annot,sigma_vec,
            n_it, control_timecourse = ctrl_tk, ctrl_name, fit_backg = TRUE)
      })
      print(paste("Consumed time bg clus fit: ",round(tm_f["elapsed"]/60,2)," min",sep=""))
      
      # 6.  Fit impulse model to the genes of data simulated under H0.
      tm_fg <- system.time({
        impulse_fit_genes_bg <- impulse_fit(Y_star_s, data_annot, sigma_vec, 1, ctrl_tk, ctrl_name,
            clustering_results_background, impulse_fit_clusters_bg, fit_backg = TRUE )
      })
      #      save(impulse_fit_genes_bg, file=file.path(getwd(),paste("impulse_fit_genes_bg_", round(runif(1)*10000),".RData",sep="")))
      print(paste("Consumed time bg gene fit: ",round(tm_fg["elapsed"]/60,2)," min",sep=""))
      
      #      if(is.null(imp_fit$control_model)){
      #        control_SSE = impulse_fit_genes_bg[[5]][rownames(Y_star_s),"SSE"]
      #      } else {
      #        control_SSE = apply((data_list[[3]]-rowMeans(data_list)),1,function(x){sum(x^2)})
      #      }
      ### calculate background
      #      background = (impulse_fit_genes_bg[[1]][rownames(Y_star_s),"SSE"] - (impulse_fit_genes_bg[[3]][rownames(Y_star_s),"SSE"] +
      #          control_SSE)) / (impulse_fit_genes_bg[[3]][rownames(Y_star_s),"SSE"] + control_SSE) 
      
      # 7.  Compute sum of squares for H0 and H1.
      RESD_1_s_b = cbind((Y_star_s[rownames(data_list[[1]]),colnames(data_list[[2]])] -
          impulse_fit_genes_bg[[4]][rownames(data_list[[1]]),
          as.character(data_annot[colnames(data_list[[2]]),"Time"])]),
          (Y_star_s[rownames(data_list[[1]]),colnames(data_list[[3]])] -
          impulse_fit_genes_bg[[6]][rownames(data_list[[1]]),
          as.character(data_annot[colnames(data_list[[3]]),"Time"])]))
      
      colnames(RESD_1_s_b) <- c(colnames(data_list[[2]]), colnames(data_list[[3]]))
      RESD_1_s_b = RESD_1_s_b[,colnames(data_list[[1]])]
      RESD_0_s_b = Y_star_s[rownames(data_list[[1]]),colnames(data_list[[1]])] -
          impulse_fit_genes_bg[[2]][rownames(data_list[[1]]),
          as.character(data_annot[colnames(data_list[[1]]),"Time"])]
      
      SS_1_s_b = apply(RESD_1_s_b,1,function(x){sum(x^2)})
      SS_0_s_b = apply(RESD_0_s_b,1,function(x){sum(x^2)})
      
      ### ---> if no control data is present
      ### DAVID: our case
    } else {
      # 2.  Calculate residuals for random samples.
      y_i_null_s = rowMeans(data_tab_rand)
      calc_impulse_comp_data_tab_rand <- (imp_fit[[2]][rownames(data_tab),])[i_s,as.character(data_annot[,"Time"])]
      #      calc_impulse_comp_data_tab_rand <- t(apply((imp_fit$impulse_parameters_case[rownames(data_tab),])[i_s,],1,
      #           function(x){calc_impulse_comp(x[1:NPARAM],timep)}))
      RESD_1_s = data_tab_rand - calc_impulse_comp_data_tab_rand
      
      # 3.  Draw random residuals and add them to the null model:
      #     Sample n (timepoints) residuals with replacement from 
      #     n residuals of gene i to add to null model of gene i.
      #     The null model data + residuals (data simulated under
      #     H0) are in Y_star_s
      ### DAVID: random residuals are drawn and added here: Change 3
      ### Fitting following after with funcitons defined above- dont have to change
      Y_star_s = t(y_i_null_s + apply(RESD_1_s,1,function(x){x[sample(1:ncol(RESD_1_s),
          ncol(data_tab_rand),replace=TRUE)]}))
      colnames(Y_star_s) <- colnames(data_tab_rand)
      
      # 4.  Cluster data simulated under H0. Number of clusters as
      #     determined on observed data.
      tm_c <- system.time({
          clustering_results_background <- cluster_genes_for_impulse(Y_star_s, data_annot,
              ctrl_tk, ctrl_name, plot_clusters = TRUE, no_of_clusters = n_clust, n_genes = nrow(data_tab))
      })
      print(paste("Consumed time bg clustering: ",round(tm_c["elapsed"]/60,2)," min",sep=""))
      
      # 5. Fit impulse model to the clusters of data simulated under H0.
      tm_f <- system.time({
        impulse_fit_clusters_bg <- impulse_fit(clustering_results_background, data_annot, sigma_vec, n_it,
                                               control_timecourse = ctrl_tk, ctrl_name, fit_backg = TRUE)
      })
      print(paste("Consumed time bg clust fit: ",round(tm_f["elapsed"]/60,2)," min",sep=""))
      
      # 6.  Fit impulse model to the genes of data simulated under H0.
      tm_fg <- system.time({
        impulse_fit_genes_bg <- impulse_fit(Y_star_s, data_annot, sigma_vec, 1, ctrl_tk, ctrl_name,
            clustering_results_background, impulse_fit_clusters_bg, fit_backg = TRUE )
      })
      #      save(impulse_fit_genes_bg, file=file.path(getwd(),paste("impulse_fit_genes_bg_", round(runif(1)*10000),".RData",sep="")))
      print(paste("Consumed time bg gene fit: ",round(tm_fg["elapsed"]/60,2)," min",sep=""))
      
      ### calculate background
      #      background = (apply(Y_star_s,1,function(x){sum((x - mean(x))^2)}) -
      #          impulse_fit_genes_bg[[1]][rownames(Y_star_s),"SSE"]) /
      #          impulse_fit_genes_bg[[1]][rownames(Y_star_s),"SSE"]
      
      # 7.  Compute sum of squares for H0 and H1.
      RESD_1_s_b = Y_star_s - impulse_fit_genes_bg[[2]][rownames(Y_star_s),
          as.character(data_annot[colnames(Y_star_s),"Time"])]
      SS_1_s_b = apply(RESD_1_s_b,1,function(x){sum(x^2)})
      RESD_0_s_b = Y_star_s - rowMeans(Y_star_s)
      SS_0_s_b = apply(RESD_0_s_b,1,function(x){sum(x^2)})
    }

    #background = (SS_0_s_b - SS_1_s_b) / SS_1_s_b
    # 8.  Compute F-scores.
    SS_1_s_b = apply(RESD_1_s_b,1,function(x){sum(x^2)})
    SS_0_s_b = apply(RESD_0_s_b,1,function(x){sum(x^2)})
    background = (SS_0_s_b - SS_1_s_b) / SS_1_s_b
    return(background)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~ 2. Parallelise background fitting ~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  #   define a new environment, which is used for parallelization
  my.env <- new.env()
  assign("n_rands", n_rands, envir = my.env)
  assign("n_proc", n_proc, envir = my.env)
  assign("data_table", data_table, envir = my.env)
  assign("data_annotation", data_annotation, envir = my.env)
  assign("n_iter",n_iter,envir = my.env)
  assign("control_timecourse", control_timecourse, envir = my.env)
  assign("control_name", control_name, envir = my.env)
  assign("imp_fit_genes", imp_fit_genes, envir = my.env)
  assign("no_of_clus", no_of_clus, envir = my.env)
  assign("calc_impulse_comp", calc_impulse_comp, envir = my.env)
  assign("cluster_genes_for_impulse", cluster_genes_for_impulse, envir = my.env)
  assign("impulse_fit", impulse_fit, envir = my.env)
  assign("two_impulses_WLS", two_impulses_WLS, envir = my.env)
  environment(impulse_bg) <- my.env
  
  #' @import parallel
  #' @import boot
  #
  #  ### set number of nodes to maximum - 1
  mc <- min(detectCores() - 1, n_proc)
  assign("mc", mc, envir = my.env)
  cl <- makeCluster(mc, outfile="cluster_out_random.txt")
  clusterExport(cl=cl, varlist=c("impulse_bg","n_rands","n_proc","data_table",
                                 "data_annotation","n_iter","control_timecourse","control_name",
                                 "imp_fit_genes","mc", "no_of_clus","calc_impulse_comp",
                                 "cluster_genes_for_impulse","impulse_fit", "two_impulses_WLS"), envir = my.env)
  junk <- clusterEvalQ(cl,library(boot))
  clusterSetRNGStream(cl,123)
  res <- clusterEvalQ(cl,impulse_bg(data_table, data_annotation, n_iter,
      imp_fit_genes, control_timecourse, control_name, no_of_clus, n_rands, mc))
  environment(res) <- my.env
  res2 <- do.call(c,res)
  stopCluster(cl)
  return(res2)
  
  ## if not paralellized
  #  res <- impulse_bg(data_table, data_annotation, n_iter, imp_fit_genes,
  #       control_timecourse, control_name, no_of_clus, n_rands, mc)
  #  return(res)
}
