#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Impulse model fit    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Fits impulse model to a timecourse dataset
### Function is split into 3 parts
###   Check data consistency
###   1. [Subfunction impulse_fit_gene_wise] Fit impulse model to a single gene
###   2. [Subfunction impulse_fit_matrix] Fit impulse model to matrix of genes.
###       Calls impulse_fit_gene_wise, calc_impulse_comp
###   3. ["Script-body"] Prepare data and fit model by calling impulse_fit_matrix.

### Developmental note:
### Three-way separation was chosen to use 'apply' instead of for-loops

# INPUT:
#   data_input: Two possibilities, depending on cluster or gene-wise fitting:
#       1.  Cluster. cluster_results as below.
#       2.  Gene-wise. (Numeric 3D array genes x samples x replicates).
#           Contains expression values or similar locus-specific read-outs.
#   data_annotation: (Table samples x 2[time and condition]) 
#       Co-variables for the samples including condition and time points.
#       Time points must be numeric.
#   weight_mat:
#   n_iter: (scalar) [Defaul 100] Number of iterations, which are performed 
#       to fit the impulse model to the clusters.
#   control_timecourse: (bool) [Default FALSE] Control time timecourse is 
#       part of the data set (TRUE) or not (FALSE).
#   control_name: (str) name of the control condition in annotation_table.
#   cluster_results: (list ["kmeans_clus","cluster_means", "n_pre_clus",
#       "fine_clus"] x number of runs [combined, case, control])
#       kmeans_clus: (vector number of genes) [1, number of clusters] 
#           Indicates assignment (value of function Kmeans) 
#       cluster_means: (matrix number of clusters x gene dimensions 
#           [number of samples per run, i.e. in case])
#           Centroids of fine (2nd) clustering.
#       n_pre_clus: (scalar) Number of clusters used in pre-clustering
#       n_fine_clusters: (scalar) Number of clusters used in fine-clustering
#       NOTE: This variable is only passed in gene-wise fitting, the results
#       from clustering are passed as data_input in cluster fitting. This is
#       done because cluster results serve as data to be fitted in cluster
#       fitting, while they are only auxillary in gene-wise fitting.
#   start_values:
#   fit_backg: (bool) [Defaul FALSE]
#   n_process: (scalar) [Default 4] number of processes, which can be used on 
#       the machine to run the background calculation in parallel
# OUTPUT:
#   results: (list runs*2 ["impulse_parameters","impulse_fits"])
#       impulse_parameters: (matrix genes x (NPARAM+1)*runs) +1 is value of 
#           objective of optimisation (i.e. sum of squares or weighted SS)
#       impulse_fits: (matrix genes x timepoints) model values for gene data

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#impulse_fit <- function(data_input, data_annotation, weight_mat, n_iter = 100,
#    control_timecourse = FALSE, control_name = NULL, cluster_results = NULL,
#    start_values = NULL, fit_backg = FALSE, n_proc = 4){
impulse_fit <- function(data_input, data_annotation, weight_mat, n_iter = 100,
  control_timecourse = FALSE, control_name = NULL, fit_backg = FALSE, n_proc = 4,
  dispersion_vector=NULL){
  
  ### Check dataset for consistency
  if(ncol(data_annotation) != 2 ||
      FALSE %in% (colnames(data_annotation) == c("Time","Condition"))){
    stop("Please use function 'annotation_preparation' to prepare
            the annotation table")
  }
  if(control_timecourse == TRUE & is.null(control_name)){
    stop("Please specify the name of the control in the column 'Condition'")
  }
  #if(control_timecourse == TRUE & (is.null(start_values) == FALSE)){
  #  start_values <- start_values[c(1,3,5)]
  #}
  #if(control_timecourse == FALSE & (is.null(control_name) == FALSE)){
  #  start_values <- start_values[c(1,3)]
  #}
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~ 1. Fit impulse model to gene ~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  ### Fits an impulse model to a single gene
  
  # INPUT:
  #   expression_values: (numeric 2D array timepoints x replicates) 
  #   data_annotation: (Table samples x 2[time and condition]) 
  #       Co-variables for the samples including condition and time points.
  #       Time points must be numeric.
  #   weight_vector:
  #   NPARAM: (scalar) Number of model parameters
  #   n_iter: (scalar) [Default 100] Number of iterations, which are performed 
  #       to fit the impulse model to the clusters.
  #   fit_to_clusters: (bool) [Default FALSE]
  #   start_val: 
  #   fit_backg: (bool) [Default FALSE]
  #   n_process: (scalar) [Default 4] number of processes, which can be used on 
  #       the machine to run the background calculation in parallel
  # OUTPUT:
  #   pvec_and_objective: (vector NPARAM+1) +1 is value of objective of optimisation 
  #       (i.e. sum of squares or weighted SS)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  impulse_fit_gene_wise <- function(expression_values, timepoints, weight_vector, 
    NPARAM, n_iters = 100, fit_bg = FALSE, dispersion_estimate=NULL,...){
    
    objective_fun <- cost_fun_logl_comp
    optim_method <- "nlminb"
    # DEVELOPMENTAL NOTE: cannot use optim if optim is used within optim?
    # set alpha with DESeq input and get rid of optim in cost function
    
    # If handed a list, i.e. not replicates
    if (is.vector(expression_values)==TRUE){
      # Transform into 2D array
      expression_values <- array(expression_values,c(length(expression_values),1))
      # Do intitial guesses based on centroids
      expression_means <- expression_values
    } else {
      # Do initial guesses based on mean over replicates:
      expression_means <- apply(expression_values,1,mean)
    }
    
    num_points <- length(expression_means)
    max_value = max(expression_values)
    min_value = min(expression_values)
    max_mean <- max(expression_means)
    min_mean <- min(expression_means)
    max_mean_ind <- match(max_mean,expression_means)
    min_mean_ind <- match(min_mean,expression_means)
    max_middle_mean <- max(expression_means[2:(num_points-1)])
    min_middle_mean <- min(expression_means[2:(num_points-1)])
    max_middle_val <- max(expression_values[2:(num_points-1),])
    min_middle_val <- min(expression_values[2:(num_points-1),])
    # +1 to push indicices up from middle stretch to entire window (first is omitted here)
    max_middle_mean_ind <- match(max_middle_mean,expression_means[2:(num_points-1)]) + 1
    min_middle_mean_ind <- match(min_middle_mean,expression_means[2:(num_points-1)]) + 1
    # Gradients between neighbouring points
    neighbour_grad <- unlist( lapply(c(1:(num_points-1)),function(x){
      (expression_means[x+1]-expression_means[x])/(timepoints[x+1]-timepoints[x])}) )
    if(max_mean < 1.0001){max_mean <- 1.0001}
    if(min_mean < 1.0001){min_mean <- 1.0001}
    if(max_value < 1.0001){max_value <- 1.0001}
    if(min_value < 1.0001){min_value <- 1.0001}
    if(max_middle_mean < 1.0001){max_middle_mean <- 1.0001}
    if(min_middle_mean < 1.0001){min_middle_mean <- 1.0001}
    
    tmm2 <- system.time({
      
      # NOT working in log 2 space
      # 1. Peak:
      #if(max_middle_mean > 2*expression_means[1] &&
      #    max_middle_mean > 2*expression_means[num_points] &&
      #    ((expression_means[1] > 10 && expression_means[num_points] > 10) || 
      #    (max_middle_mean > 4*expression_means[1] ||
      #    max_middle_mean > 4*expression_means[num_points]))){
      # Beta: Has to be negative
      # Theta1: Low
      # Theta2: High
      # Theta3: Low
      # t1: Around first observed inflexion point
      # t2: Around second observed inflexion point
      #guess1 <- c(1,expression_means[1],max_mean,expression_means[num_points],timepoints[(max_middle_mean_ind-1)],timepoints[max_middle_mean_ind])
      lower_inflexion_point_ind <- match(max(neighbour_grad[1:(max_middle_mean_ind-1)]), neighbour_grad[1:(max_middle_mean_ind-1)])
      upper_inflexion_point_ind <- max_middle_mean_ind - 1 + match(min(neighbour_grad[max_middle_mean_ind:length(neighbour_grad)]), neighbour_grad[max_middle_mean_ind:length(neighbour_grad)])
      upper_h0 <- min( max(expression_values[1:(max_middle_mean_ind-1),]), max_middle_mean )
      upper_h2 <- min( max(expression_values[(max_middle_mean_ind+1):num_points,]), max_middle_mean )
      guess1 <- c(1,expression_means[1],max_middle_mean,expression_means[num_points],
        (timepoints[lower_inflexion_point_ind]+timepoints[lower_inflexion_point_ind+1])/2,
        (timepoints[upper_inflexion_point_ind]+timepoints[upper_inflexion_point_ind+1])/2 )
      lower_b <- c(0.01,1.0001,max(expression_means[1],expression_means[num_points]),1.0001,(timepoints[1]+timepoints[2])/3,timepoints[max_middle_mean_ind])
      upper_b <- c(100,upper_h0,max_middle_val,upper_h2,timepoints[max_middle_mean_ind],timepoints[num_points])
      # check that guess fullfills bounds
      guess1[guess1 <= lower_b] <- lower_b[guess1 <= lower_b] + 0.0001
      guess1[guess1 >= upper_b] <- upper_b[guess1 >= upper_b] - 0.0001
      # check validity of bounds
      upper_b[upper_b <= lower_b] <- lower_b[upper_b <= lower_b] + 1
      if("optim" %in% optim_method){
        fit_peak1 <- unlist(optim(guess1, objective_fun, x_vec=timepoints,
          y_mat=expression_values, weight_vec=weight_vector, disp_est=dispersion_estimate,
          method="L-BFGS-B",lower=lower_b, upper=upper_b)[c("par","value")])
      }
      if("nlminb" %in% optim_method){
        fit_peak2 <- unlist(nlminb(guess1, objective_fun, x_vec=timepoints,
          y_mat=expression_values, weight_vec=weight_vector, disp_est=dispersion_estimate,
          lower=lower_b, upper=upper_b)[c("par","objective")])
      }
      # 2. Valley
      #} else if(min_middle_mean < expression_means[1]/2 &&
      #    min_middle_mean < expression_means[num_points]/2 &&
      #    expression_means[1] > 10 && expression_means[num_points] > 10){
      # Beta: Has to be negative
      # Theta1: High
      # Theta2: Low
      # Theta3: High
      # t1: Around first observed inflexion point
      # t2: Around second observed inflexion point
      #guess1 <- c(1,expression_means[1],min_mean,expression_means[num_points],timepoints[min_middle_mean_ind-1],timepoints[min_middle_mean_ind+1])
      lower_inflexion_point_ind <- match(min(neighbour_grad[1:(min_middle_mean_ind-1)]), neighbour_grad[1:(min_middle_mean_ind-1)])
      upper_inflexion_point_ind <- min_middle_mean_ind - 1 + match(max(neighbour_grad[min_middle_mean_ind:(num_points-1)]), neighbour_grad[min_middle_mean_ind:(num_points-1)])
      upper_h0 <- max(expression_values[1:(min_middle_mean_ind-1),])
      upper_h2 <- max(expression_values[(min_middle_mean_ind+1):num_points,])
      guess1 <- c(-1,expression_means[1],min_middle_mean,expression_means[num_points],
        (timepoints[lower_inflexion_point_ind]+timepoints[lower_inflexion_point_ind+1])/2,
        (timepoints[upper_inflexion_point_ind]+timepoints[upper_inflexion_point_ind+1])/2 )
      lower_b <- c(-0.01,min_middle_mean,0.0001,min_middle_mean,(timepoints[1]+timepoints[2])/3,timepoints[min_middle_mean_ind])
      upper_b <- c(-100,upper_h0,min(expression_means[1],expression_means[num_points]),upper_h2,timepoints[min_middle_mean_ind],timepoints[num_points])          
      # check that guess fullfills bounds
      guess1[guess1 <= lower_b] <- lower_b[guess1 <= lower_b] + 0.0001
      guess1[guess1 >= upper_b] <- upper_b[guess1 >= upper_b] - 0.0001
      # check validity of bounds
      upper_b[upper_b <= lower_b] <- lower_b[upper_b <= lower_b] + 1
      if("optim" %in% optim_method){
        fit_valley1 <- unlist(optim(guess1, objective_fun, x_vec=timepoints,
          y_mat=expression_values, weight_vec=weight_vector, disp_est=dispersion_estimate,
          method="L-BFGS-B",lower=lower_b, upper=upper_b)[c("par","value")])
      }
      if("nlminb" %in% optim_method){
        fit_valley2 <- unlist(nlminb(guess1, objective_fun, x_vec=timepoints,
          y_mat=expression_values, weight_vec=weight_vector, disp_est=dispersion_estimate,
          lower=lower_b, upper=upper_b)[c("par","objective")])
      }
      # 3. Up
      #} else if(expression_means[1] < expression_means[num_points]){
      #guess1 <- c(10,min_mean,max_mean,min_mean,timepoints[2],mx_tm+timepoints[2])
      # Beta: Has to be positive
      # Theta1: Low
      # Theta2: High
      # Theta3: Arbitrary (outside reach - defined through t2)
      # t1: Around observed inflexion point
      # t2: Should be well outside reach so that model doesn't use the 2nd sigmoid.
      inflexion_point_ind <- match(max(neighbour_grad), neighbour_grad)
      guess1 <- c(1,min_mean,max_mean,min_mean,
        (timepoints[inflexion_point_ind]+timepoints[inflexion_point_ind+1])/2,
        2.5*timepoints[num_points])
      lower_b <- c(0.01,0.0001,max_middle_mean,0.0001,(timepoints[1]+timepoints[2])/3,2*timepoints[num_points])
      upper_b <- c(100,max_middle_mean,max_value,max_middle_mean,timepoints[num_points],3*timepoints[num_points])
      # check that guess fullfills bounds
      guess1[guess1 <= lower_b] <- lower_b[guess1 <= lower_b] + 0.0001
      guess1[guess1 >= upper_b] <- upper_b[guess1 >= upper_b] - 0.0001
      # check validity of bounds
      upper_b[upper_b <= lower_b] <- lower_b[upper_b <= lower_b] + 1
      if("optim" %in% optim_method){
        fit_up1 <- unlist(optim(guess1, objective_fun, x_vec=timepoints,
          y_mat=expression_values, weight_vec=weight_vector, disp_est=dispersion_estimate,
          method="L-BFGS-B",lower=lower_b, upper=upper_b)[c("par","value")])
      }
      if("nlminb" %in% optim_method){
        fit_up2 <- unlist(nlminb(guess1, objective_fun, x_vec=timepoints,
          y_mat=expression_values, weight_vec=weight_vector, disp_est=dispersion_estimate,
          lower=lower_b, upper=upper_b)[c("par","objective")])
      }
      # 4. Down
      #} else if(expression_means[1] >= expression_means[num_points]){
      #guess1 <- c(10,max_mean,min_mean,max_mean,timepoints[2],mx_tm+timepoints[2])
      # Beta: Has to be negative
      # Theta1: High
      # Theta2: Low
      # Theta3: Arbitrary (outside reach - defined through t2)
      # t1: Around observed inflexion point
      # t2: Should be well outside reach so that model doesn't use the 2nd sigmoid.
      inflexion_point_ind <- match(min(neighbour_grad), neighbour_grad)
      guess1 <- c(-1,max_mean,min_mean,max_mean,
        (timepoints[inflexion_point_ind]+timepoints[inflexion_point_ind+1])/2,
        2.5*timepoints[num_points])
      lower_b <- c(-0.01,min_middle_mean,0.0001,0.0001,(timepoints[1]+timepoints[2])/3,2*timepoints[num_points])
      upper_b <- c(-100,max_value,max_middle_mean,max_value,timepoints[num_points],3*timepoints[num_points])
      # check that guess fullfills bounds
      guess1[guess1 <= lower_b] <- lower_b[guess1 <= lower_b] + 0.0001
      guess1[guess1 >= upper_b] <- upper_b[guess1 >= upper_b] - 0.0001
      # check validity of bounds
      upper_b[upper_b <= lower_b] <- lower_b[upper_b <= lower_b] + 1
      if("optim" %in% optim_method){
        fit_down1 <- unlist(optim(guess1, objective_fun, x_vec=timepoints,
          y_mat=expression_values, weight_vec=weight_vector, disp_est=dispersion_estimate,
          method="L-BFGS-B",lower=lower_b, upper=upper_b)[c("par","value")])
      }
      if("nlminb" %in% optim_method){
        fit_down2 <- unlist(nlminb(guess1, objective_fun, x_vec=timepoints,
          y_mat=expression_values, weight_vec=weight_vector, disp_est=dispersion_estimate,
          lower=lower_b, upper=upper_b)[c("par","objective")])
      }
      
      if(("optim" %in% optim_method) && ("nlminb" %in% optim_method)){
        fmin_outs <- cbind(fit_peak1,fit_peak2,fit_valley1,fit_valley2,fit_up1,fit_up2,fit_down1,fit_down2)
      }
      if(("optim" %in% optim_method) && !("nlminb" %in% optim_method)){
        fmin_outs <- cbind(fit_peak1,fit_valley1,fit_up1,fit_down1)
      }
      if(!("optim" %in% optim_method) && ("nlminb" %in% optim_method)){
        fmin_outs <- cbind(fit_peak2,fit_valley2,fit_up2,fit_down2)
      } 
      
      # Name row containing value of objective as value 
      # This name differs between optimisation functions
      rownames(fmin_outs) <- c(rownames(fmin_outs[1:(nrow(fmin_outs)-1),]),"value")
    })
    
    if(is.null(dim(fmin_outs[ ,fmin_outs["value",] == min(fmin_outs["value",])])) == TRUE){
      pvec_and_objective = fmin_outs[ ,fmin_outs["value",] == min(fmin_outs["value",])]
    } else {
      # If two or more randomization have the same value of the objective,
      # choose the first one
      pvec_and_objective = fmin_outs[ ,fmin_outs["value",] ==
          min(fmin_outs["value",])][,1]
    }
    
    return(pvec_and_objective)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~ 2. Fit impulse model to matrix ~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  ### Fits an impulse model to all genes of a dataset
  
  # INPUT:
  #   data_arr: (Numeric 3D array genes x samples x replicates).
  #       Contains expression values or similar locus-specific read-outs.
  #   timepoints: (Table samples x 2 [time and condition]) 
  #       Co-variables for the samples including condition and time points.
  #       Time points must be numeric.
  #   weight_mat:
  #   NPARAM: 
  #   n_it: (scalar) [Defaul 100] Number of iterations, which are performed 
  #       to fit the impulse model to the clusters.
  #   ctrl_tc
  #   ctrl: (bool)
  #   fit_to_clus: (bool) [Default FALSE]
  #   start_val: 
  #   fit_bg: (bool) [Default FALSE]
  #   n_process: (scalar) [Default 4] number of processes, which can be used on 
  #       the machine to run the background calculation in parallel
  #   ...
  # OUTPUT:
  #   res: (list 2 ["impulse_parameters","impulse_fits"])
  #       impulse_parameters: (matrix genes x (NPARAM+1)*runs) +1 is value of 
  #           objective of optimisation (i.e. sum of squares or weighted SS)
  #       impulse_fits: (matrix genes x timepoints) model values for gene data
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  impulse_fit_matrix <- function(data_arr, timepoints, weight_mat, n_it = 100,
    ctrl_tc = FALSE, ctrl = NULL, fit_bg = FALSE, n_process = 4, dispersion_vec=NULL, ...){
    
    NPARAM = 6 # Number of model parameters
    
    mc <- min(detectCores() - 1, n_process)
    
    # Use parallelisation if number of genes/centroids to fit is large
    if(fit_bg == FALSE && nrow(data_arr) > max(2*mc,10)){
      # Use parallisation
      # Define partitioning of genes onto nodes: ind_list
      ind_list = list()
      bord = floor(nrow(data_arr)/mc)
      for (i in 1:mc){
        ind_list[[i]] <- ((i-1)*bord+1):(i*bord)
      }
      if(mc*bord < nrow(data_arr)){
        # Add remaining genes to last node
        ind_list[[mc]] <-  c(ind_list[[mc]],(mc*bord+1):nrow(data_arr))
      }
      
      cl <- makeCluster(mc, outfile = "clus_out_impulse_fit.txt")
      my.env <- new.env()
      assign("data_arr", data_arr, envir = my.env)
      assign("weight_mat", weight_mat, envir = my.env)
      assign("calc_impulse_comp", calc_impulse_comp, envir = my.env)
      assign("n_it", n_it, envir = my.env)
      #assign("cluster_genes_for_impulse", cluster_genes_for_impulse,
      #  envir = my.env)
      assign("impulse_fit", impulse_fit, envir = my.env)
      assign("cost_fun_logl_comp", cost_fun_logl_comp, envir = my.env)
      assign("cost_fun_WLS_comp", cost_fun_WLS_comp, envir = my.env)
      assign("cost_fun_OLS_comp", cost_fun_OLS_comp, envir = my.env)
      assign("impulse_fit_gene_wise",impulse_fit_gene_wise, envir = my.env)
      assign("impulse_fit_matrix",impulse_fit_matrix, envir = my.env)
      assign("timepoints", timepoints, envir = my.env)
      #assign("fit_to_clus", fit_to_clus, envir = my.env)
      #assign("start_val", start_val, envir = my.env)
      assign("fit_bg", fit_bg, envir = my.env)
      assign("NPARAM", NPARAM, envir = my.env)
      assign("dispersion_vec", dispersion_vec, envir = my.env)
      
      #clusterExport(cl=cl, varlist=c("data_arr","weight_mat",
      #  "calc_impulse_comp","n_it","cluster_genes_for_impulse","impulse_fit",
      #  "cost_fun_WLS_comp", "cost_fun_OLS_comp", "impulse_fit_gene_wise",
      #  "impulse_fit_matrix", "timepoints","fit_to_clus","start_val",
      #  "fit_bg","NPARAM"), envir = my.env)
      clusterExport(cl=cl, varlist=c("data_arr","weight_mat",
        "calc_impulse_comp","n_it","impulse_fit","cost_fun_logl_comp",
        "cost_fun_WLS_comp", "cost_fun_OLS_comp", "impulse_fit_gene_wise",
        "impulse_fit_matrix", "timepoints", "fit_bg","NPARAM","dispersion_vec"), envir = my.env)
      
      junk <- clusterEvalQ(cl,library(MASS))
      # Fit impulse model to each gene of matrix and get impulse parameters:
      # clusterApply runs the function impulse_fit_gene_wise
      # The input data are distributed to nodes by ind_list partitioning
      resi <- clusterApply(cl, 1:length(ind_list), function(z){t(sapply(ind_list[[z]],
        function(x){impulse_fit_gene_wise(data_arr[x,,],timepoints,weight_mat[x,],
          NPARAM,n_it, fit_bg, dispersion_vec[x])}))})
      # Give output rownames again, which are lost above
      for(i in 1:length(ind_list)){
        rownames(resi[[i]]) <- rownames(data_arr[ind_list[[i]],,])
      }
      # Concatenate the output objects of each node
      ### DAVID bug fix: = to <-
      resmat <- do.call(rbind,resi)
      stopCluster(cl)
      # Parallelisation ends here
      
      # Do not use parallelisation if number of genes/centroids to fit is small
    } else {
      # Fit impulse model to each gene of matrix and get impulse parameters
      ind_list_all <- 1:nrow(data_arr)
      resmat <- lapply(ind_list_all,function(x){
        impulse_fit_gene_wise(data_arr[x,,],timepoints,weight_mat[x,],
          NPARAM,n_it,fit_bg, dispersion_vec[x])})
      resmat_temp <- array(0,c(length(resmat),length(resmat[[1]])))
      for (i in 1:length(resmat)){
        resmat_temp[i,] <- resmat[[i]]
      }
      resmat <- resmat_temp
    }   
    
    # Use obtained impulse parameters to calculate impulse fit values
    colnames(resmat) <- c("beta","h0","h1","h2","t1","t2","objective")
    if(nrow(resmat) == 1){      # if matrix contains only one gene
      resmat2 <- as.data.frame(t(calc_impulse_comp(resmat[,1:NPARAM],
        unique(sort(timepoints)))),row.names=rownames(resmat))
      colnames(resmat2) = unique(sort(timepoints))
    } else {                    # if matrix contains > 1 genes
      resmat2 <- t(apply(resmat[,1:NPARAM],1,function(x){calc_impulse_comp(x,
        unique(sort(timepoints)))}))
      colnames(resmat2) = unique(sort(timepoints))
    } 
    
    # Report results.
    res <- list()
    rownames(resmat) <- rownames(data_arr)
    rownames(resmat2) <- rownames(data_arr)
    res[[1]] <- resmat    # "beta","h0","h1","h2","t1","t2","value" [x3 if with control]
    res[[2]] <- resmat2   # values of impulse model on input genes
    names(res) <- c("impulse_parameters","impulse_fits")
    return(res)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~ 3. Prepare data for impulse model fit ~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  ### assume only one timecourse as default and perform only 1 run
  g_names = rownames(data_input)
  runs = 1
  label = "case"
  results = list()
  if(control_timecourse == TRUE){
    runs = 3
    label = c("combined","case","control")
  }
  
  # data_input is defined in the following as:
  #   1) A list of three data 3D arrays combined, case and control
  #   2) A list of a single data 3D array if only have case

  # Perform 3 runs (combined, case and control) if control timecourse is present
  # datX contains expression data of genes associates with cluster
  # datX_n contains expression data of genes NOT in cluster (low variation)
  # DAVID dat are defined as data.frames, which i cannot maintain
  # it doesn seem to have relevance downstream though, change to 3D array
  if(control_timecourse == TRUE){
    # Note on variables used below:
    # cluster_results[[1,5,9]] are kmeans_clus for [combined, case, control]
    
    # Combined: Expression values of both case and control (all columns).
    dat1 = data_input
    dat1_n = data_input[!(rownames(data_input) %in% rownames(dat1)),,] 
    dat1_n = dat1_n[,colnames(dat1),]
    
    # Case: Expression values of case.
    dat2 = data_input[,!(data_annotation$Condition %in% control_name),]
    dat2_n = data_input[!(rownames(data_input) %in% rownames(dat2)),,]
    dat2_n = dat2_n[,colnames(dat2),]
    
    # Control: Expression values of control.  
    dat3 = data_input[,(data_annotation$Condition %in% control_name),]
    dat3_n = data_input[!(rownames(data_input) %in% rownames(dat3)),,]
    dat3_n = dat3_n[,colnames(dat3),]
    
    # Formatting into 3D breaks if only single genes are selected:
    # (run 1)
    # If exactly one gene was excluded:
    if(dim(dat1)[1]==dim(data_input)[1]-1){ 
      dat1_n_temp <- array(NA,c(1,dim(data_input)[2],dim(data_input)[3]))
      for(z in 1:dim(data_input)[3]){
        dat1_n_temp[1,,z] <- dat1_n[,z]
      }
      dat1_n <- dat1_n_temp
      rownames(dat1_n) <- rownames(data_input)[!(rownames(data_input) %in% rownames(dat1))]
      colnames(dat1_n) <- colnames(data_input)
    }
    # If exactly one gene was included:
    if(dim(dat1)[1]==1){ 
      dat1_temp <- array(NA,c(1,dim(data_input)[2],dim(data_input)[3]))
      for(z in 1:dim(data_input)[3]){
        dat1_temp[1,,z] <- dat1[,z]
      }
      dat1 <- dat1_temp
      rownames(dat1) <- rownames(data_input)[(rownames(data_input) %in% rownames(dat1))]
      colnames(dat1) <- colnames(data_input)
    }
    # (run 2)
    # If exactly one gene was excluded:
    if(dim(dat2)[1]==dim(data_input)[1]-1){ 
      dat2_n_temp <- array(NA,c(1,dim(data_input)[2],dim(data_input)[3]))
      for(z in 1:dim(data_input)[3]){
        dat2_n_temp[1,,z] <- dat2_n[,z]
      }
      dat2_n <- dat2_n_temp
      rownames(dat2_n) <- rownames(data_input)[!(rownames(data_input) %in% rownames(dat2))]
      colnames(dat2_n) <- colnames(data_input[,!(data_annotation$Condition %in% control_name),])
    }
    # If exactly one gene was included:
    if(dim(dat2)[1]==1){ 
      dat2_temp <- array(NA,c(1,dim(data_input)[2],dim(data_input)[3]))
      for(z in 1:dim(data_input)[3]){
        dat2_temp[1,,z] <- dat2[,z]
      }
      dat2 <- dat2_temp
      rownames(dat2) <- rownames(data_input)[(rownames(data_input) %in% rownames(dat2))]
      colnames(dat2) <- colnames(data_input[,!(data_annotation$Condition %in% control_name),])
    }
    # (run 3)
    # If exactly one gene was excluded:
    if(dim(dat2)[1]==dim(data_input)[1]-1){ 
      dat3_n_temp <- array(NA,c(1,dim(data_input)[2],dim(data_input)[3]))
      for(z in 1:dim(data_input)[3]){
        dat3_n_temp[1,,z] <- dat3_n[,z]
      }
      dat3_n <- dat3_n_temp
      rownames(dat3_n) <- rownames(data_input)[!(rownames(data_input) %in% rownames(dat3))]
      colnames(dat3_n) <- colnames(data_input[,(data_annotation$Condition %in% control_name),])
    }
    # If exactly one gene was included:
    if(dim(dat3)[1]==1){ 
      dat3_temp <- array(NA,c(1,dim(data_input)[2],dim(data_input)[3]))
      for(z in 1:dim(data_input)[3]){
        dat3_temp[1,,z] <- dat3[,z]
      }
      dat3 <- dat3_temp
      rownames(dat3) <- rownames(data_input)[(rownames(data_input) %in% rownames(dat3))]
      colnames(dat3) <- colnames(data_input[,(data_annotation$Condition %in% control_name),])
    }
    
    dat_ind <- c("dat1_n","dat2_n", "dat3_n")
    data_input <- list(dat1, dat2, dat3)
    
  } else if(control_timecourse == FALSE){
    dat1 = data_input[,,]
    dat1_n = data_input[!(rownames(data_input) %in% rownames(dat1)),,]
    # Formatting into 3D breaks if only single gene is selected:
    # If exactly one gene was excluded:
    if(dim(dat1)[1]==dim(data_input)[1]-1){ 
      dat1_n_temp <- array(NA,c(1,dim(data_input)[2],dim(data_input)[3]))
      for(z in 1:dim(data_input)[3]){
        dat1_n_temp[1,,z] <- dat1_n[,z]
      }
      dat1_n <- dat1_n_temp
      rownames(dat1_n) <- rownames(data_input)[!(rownames(data_input) %in% rownames(dat1))]
      colnames(dat1_n) <- colnames(data_input)
    }
    # If exactly one gene was included:
    if(dim(dat1)[1]==1){ 
      dat1_temp <- array(NA,c(1,dim(data_input)[2],dim(data_input)[3]))
      for(z in 1:dim(data_input)[3]){
        dat1_temp[1,,z] <- dat1[,z]
      }
      dat1 <- dat1_temp
      rownames(dat1) <- rownames(data_input)[(rownames(data_input) %in% rownames(dat1))]
      colnames(dat1) <- colnames(data_input)
    }
    
    dat_ind <- c("dat1_n")
    data_input <- list(dat1)
  }
  
  # Fitting for different runs
  for (c_runs in 1:runs){
    imp_res <- NULL
    
    imp_res <- impulse_fit_matrix(data_input[[c_runs]],
      as.numeric(as.character(data_annotation[colnames(data_input[[c_runs]]),"Time"])), 
      weight_mat = weight_mat,n_it = n_iter, ctrl_tc = control_timecourse, 
      ctrl = control_name, fit_bg = fit_backg, n_process = n_proc, dispersion_vec=dispersion_vector)
    
    # Use 'get' to pull set of not clustered genes in current run (from current environment)
    # If there is >0 genes not clustered in current run
    if(nrow(get(dat_ind[c_runs])) != 0){
      # Report results obtained so far
      tump1 <- imp_res[[1]]
      tump2 <- imp_res[[2]]
      
      # Set model value for non clustered genes to mean, this is the null hypothesis
      tmpp1 <- t(apply(get(dat_ind[c_runs]),1,function(x){rep(mean(x),ncol(tump2))}))
      colnames(tmpp1) <- colnames(tump2)
      tump2 <- rbind(tump2, tmpp1)
      
      # Fill parameter fits of non clustered genes with NA and sum of squares
      # position with squared deviation to mean (H0)
      # DAVID: leave as NA for now, do these still exist?
      datn_WSS <- rep(NA,length(rownames(get(dat_ind[c_runs]))))
      names(datn_WSS) <- rownames(get(dat_ind[c_runs]))
      #for(iExclGenes in rownames(get(dat_ind[c_runs]))){
      #  datn_WSS[[iExclGenes]] <-  sum( apply(get(dat_ind[c_runs])[iExclGenes,,],2,function(y_vec){
      #    t((tmpp1[iExclGenes,] - y_vec)^2) %*% (1/weight_mat[iExclGenes,]^2)}))
      #}
      tmpp2 <-  cbind( matrix(NA,nrow(get(dat_ind[c_runs])),ncol(tump1)-1), datn_WSS)
      
      rownames(tmpp2) <- rownames(get(dat_ind[c_runs]))
      tump1 <- rbind(tump1, tmpp2)
      
      # Report results previously obtained through fitting and results added for
      # non-clustered genes.
      tempi1 <- tump1
      tempi2 <- tump2
      # If all genes were clustered and fitted
    } else {
      # Simply report results, all genes are covered.
      tempi1 <- imp_res[[1]]
      tempi2 <- imp_res[[2]]
    }
    
    imp_res$impulse_parameters <- tempi1[g_names,]
    imp_res$impulse_fits <- tempi2[g_names,]
    
    names(imp_res) <- paste(c("impulse_parameters","impulse_fits"),label[c_runs],sep="_")
    results[[names(imp_res)[1]]] <- imp_res[[1]]
    results[[names(imp_res)[2]]] <- imp_res[[2]]
  }
  
  return(results)
}

