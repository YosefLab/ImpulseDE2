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

impulse_fit <- function(data_input, data_annotation,
  control_timecourse = FALSE, control_name = NULL, n_proc = 4,
  dispersion_vector=NULL,NPARAM=6){
  
  # Check dataset for consistency
  if(ncol(data_annotation) != 2 ||
      FALSE %in% (colnames(data_annotation) == c("Time","Condition"))){
    stop("Please use function 'annotation_preparation' to prepare
            the annotation table")
  }
  if(control_timecourse == TRUE & is.null(control_name)){
    stop("Please specify the name of the control in the column 'Condition'")
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~ 1. Fit impulse model to gene ~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  ### Fits an impulse model to a single gene. Fitted are 1. peak, 2. valley,
  ### 3. up, 4. down and 5. mean model, the best fit is kept.
  
  # INPUT:
  #   expression_values: (numeric 2D array timepoints x replicates) 
  #   data_annotation: (Table samples x 2[time and condition]) 
  #       Co-variables for the samples including condition and time points.
  #       Time points must be numeric.
  #   NPARAM: (scalar) Number of model parameters
  #   n_iter: (scalar) [Default 100] Number of iterations, which are performed 
  #       to fit the impulse model to the clusters.
  #   n_process: (scalar) [Default 4] number of processes, which can be used on 
  #       the machine to run the background calculation in parallel
  # OUTPUT:
  #   pvec_and_objective: (vector NPARAM+1) +1 is value of objective of optimisation 
  #       (i.e. sum of squares or weighted SS)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  impulse_fit_gene_wise <- function(expression_values, timepts, 
    NPARAM=6, MAXIT=1000, dispersion_estimate=NULL,...){
    
    optim_method <- "optim"
    #optim_method <- "nlminb"
    #optim_method <- c("optim","nlminb")
    
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
    max_middle_mean <- max(expression_means[2:(num_points-1)])
    min_middle_mean <- min(expression_means[2:(num_points-1)])
    max_middle_val <- max(expression_values[2:(num_points-1),])
    min_middle_val <- min(expression_values[2:(num_points-1),])
    # +1 to push indicices up from middle stretch to entire window (first is omitted here)
    max_middle_mean_ind <- match(max_middle_mean,expression_means[2:(num_points-1)]) + 1
    min_middle_mean_ind <- match(min_middle_mean,expression_means[2:(num_points-1)]) + 1
    # Gradients between neighbouring points
    neighbour_grad <- unlist( lapply(c(1:(num_points-1)),function(x){
      (expression_means[x+1]-expression_means[x])/(timepts[x+1]-timepts[x])}) )
    
    tmm2 <- system.time({
      
      # 1. Alternative model: Peak
      # Beta: Has to be negative
      # Theta1: Low
      # Theta2: High
      # Theta3: Low
      # t1: Around first observed inflexion point
      # t2: Around second observed inflexion point
      lower_inflexion_point_ind <- match(max(neighbour_grad[1:(max_middle_mean_ind-1)]), neighbour_grad[1:(max_middle_mean_ind-1)])
      upper_inflexion_point_ind <- max_middle_mean_ind - 1 + match(min(neighbour_grad[max_middle_mean_ind:length(neighbour_grad)]), neighbour_grad[max_middle_mean_ind:length(neighbour_grad)])
      par_guess <- c(1,log(expression_means[1]+1),log(max_middle_mean+1),log(expression_means[num_points]+1),
        (timepts[lower_inflexion_point_ind]+timepts[lower_inflexion_point_ind+1])/2,
        (timepts[upper_inflexion_point_ind]+timepts[upper_inflexion_point_ind+1])/2)
      #par_guess <- c(1,1,2,1,timepts[2],timepts[3])
      
      if("optim" %in% optim_method){
        fit_peak1 <- unlist( optim(par=par_guess, fn=cost_fun_logl_comp, x_vec=timepts,
          y_mat=expression_values, disp_est=dispersion_estimate,
          method="BFGS", control=list(maxit=MAXIT,fnscale=-1))[c("par","value","convergence")] )
      }
      if("nlminb" %in% optim_method){
        stop("switched optimisation to maximisation which is not used here")
        fit_peak2 <- unlist(nlminb(start=par_guess, objective=cost_fun_logl_comp, x_vec=timepts,
          y_mat=expression_values, disp_est=dispersion_estimate,
          lower=lower_b, upper=upper_b)[c("par","objective")])
      }
      
      # 2. Alternative model: Valley
      # Beta: Has to be negative
      # Theta1: High
      # Theta2: Low
      # Theta3: High
      # t1: Around first observed inflexion point
      # t2: Around second observed inflexion point
      lower_inflexion_point_ind <- match(min(neighbour_grad[1:(min_middle_mean_ind-1)]), neighbour_grad[1:(min_middle_mean_ind-1)])
      upper_inflexion_point_ind <- min_middle_mean_ind - 1 + match(max(neighbour_grad[min_middle_mean_ind:(num_points-1)]), neighbour_grad[min_middle_mean_ind:(num_points-1)])
      par_guess <- c(1,log(expression_means[1]+1),log(min_middle_mean+1),log(expression_means[num_points]+1),
        (timepts[lower_inflexion_point_ind]+timepts[lower_inflexion_point_ind+1])/2,
        (timepts[upper_inflexion_point_ind]+timepts[upper_inflexion_point_ind+1])/2 )
      #par_guess <- c(1,1,2,1,timepts[2],timepts[3])
      
      if("optim" %in% optim_method){
        fit_valley1 <- unlist( optim(par=par_guess, fn=cost_fun_logl_comp, x_vec=timepts,
          y_mat=expression_values, disp_est=dispersion_estimate,
          method="BFGS", control=list(maxit=MAXIT,fnscale=-1))[c("par","value","convergence")] )
      }
      if("nlminb" %in% optim_method){
        stop("switched optimisation to maximisation which is not used here")
        fit_valley2 <- unlist(nlminb(start=par_guess, objective=cost_fun_logl_comp, x_vec=timepts,
          y_mat=expression_values, disp_est=dispersion_estimate,
          lower=lower_b, upper=upper_b)[c("par","objective")])
      }
      
      # 3. Null model: Single mean
      # Parameter estimate: Overall mean
      mu_guess <- log(mean(expression_values))
      fit_mean <- unlist(optim(par=mu_guess, fn=cost_fun_logl_meanfit_comp,
        y_mat=expression_values, disp_est=dispersion_estimate,
        method="BFGS", control=list(maxit=MAXIT,fnscale=-1))[c("par","value","convergence")])
      
      # Summarise results
      if(("optim" %in% optim_method) && ("nlminb" %in% optim_method)){
        fmax_outs <- cbind(fit_peak1,fit_peak2,fit_valley1,fit_valley2)
      }
      if(("optim" %in% optim_method) && !("nlminb" %in% optim_method)){
        fmax_outs <- cbind(fit_peak1,fit_valley1)
      }
      if(!("optim" %in% optim_method) && ("nlminb" %in% optim_method)){
        fmax_outs <- cbind(fit_peak2,fit_valley2)
      } 
      
      # Name row containing value of objective as value 
      # This name differs between optimisation functions
      # DAVID: take out when take out nlminb
      rownames(fmax_outs) <- c(rownames(fmax_outs[1:(nrow(fmax_outs)-1),]),"value")
    })
    
    # Select best fit and report fit type
    # Report mean fit objective value as null hypothesis, too.
    # match() selects first hit if maximum occurs multiple times
    best_fit <- match(max(fmax_outs["value",]),fmax_outs["value",])
    fit_and_null <- c(fmax_outs[,best_fit],fit_mean)
    names(fit_and_null) <- c("beta","h0","h1","h2","t1","t2",
      "logL_H1","converge_H1","mu","logL_H0","converge_H0")
    
    return(fit_and_null)
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
  #   NPARAM: 
  #   n_it: (scalar) [Defaul 100] Number of iterations, which are performed 
  #       to fit the impulse model to the clusters.
  #   ctrl_tc
  #   ctrl: (bool)
  #   n_process: (scalar) [Default 4] number of processes, which can be used on 
  #       the machine to run the background calculation in parallel
  #   ...
  # OUTPUT:
  #   res: (list 2 ["impulse_parameters","impulse_fits"])
  #       impulse_parameters: (matrix genes x (NPARAM+1)*runs) +1 is value of 
  #           objective of optimisation (i.e. sum of squares or weighted SS)
  #       impulse_fits: (matrix genes x timepoints) model values for gene data
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  impulse_fit_matrix <- function(data_arr, timepoints,
    ctrl_tc = FALSE, ctrl = NULL, n_process = 4, dispersion_vec=NULL,NPARAM=6){
    
    MAXIT <- 1000 # Number of maximum iterations for MLE fitting
    mc <- min(detectCores() - 1, n_process)
    
    # Use parallelisation if number of genes/centroids to fit is large
    if(nrow(data_arr) > max(2*mc,10)){
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
      assign("calc_impulse_comp", calc_impulse_comp, envir = my.env)
      assign("impulse_fit", impulse_fit, envir = my.env)
      assign("cost_fun_logl_comp", cost_fun_logl_comp, envir = my.env)
      assign("cost_fun_logl_meanfit_comp", cost_fun_logl_meanfit_comp, envir = my.env)
      assign("impulse_fit_gene_wise",impulse_fit_gene_wise, envir = my.env)
      assign("impulse_fit_matrix",impulse_fit_matrix, envir = my.env)
      assign("timepoints", timepoints, envir = my.env)
      assign("NPARAM", NPARAM, envir = my.env)
      assign("MAXIT", MAXIT, envir = my.env)
      assign("dispersion_vec", dispersion_vec, envir = my.env)
      
      clusterExport(cl=cl, varlist=c("data_arr","dispersion_vec",
        "calc_impulse_comp","impulse_fit","impulse_fit_gene_wise",
        "cost_fun_logl_comp","cost_fun_logl_meanfit_comp",
        "impulse_fit_matrix", "timepoints",
        "MAXIT","NPARAM"), envir = my.env)
      
      # Fit impulse model to each gene of matrix and get impulse parameters:
      # clusterApply runs the function impulse_fit_gene_wise
      # The input data are distributed to nodes by ind_list partitioning
      lsmatFits <- clusterApply(cl, 1:length(ind_list), function(z){t(sapply(ind_list[[z]],
        function(x){impulse_fit_gene_wise(expression_values=data_arr[x,,],
          timepts=timepoints,NPARAM=NPARAM,MAXIT=MAXIT,
          dispersion_estimate=dispersion_vec[x])}))})
      # Give output rownames again, which are lost above
      for(i in 1:length(ind_list)){
        rownames(lsmatFits[[i]]) <- rownames(data_arr[ind_list[[i]],,])
      }
      # Concatenate the output objects of each node
      matFits <- do.call(rbind,lsmatFits)
      stopCluster(cl)
      # Parallelisation ends here
      
      # Do not use parallelisation if number of genes/centroids to fit is small
    } else {
      # Fit impulse model to each gene of matrix and get impulse parameters
      ind_list_all <- 1:nrow(data_arr)
      matFits <- lapply(ind_list_all,function(x){
        impulse_fit_gene_wise(expression_values=data_arr[x,,],
          timepts=timepoints,NPARAM=NPARAM,MAXIT=MAXIT,
          dispersion_estimate=dispersion_vec[x])})
      matFits_temp <- array(NA,c(length(matFits),length(matFits[[1]])))
      for (i in 1:length(matFits)){
        matFits_temp[i,] <- matFits[[i]]
      }
      matFits <- matFits_temp
    }   
    colnames(matFits) <- c("beta","h0","h1","h2","t1","t2","logL_H1",
      "converge_H1","mu","logL_H0","converge_H0")
    rownames(matFits) <- rownames(data_arr)
    
    ### DAVID: can reduce this if clause?
    # Use obtained impulse parameters to calculate impulse fit values   
    if(nrow(matFits) == 1){
      # if matrix contains only one gene
      matImpulseValues <- as.data.frame(t(calc_impulse_comp(matFits[,1:NPARAM],
        unique(sort(timepoints)))),row.names=rownames(matFits))
    } else {                    
      # if matrix contains > 1 genes
      matImpulseValues <- t(apply(matFits[,1:NPARAM],1,function(x){calc_impulse_comp(x,
        unique(sort(timepoints)))}))
    } 
    colnames(matImpulseValues) = unique(sort(timepoints))
    rownames(matImpulseValues) <- rownames(data_arr)
    
    # Report results.
    res <- list()
    res[[1]] <- matFits   # "beta","h0","h1","h2","t1","t2","logL_H1",
                          # "converge_H1","mu","logL_H0","converge_H0" [x3 if with control]
    res[[2]] <- matImpulseValues   # values of impulse model on input timepoints and given genes
    names(res) <- c("impulse_parameters","impulse_fits")
    return(res)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~ 3. Prepare data for impulse model fit ~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  #g_names = rownames(data_input)
  results = list()
  
  if(control_timecourse == TRUE){
    runs = 3
    label = c("combined","case","control")
  } else {
    runs = 1
    label = "case"
  }
  
  # data_input is defined in the following as:
  #   1) A list of three data 3D arrays combined, case and control
  #   2) A list of a single data 3D array if only have case
  
  # Perform 3 runs (combined, case and control) if control timecourse is present
  if(control_timecourse == TRUE){    
    # Combined: Expression values of both case and control (all columns).
    dat1 = data_input
    # Case: Expression values of case.
    dat2 = data_input[,!(data_annotation$Condition %in% control_name),]
    # Control: Expression values of control.  
    dat3 = data_input[,(data_annotation$Condition %in% control_name),]
    
    # Formatting into 3D breaks if single gene is selected:
    # If exactly one gene is tested:
    if(dim(dat1)[1]==1){ 
      dat1_temp <- array(NA,c(1,dim(data_input)[2],dim(data_input)[3]))
      for(z in 1:dim(data_input)[3]){
        dat1_temp[1,,z] <- dat1[,z]
        dat2_temp[1,,z] <- dat2[,z]
        dat3_temp[1,,z] <- dat3[,z]
      }
      dat1 <- dat1_temp
      dat2 <- dat2_temp
      dat3 <- dat3_temp
      rownames(dat1) <- rownames(data_input)
      colnames(dat1) <- colnames(data_input)
      rownames(dat2) <- rownames(data_input)
      colnames(dat2) <- colnames(data_input[,!(data_annotation$Condition %in% control_name),])
      rownames(dat3) <- rownames(data_input)
      colnames(dat3) <- colnames(data_input[,(data_annotation$Condition %in% control_name),])
    }
    data_input <- list(dat1, dat2, dat3)
    
  } else if(control_timecourse == FALSE){
    dat1 = data_input[,,]
    # Formatting into 3D breaks if only single gene is selected:
    # If exactly one gene is tested:
    if(dim(dat1)[1]==1){ 
      dat1_temp <- array(NA,c(1,dim(data_input)[2],dim(data_input)[3]))
      for(z in 1:dim(data_input)[3]){
        dat1_temp[1,,z] <- dat1[,z]
      }
      dat1 <- dat1_temp
      rownames(dat1) <- rownames(data_input)
      colnames(dat1) <- colnames(data_input)
    }
    data_input <- list(dat1)
  }
  
  # Fitting for different runs
  for (c_runs in 1:runs){
    imp_res <- impulse_fit_matrix(data_arr=data_input[[c_runs]],
      timepoints=as.numeric(as.character(
        data_annotation[colnames(data_input[[c_runs]]),"Time"])), 
      ctrl_tc = control_timecourse,ctrl = control_name,
      n_process = n_proc, dispersion_vec=dispersion_vector,NPARAM=NPARAM)
    
    # DAVID: do i still need this?
    #imp_res$impulse_parameters <- (imp_res$impulse_parameters)[g_names,]
    #imp_res$impulse_fits <- (imp_res$impulse_fits)[g_names,]
    
    names(imp_res) <- paste(c("impulse_parameters","impulse_fits"),label[c_runs],sep="_")
    results[[names(imp_res)[1]]] <- imp_res[[1]]
    results[[names(imp_res)[2]]] <- imp_res[[2]]
  }
  
  return(results)
}

