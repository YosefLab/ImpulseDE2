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
  
  impulse_fit_gene_wise <- function(expression_values, timepts, 
    NPARAM=6, MAXIT=10000, dispersion_estimate=NULL,...){
    
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
      (expression_means[x+1]-expression_means[x])/(timepts[x+1]-timepts[x])}) )
    if(max_mean < 1.001){max_mean <- 1.001}
    if(min_mean < 0.001){min_mean <- 0.0001}
    if(max_value < 1.001){max_value <- 1.001}
    if(min_value < 0.001){min_value <- 0.001}
    if(max_middle_mean < 1.001){max_middle_mean <- 1.001}
    if(min_middle_mean < 0.001){min_middle_mean <- 0.001}
    
    tmm2 <- system.time({
      
      # 1. Peak:
      # Beta: Has to be negative
      # Theta1: Low
      # Theta2: High
      # Theta3: Low
      # t1: Around first observed inflexion point
      # t2: Around second observed inflexion point
      lower_inflexion_point_ind <- match(max(neighbour_grad[1:(max_middle_mean_ind-1)]), neighbour_grad[1:(max_middle_mean_ind-1)])
      upper_inflexion_point_ind <- max_middle_mean_ind - 1 + match(min(neighbour_grad[max_middle_mean_ind:length(neighbour_grad)]), neighbour_grad[max_middle_mean_ind:length(neighbour_grad)])
      #upper_h0 <- min( max(expression_values[1:(max_middle_mean_ind-1),]), max_middle_mean )
      #upper_h2 <- min( max(expression_values[(max_middle_mean_ind+1):num_points,]), max_middle_mean )
      par_guess <- c(1,log(expression_means[1]+1),log(max_middle_mean+1),log(expression_means[num_points]+1),
        (timepts[lower_inflexion_point_ind]+timepts[lower_inflexion_point_ind+1])/2,
        (timepts[upper_inflexion_point_ind]+timepts[upper_inflexion_point_ind+1])/2)
      #lower_b <- c(0.01,0.001,max(expression_means[1],expression_means[num_points]),0.001,(timepts[1]+timepts[2])/3,timepts[max_middle_mean_ind])
      #upper_b <- c(100,upper_h0,max_middle_val,upper_h2,timepts[max_middle_mean_ind],timepts[num_points])
      # check validity of bounds
      #upper_b[upper_b <= lower_b] <- lower_b[upper_b <= lower_b] + 1
      # check that guess fullfills bounds
      #par_guess[par_guess <= lower_b] <- lower_b[par_guess <= lower_b] + 0.01
      #par_guess[par_guess >= upper_b] <- upper_b[par_guess >= upper_b] - 0.01
      if("optim" %in% optim_method){
        fit_peak1 <- unlist( optim(par=par_guess, fn=cost_fun_logl_comp, x_vec=timepts,
          y_mat=expression_values, disp_est=dispersion_estimate,
          method="BFGS", control=list(maxit=MAXIT))[c("par","convergence","value")] )
        #method="L-BFGS-B",lower=lower_b, upper=upper_b)[c("par","value")])
      }
      if("nlminb" %in% optim_method){
        fit_peak2 <- unlist(nlminb(start=par_guess, objective=cost_fun_logl_comp, x_vec=timepts,
          y_mat=expression_values, disp_est=dispersion_estimate,
          lower=lower_b, upper=upper_b)[c("par","objective")])
      }
      # 2. Valley
      # Beta: Has to be negative
      # Theta1: High
      # Theta2: Low
      # Theta3: High
      # t1: Around first observed inflexion point
      # t2: Around second observed inflexion point
      lower_inflexion_point_ind <- match(min(neighbour_grad[1:(min_middle_mean_ind-1)]), neighbour_grad[1:(min_middle_mean_ind-1)])
      upper_inflexion_point_ind <- min_middle_mean_ind - 1 + match(max(neighbour_grad[min_middle_mean_ind:(num_points-1)]), neighbour_grad[min_middle_mean_ind:(num_points-1)])
      #upper_h0 <- max(expression_values[1:(min_middle_mean_ind-1),])
      #upper_h2 <- max(expression_values[(min_middle_mean_ind+1):num_points,])
      par_guess <- c(1,log(expression_means[1]+1),log(min_middle_mean+1),log(expression_means[num_points]+1),
        (timepts[lower_inflexion_point_ind]+timepts[lower_inflexion_point_ind+1])/2,
        (timepts[upper_inflexion_point_ind]+timepts[upper_inflexion_point_ind+1])/2 )
      #lower_b <- c(-0.01,min_middle_mean,0.001,min_middle_mean,(timepts[1]+timepts[2])/3,timepts[min_middle_mean_ind])
      #upper_b <- c(-100,upper_h0,min(expression_means[1],expression_means[num_points]),upper_h2,timepts[min_middle_mean_ind],timepts[num_points])
      # check validity of bounds
      #upper_b[upper_b <= lower_b] <- lower_b[upper_b <= lower_b] + 1      
      # check that guess fullfills bounds
      #par_guess[par_guess <= lower_b] <- lower_b[par_guess <= lower_b] + 0.01
      #par_guess[par_guess >= upper_b] <- upper_b[par_guess >= upper_b] - 0.01
      if("optim" %in% optim_method){
        fit_valley1 <- unlist( optim(par=par_guess, fn=cost_fun_logl_comp, x_vec=timepts,
          y_mat=expression_values, disp_est=dispersion_estimate,
          method="BFGS", control=list(maxit=MAXIT))[c("par","convergence","value")] )
        #method="L-BFGS-B",lower=lower_b, upper=upper_b)[c("par","value")])
      }
      if("nlminb" %in% optim_method){
        fit_valley2 <- unlist(nlminb(start=par_guess, objective=cost_fun_logl_comp, x_vec=timepts,
          y_mat=expression_values, disp_est=dispersion_estimate,
          lower=lower_b, upper=upper_b)[c("par","objective")])
      }
      # DAVID take out these initialisations to save time
      if(FALSE){
        # 3. Up
        # Beta: Has to be positive
        # Theta1: Low
        # Theta2: High
        # Theta3: Arbitrary (outside reach - defined through t2)
        # t1: Around observed inflexion point
        # t2: Should be well outside reach so that model doesn't use the 2nd sigmoid.
        inflexion_point_ind <- match(max(neighbour_grad), neighbour_grad)
        par_guess <- c(1,min_mean,max_mean,min_mean,
          (timepts[inflexion_point_ind]+timepts[inflexion_point_ind+1])/2,
          2.5*timepts[num_points])
        #lower_b <- c(0.01,0.001,max_middle_mean,0.001,(timepts[1]+timepts[2])/3,2*timepts[num_points])
        #upper_b <- c(100,max_middle_mean,max_value,max_middle_mean,timepts[num_points],3*timepts[num_points])
        # check validity of bounds
        #upper_b[upper_b <= lower_b] <- lower_b[upper_b <= lower_b] + 1      
        # check that guess fullfills bounds
        #par_guess[par_guess <= lower_b] <- lower_b[par_guess <= lower_b] + 0.01
        #par_guess[par_guess >= upper_b] <- upper_b[par_guess >= upper_b] - 0.01
        if("optim" %in% optim_method){
          fit_up1 <- unlist(optim(par=par_guess, fn=cost_fun_logl_comp, x_vec=timepts,
            y_mat=expression_values, disp_est=dispersion_estimate,
            method="L-BFGS-B")[c("par","value")])
          #method="L-BFGS-B",lower=lower_b, upper=upper_b)[c("par","value")])
        }
        if("nlminb" %in% optim_method){
          fit_up2 <- unlist(nlminb(start=par_guess, objective=cost_fun_logl_comp, x_vec=timepts,
            y_mat=expression_values, disp_est=dispersion_estimate,
            lower=lower_b, upper=upper_b)[c("par","objective")])
        }
        # 4. Down
        # Beta: Has to be negative
        # Theta1: High
        # Theta2: Low
        # Theta3: Arbitrary (outside reach - defined through t2)
        # t1: Around observed inflexion point
        # t2: Should be well outside reach so that model doesn't use the 2nd sigmoid.
        inflexion_point_ind <- match(min(neighbour_grad), neighbour_grad)
        par_guess <- c(-1,max_mean,min_mean,max_mean,
          (timepts[inflexion_point_ind]+timepts[inflexion_point_ind+1])/2,
          2.5*timepts[num_points])
        #lower_b <- c(-0.01,min_middle_mean,0.001,0.001,(timepts[1]+timepts[2])/3,2*timepts[num_points])
        #upper_b <- c(-100,max_value,max_middle_mean,max_value,timepts[num_points],3*timepts[num_points])
        # check validity of bounds
        #upper_b[upper_b <= lower_b] <- lower_b[upper_b <= lower_b] + 1
        # check that guess fullfills bounds
        #par_guess[par_guess <= lower_b] <- lower_b[par_guess <= lower_b] + 0.01
        #par_guess[par_guess >= upper_b] <- upper_b[par_guess >= upper_b] - 0.01
        if("optim" %in% optim_method){
          fit_down1 <- unlist(optim(par=par_guess, fn=cost_fun_logl_comp, x_vec=timepts,
            y_mat=expression_values, disp_est=dispersion_estimate,
            method="L-BFGS-B")[c("par","value")])
          #method="L-BFGS-B",lower=lower_b, upper=upper_b)[c("par","value")])
        }
        if("nlminb" %in% optim_method){
          fit_down2 <- unlist(nlminb(start=par_guess, objective=cost_fun_logl_comp, x_vec=timepts,
            y_mat=expression_values, disp_est=dispersion_estimate,
            lower=lower_b, upper=upper_b)[c("par","objective")])
        }
      }
      # 5. Mean
      # Parameter estimates: Only h1 is in modeled time interval and represents mean
      mu_guess <- log(mean(expression_values))
      #lower_b <- c(min_value-1)
      #upper_b <- c(max_value+1)
      # check validity of bounds
      #upper_b[upper_b <= lower_b] <- lower_b[upper_b <= lower_b] + 1
      # check that guess fullfills bounds
      #if(mu_guess <= lower_b){mu_guess <- lower_b + 0.01}
      #if(mu_guess >= upper_b){mu_guess <- upper_b - 0.01}
      fit_mean <- unlist(optim(par=mu_guess, fn=cost_fun_logl_meanfit_comp,
        y_mat=expression_values, disp_est=dispersion_estimate,
        method="L-BFGS-B", control=list(maxit=MAXIT))[c("par","value","convergence")])
      #method="L-BFGS-B",lower=lower_b, upper=upper_b)[c("par","value")])
      # fit_mean is object of 1D optimisation, fill up missing parameters
      # to conform with impulse model reporting scheme
      # This impulse model corresponds to a peak which is almost 
      # constant in the modelled time range and only starts peaking 
      # much later (3*tmax in modelled range)
      #fit_mean <- c(100,fit_mean[1],fit_mean[1],fit_mean[1],
      #  timepts[1],timepts[num_points],fit_mean[2])
      
      # Summarise results
      if(("optim" %in% optim_method) && ("nlminb" %in% optim_method)){
        fmin_outs <- cbind(fit_peak1,fit_peak2,fit_valley1,fit_valley2,fit_up1,fit_up2,
          fit_down1,fit_down2,fit_mean)
      }
      if(("optim" %in% optim_method) && !("nlminb" %in% optim_method)){
        #fmin_outs <- cbind(fit_peak1,fit_valley1,fit_up1,fit_down1,fit_mean)
        fmin_outs <- cbind(fit_peak1,fit_valley1)
      }
      if(!("optim" %in% optim_method) && ("nlminb" %in% optim_method)){
        fmin_outs <- cbind(fit_peak2,fit_valley2,fit_up2,fit_down2,fit_mean)
      } 
      
      # Name row containing value of objective as value 
      # This name differs between optimisation functions
      # DAVID: take out when take out nlminb
      rownames(fmin_outs) <- c(rownames(fmin_outs[1:(nrow(fmin_outs)-1),]),"value")
    })
    
    # Select best fit and report fit type
    # Report mean fit objective value as null hypothesis, too.
    # match() selects first hit if minimum occurs multiple times
    best_fit <- match(min(fmin_outs["value",]),fmin_outs["value",])
    fit_and_null <- c(fmin_outs[,best_fit],fit_mean[2])
    names(fit_and_null) <- c("beta","h0","h1","h2","t1","t2","convergence","objective","nullfit")
    
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
  
  impulse_fit_matrix <- function(data_arr, timepoints,
    ctrl_tc = FALSE, ctrl = NULL, n_process = 4, dispersion_vec=NULL){
    
    NPARAM = 6 # Number of model parameters
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
      
      junk <- clusterEvalQ(cl,library(MASS))
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
    
    ### DAVID: can reduce this if clause?
    # Use obtained impulse parameters to calculate impulse fit values
    colnames(matFits) <- c("beta","h0","h1","h2","t1","t2","convergence","objective","nullfit")
    if(nrow(matFits) == 1){      # if matrix contains only one gene
      matImpulseValues <- as.data.frame(t(calc_impulse_comp(matFits[,1:NPARAM],
        unique(sort(timepoints)))),row.names=rownames(matFits))
      colnames(matImpulseValues) = unique(sort(timepoints))
    } else {                    # if matrix contains > 1 genes
      matImpulseValues <- t(apply(matFits[,1:NPARAM],1,function(x){calc_impulse_comp(x,
        unique(sort(timepoints)))}))
      colnames(matImpulseValues) = unique(sort(timepoints))
    } 
    
    # Report results.
    res <- list()
    rownames(matFits) <- rownames(data_arr)
    rownames(matImpulseValues) <- rownames(data_arr)
    res[[1]] <- matFits    # "beta","h0","h1","h2","t1","t2","objective","nullfit" [x3 if with control]
    res[[2]] <- matImpulseValues   # values of impulse model on input timepoints and given genes
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
      n_process = n_proc, dispersion_vec=dispersion_vector)
    
    # DAVID: do i still need this?
    imp_res$impulse_parameters <- (imp_res$impulse_parameters)[g_names,]
    imp_res$impulse_fits <- (imp_res$impulse_fits)[g_names,]
    
    names(imp_res) <- paste(c("impulse_parameters","impulse_fits"),label[c_runs],sep="_")
    results[[names(imp_res)[1]]] <- imp_res[[1]]
    results[[names(imp_res)[2]]] <- imp_res[[2]]
  }
  
  return(results)
}

