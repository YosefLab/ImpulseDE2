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
#   sigma_mat:
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

impulse_fit <- function(data_input, data_annotation, sigma_mat, n_iter = 100,
    control_timecourse = FALSE, control_name = NULL, cluster_results = NULL,
    start_values = NULL, fit_backg = FALSE, n_proc = 4){
  
  ### Check dataset for consistency
  if(ncol(data_annotation) != 2 ||
      FALSE %in% (colnames(data_annotation) == c("Time","Condition"))){
        stop("Please use function 'annotation_preparation' to prepare
            the annotation table")
  }
  #  if(length(summary(as.factor(data_annotation$Condition)))>=2 &
  #     control_timecourse == FALSE){
  #     stop("Too many conditions although no control timecourse is specified.
  #         Please use function 'annotation_preparation' to prepare the annotation
  #         table or set variables 'control_timecourse' and 'control_name'")
  #  }
  if(control_timecourse == TRUE & is.null(control_name)){
    stop("Please specify the name of the control in the column 'Condition'")
  }
  if(control_timecourse == TRUE & (is.null(start_values) == FALSE)){
    start_values <- start_values[c(1,3,5)]
  }
  if(control_timecourse == FALSE & (is.null(control_name) == FALSE)){
    start_values <- start_values[c(1,3)]
  }
   
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~ 1. Fit impulse model to gene ~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  ### Fits an impulse model to a single gene
  
  # INPUT:
  #   expression_values: (numeric 2D array timepoints x replicates) 
  #   data_annotation: (Table samples x 2[time and condition]) 
  #       Co-variables for the samples including condition and time points.
  #       Time points must be numeric.
  #   sigma_vec:
  #   NPARAM: (scalar) Number of model parameters
  #   n_iter: (scalar) [Default 100] Number of iterations, which are performed 
  #       to fit the impulse model to the clusters.
  #   fit_to_clusters: (bool) [Default FALSE]
  #   start_val: 
  #   fit_backg: (bool) [Default FALSE]
  #   n_process: (scalar) [Default 4] number of processes, which can be used on 
  #       the machine to run the background calculation in parallel
  # OUTPUT:
  #   pvec_and_SSE: (vector NPARAM+1) +1 is value of objective of optimisation 
  #       (i.e. sum of squares or weighted SS)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  impulse_fit_gene_wise <- function(expression_values, timepoints, sigma_vec, 
      NPARAM, n_iters = 100, fit_to_clusters = FALSE, start_val = NULL, 
      fit_bg = FALSE, ...){
    
    # Initial parameter guesses:
    # beta: slope
    # h0, h2: first and last value 
    # h1: peak value
    # t1, t2: peak time difference to ends
    # Do initial guesses based on centroids
    
    # If handed a list, i.e. not replicates, e.g. fitting centroids
    if (is.vector(expression_values)==TRUE){
      # Transform into 2D array
      expression_values <- array(expression_values,c(length(expression_values),1))
      # Do intitial guesses based on centroids
      expression_means <- expression_values
    } else {
      # Do initial guesses based on mean over replicates:
      expression_means <- apply(expression_values,1,mean)
    }
    mat1 = c(0,abs(expression_means[2:length(expression_means)] -
        expression_means[1:(length(expression_means)-1)]) *
        expression_means[2:length(expression_means)])
    mn_beta = 0
    mx_beta = 10
    middle_ind = which(timepoints > min(timepoints) & 
        timepoints < max(timepoints))
    min_ind = which(timepoints == min(timepoints))[1]
    max_ind = which(timepoints == max(timepoints))[1]
    mx_tm = max(timepoints)
    mn_tm = min(timepoints)
    tmp_ind = which(mat1 == max(mat1))[1]
    peaktime = timepoints[tmp_ind]
    peak = min(which(timepoints == peaktime))
    beta1 = abs(log(abs((expression_means[peak] - 
        expression_means[min_ind])/(peaktime-
        timepoints[min_ind]))));
    ### DAVID new param
    num_points <- length(expression_means)
    max_value = max(expression_values)
    min_value = min(expression_values)
    max_mean <- max(expression_means)
    min_mean <- min(expression_means)
    max_mean_ind <- match(max_mean,expression_means)
    min_mean_ind <- match(min_mean,expression_means)
    max_middle_mean <- max(expression_means[2:num_points-1])
    min_middle_mean <- min(expression_means[2:num_points-1])
    # +1 to push indicices up from middle stretch to entire window (first is omitted here)
    max_middle_mean_ind <- match(max_middle_mean,expression_means[2:num_points-1]) + 1
    min_middle_mean_ind <- match(min_middle_mean,expression_means[2:num_points-1]) + 1
    # Gradients between neighbouring points
    neighbour_grad <- unlist( lapply(c(1:(num_points-1)),function(x){
      (expression_means[x+1]-expression_means[x])/(timepoints[x+1]-timepoints[x])}) )
    if(max_mean <= 1){max_mean <- 1.0001}
    if(min_mean <= 1){min_mean <- 1.0001}
    if(max_value <= 1){max_value <- 1.0001}
    if(min_value <= 1){min_value <- 1.0001}
    if(max_middle_mean <= 1){max_middle_mean <- 1.0001}
    if(min_middle_mean <= 1){min_middle_mean <- 1.0001}
    
    ### set beta to 1 if calculcated value is infinite
    if(is.finite(beta1) == FALSE || is.na(beta1)) { beta1 = 1 }
    
    ### define start value theta based on the gene's data
    orig_theta = c(beta1,
                   expression_means[min_ind],                   # h0
                   expression_means[peak],                      # h1
                   expression_means[max_ind],                   # h2
                   (peaktime-timepoints[min_ind])/2,            # t1
                   peaktime+(timepoints[max_ind] - peaktime)/2) # t2
    names(orig_theta) = c("beta","h0","h1","h2","t1","t2")
    theta = orig_theta
    
    ### replace theta estimates of 0 by small value, because optimization
    ### function has problems with 0 as input
    rand_num <- runif(length(which(theta == 0)),0,0.0001)
    theta[names(which(theta == 0))] <- rand_num
    names(theta) = c("beta","h0","h1","h2","t1","t2")
    
    #####################
    # 1.  Fit to clusters
    if(fit_to_clusters == TRUE & is.null(start_val)){
      
      ### DAVID what is this, i ll ignore it, handing means above so should be fine
      if(fit_bg == FALSE){
        thetas <- cbind(as.numeric(theta),
          apply(matrix(rep(runif((n_iters-1)*NPARAM)),NPARAM,(n_iters-1)),
            2, function(x){c(mn_beta + (mx_beta-mn_beta)*x[1],
            min(expression_means) + (max(expression_means) - 
              min(expression_means))*x[2:4],
            mn_tm + (mx_tm-mn_tm)*x[5],
            (mn_tm + (mx_tm-mn_tm)*x[5]) + (mx_tm-(mn_tm + 
              (mx_tm-mn_tm)*x[5]))*x[6])}))
      } else {
        thetas <- apply(matrix(rep(runif((n_iters)*NPARAM)),NPARAM,(n_iters)),
            2, function(x){c(mn_beta + (mx_beta-mn_beta)*x[1],
            min(expression_means) + (max(expression_means) -
              min(expression_means))*x[2:4],
            mn_tm + (mx_tm-mn_tm)*x[5],
            (mn_tm + (mx_tm-mn_tm)*x[5]) + (mx_tm-(mn_tm + 
                (mx_tm-mn_tm)*x[5]))*x[6])})
      }
      # Use ordinary least squares on fitting centroids: Variance information
      # does not apply here.
      tmm1 <- system.time({
        fmin_outs <- cbind(
          # fitting via quasi-Newton method (optim, "BFGS")
          apply(thetas[,1:(floor(n_iters/2))],
            2,function(x){unlist(optim(x,two_impulses_OLS_comp, x_vec = timepoints,
            y_mat = expression_values, method="L-BFGS-B", 
              lower=c(-100,0.0001,0.0001,0.0001,0.0001,0.0001), upper=c(100,max_value,max_value,max_value,mx_tm,2*mx_tm))[c("par","value")])}),
          # fitting via PORT routines (nlminb)
          apply(thetas[, c(1,(floor(n_iters/2)+1):n_iters)],
            2,function(x){unlist(nlminb(x, two_impulses_OLS_comp, x_vec = timepoints,
            y_mat = expression_values, 
              lower=c(-100,0.0001,0.0001,0.0001,0.0001,0.0001), upper=c(100,max_value,max_value,max_value,mx_tm,2*mx_tm))[c("par","objective")])})
        )
      })
      
      
      ### return results of 3 best fits, which will serve as start values for
      ### the fits to the genes
      pvec_and_SSE = as.vector(fmin_outs[,order(fmin_outs["value",])][,1:3])
      
    #################
    # 2. Fit to genes
    } else if(fit_to_clusters == FALSE & is.null(start_val) == FALSE){
      tmm2 <- system.time({
        
        # Toll has the parameter intialisations
        ### use first 3 best hits from fit to clusters plus initial guess 'theta'
        ### as start values
        if(fit_bg == TRUE){
          toll <- matrix(start_val,7,3)[1:NPARAM,]
        } else {
          toll <- cbind(theta,matrix(start_val,7,3)[1:NPARAM,][,1:3])
        }
        
        # NOT working in log 2 space
        # 1. Peak:
        if(max_middle_mean > 2*expression_means[1] &&
            max_middle_mean > 2*expression_means[num_points] &&
            ((expression_means[1] > 10 && expression_means[num_points] > 10) || 
            (max_middle_mean > 4*expression_means[1] ||
            max_middle_mean > 4*expression_means[num_points]))){
          #guess1 <- c(1,expression_means[1],max_mean,expression_means[num_points],timepoints[max_middle_mean_ind-1],timepoints[max_middle_mean_ind])
          guess1 <- c(max(neighbour_grad),expression_means[1],max_middle_mean,expression_means[num_points],
            timepoints[match(max(neighbour_grad[1:max_middle_mean_ind-1]), neighbour_grad[1:max_middle_mean_ind-1])],
            timepoints[match(max(neighbour_grad[max_middle_mean_ind:num_points-1]), neighbour_grad[max_middle_mean_ind:num_points-1])])
          lower_b <- c(0.0001,0.0001,1.001,max(expression_means[1],expression_means[num_points]),0.0001,0.0001)
          upper_b <- c(1000,max_middle_mean,max_value,max_middle_mean,timepoints[num_points],timepoints[num_points])
        # 2. Valley
        } else if(min_middle_mean < expression_means[1]/2 &&
            min_middle_mean < expression_means[num_points]/2 &&
            expression_means[1] > 10 && expression_means[num_points] > 10){
          #guess1 <- c(1,expression_means[1],min_mean,expression_means[num_points],timepoints[min_middle_mean_ind-1],timepoints[min_middle_mean_ind+1])
          guess1 <- c(min(neighbour_grad),expression_means[1],min_middle_mean,expression_means[num_points],
            timepoints[match(min(neighbour_grad[1:min_middle_mean_ind-1]), neighbour_grad[1:min_middle_mean_ind-1])],
            timepoints[match(min(neighbour_grad[min_middle_mean_ind:num_points-1]), neighbour_grad[min_middle_mean_ind:num_points-1])])
          lower_b <- c(-0.0001,min_middle_mean,1.0001,min_middle_mean,0.0001,0.0001)
          upper_b <- c(-1000,max_value,min(expression_means[1],expression_means[num_points]),max_value,timepoints[num_points],timepoints[num_points])
        # 3. Up
        } else if(expression_means[1] < expression_means[num_points]){
          #guess1 <- c(10,min_mean,max_mean,min_mean,timepoints[2],mx_tm+timepoints[2])
          guess1 <- c(max(neighbour_grad),min_mean,max_mean,min_mean,
            timepoints[match(max(neighbour_grad[1:num_points-1]), neighbour_grad[1:num_points-1])],
            mx_tm+timepoints[num_points-1])
          lower_b <- c(0.0001,0.0001,min_middle_mean,0.0001,0.0001,timepoints[num_points])
          upper_b <- c(1000,max_middle_mean,max_value,max_middle_mean,timepoints[num_points],2*timepoints[num_points])
        # 4. Down
        } else if(expression_means[1] >= expression_means[num_points]){
          #guess1 <- c(10,max_mean,min_mean,max_mean,timepoints[2],mx_tm+timepoints[2])
          guess1 <- c(min(neighbour_grad),max_mean,min_mean,max_mean,
            timepoints[match(min(neighbour_grad[1:num_points-1]), neighbour_grad[1:num_points-1])+1],
            mx_tm+timepoints[num_points-1])
          lower_b <- c(-0.0001,min_middle_mean,1.0001,min_middle_mean,0.0001,timepoints[num_points])
          upper_b <- c(-1000,max_value,max_middle_mean,max_value,timepoints[num_points],2*timepoints[num_points])
        }
        print(guess1)
        # check that lower time bound is not 0
        if(guess1[[5]] < 0.001){guess1[5]<- 0.001}
        toll <- cbind(toll,guess1)
        
        # Use weighted least squares when fitting individual genes.
        fmin_outs <- cbind(
          # fitting via quasi-Newton method (optim, "BFGS")
          ### DAVID: took out as.numeric for y_mat here
          apply(toll,
            2,function(x){unlist(optim(x, two_impulses_WLS_comp, x_vec=timepoints,
            y_mat=expression_values, sigma_vec=sigma_vec, method="L-BFGS-B",
              lower=c(-100,0.0001,0.0001,0.0001,0.0001,0.0001), upper=c(100,max_value,max_value,max_value,mx_tm,2*mx_tm))[c("par","value")])}),
              #lower=lower_b, upper=upper_b)[c("par","value")])}),
          # again fitting via PORT routines (nlminb)
          apply(toll,
            2,function(x){unlist(nlminb(x, two_impulses_WLS_comp, x_vec=timepoints,
            y_mat=expression_values, sigma_vec=sigma_vec,
              lower=c(-100,0.0001,0.0001,0.0001,0.0001,0.0001), upper=c(100,max_value,max_value,max_value,mx_tm,2*mx_tm))[c("par","objective")])})
              #lower=lower_b, upper=upper_b)[c("par","objective")])})
        )
      })
      if(is.null(dim(fmin_outs[ ,fmin_outs["value",] == min(fmin_outs["value",])])) == TRUE){
        pvec_and_SSE = fmin_outs[ ,fmin_outs["value",] == min(fmin_outs["value",])]
      } else {
        ### if two or more randomization have the same minimum Sum of
        ### Squared Errors (SSE), choose the first one
        pvec_and_SSE = fmin_outs[ ,fmin_outs["value",] ==
            min(fmin_outs["value",])][,1]
      }
    }
    
    return(pvec_and_SSE)
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
  #   sigma_mat:
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
  
  impulse_fit_matrix <- function(data_arr, timepoints, sigma_mat, n_it = 100,
    ctrl_tc = FALSE, ctrl = NULL, fit_to_clus = FALSE, start_val = NULL,
    fit_bg = FALSE, n_process = 4, ...){
    
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
      assign("sigma_mat", sigma_mat, envir = my.env)
      assign("calc_impulse_comp", calc_impulse_comp, envir = my.env)
      assign("n_it", n_it, envir = my.env)
      assign("cluster_genes_for_impulse", cluster_genes_for_impulse,
             envir = my.env)
      assign("impulse_fit", impulse_fit, envir = my.env)
      assign("two_impulses_WLS_comp", two_impulses_WLS_comp, envir = my.env)
      assign("two_impulses_OLS_comp", two_impulses_OLS_comp, envir = my.env)
      assign("impulse_fit_gene_wise",impulse_fit_gene_wise, envir = my.env)
      assign("impulse_fit_matrix",impulse_fit_matrix, envir = my.env)
      assign("timepoints", timepoints, envir = my.env)
      assign("fit_to_clus", fit_to_clus, envir = my.env)
      assign("start_val", start_val, envir = my.env)
      assign("fit_bg", fit_bg, envir = my.env)
      assign("NPARAM", NPARAM, envir = my.env)
      
      clusterExport(cl=cl, varlist=c("data_arr","sigma_mat",
        "calc_impulse_comp","n_it","cluster_genes_for_impulse","impulse_fit",
        "two_impulses_WLS_comp", "two_impulses_OLS_comp", "impulse_fit_gene_wise",
        "impulse_fit_matrix", "timepoints","fit_to_clus","start_val",
        "fit_bg","NPARAM"), envir = my.env)
      
      # Fit impulse model to each gene of matrix and get impulse parameters:
      # clusterApply runs the function impulse_fit_gene_wise
      # The input data are distributed to nodes by ind_list partitioning
      resi <- clusterApply(cl, 1:length(ind_list), function(z){t(sapply(ind_list[[z]],
                function(x){impulse_fit_gene_wise(data_arr[x,,],timepoints,sigma_mat[x,],
                NPARAM,n_it,fit_to_clus, start_val, fit_bg)}))})
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
      resmat <- lapply(ind_list_all,function(z){
          impulse_fit_gene_wise(data_arr[z,,],timepoints,sigma_mat[z,],
          NPARAM,n_it,fit_to_clus, start_val, fit_bg)})
      resmat_temp <- array(0,c(length(resmat),length(resmat[[1]])))
      for (i in 1:length(resmat)){
        resmat_temp[i,] <- resmat[[i]]
      }
      resmat <- resmat_temp
    }   
    
    # Use obtained impulse parameters to calculate impulse fit values
    ############################
    # 1.  Value for fit to genes
    if(fit_to_clus == FALSE){
      colnames(resmat) <- c("beta","h0","h1","h2","t1","t2","SSE")
      if(nrow(resmat) == 1){      # if matrix contains only one gene
        resmat2 <- as.data.frame(t(calc_impulse_comp(resmat[,1:NPARAM],
          unique(sort(timepoints)))),row.names=rownames(resmat))
        colnames(resmat2) = unique(sort(timepoints))
      } else {                    # if matrix contains > 1 genes
        resmat2 <- t(apply(resmat[,1:NPARAM],1,function(x){calc_impulse_comp(x,
          unique(sort(timepoints)))}))
        colnames(resmat2) = unique(sort(timepoints))
      }
    
    ##############################
    # 2. Value for fit to clusters
      ### more complex because results from 3 best fits need to be saved and later
      ### used for the fit to the genes
    } else {
      colnames(resmat) <- c(paste(c("beta","h0","h1","h2","t1","t2","SSE"),1,sep="_"),
                            paste(c("beta","h0","h1","h2","t1","t2","SSE"),2,sep="_"),
                            paste(c("beta","h0","h1","h2","t1","t2","SSE"),3,sep="_"))
      resmat2 <- t(apply(resmat,
        1,function(x){apply(matrix(x,7,3)[1:NPARAM,],
        2,function(y){calc_impulse_comp(y,unique(sort(timepoints)))})}))
      head(resmat2)
      colnames(resmat2) <- c(paste(unique(sort(timepoints)),1,sep="_"),
                             paste(unique(sort(timepoints)),2,sep="_"),
                             paste(unique(sort(timepoints)),3,sep="_"))
    }
    
    # Report results.
    res <- list()
    rownames(resmat) <- rownames(data_arr)
    rownames(resmat2) <- rownames(data_arr)
    res[[1]] <- resmat    # "beta","h0","h1","h2","t1","t2","SSE" [x3 if with control]
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
  # The following if - else if statement distinguishes whether impulse_fit
  # is called to fit the model to a) clusters or b) genes The difference
  # is that fitting to clusters is done as the first fitting stage and
  # therefore does not require starting values. In contrast to that,
  # fitting to genes assumes that genes are associated with clusters, 
  # which have have been fitted before. The cluster parameter sets are 
  # used to initialise the gene-wise models.
  
  #################
  # 1) Fit to genes
  if(is.null(cluster_results) == FALSE & is.null(start_values) == FALSE){
    file_ext = "genes"
    
    # Perform 3 runs (combined, case and control) if control timecourse is present
    # datX contains expression data of genes associates with cluster
    # datX_n contains expression data of genes NOT in cluster (low variation)
    # DAVID dat are defined as data.frames, which i cannot maintain
    # it doesn seem to have relevance downstream though, change to 3D array
    if(control_timecourse == TRUE){
      # Note on variables used below:
      # cluster_results[[1,5,9]] are kmeans_clus for [combined, case, control]
      
      # Combined: Expression values of both case and control (all columns).
      # Clustering based on combined.
      dat1 = data_input[names(cluster_results[[1]]),,]
      dat1_n = data_input[!(rownames(data_input) %in% rownames(dat1)),,] 
      dat1_n = dat1_n[,colnames(dat1),]
      
      # Case: Expression values of case. Clustering based on case.
      dat2 = data_input[,!(data_annotation$Condition %in% control_name),]
      dat2 = dat2[names(cluster_results[[5]]),,]
      dat2_n = data_input[!(rownames(data_input) %in% rownames(dat2)),,]
      dat2_n = dat2_n[,colnames(dat2),]
      
      # Control: Expression values of control. Clustering based on control. 
      dat3 = data_input[,(data_annotation$Condition %in% control_name),]
      dat3 = dat3[names(cluster_results[[9]]),]
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
      dat1 = data_input[names(cluster_results[[1]]),,]
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
    
  ####################  
  # 2) Fit to clusters
    ### DAVID unaffected by 3D
  }  else if(is.null(cluster_results) & is.null(start_values)) {
    file_ext = "clusters"

    # Input are the centroids of the clusters
    # 2 is combined, 6 is case and 10 is control
    # This selects 1 set of centroid if no control is present
    # (runs==1) or 3 sets of centroids if control is present
    # (runs==3).
    data_input <- data_input[c(2,6,10)[1:runs]]
  }
  
  # Fitting for different runs
  for (c_runs in 1:runs){
    imp_res <- NULL
    
    ####################  
    # 1) Fit to genes
    if(is.null(cluster_results) == FALSE & is.null(start_values) == FALSE){
      
      # Split genes into the clusters and fit impulse model to the genes of a cluster:
      # Split expression array of current run into list of arrays for each cluster:
      cluster_res_list <- list()
      for (iClusters in unique(cluster_results[[(c_runs-1)*4 +1]])){
          cluster_res_list[[iClusters]] <- data_input[[c_runs]][
              cluster_results[[(c_runs-1)*4 +1]] %in% iClusters,,]  
      }
      # Get list of cluster indices with numeric index as character string as name of
      # each list element.
      ind <- split(1:max(cluster_results[[(c_runs-1)*4 +1]]), 1:max(cluster_results[[(c_runs-1)*4 +1]]))
      
      # Fit model to each cluster matrix
      imp_res_list <- lapply(ind, function(x){impulse_fit_matrix(cluster_res_list[[x]],
          as.numeric(as.character(data_annotation[colnames(data_input[[c_runs]]),"Time"])), 
          sigma_mat = sigma_mat,n_it = n_iter, ctrl_tc = control_timecourse, 
          ctrl = control_name, fit_to_clus = FALSE,start_val = start_values[[c_runs]][x,], 
          fit_bg = fit_backg, n_process = n_proc)})
      
      # Use 'get' to pull set of not clustered genes in current run (from current environment)
      # If there is >0 genes not clustered in current run
      if(nrow(get(dat_ind[c_runs])) != 0){
        # Report results obtained so far
        tump1  <- do.call(rbind,lapply(1:(length(imp_res_list)), function(x) imp_res_list[[x]][[1]]))
        tump2  <- do.call(rbind,lapply(1:(length(imp_res_list)), function(x) imp_res_list[[x]][[2]]))
        
        # Set model value for non clustered genes to mean, this is the null hypothesis
        tmpp1 <- t(apply(get(dat_ind[c_runs]),1,function(x){rep(mean(x),ncol(tump2))}))
        colnames(tmpp1) <- colnames(tump2)
        tump2 <- rbind(tump2, tmpp1)
        
        # Fill parameter fits of non clustered genes with NA and sum of squares
        # position with squared deviation to mean (H0)
        # DAVID: replace with weighted SS
        datn_WSS <- rep(NA,length(rownames(get(dat_ind[c_runs]))))
        names(datn_WSS) <- rownames(get(dat_ind[c_runs]))
        for(iExclGenes in rownames(get(dat_ind[c_runs]))){
            datn_WSS[[iExclGenes]] <-  sum( apply(get(dat_ind[c_runs])[iExclGenes,,],2,function(y_vec){
              t((tmpp1[iExclGenes,] - y_vec)^2) %*% (1/sigma_mat[iExclGenes,]^2)}))
        }
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
        tempi1  <- do.call(rbind,lapply(1:(length(imp_res_list)), function (x) imp_res_list[[x]][[1]]))
        tempi2  <- do.call(rbind,lapply(1:(length(imp_res_list)), function(x) imp_res_list[[x]][[2]]))
      }

      imp_res$impulse_parameters <- tempi1[g_names,]
      imp_res$impulse_fits <- tempi2[g_names,]
      
    ####################  
    # 2) Fit to clusters
    } else if(is.null(cluster_results) & is.null(start_values)) {
      # Transform input matrix into 3D array
      centroid_3Darray <- array(0,c(nrow(data_input[[c_runs]]),ncol(data_input[[c_runs]]),1))
      centroid_3Darray[,,1] <- data_input[[c_runs]]
      rownames(centroid_3Darray) <- rownames(data_input[[c_runs]])
      colnames(centroid_3Darray) <- colnames(data_input[[c_runs]])
      # Fit impulse model to the clustering centroids:
      imp_res <-  impulse_fit_matrix(centroid_3Darray,
        as.numeric(as.character(data_annotation[colnames(data_input[[c_runs]]),"Time"])),
        sigma_mat = sigma_mat,n_it = n_iter, ctrl_tc = control_timecourse, 
        ctrl = control_name, fit_to_clus = TRUE, fit_bg = fit_backg, n_process = n_proc)
    }
    names(imp_res) <- paste(c("impulse_parameters","impulse_fits"),label[c_runs],sep="_")
    results[[names(imp_res)[1]]] <- imp_res[[1]]
    results[[names(imp_res)[2]]] <- imp_res[[2]]
  }
  
  return(results)
}

