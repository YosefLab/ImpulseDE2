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

### Developmentat note:
### Three-way separation was chosen to use 'apply' instead of for-loops

# INPUT:
#   data_input: (Numeric matrix genes x samples) Expression values 
#   data_annotation: (Table samples x 2[time and condition]) 
#       Co-variables for the samples including condition and time points.
#       Time points must be numeric.
#   sigma_vec
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
#   start_values:
#   fit_backg: (bool) [Defaul FALSE]
#   n_process: (scalar) [Default 4] number of processes, which can be used on 
#       the machine to run the background calculation in parallel
# OUTPUT:
#   results: (list runs*2 ["impulse_parameters","impulse_fits"])
#       impulse_parameters: (matrix genes x (NPARAM+1)*runs) +1 is value of 
#           objective of optimisation (i.e. sum of squares or weighted SS)
#       impulse_fits: (matrix genes x runs) model values for gene data

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

impulse_fit <- function(data_input, data_annotation, sigma_vec, n_iter = 100,
    control_timecourse = FALSE, control_name = NULL, cluster_results = NULL,
    start_values = NULL, fit_backg = FALSE, n_proc = 4, ...){
  
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
  #   vec: (vector) 
  #   data_annotation: (Table samples x 2[time and condition]) 
  #       Co-variables for the samples including condition and time points.
  #       Time points must be numeric.
  #   sigma_vec:
  #   NPARAM: (scalar) Number of model parameters
  #   n_iter: (scalar) [Defaul 100] Number of iterations, which are performed 
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
  
  impulse_fit_gene_wise <- function(vec, timepoints, sigma_vec, NPARAM, n_iters = 100,
         fit_to_clusters = FALSE, start_val = NULL, fit_bg = FALSE, ...){
    
    # Initial parameter guesses:
    # beta: slope
    # h0, h2: first and last value 
    # h1: peak value
    # t1, t2: peak time difference to ends
    mat1 = c(0,abs(vec[2:length(vec)] - vec[1:(length(vec)-1)]) *
               vec[2:length(vec)])
    mn_beta = 0
    mx_beta = 10
    middle_ind = which(timepoints > min(timepoints) & timepoints < max(timepoints))
    min_ind = which(timepoints == min(timepoints))[1]
    max_ind = which(timepoints == max(timepoints))[1]
    mx_tm = max(timepoints)
    mn_tm = min(timepoints)
    tmp_ind = which(mat1 == max(mat1))[1]
    peaktime = timepoints[tmp_ind]
    peak = min(which(timepoints == peaktime))
    beta1 = abs(log(abs((vec[peak]-vec[min_ind])/(peaktime-
               timepoints[min_ind]))));
    
    ### set beta to 1 if calculcated value is infinite
    if(is.finite(beta1) == FALSE || is.na(beta1)) { beta1 = 1 }
    
    ### define start value theta based on the gene's data
    orig_theta = c(beta1,
                   vec[min_ind],                    # h0
                   vec[peak],                       # h1
                   vec[max_ind],                    # h2
                   (peaktime-timepoints[min_ind])/2,            # t1
                   peaktime+(timepoints[max_ind] - peaktime)/2) # t2
    names(orig_theta) = c("beta","h0","h1","h2","t1","t2")
    theta = orig_theta
    
    ### replace theta estimates of 0 by small value, because optimization
    ### function has problems with 0 as input
    rand_num <- runif(length(which(theta == 0)),0,0.0001)
    theta[names(which(theta == 0))] <- rand_num
    names(theta) = c("beta","h0","h1","h2","t1","t2")
    
    ### ---> fit model to clusters
    ### Note that this splits into fit to cluster = TRUE and FALSE
    if(fit_to_clusters == TRUE & is.null(start_val)){
      
      ### DAVID what is this, i ll ignore it, handing means above so should be fine
      if(fit_bg == FALSE){
        thetas <- cbind(as.numeric(theta),
          apply(matrix(rep(runif((n_iters-1)*NPARAM)),NPARAM,(n_iters-1)),
            2, function(x){c(mn_beta + (mx_beta-mn_beta)*x[1],
                      min(vec) + (max(vec)-min(vec))*x[2:4],
                      mn_tm + (mx_tm-mn_tm)*x[5],
                      (mn_tm + (mx_tm-mn_tm)*x[5]) + (mx_tm-(mn_tm + (mx_tm-mn_tm)*x[5]))*x[6])}))
      } else {
        thetas <- apply(matrix(rep(runif((n_iters)*NPARAM)),NPARAM,(n_iters)),
            2, function(x){c(mn_beta + (mx_beta-mn_beta)*x[1],
                      min(vec) + (max(vec)-min(vec))*x[2:4],
                      mn_tm + (mx_tm-mn_tm)*x[5],
                      (mn_tm + (mx_tm-mn_tm)*x[5]) + (mx_tm-(mn_tm + (mx_tm-mn_tm)*x[5]))*x[6])})
      }
      tmm1 <- system.time({
        fmin_outs <- cbind(
          # fitting via quasi-Newton method (optim, "BFGS")
          apply(thetas[,1:(floor(n_iters/2))],
            2,function(x){unlist(optim(x,two_impulses_WLS, x_vec = timepoints,
            y_vec = vec, sigma_vec = sigma_vec, method="BFGS")[c("par","value")])}),
          # fitting via PORT routines (nlminb)
          apply(thetas[, c(1,(floor(n_iters/2)+1):n_iters)],
            2,function(x){unlist(nlminb(x, two_impulses_WLS, x_vec = timepoints,
            y_vec = vec, sigma_vec = sigma_vec)[c("par","objective")])})
        )
      })
      
      
      ### return results of 3 best fits, which will serve as start values for
      ### the fits to the genes
      pvec_and_SSE = as.vector(fmin_outs[,order(fmin_outs["value",])][,1:3])
      
      ### ---> fit model to genes
    } else if(fit_to_clusters == FALSE & is.null(start_val) == FALSE){
      tmm2 <- system.time({
        
        ### use first 3 best hits from fit to clusters plus initial guess 'theta'
        ### as start values
        if(fit_bg == TRUE){
          toll <- matrix(start_val,7,3)[1:NPARAM,]
        } else {
          toll <- cbind(theta,matrix(start_val,7,3)[1:NPARAM,][,1:3])
        }
        fmin_outs <- cbind(
          # fitting via quasi-Newton method (optim, "BFGS")
          ### toll has the parameter intialisation in it
          apply(toll,
            2,function(x){unlist(optim(x, two_impulses_WLS, x_vec = timepoints,
            y_vec = as.numeric(vec),method="BFGS" )[c("par","value")])}),
          # again fitting via PORT routines (nlminb)
          apply(toll,
            2,function(x){unlist(nlminb(x, two_impulses_WLS, x_vec = timepoints,
            y_vec = as.numeric(vec) )[c("par","objective")])})
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
  #   data_tab: (Numeric matrix genes x samples) Expression values of genes in
  #       cluster or entire matrix, depending on which fitting mode runs.
  #   timepoints: (Table samples x 2 [time and condition]) 
  #       Co-variables for the samples including condition and time points.
  #       Time points must be numeric.
  #   sigma_vec:
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
  #       impulse_fits: (matrix genes x runs) model values for gene data
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  impulse_fit_matrix <- function(data_tab, timepoints, sigma_vec, n_it = 100,
    ctrl_tc = FALSE, ctrl = NULL, fit_to_clus = FALSE, start_val = NULL,
    fit_bg = FALSE, n_process = 4, ...){
    
    NPARAM = 6 # Number of model parameters
    mat1 = NULL
    res <- list()
    for (j in 2:ncol(data_tab)){
      mat1 = cbind(mat1, abs(data_tab[,j] - data_tab[,j-1]) * data_tab[,j])
    }
    
    mc <- min(detectCores() - 1, n_process)
    # Don't split if only have few genes
    if(fit_bg == FALSE && nrow(data_tab) > max(2*mc,10)){
      # Use parallisation
      # Define partitioning of genes onto nodes: ind_list
      ind_list = list()
      bord = floor(nrow(data_tab)/mc)
      for (i in 1:mc){
        ind_list[[i]] <- ((i-1)*bord+1):(i*bord)
      }
      if(mc*bord < nrow(data_tab)){
        # Add remaining genes to last node
        ind_list[[mc]] <-  c(ind_list[[mc]],(mc*bord+1):nrow(data_tab))
      }
      
      cl <- makeCluster(mc, outfile = "clus_out2.txt")
      my.env <- new.env()
      assign("data_tab", data_tab, envir = my.env)
      assign("calc_impulse_comp", calc_impulse_comp, envir = my.env)
      assign("n_it", calc_impulse_comp, envir = my.env)
      assign("cluster_genes_for_impulse", cluster_genes_for_impulse,
             envir = my.env)
      assign("impulse_fit", impulse_fit, envir = my.env)
      assign("two_impulses_WLS", two_impulses_WLS, envir = my.env)
      assign("impulse_fit_gene_wise",impulse_fit_gene_wise, envir = my.env)
      assign("impulse_fit_matrix",impulse_fit_matrix, envir = my.env)
      assign("timepoints", timepoints, envir = my.env)
      assign("fit_to_clus", fit_to_clus, envir = my.env)
      assign("start_val", start_val, envir = my.env)
      assign("fit_bg", fit_bg, envir = my.env)
      assign("NPARAM", NPARAM, envir = my.env)
      
      clusterExport(cl=cl, varlist=c("data_tab",
        "calc_impulse_comp","n_it","cluster_genes_for_impulse","impulse_fit",
        "two_impulses_WLS", "impulse_fit_gene_wise",
        "impulse_fit_matrix", "timepoints","fit_to_clus","start_val",
        "fit_bg","NPARAM"), envir = my.env)
      
      # Fit impulse model to each gene of matrix and get impulse parameters:
      # clusterApply runs the function impulse_fit_gene_wise.
      # The input data are distributed to nodes by ind_list partitioning
      resi <- clusterApply(cl, ind_list, function(z){t(apply(data_tab[z,],1,
                function(x){impulse_fit_gene_wise(x,timepoints,sigma_vec,
                NPARAM,n_it,fit_to_clus, start_val, fit_bg)}))})
      # Concatenate the output objects of each node
      resmat = do.call(rbind,resi)
      resmat <- resmat[rownames(data_tab),]
      stopCluster(cl)
      # Parallelisation ends here
    } else {
      # Do not use parallelisation
      # fit impulse model to each gene of matrix and get impulse parameters
      resmat <- t(apply(data_tab,1,function(x){impulse_fit_gene_wise(x,
                  timepoints,sigma_vec,NPARAM,n_it,fit_to_clus, start_val, fit_bg)}))
    }   
    
    ### use obtained impulse parameters to calculate impulse fit values
    ### ---> if fit to genes
    if(fit_to_clus == FALSE){
      colnames(resmat) <- c("beta","h0","h1","h2","t1","t2","SSE")
      if(nrow(resmat) == 1){      # if matrix contains only one gene
        resmat2 <- as.data.frame(t(calc_impulse_comp(resmat[,1:NPARAM],
          unique(sort(timepoints)))),
          row.names=rownames(resmat))
        colnames(resmat2) = unique(sort(timepoints))
      } else {                    # if matrix contains > 1 genes
        resmat2 <- t(apply(resmat[,1:NPARAM],1,function(x){calc_impulse_comp(x,
          unique(sort(timepoints)))}))
        colnames(resmat2) = unique(sort(timepoints))
      }
      
      ### ---> if fit to clusters
      ### more complex because results from 3 best fits need to be saved and later
      ### used for the fit to the genes
    } else {
      colnames(resmat) <- c(paste(c("beta","h0","h1","h2","t1","t2","SSE"),1,sep="_"),
                            paste(c("beta","h0","h1","h2","t1","t2","SSE"),2,sep="_"),
                            paste(c("beta","h0","h1","h2","t1","t2","SSE"),3,sep="_"))
      resmat2 <- t(apply(resmat,
        1,function(x){apply(matrix(x,7,3)[1:NPARAM,],
        2,function(y){calc_impulse_comp(y,unique(sort(timepoints)))})}))
      colnames(resmat2) <- c(paste(unique(sort(timepoints)),1,sep="_"),
                             paste(unique(sort(timepoints)),2,sep="_"),
                             paste(unique(sort(timepoints)),3,sep="_"))
    }
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
  #   a) A list of three data matrices combined, case and control
  #   b) A list of a single data matrix if only have case
  # data_input takes data from 
  ### ---> fit to genes
  if(is.null(cluster_results) == FALSE & is.null(start_values) == FALSE){
    file_ext = "genes"
    c1 = "genes"
    
    ### perform 3 runs if control timecourse data is present
    if(control_timecourse == TRUE){
      ### generate a data list containing the combined, case and control data if
      ### control timecourse is present
      # Note on variables used below:
      # cluster_results[[1,5,9]] are kmeans_clus for [combined, case, control]
      
      # Combined: Expression values of both case and control (all columns).
      # Clustering based on combined.
      dat1 = as.data.frame(data_input)[names(cluster_results[[1]]),]
      dat1_n = as.data.frame(data_input)[!(rownames(data_input) %in% rownames(dat1)),]
      dat1_n = dat1_n[,colnames(dat1)]
      
      # Case: Expression values of case. Clustering based on case.
      dat2 = as.data.frame(data_input)[,!(data_annotation$Condition %in% control_name)]
      dat2 = dat2[names(cluster_results[[5]]),]
      dat2_n = as.data.frame(data_input)[!(rownames(data_input) %in% rownames(dat2)),]
      dat2_n = dat2_n[,colnames(dat2)]
      
      # Control: Expression values of control. Clustering based on control. 
      dat3 = as.data.frame(data_input)[,(data_annotation$Condition %in% control_name)]
      dat3 = dat3[names(cluster_results[[9]]),]
      dat3_n = as.data.frame(data_input)[!(rownames(data_input) %in% rownames(dat3)),]
      dat3_n = dat3_n[,colnames(dat3)]
      
      dat_ind <- c("dat1_n","dat2_n", "dat3_n")
      data_input <- list(dat1, dat2, dat3)
      
    } else if(control_timecourse == FALSE){
      dat1 = as.data.frame(data_input)[names(cluster_results[[1]]),]
      dat1_n = as.data.frame(data_input)[!(rownames(data_input) %in% rownames(dat1)),]
      
      dat_ind <- c("dat1_n")
      data_input <- list(dat1)
    }
    
    
    ### ---> fit to clusters
    ### DAVID BUG: missing option for no control ?? 
  }  else if(is.null(cluster_results) & is.null(start_values)) {
    file_ext = "clusters"
    c1 = "cluster"
    ### Input are the centroids of the clusters
    ### 2 is combined, 6 is case and 10 is control
    ### if runs == 1 or runs == 2, then 6 and 10 or just 10 are empty
    data_input <- data_input[c(2,6,10)[1:runs]]
  }
  
  ### fitting for different runs
  for (c_runs in 1:runs){
    imp_res <- NULL
    
    ### ---> fit to genes
    if(is.null(cluster_results) == FALSE & is.null(start_values) == FALSE){
      
      ### split genes into the clusters and fit impulse model to the genes of a cluster
      # Split expression matrix of current run into matrices for each cluster
      cluster_res_list <- split(data_input[[c_runs]],cluster_results[[(c_runs-1)*4 +1]])
      ind <- split(1:max(cluster_results[[(c_runs-1)*4 +1]]), 1:max(cluster_results[[(c_runs-1)*4 +1]]))
      
      # Fit model to each cluster matrix
      imp_res_list <- lapply(ind, function(x){impulse_fit_matrix(cluster_res_list[[x]],
         as.numeric(as.character(data_annotation[colnames(data_input[[c_runs]]),"Time"])), 
         sigma_vec = sigma_vec,
         n_it = n_iter, ctrl_tc = control_timecourse, ctrl = control_name, fit_to_clus = FALSE,
         start_val = start_values[[c_runs]][x,], fit_bg = fit_backg, n_process = n_proc)})
      
      if(nrow(get(dat_ind[c_runs])) != 0){
        tump1  <- do.call(rbind,lapply(1:(length(imp_res_list)), function(x) imp_res_list[[x]][[1]]))
        tump2  <- do.call(rbind,lapply(1:(length(imp_res_list)), function(x) imp_res_list[[x]][[2]]))
        
        tmpp1 <- t(apply(get(dat_ind[c_runs]),1,function(x){rep(mean(x),ncol(tump2))}))
        colnames(tmpp1) <- colnames(tump2)
        tump2 <- rbind(tump2, tmpp1)
        tmpp2 <-  cbind( matrix(NA,nrow(get(dat_ind[c_runs])),ncol(tump1)-1),
                         apply(cbind(get(dat_ind[c_runs]), tmpp1[,1]),
                            1,function(x){sum((x[1:(length(x)-1)] - x[length(x)])^2)}))
        rownames(tmpp2) <- rownames(get(dat_ind[c_runs]))
        tump1 <- rbind(tump1, tmpp2)
        
        tempi1 <- tump1
        tempi2 <- tump2
      } else {
        tempi1  <- do.call(rbind,lapply(1:(length(imp_res_list)), function (x) imp_res_list[[x]][[1]]))
        tempi2  <- do.call(rbind,lapply(1:(length(imp_res_list)), function(x) imp_res_list[[x]][[2]]))
      }
      imp_res$impulse_parameters <- tempi1[g_names,]
      imp_res$impulse_fits <- tempi2[g_names,]
      
      ### ---> fit to clusters
    } else if(is.null(cluster_results) & is.null(start_values)) {
      ### fit impulse model to the mean of the cluster
      # Fit model to entire data matrix, plain running of impulse_fit_matrix
      imp_res <-  impulse_fit_matrix(data_input[[c_runs]],
        as.numeric(as.character(data_annotation[colnames(data_input[[c_runs]]),"Time"])),
        sigma_vec = sigma_vec,
        n_it = n_iter, ctrl_tc = control_timecourse, ctrl = control_name, 
        fit_to_clus = TRUE, fit_bg = fit_backg, n_process = n_proc)
    }
    names(imp_res) <- paste(c("impulse_parameters","impulse_fits"),label[c_runs],sep="_")
    results[[names(imp_res)[1]]] <- imp_res[[1]]
    results[[names(imp_res)[2]]] <- imp_res[[2]]
  }
  
  ### Name output files
  ### DAVID: c2 is not used anywhere - take out this line
  if(fit_backg == TRUE){ c2 = "bg" } else {c2 = "data"}
  
  return(results)
  }
