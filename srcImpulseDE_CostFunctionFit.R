#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++     Cost Functions    +++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# Ordinary least squares for 1D array (genes):
#   cost_fun_OLS
# Weighted least squares for 2D array (genes x replicates):
#   cost_fun_WLS

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Cost function for parametric model fit: Ordinary least squares
### This cost function is used for fitting cluster centroids. (1 replicate)

# INPUT:
#   theta:  (vector number of parameters [6]) Impulse model parameters
#   x_vec:  (vector number of timepoints) Time-points at which gene was sampled
#   y_mat:  (2D array timepoints x replicates) Observed expression values for 
#     given gene.
# OUTPUT:
#   f_SS: (scalar) Value of cost function for given gene

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

cost_fun_OLS <- function(theta,x_vec,y_mat,...){
  # Model parameters
  beta1 = theta[1]
  h0 = theta[2]
  h1 = theta[3]
  h2 = theta[4]
  t1 = theta[5]
  t2 = theta[6]
  
  # Cost function: ordinary least squares
  f_SS = sum(apply(y_mat,2,function(y_vec){
      sum((unlist(lapply(x_vec, function(x) {(1/h1) * (h0 + (h1-h0) *
      (1/(1+exp(-beta1*(x-t1))))) * (h2 + (h1-h2) *
      (1/(1+exp(beta1*(x-t2)))))})) - y_vec)^2)}))
  
  return(f_SS)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Cost function for parametric model fit: weighted least squares
### This cost function is used for fitting individual genes with replicates

# INPUT:
#   theta:  (vector number of parameters [6]) Impulse model parameters
#   x_vec:  (vector number of timepoints) Time-points at which gene was sampled
#   y_mat:  (2D array timepoints x replicates) Observed expression values for 
#     given gene.
#   weight_vec: (vector number of timepoints) Weights of data points for objective.
# OUTPUT:
#   f_WSS: (scalar) Value of cost function for given gene

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

cost_fun_WLS <- function(theta,x_vec,y_mat,weight_vec,...){
  
  # Cost function: weighted least squares
  impulse_value = calc_impulse_comp(theta,x_vec)

  f_WSS = sum(apply(y_mat,2,function(y_vec){
    sum(((impulse_value - y_vec)^2) * (1/weight_vec^2))}))
  # Take out weighting of SS:
  #f_WSS = sum(apply(y_mat,2,function(y_vec){
  #  sum(((impulse_value - y_vec)^2))}))
  
  return(f_WSS)
}

# objective is log likelihood
cost_fun_logl <- function(theta,x_vec,y_mat,disp_est,...){
  # weight_vec not used
  
  # Compute function value - mean
  impulse_value = calc_impulse_comp(theta,x_vec)
  
  # Get the maximum likelihood negative binomial model
  # at each time point by maximising over the dispersion
  # coefficient. This routine returns the objective:
  # log likelihood of the model (at the MLE).
  start <- list()
  start$size <- 0.1
  lower <- list()
  lower$size <- 0.0001
  
  # With estimation of dispersion
  #nb_estimation <- lapply(c(1:length(x_vec)), function(tp){
  #  nb_estimate <- try( fitdistr(y_mat[tp,], "Negative Binomial", method="L-BFGS-B", mu=impulse_value[tp], start=start, lower=lower), silent=TRUE );
  #  if(is.na(nb_estimate[2])){nb_estimate_fake <- list(); nb_estimate_fake$loglik <- c(-99); return(nb_estimate_fake)}
  #  else {return(nb_estimate)} } )
  
  # Get likelihood of model: Add up log likelihoods of timepoints
  #logl_impulse <- sum( unlist(lapply(nb_estimation,function(x){x$loglik})) )
  
  # Without dispersion estimation, take given
  logl_tps <- unlist( lapply(c(1:length(x_vec)), function(tp){
    prod(dnbinom(y_mat[tp,], mu=impulse_value[tp], size=disp_est))}) )
  logl_impulse <- log(prod( logl_tps ))
  
  # Maximise log likelihood
  return(-logl_impulse)
}