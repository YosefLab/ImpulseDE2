#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++     Cost Functions    +++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Cost function for parametric model fit: log likelihood
### This cost function is used for fitting individual genes with replicates

# INPUT:
#   theta:  (vector number of parameters [6]) Impulse model parameters
#   x_vec:  (vector number of timepoints) Time-points at which gene was sampled
#   y_mat:  (2D array timepoints x replicates) Observed expression values for 
#     given gene.
#   disp_est: (scalar) Dispersion estimate for given gene.
# OUTPUT:
#   logl_impulse: (scalar) Value of cost function (likelihood) for given gene

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# Cost function for fitting impulse model
cost_fun_logl <- function(theta,x_vec,y_mat,disp_est){  
  # Compute function value - mean
  impulse_value = calc_impulse_comp(theta,x_vec)
  
  # Without dispersion estimation, take given
  logl_impulse <- sum(unlist( lapply(c(1:length(x_vec)), function(tp){
    sum(dnbinom(y_mat[tp,], mu=impulse_value[tp], size=disp_est, log=TRUE))}) ))
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(logl_impulse)
}

# Cost function for optim for mean fitting without impulse model
# Only the mean (h0) is fitted in this case
cost_fun_logl_meanfit <- function(mu_est,y_mat,disp_est){
  # Compute log likelihood assuming constant mean of negative binomial
  logl_impulse <- sum (dnbinom(y_mat, mu=exp(mu_est), size=disp_est, log=TRUE))
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(logl_impulse)
}
