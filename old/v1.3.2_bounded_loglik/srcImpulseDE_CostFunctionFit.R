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
#   logl_impulse: (scalar) Value of cost function for given gene

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

cost_fun_logl <- function(theta,x_vec,y_mat,disp_est){  
  # Compute function value - mean
  impulse_value = calc_impulse_comp(theta,x_vec)
  
  # Without dispersion estimation, take given
  logl_impulse <- sum(unlist( lapply(c(1:length(x_vec)), function(tp){
    sum(dnbinom(y_mat[tp,], mu=impulse_value[tp], size=disp_est, log=TRUE))}) ))
  
  # log(10^(-323))=-743.7469 is the smallest number not to give an
  # error in log. -800 is smaller than that - sets the boundary
  # for optimisation to move away from.
  #if(is.na(logl_impulse) | !is.finite(logl_impulse)){logl_impulse <- -800}
  
  # Maximise log likelihood, optim minimises objective -> take negative
  return(-logl_impulse)
}

# Cost function for optim for mean fitting without impulse model
# Only the mean (h0) is fitted in this case
cost_fun_logl_meanfit <- function(mu_est,y_mat,disp_est){
  # Compute log likelihood assuming constant mean of negative binomial
  logl_impulse <- sum (dnbinom(y_mat, mu=exp(mu_est), size=disp_est, log=TRUE))
  
  # log(10^(-323))=-743.7469 is the smallest number not to give an
  # error in log. -800 is smaller than that - sets the boundary
  # for optimisation to move away from.
  #if(is.na(logl_impulse) | !is.finite(logl_impulse)){logl_impulse <- -800}
  
  # Maximise log likelihood, optim minimises objective -> take negative
  return(-logl_impulse)
}

# To be deprecated: Fit mean as impulse model, below: fit mean as constant
# Cost function for optim for mean fitting:
# Only mu_est (the mean) is fitted in this case
cost_fun_logl_meanimpulsefit <- function(mu_est,x_vec,y_mat,disp_est){  
  
  # Construct theta vector:
  # Beta: high to keep slope in modelled range small
  # h0: represents mean
  # h1,h2: close to h0, lie outside modelled data range
  # t1,t2: outside modelled range
  theta <- c(100,mu_est,mu_est+1,mu_est,
    3*x_vec[length(x_vec)],4*x_vec[length(x_vec)])
  # Compute function value - mean
  impulse_value = calc_impulse_comp(theta,x_vec)
  
  # Without dispersion estimation, take given
  lik_tps <- unlist( lapply(c(1:length(x_vec)), function(tp){
    prod(dnbinom(y_mat[tp,], mu=impulse_value[tp], size=disp_est))}) )
  logl_impulse <- log(prod( lik_tps ))
  
  # log(10^(-323))=-743.7469 is the smallest number not to give an
  # error in log. -800 is smaller than that - sets the boundary
  # for optimisation to move away from.
  if(is.na(logl_impulse) || !is.finite(logl_impulse)){logl_impulse <- -800}
  
  # Maximise log likelihood, optim minimises objective -> take negative
  return(-logl_impulse)
}
