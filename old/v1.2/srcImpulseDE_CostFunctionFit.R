#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++     Cost Functions    +++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# Ordinary least squares for 1D array (genes):
#   two_impulses_OLS
# Weighted least squares for 2D array (genes x replicates):
#   two_impulses_WLS

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

two_impulses_OLS <- function(theta,x_vec,y_mat){
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
#   sigma_vec: (vector number of timepoints) Standard deviation of 
#     observed expression values for given gene.
# OUTPUT:
#   f_WSS: (scalar) Value of cost function for given gene

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

two_impulses_WLS <- function(theta,x_vec,y_mat,sigma_vec){
  
  # Cost function: weighted least squares
  impulse_value = calc_impulse_comp(theta,x_vec)

  f_WSS = sum(apply(y_mat,2,function(y_vec){
    sum(((impulse_value - y_vec)^2) * (1/sigma_vec^2))}))
  # Take out weighting of SS:
  #f_WSS = sum(apply(y_mat,2,function(y_vec){
  #  sum(((impulse_value - y_vec)^2))}))
  
  return(f_WSS)
}