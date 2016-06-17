#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++     Calculation of input for Optimization    ++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Cost function for parametric model fit: Ordinary least squares (two_impulses)
### or weighted least squares (two_impulses_WLS)

# INPUT:
#   theta: (vector number of parameters [6]) Impulse model parameters
#   x_vec: (vec number of timepoints) Time-points at which gene was sampled
#   y-vec: (vec number of timepoints) Observed expression values of given gene for time points
# OUTPUT:
#   f_SS/f_WSS: (scalar) Value of cost function for given gene

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# Objective for ordinary least sqaures: sum of squares
two_impulses <- function(theta,x_vec,y_vec){
  # Model parameters
  beta1 = theta[1]
  h0 = theta[2]
  h1 = theta[3]
  h2 = theta[4]
  t1 = theta[5]
  t2 = theta[6]
  
  # Cost function: ordinary least squares
  f_SS = sum((unlist(lapply(x_vec, function(x) {(1/h1) * (h0 + (h1-h0) *
          (1/(1+exp(-beta1*(x-t1))))) * (h2 + (h1-h2) *
          (1/(1+exp(beta1*(x-t2)))))})) - y_vec)^2)
  
  return(f_SS)
}

# Objective for weighted least sqaures: weighted sum of squares
two_impulses_WLS <- function(theta,x_vec,y_vec,sigma_vec){
  # Model parameters
  beta1 = theta[1]
  h0 = theta[2]
  h1 = theta[3]
  h2 = theta[4]
  t1 = theta[5]
  t2 = theta[6]
  
  # Cost function: weighted least squares
  f_SS = sum((unlist(lapply(x_vec, function(x) {(1/h1) * (h0 + (h1-h0) *
          (1/(1+exp(-beta1*(x-t1))))) * (h2 + (h1-h2) *
          (1/(1+exp(beta1*(x-t2)))))})) - y_vec)^2)
  # Weight sum of squares by variance
  f_WSS = f_SS/sigma_vec^2
  
  return(f_WSS)
}