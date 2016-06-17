#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++     Impulse model value prediction    ++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Compute value of impulse function given parameters.

# INPUT:
#   theta:  (vector number of parameters [6]) Impulse model parameters
#   timepoints: (vec number of timepoints) Time-points at which gene was sampled
# OUTPUT:
#   y_vec:  (vec number of timepoints) Model expression values of given gene 
#           for time points

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Impulse model value prediction
#'
#' Calculates impulse model values for given timepoints and predicted
#' impulse parameters.
#' @aliases calc_impulse
#' @param theta numerical vector of impulse parameters with the order
#' beta, h0, h1, h2, t1, t2.
#' @param timepoints numercial vector of time point(s).
#' @return The predicted impulse model values for the given time point(s).
#' @export
calc_impulse <- function(theta,timepoints){
  beta1 = theta[1]
  h0 = theta[2]
  h1 = theta[3]
  h2 = theta[4]
  t1 = theta[5]
  t2 = theta[6]
  
  # DAVID replaced this by construct usde in src_Impulse_CostFuntionFit.R
  #res = NULL
  #for (i in 1:length(timepoints)){
  #  res[i] = (1/h1) * (h0 + (h1-h0) * (1/(1+exp(-beta1*(timepoints[i]-t1))))) *
  #    (h2 + (h1-h2) * (1/(1+exp(beta1*(timepoints[i]-t2)))))
  #}
  #res = unlist(res)
  #  names(res) <- timepoints
  y_vec = unlist(lapply(timepoints, function(x) {(1/h1) * 
      (h0 + (h1-h0) * (1/(1+exp(-beta1*(x-t1))))) *
      (h2 + (h1-h2) * (1/(1+exp(beta1*(x-t2)))))}))
  
  return(y_vec)
}