#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++     Cost Functions    +++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Cost function for parametric model fit: log likelihood
### This cost function is used for fitting individual genes with replicates

# INPUT:
#   vecTheta:  (vector number of parameters [6]) Impulse model parameters
#   vecX:  (vector number of timepoints) Time-points at which gene was sampled
#   matY:  (2D array timepoints x replicates) Observed expression values for 
#     given gene.
#   scaDispEst: (scalar) Dispersion estimate for given gene.
# OUTPUT:
#   scaLogLik: (scalar) Value of cost function (likelihood) for given gene

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# Cost function for fitting impulse model
evalLogLikImpulse <- function(vecTheta,vecX,matY,scaDispEst){  
  # Compute function value - mean
  vecImpulseValue = calcImpulse_comp(vecTheta,vecX)
  
  # Without dispersion estimation, take given
  scaLogLikImpulse <- sum(unlist( lapply(c(1:length(vecX)), function(tp){
    sum(dnbinom(matY[tp,], mu=vecImpulseValue[tp], size=scaDispEst, log=TRUE))}) ))
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}

# Cost function for optim for mean fitting without impulse model
# Only the mean (h0) is fitted in this case
evalLogLikMean <- function(scaMuEst,matY,scaDispEst){
  # Compute log likelihood assuming constant mean of negative binomial
  scaLogLikMean <- sum (dnbinom(matY, mu=exp(scaMuEst), size=scaDispEst, log=TRUE))
  
  # Maximise log likelihood: Return likelihood as value to optimisation routine
  return(scaLogLik)
}
