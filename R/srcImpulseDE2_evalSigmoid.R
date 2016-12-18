#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++     Impulse model value prediction    ++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute value of impulse function given parameters.
#' 
#' Amplitude parameters are assumed to be in log space.
#' 
#' @aliases calcImpulse_comp
#' 
#' @seealso Called by \code{evalLogLikImpulse},\code{evalLogLikMean}, 
#'    \code{plotDEGenes}.
#' 
#' @param vecSigmoidParam (vector number of impulse model parameters)
#'  { beta, h0, h1, h2, t1, t2 }
#'  Vector of impulse parameters.
#' @param vecTimepoints (vector number vecTimepoints) 
#'    Observed vecTimepoints, numeric.
#' 
#' @return vecImpulseValue (vec number of vecTimepoints) 
#'    Model expression values of given gene for time points
#' @export

evalSigmoid <- function(vecSigmoidParam,
                        vecTimepoints){
  
  # beta is vecSigmoidParam[1]
  # h0 is vecSigmoidParam[2]
  # h1 is vecSigmoidParam[3]
  # t1 is vecSigmoidParam[4]
  
  vecImpulseValue <- sapply(vecTimepoints, function(t){
    (1/vecSigmoidParam[3]) * 
      (vecSigmoidParam[2] + (vecSigmoidParam[3]-vecSigmoidParam[2])*
         (1/(1+exp(-vecSigmoidParam[1]*(t-vecSigmoidParam[5]))))) *
      (vecSigmoidParam[4] + (vecSigmoidParam[3]-vecSigmoidParam[4])*
         (1/(1+exp(vecSigmoidParam[1]*(t-vecSigmoidParam[6])))))
  })
  
  # Catch lower bound on mu space
  vecImpulseValue[vecImpulseValue < 10^(-10)] <- 10^(-10)
  
  return(vecImpulseValue)
}