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
#' @param vecImpulseParam (vector number of impulse model parameters)
#'  { beta, h0, h1, h2, t1, t2 }
#'  Vector of impulse parameters.
#' @param vecTimepoints (vector number vecTimepoints) 
#'    Observed vecTimepoints, numeric.
#' 
#' @return vecImpulseValue (vec number of vecTimepoints) 
#'    Model expression values of given gene for time points
#' @export

calcImpulse <- function(vecImpulseParam,
                        vecTimepoints){
  
  # beta is vecImpulseParam[1]
  # h0 is vecImpulseParam[2]
  # h1 is vecImpulseParam[3]
  # h2 is vecImpulseParam[4]
  # t1 is vecTheta[5]
  # t2 is vecTheta[6]
  
  vecImpulseValue <- sapply(vecTimepoints, function(t){
    (1/vecImpulseParam[3]) * 
      (vecImpulseParam[2] + (vecImpulseParam[3]-vecImpulseParam[2])*
         (1/(1+exp(-vecImpulseParam[1]*(t-vecImpulseParam[5]))))) *
      (vecImpulseParam[4] + (vecImpulseParam[3]-vecImpulseParam[4])*
         (1/(1+exp(vecImpulseParam[1]*(t-vecImpulseParam[6])))))
  })
  
  # Catch lower bound on mu space
  vecImpulseValue[vecImpulseValue < 10^(-10)] <- 10^(-10)
  
  return(vecImpulseValue)
}