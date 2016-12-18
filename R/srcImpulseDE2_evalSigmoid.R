#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++     Impulse model value prediction    ++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute value of sigmoid function given parameters.
#' 
#' Compute value of sigmoid function given parameters.
#' 
#' @aliases evalSigmoid_comp
#' 
#' @seealso 
#' 
#' @param vecSigmoidParam (numeric vector number of sigmoid model parameters)
#'  { beta, h0, h1, t1 }
#'  Vector of impulse parameters.
#' @param vecTimepoints (numeric vector number vecTimepoints) 
#'    Time points to be evaluated.
#' 
#' @return vecSigmoidValue (numeric vector length of vecTimepoints) 
#'    Model expression values of given gene for time points
#' @export

evalSigmoid <- function(vecSigmoidParam,
                        vecTimepoints){
  
  # beta is vecSigmoidParam[1]
  # h0 is vecSigmoidParam[2]
  # h1 is vecSigmoidParam[3]
  # t1 is vecSigmoidParam[4]
  
  vecSigmoidValue <- sapply(vecTimepoints, function(t){
    vecSigmoidParam[2] + (vecSigmoidParam[3]-vecSigmoidParam[2])*
      (1/(1+exp(-vecSigmoidParam[1]*(t-vecSigmoidParam[4]))))
  })
  
  # Catch lower bound on mu space
  vecSigmoidValue[vecSigmoidValue < 10^(-10)] <- 10^(-10)
  
  return(vecSigmoidValue)
}