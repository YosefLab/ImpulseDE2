### Compute value of sigmoidal model

#' Compute value of sigmoidal model given parameters.
#' 
#' Compute value of sigmoidal model given parameters.
#' Enforces lower bound on value of function to avoid numerical
#' errors during model fitting.
#' 
#' @seealso Compiled version: \link{evalSigmoid_comp}
#' 
#' @param vecSigmoidParam (numeric vector number of sigmoid model parameters)
#' \{beta, h0, h1, t1\}
#' Vector of sigmoidal model parameters.
#' @param vecTimepoints (numeric vector length number of time points) 
#' Time points to be evaluated.
#' 
#' @return vecSigmoidValue (numeric vector length of vecTimepoints) 
#' Model values for given time points.
#' 
#' @author David Sebastian Fischer
evalSigmoid <- function(vecSigmoidParam, vecTimepoints) {
    
    # beta is vecSigmoidParam[1] h0 is vecSigmoidParam[2] h1 is
    # vecSigmoidParam[3] t1 is vecSigmoidParam[4]
    
    vecSigmoidValue <- sapply(vecTimepoints, function(t) {
        vecSigmoidParam[2] + (vecSigmoidParam[3] - vecSigmoidParam[2]) * 
            (1/(1 + exp(-vecSigmoidParam[1] * (t - vecSigmoidParam[4]))))
    })
    vecSigmoidValue[vecSigmoidValue < 10^(-10)] <- 10^(-10)
    
    return(vecSigmoidValue)
}

#' Compiled function: evalSigmoid
#' 
#' Pre-compile heavily used functions.
#' Refer to \link{evalSigmoid}.
#' 
#' @param vecSigmoidParam (numeric vector number of sigmoid model parameters)
#' \{beta, h0, h1, t1\}
#' Vector of sigmoidal model parameters.
#' @param vecTimepoints (numeric vector length number of time points) 
#' Time points to be evaluated.
#' 
#' @return vecSigmoidValue (numeric vector length of vecTimepoints) 
#' Model values for given time points.
#' 
#' @author David Sebastian Fischer
evalSigmoid_comp <- cmpfun(evalSigmoid)
