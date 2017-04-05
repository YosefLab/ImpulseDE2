### Compute value of impulse model

#' Compute value of impulse function given parameters.
#' 
#' Compute value of impulse function given parameters.
#' Enforces lower bound on value of function to avoid numerical
#' errors during model fitting.
#' 
#' @seealso Compiled version: \link{evalImpulse_comp}
#' 
#' @param vecImpulseParam (numeric vector number of impulse model parameters)
#' \{beta, h0, h1, h2, t1, t2\}
#' Vector of impulse model parameters.
#' @param vecTimepoints (numeric vector length number of time points) 
#' Time points to be evaluated.
#' 
#' @return vecImpulseValue (vec number of vecTimepoints) 
#' Model values for given time points.
#' 
#' @author David Sebastian Fischer
evalImpulse <- function(vecImpulseParam, vecTimepoints) {
    
    # beta is vecImpulseParam[1] h0 is vecImpulseParam[2] h1 is
    # vecImpulseParam[3] h2 is vecImpulseParam[4] t1 is vecImpulseParam[5]
    # t2 is vecImpulseParam[6]
    
    vecImpulseValue <- sapply(vecTimepoints, function(t) {
        (1/vecImpulseParam[3]) * 
            (vecImpulseParam[2] + (vecImpulseParam[3] - vecImpulseParam[2]) *
                 (1/(1 + exp(-vecImpulseParam[1] * (t - vecImpulseParam[5]))))) *
            (vecImpulseParam[4] + (vecImpulseParam[3] - vecImpulseParam[4]) *
                 (1/(1 + exp(vecImpulseParam[1] * (t - vecImpulseParam[6])))))
    })
    vecImpulseValue[vecImpulseValue < 10^(-10)] <- 10^(-10)
    
    return(vecImpulseValue)
}

#' Compiled function: evalImpulse
#' 
#' Pre-compile heavily used functions.
#' Refer to \link{evalImpulse}.
#' 
#' @param vecImpulseParam (numeric vector number of impulse model parameters)
#' \{beta, h0, h1, h2, t1, t2\}
#' Vector of impulse model parameters.
#' @param vecTimepoints (numeric vector length number of time points) 
#' Time points to be evaluated.
#' 
#' @return vecImpulseValue (vec number of vecTimepoints) 
#' Model values for given time points.
#' 
#' @author David Sebastian Fischer
evalImpulse_comp <- cmpfun(evalImpulse)
