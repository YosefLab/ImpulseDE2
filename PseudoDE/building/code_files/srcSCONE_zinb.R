#' Parameter estimation of zero-inflated negative binomial model
#'
#' This function implements an EM algorithm to estimate the parameters of a zero-inflated negative binomial model
#'
#' This function implements an expectation-maximization algorithm for a zero-inflated
#' model, with a constant gene-specific mean parameter, gene-specific dispersion and
#' a mixing probability that depends solely on the mean.
#'
#' @param Y matrix. An integer data matrix (genes in rows, cells in columns)
#' @param maxiter numeric. The maximum number of iterations.
#' @param verbose logical. Whether or not to print the value of the likelihood at each value.
#'
#' @importFrom MASS glm.nb
#' @export
#'
#' @return a list with the following elements:
#' \itemize{
#' \item{mu}{the mean of the negative binomial component}
#' \item{theta}{the dispersion parameter of the negative binomial component}
#' \item{pi}{the mixing probability}
#' \item{coefs}{the coefficients of the logistic regression}
#' \item{p_z}{the probability P(Z=0|Y), useful as weights in weighted principal components analysis}
#' \item{expected_value}{the expected value E[Y]}
#' \item{variance}{the variance Var(Y)}
#' \item{loglik}{the log-likelihood}
#' \item{convergence}{0 if the algorithm converged and 1 if maxiter was reached}
#' }
estimate_zinb <- function(Y, vecClusterAssign, maxiter=10, verbose=FALSE) {
  
  if(!all( .isWholeNumber(Y) )){
    stop("Expression matrix contains non-integer values.")
  }
  boolOnePiPerCell <- TRUE
  boolPiByCell <- TRUE
  
  n <- ncol(Y)
  J <- nrow(Y)
  
  # 0. Initial estimate of mu, pi and z
  # Add pseudocounts
  print("Initialise NB parameters")
  vecClusters <- unique(vecClusterAssign)
  vecindClusterAssing <- match(vecClusterAssign, vecClusters)
  if(length(vecClusters)>1){
    mu0 <- lapply(vecClusters ,function(cluster){rowMeans(Y[,vecClusterAssign==cluster]+1)})
    mu0 <- lapply(vecClusters ,function(cluster){apply(Y[,vecClusterAssign==cluster],1,
      function(gene){mean(gene[gene>0])}
    )})
    mu0 <- do.call(cbind, mu0)
  } else {
    mu0 <- rowMeans(Y)
    mu0 <- matrix(mu0,nrow=length(mu0),ncol=1)
  }
  muhat0 <- mu0[,vecindClusterAssing]
  thetahat <- array(1,c(J))
  matthetahat <- matrix(thetahat, nrow=J, ncol=n, byrow=FALSE)
  
  if(dim(mu0)[2]>1){
    mu0_nonintersection <- log(mu0[,2:dim(mu0)[2]]) - log(mu0[,1])
    coefs_mu <- cbind( log(mu0[,1]), mu0_nonintersection )
  } else {
    coefs_mu <- log(mu0[,1])
  }
  
  if(boolOnePiPerCell){
    
    # for MLE
    matboolNotZeroObserved <- !is.na(Y) & (Y>0)
    matboolZero <- Y==0
    vecParamGuess <- apply(Y==0, 2, function(cell){sum(cell)/length(cell)})
    
    print("Fitting dropout rates")
    fit_pi <- bplapply(seq(1,n), function(i) {
      # Find MLE through numerical optimisation
      fit <- unlist( optim(
        par=vecParamGuess[i],
        fn=evalLogLikHurdleDrop_comp,
        vecY=Y[,i],
        vecMeanEst=muhat0[,i],
        vecDispEst=thetahat,
        vecboolNotZeroObserved=matboolNotZeroObserved[,i], 
        vecboolZero=matboolZero[,i],
        method="Brent",
        lower=0,upper=1,
        control=list(maxit=1000,fnscale=-1)
      )[c("par","value","convergence")] )
      if(fit[[3]]){
        print(paste0("Not converged MLE pi: ",i))
      }
      return(list(fitted=fit[[1]]))
    })
    pihat <- matrix( sapply(fit_pi, function(x) x$fitted), nrow=J, ncol=n, byrow=TRUE)
  } else {
    pi0 <- sapply(seq_len(n), function(i) {
      y <- as.numeric(Y[,i]<=0)
      fit <- suppressWarnings( glm(y~log(mu0[,vecindClusterAssing[i]]), family = binomial, control = list(maxit = 100)) )
      if(!fit$converged){print("Not converged in glm pi0")}
      return(fitted.values(fit))
    })
    
    #coefs_mu <- log(mu0)
    fit_pi <- bplapply(seq_len(n), function(i) {
      fit <- suppressWarnings( glm(zhat[,i]~muhat0[,i], family = binomial(link = logit), control = list(maxit = 100)) )
      if(!fit$converged){print("Not converged in glm fit_pi")}
      return(list(coefs=coefficients(fit), fitted=fitted.values(fit)))
    })
    coefs_pi <- sapply(fit_pi, function(x) x$coefs)
    
    pihat <- lapply(fit_pi, function(x) x$fitted)
    pihat <- do.call(cbind, pihat)
  }
  zhat <- pihat/(pihat + (1 - pihat) * dnbinom(0, size = 1, mu = muhat0))
  zhat[Y>0] <- 0
  
  #X <- array(0,c(dim(coefs_mu)[2],n))
  #X[1,] <- 1
  #if(dim(coefs_mu)[2] > 1){
  #  for(i in seq(1,dim(coefs_mu)[2])){ X[i,vecindClusterAssing==i] <- 1 }
  #} 
  #W <- array(1,c(J,dim(coefs_pi)[1]))
  
  linkobj <- binomial()
  
  #ll_new_ref <- loglik_small( beta=coefs_mu, alpha=coefs_pi, theta=thetahat, Y, Y>0, X, W, 0, 0, linkobj)
  loglik0 <- sum(log( pihat[Y==0] + (1 - pihat[Y==0]) * dnbinom(0, mu=muhat0[Y==0], size=matthetahat[Y==0], log=FALSE) ))
  loglik1 <- sum( log(1 - pihat[Y>0]) + dnbinom(Y[Y>0], mu=muhat0[Y>0], size=matthetahat[Y>0], log=TRUE) )
  ll_new <- loglik0 + loglik1
  ll_old <- ll_new
  
  if(verbose) {
    print(ll_new)
  }
  
  ## EM iteration
  iter <- 0
  while (abs((ll_old - ll_new)/ll_old) > 1e-4 & iter<maxiter | iter==0) {
    print(paste0("Iteration ", iter+1))
    ll_old <- ll_new
    print("Fitting means")
    fit_mu <- bplapply(seq_len(J), function(i) {
      fit <- glm.nb(Y[i,] ~ vecClusterAssign, weights = (1 - zhat[i,]), init.theta = thetahat[i], start=coefs_mu[i,])
      if(!fit$converged){
        print(paste0("Not converged in glm mu: ",i))
      }
      return(list(fitted=fit$fitted.values, coefs=fit$coefficients, theta=fit$theta))
      #return(list(coefs=coefficients(fit), theta=fit$theta))
    })
    
    coefs_mu <- lapply(fit_mu, function(x) x$coefs)
    coefs_mu <- do.call(rbind, coefs_mu)
    
    #save(coefs_mu,file=file.path(getwd(),"PseudoDE_coefs_mu.RData"))
    #save(fit_mu,file=file.path(getwd(),"PseudoDE_fit_mu.RData"))
    #save(Y,file=file.path(getwd(),"PseudoDE_Y.RData"))
    #save(zhat,file=file.path(getwd(),"PseudoDE_zhat.RData"))
    #save(thetahat,file=file.path(getwd(),"PseudoDE_thetahat.RData"))
    
    #muhat <- exp(coefs_mu)
    muhat <- lapply(fit_mu, function(x) x$fitted)
    muhat <- do.call(rbind, muhat)
    
    thetahat <- unlist(lapply(fit_mu, function(x) x$theta))
    matthetahat <- matrix(thetahat, nrow=J, ncol=n, byrow=FALSE)
    
    # for MLE
    vecParamGuess <- pihat[1,]
    
    # Only estimate on highly expressed genes
    vecboolHighMean <- rowMeans(Y)>=20 
    
    print("Fitting dropout rates")
    if(!boolPiByCell){
      # Find MLE through numerical optimisation
      fit_pi <- unlist( optim(
        par=vecParamGuess,
        fn=evalLogLikHurdleDrop_comp,
        vecY=Y[vecboolHighMean,],
        vecMeanEst=muhat[vecboolHighMean,],
        vecDispEst=thetahat[vecboolHighMean],
        vecboolNotZeroObserved=matboolNotZeroObserved[vecboolHighMean,], 
        vecboolZero=matboolZero[vecboolHighMean,],
        method="Brent",
        lower=0,upper=1,
        control=list(maxit=1000,fnscale=-1)
      )[c("par","value","convergence")] )
      if(fit[[3]]){
        print(paste0("Not converged MLE pi: ",i))
      }
    } else {
      fit_pi <- bplapply(seq(1,n), function(i) {
        #fit <- suppressWarnings(glm(zhat[,i]~log(muhat[,vecindClusterAssing[i]]), family = binomial(link = logit), start=coefs_pi[,i]))
        if(boolOnePiPerCell){
          # Find MLE through numerical optimisation
          vecParamGuessesLocal <- list(0,vecParamGuess[i],1)
          fits <- lapply(vecParamGuessesLocal, function(x){unlist( 
            optim(
              par=x,
              fn=evalLogLikHurdleDrop_comp,
              vecY=Y[vecboolHighMean,i],
              vecMeanEst=muhat[vecboolHighMean,i],
              vecDispEst=thetahat[vecboolHighMean],
              vecboolNotZeroObserved=matboolNotZeroObserved[vecboolHighMean,i], 
              vecboolZero=matboolZero[vecboolHighMean,i],
              method="Brent",
              lower=0,upper=1,
              control=list(maxit=1000,fnscale=-1)
            )[c("par","value","convergence")] ) 
          })
          vecFitValues <- unlist(lapply(fits, function(x) x["value"]))
          fit <- fits[[match(min(vecFitValues),vecFitValues)]]
          if(fit[[3]]){
            print(paste0("Not converged MLE pi: ",i))
          }
          return(list(fitted=fit[[1]]))
        } else {
          fit <- suppressWarnings( glm(zhat[,i]~log(muhat[,i]), family = binomial(link = logit), start=coefs_pi[,i],
            control = list(maxit = 100)) )
          if(!fit$converged){
            print(paste0("Not converged in glm pi: ",i))
            fit <- suppressWarnings( glm(zhat[,i]~log(muhat[,i]), family = binomial(link = logit), start=coefs_pi[,i],
              control = list(maxit = 10000)) )
            if(!fit$converged){ print("Still Not converged in glm pi") }
          }
          return(list(coefs=coefficients(fit), fitted=fitted.values(fit)))
        }
      })
      #save(fit_pi,file=file.path(getwd(),"PseudoDE_fit_pi.RData"))
      if(boolOnePiPerCell){
        pihat <- matrix( sapply(fit_pi, function(x) x$fitted), nrow=J, ncol=n, byrow=TRUE)
      } else {
        coefs_pi <- sapply(fit_pi, function(x) x$par)
        pihat <- lapply(fit_pi, function(x) x$fitted)
        pihat <- do.call(cbind, pihat)
      }
    }
    
    #zhat <- pihat/(pihat + (1 - pihat) * dnbinom(0, size = matrix(thetahat, nrow=J, ncol=n), mu = muhat[,vecindClusterAssing]))
    zhat <- pihat/(pihat + (1 - pihat) * dnbinom(0, size = matthetahat, mu = muhat))
    zhat[Y>0] <- 0
    
    print("Compute Loglik")
    #W <- model.matrix(~log(coefs_mu))
    #W <- model.matrix(~coefs_mu)
    #ll_new_ref <- loglik_small(c(coefs_mu, coefs_pi[1,], coefs_pi[2,], log(thetahat)), Y, Y>0, X, W, J, n*2, 0, 0, linkobj)
    loglik0 <- sum(log( pihat[Y==0] + (1 - pihat[Y==0]) * dnbinom(Y[Y==0], mu=muhat[Y==0], size=matthetahat[Y==0], log=FALSE) ))
    loglik1 <- sum(log( (1 - pihat[Y>0]) * dnbinom(Y[Y>0], mu=muhat[Y>0], size=matthetahat[Y>0], log=FALSE) ))
    ll_new <- loglik0 + loglik1
    
    if(verbose) {
      if(FALSE){
        print(cbind( zeros=unlist(sapply(seq(1,J),function(i){
          sum(log( pihat[i,Y[i,]==0] + (1 - pihat[i,Y[i,]==0]) * dnbinom(Y[i,Y[i,]==0], mu=muhat[i,Y[i,]==0], size=matthetahat[i,Y[i,]==0], log=FALSE) ))
        })),
          nonzeros=unlist(sapply(seq(1,J),function(i){
            sum(log( (1 - pihat[i,Y[i,]>0]) * dnbinom(Y[i,Y[i,]>0], mu=muhat[i,Y[i,]>0], size=matthetahat[i,Y[i,]>0], log=FALSE) ))
          })),
          mean=rowMeans(Y)
        ))
      }
      print(ll_new)
    }
    
    iter <- iter + 1
  }
  
  convergence <- 0
  if(iter == maxiter) {
    convergence <- 1
  }
  
  #p_z <- 1 - (Y == 0) * pihat / (pihat + (1 - pihat) * (1 + muhat[,vecindClusterAssing] / matthetahat[,vecindClusterAssing])^(-matthetahat[,vecindClusterAssing]))
  #eval <- (1 - pihat) * muhat[,vecindClusterAssing]
  #variance <- (1 - pihat) * muhat[,vecindClusterAssing] * (1 + muhat[,vecindClusterAssing] * (1/matthetahat[,vecindClusterAssing] + pihat))
  p_z <- 1 - (Y == 0) * pihat / (pihat + (1 - pihat) * (1 + muhat / matthetahat)^(-matthetahat))
  eval <- (1 - pihat) * muhat
  variance <- (1 - pihat) * muhat * (1 + muhat * (1/matthetahat + pihat))
  
  return(list(mu=muhat, theta=thetahat, pi=pihat,
    p_z=p_z, expected_value=eval, variance=variance,
    loglik=ll_new, convergence=convergence))
}

#' Log-likelihood function of the zero-inflated negative binomial model
#'
#' This function computes the log-likelihood of a standard regression model
#'
#' This is a (hopefully) memory-efficient implementation of the log-likelihood of a
#' zero-inflated negative binomial regression model.
#' In this attempt, the design matrices don't have n*J rows, but n and J, respectively.
#' The computation is a bit slower, but the memory usage should be much smaller for
#' large J and n.
#'
#' @param parms a vector of parameters: should contain the values of beta, followed by those of alpha, followed by the log(1/phi)
#' @param Y the data matrix (genes in rows, cells in columns)
#' @param Y1 a logical indicator of Y>0
#' @param X the design matrix for the regression on mu (n x k_X)
#' @param W the design matrix for the regression on pi (J x k_W)
#' @param kx the number of beta parameters
#' @param kw the number of alpha parameters
#' @param offsetx the offset for the regression on X
#' @param offsetw the offset for the regression on W
#' @param linkobj the link function object for the regression on pi (typically the result of binomial())
loglik_small <- function(beta, alpha, theta, Y, Y1, X, W, offsetx, offsetw, linkobj) {
  
  J <- nrow(Y)
  n <- ncol(Y)
  
  mu <- exp(beta %*% X + offsetx)
  
  eta <- W %*% alpha + offsetw
  pi <- logistic(eta)
  
  theta <- matrix(exp(theta), nrow=J, ncol=n, byrow=FALSE)
  
  loglik0 <- log(pi + exp(log(1 - pi) + suppressWarnings(dnbinom(0, size = theta, mu = mu, log = TRUE))))
  loglik1 <- log(1 - pi) + suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))
  
  return(sum(loglik0[which(Y==0)]) + sum(loglik1[which(Y>0)]))
}

#' Imputation of zero counts based on zero-inflated negative binomial model
#'
#' This function is used to impute the data, by weighting the observed zeros by their
#' probability of coming from the zero-inflation part of the distribution P(Z=1|Y).
#'
#' @details The imputation is carried out with the following formula:
#' y_{ij}* =  y_{ij} * Pr(Z_{ij} = 0 | Y_{ij}=y_{ij}) + mu_{ij} * Pr(Z_{ij} = 1 | Y_{ij} = y_{ij}).
#' Note that for y_{ij} > 0, Pr(Z_{ij} = 0 | Y_{ij}=y_{ij}) = 1 and hence the data are not imputed.
#'
#' @export
#'
#' @param expression the data matrix (genes in rows, cells in columns)
#'
#' @return the imputed expression matrix.
impute_zinb <- function(expression) {
  pars <- estimate_zinb(expression)
  w <- pars$p_z
  imputed <- expression * w + pars$mu * (1 - w)
  return(imputed)
}

logit <- binomial()$linkfun
logistic <- binomial()$linkinv

.isWholeNumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}
