#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++    computePval   +++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute p-values for model fit
#' 
#' Compute p-value of differential expression based on chi-squared distribution
#' of deviance of liklihoods and report summary of statistics.
#' 
#' @seealso Called by \code{runImpulseDE2}.
#' 
#' @param ImpulseDE2_arr3DCountData: (3D array genes x samples x replicates)
#'    Count data: \code{arr2DCountData} reshaped into a 3D array. 
#'    For internal use.
#' @param vecDispersions (vector number of genes) Inverse of gene-wise 
#'    negative binomial dispersion coefficients computed by DESeq2.
#' @param dfAnnotationRed (data frame) Reduced version of 
#'    \code{dfAnnotationFull}. Lists co-variables of samples: 
#'    Sample, Condition, Time. Time must be numeric. For internal use.
#' @param lsImpulseFits (list length 2 or 6) List of matrices which
#'    contain parameter fits and model values for given time course for the
#'    case condition (and control and combined if control is present).
#'    Each parameter matrix is called parameter_'condition' and has the form
#'    (genes x [beta, h0, h1, h2, t1, t2, logL_H1, converge_H1, mu, logL_H0, 
#'    converge_H0]) where beta to t2 are parameters of the impulse
#'    model, mu is the single parameter of the mean model, logL are
#'    log likelihoods of full (H1) and reduced model (H0) respectively, converge
#'    is convergence status of numerical optimisation of model fitting by
#'    \code{optim} from \code{stats} of either model. Each value matrix is called
#'    value_'condition' and has the form (genes x time points) and contains the
#'    counts predicted by the impulse model at the observed time points.
#' @param strCaseName (str) Name of the case condition in \code{dfAnnotationRedFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#'    \code{dfAnnotationRedFull}.
#' @param NPARAM (scalar) [Default 6] Number of parameters of impulse model.
#' 
#' @return dfDEAnalysis (data frame genes x fitting characteristics) 
#'    Summary of fitting procedure for each gene.
#' @export

computePval <- function(arr3DCountData,vecDispersions,
  dfAnnotationRed, lsImpulseFits,
  strCaseName=NULL, strControlName=NULL, strMode="batch",
  NPARAM=6){
  
  if(is.null(strControlName)){
    # Without control data:
    # Take Values from the only fitting performed: case
    # The full model is the impulse fit.
    vecLogLikFull <- lsImpulseFits$parameters_case[,"logL_H1"]
    # The reduced model is the mean fit.
    vecLogLikRed <- lsImpulseFits$parameters_case[,"logL_H0"]
    # Mean inferred expression:
    vecMu <- lsImpulseFits$parameters_case[,"mu"]
    if(strMode=="batch"){
      # Parameters and 1 dispersion estimate
      scaDegFreedomFull <- NPARAM + 1
      # 1 dispersion estimate and overall mean estimate
      scaDegFreedomRed <- 1 + 1
    } else if(strMode=="timecourses"){
      # Parameters, 1 dispersion estimate and 
      # mean estimates for each time course + 1 overall mean estimate
      scaDegFreedomFull <- NPARAM + 1 + dim(arr3DCountData)[3]  + 1
      # 1 dispersion estimate and 1 mean estimate for each time course
      scaDegFreedomRed <- 1 + dim(arr3DCountData)[3] 
    }
  } else {
    # With control data:
    # Full model: Case and control model separate:
    # The log likelihood of the full model is the sum of the
    # log likelihoods of case and control fits.
    vecLogLikFull <- lsImpulseFits$parameters_case[,"logL_H1"] + lsImpulseFits$parameters_control[,"logL_H1"]
    # The reduced model is the combined data fit.
    vecLogLikRed <- lsImpulseFits$parameters_combined[,"logL_H1"]
    # Mean inferred expression: On combined data
    vecMu <- lsImpulseFits$parameters_combined[,"mu"]
    if(strMode=="batch"){
      # Parameters of both models (case and control) and 1 dispersion estimate
      scaDegFreedomFull <- NPARAM*2 + 1
      # Parameters of one model (combined) and 1 dispersion estimate
      scaDegFreedomRed <- NPARAM + 1
    } else if(strMode=="timecourses"){
      # Parameters of both models (case and control), 1 dispersion estimate and 
      # mean estimates for each time course + 1 overall mean estimate
      scaDegFreedomFull <- NPARAM*2 + 1 + dim(arr3DCountData)[3]  + 1
      # Parameters of one model (combined), 1 dispersion estimate 
      # and 1 mean estimate for each time course
      scaDegFreedomRed <- 1 + dim(arr3DCountData)[3] 
    }
  }
  
  # Compute difference in degrees of freedom between null model and alternative model.
  scaDeltaDegFreedom <- scaDegFreedomFull - scaDegFreedomRed
  # Compute test statistic: Deviance
  vecDeviance <- 2*(vecLogLikFull - vecLogLikRed)
  # Get p-values from Chi-square distribution (assumption about null model)
  vecPvalue <- pchisq(vecDeviance,df=scaDeltaDegFreedom,lower.tail=FALSE)
  # Multiple testing correction (Benjamini-Hochberg)
  vecPvalueBH = p.adjust(vecPvalue, method = "BH")
  
  if(is.null(strControlName)){
    # Without control data:
    dfDEAnalysis =   as.data.frame(cbind(
      "Gene" = row.names(arr3DCountData),
      "p"=as.numeric(vecPvalue),
      "adj.p"=as.numeric(vecPvalueBH),
      "loglik_full"=vecLogLikFull,
      "loglik_red"=vecLogLikRed,
      "deviance"=vecDeviance,
      "mean"=vecMu,
      "size"=vecDispersions,
      "converge_impulse"=lsImpulseFits$parameters_case[,"converge_H1"],
      "converge_mean"=lsImpulseFits$parameters_case[,"converge_H0"],
      stringsAsFactors = FALSE))
  } else {
    # With control data:
    dfDEAnalysis =   as.data.frame(cbind(
      "Gene" = row.names(arr3DCountData),
      "p"=as.numeric(vecPvalue),
      "adj.p"=as.numeric(vecPvalueBH),
      "loglik_full"=vecLogLikFull,
      "loglik_red"=vecLogLikRed,
      "deviance"=vecDeviance,
      "mean"=vecMu,
      "size"=vecDispersions,
      "converge_combined"=lsImpulseFits$impulse_parameters_combined[,"converge_H1"],
      "converge_case"=lsImpulseFits$impulse_parameters_case[,"converge_H1"],
      "converge_control"=lsImpulseFits$impulse_parameters_control[,"converge_H1"],
      "converge_mean"=lsImpulseFits$impulse_parameters_combined[,"converge_H0"],
      stringsAsFactors = FALSE))
  }
  
  # Order data frame by adjusted p-value
  dfDEAnalysis$adj.p <- as.numeric(as.character(dfDEAnalysis$adj.p))
  dfDEAnalysis = dfDEAnalysis[order(dfDEAnalysis$adj.p),]
  
  return(dfDEAnalysis)
}