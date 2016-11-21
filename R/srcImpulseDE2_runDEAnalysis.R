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
#' @param matCountDataProc: (matrix genes x samples)
#'    Count data: Reduced version of \code{matCountData}. 
#'    For internal use.
#' @param dfAnnotationProc: (Table) Processed annotation table. 
#'    Lists co-variables of samples: 
#'    Sample, Condition, Time (numeric), TimeCateg (categorial)
#'    (and LongitudinalSeries). For internal use.
#' @param vecDispersions (vector number of genes) Inverse of gene-wise 
#'    negative binomial dispersion coefficients computed by DESeq2.
#' @param lsImpulseFits (list length 2 or 6) List of matrices which
#'    contain parameter fits and model values for given time course for the
#'    case condition (and control and combined if control is present).
#'    Each parameter matrix is called parameter_'condition' and has the form
#'    (genes x [beta, h0, h1, h2, t1, t2, logL_H1, converge_H1, mu, logL_H0, 
#'    ) where beta to t2 are parameters of the impulse
#'    model, mu is the single parameter of the mean model, logL are
#'    log likelihoods of full (H1) and reduced model (H0) respectively, converge
#'    is convergence status of numerical optimisation of model fitting by
#'    \code{optim} from \code{stats}. Each value matrix is called
#'    value_'condition' and has the form (genes x time points) and contains the
#'    counts predicted by the impulse model at the observed time points.
#' @param dfDEAnalysis (data frame genes x fitting characteristics) 
#' @param strCaseName (str) Name of the case condition in \code{dfAnnotationRedFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#'    \code{dfAnnotationRedFull}.
#' @param strMode: (str) [Default "batch"] {"batch","longitudinal","singlecell"}
#'    Mode of model fitting.
#' @param NPARAM (scalar) [Default 6] Number of parameters of impulse model.
#' 
#' @return dfDEAnalysis (data frame genes x fitting characteristics) 
#'    Summary of fitting procedure for each gene.
#' @export

runDEAnalysis <- function(matCountDataProc,
                        vecDispersions,
                        dfAnnotationProc, 
                        lsModelFits,
                        strCaseName, 
                        strControlName, 
                        strMode){
  
  if(is.null(strControlName)){
    # Without control data:
    # Take Values from the only fitting performed: case
    # The full model is the impulse fit.
    vecLogLikFull <- sapply(lsModelFits$case, function(fit) fit$lsImpulseFit$scaLL )
    # The reduced model is the mean fit.
    vecLogLikRed <- sapply(lsModelFits$case, function(fit) fit$lsConstFit$scaLL )
    # Mean inferred expression:
    vecMu <- sapply(lsModelFits$case, function(fit) fit$lsConstFit$scaMu )
    if(strMode=="singlebatch"){
      # Impulse model parameters and 1 dispersion estimate
      scaDegFreedomFull <- 6 + 1
      # 1 dispersion estimate and overall mean estimate
      scaDegFreedomRed <- 1 + 1
    } else if(strMode=="batcheffects"){
      # 6 impulse model parameters, 1 dispersion estimate and 
      # 1 scaling factor for each batch (except for the first one).
      scaDegFreedomFull <- 6 + 1 + length(unique(dfAnnotationProc$Batch)) - 1
      # 1 mean parameter, 1 dispersion parmeter and
      # 1 scaling factor for each batch (except for the first one).
      scaDegFreedomRed <- 1 + 1 + length(unique(dfAnnotationProc$Batch)) - 1
    }
    vecConvergenceImpulse <- sapply(lsModelFits$case, function(fit) fit$lsImpulseFit$scaConvergence )
    vecConvergenceConst <- sapply(lsModelFits$case, function(fit) fit$lsConstFit$scaConvergence )
  } else {
    # With control data:
    # Full model: Case and control model separate:
    # The log likelihood of the full model is the sum of the
    # log likelihoods of case and control fits.
    vecLogLikFull <- sapply(lsModelFits$case, function(fit) fit$lsImpulseFit$scaLL )+
      sapply(lsModelFits$control, function(fit) fit$lsImpulseFit$scaLL )
    # The reduced model is the combined data fit.
    vecLogLikRed <- sapply(lsModelFits$combined, function(fit) fit$lsImpulseFit$scaLL )
    # Mean inferred expression: On combined data
    vecMu <- sapply(lsModelFits$combined, function(fit) fit$lsConstFit$scaMu )
    if(strMode=="singlebatch"){
      # Parameters of both models (case and control) and 1 dispersion estimate
      scaDegFreedomFull <- 6*2 + 1
      # Parameters of one model (combined) and 1 dispersion estimate
      scaDegFreedomRed <- 6 + 1
    } else if(strMode=="batcheffects"){
      # 6 impulse model parameters for each model, 1 dispersion estimate and 
      # 1 scaling factor for each batch (except for the first one) for each model.
      scaDegFreedomFull <- 6*2 + 1 + 2*(length(unique(dfAnnotationProc$Batch)) - 1)
      # 1 mean parameter, 1 dispersion parmeter and
      # 1 scaling factor for each batch (except for the first one).
      scaDegFreedomRed <- 6 + 1 + length(unique(dfAnnotationProc$Batch)) - 1
    }
    vecConvergenceImpulseCombined <- sapply(lsModelFits$combined, function(fit) fit$lsImpulseFit$scaConvergence )
    vecConvergenceImpulseCase <- sapply(lsModelFits$case, function(fit) fit$lsImpulseFit$scaConvergence )
    vecConvergenceImpulseControl <- sapply(lsModelFits$control, function(fit) fit$lsImpulseFit$scaConvergence )
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
    dfDEAnalysis =   data.frame(
      Gene = row.names(matCountDataProc),
      p=vecPvalue,
      padj=vecPvalueBH,
      loglik_full=vecLogLikFull,
      loglik_red=vecLogLikRed,
      deviance=vecDeviance,
      mean=vecMu,
      size=vecDispersions,
      converge_impulse=vecConvergenceImpulse,
      converge_const=vecConvergenceConst,
      stringsAsFactors = FALSE)
  } else {
    # With control data:
    dfDEAnalysis =   data.frame(
      Gene = row.names(matCountDataProc),
      p=as.numeric(vecPvalue),
      padj=as.numeric(vecPvalueBH),
      loglik_full=vecLogLikFull,
      loglik_red=vecLogLikRed,
      deviance=vecDeviance,
      mean=vecMu,
      size=vecDispersions,
      converge_combined=vecConvergenceImpulseCombined,
      converge_case=vecConvergenceImpulseCase,
      converge_control=vecConvergenceImpulseControl,
      stringsAsFactors = FALSE)
  }
  
  # Order data frame by adjusted p-value
  dfDEAnalysis$padj <- as.numeric(as.character(dfDEAnalysis$padj))
  dfDEAnalysis = dfDEAnalysis[order(dfDEAnalysis$padj),]
  
  return(dfDEAnalysis)
}