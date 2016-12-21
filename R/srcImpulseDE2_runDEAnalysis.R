#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++   Perform differential expression analysis and identify transient genes  ++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Perform differential expression analysis and identification
#' of transiently activated or deactivated genes.
#' 
#' Performs model selction based on loglikelihood ratio tests.
#' The primary model selection is the differential expression analysis.
#' The secondary model selection is the selection between a sigmoidal
#' and an impulse fit for differentially expressed genes which is used
#' to define transiently activated or deactivated genes.
#' 
#' @seealso Called by \code{runImpulseDE2}.
#' 
#' @param matCountDataProc: (matrix genes x samples)
#'    Read count data.
#' @param dfAnnotationProc: (data frame samples x covariates) 
#'    {Sample, Condition, Time (numeric), TimeCateg (str)
#'    (and confounding variables if given).}
#'    Processed annotation table with covariates for each sample.
#' @parm lsModelFits: (list length number of conditions fit (1 or 3))
#'    {"case"} or {"case", "control", "combined"}
#'    This is the lsModelFits object handed to this function with additional
#'    sigmoid fit entries for every gene for the given condition.
#'    One model fitting object for each condition:
#'    In case-only DE analysis, only the condition {"case"} is fit.
#'    In case-control DE analysis, the conditions 
#'    {"case", "control","combined} are fit.
#'    Each condition entry is a list of model fits for each gene.
#'    Each gene entry is a list of model fits to the individual models:
#'    Impulse model, constant model and sigmoidal fit.
#'    Each model fit per gene is a list of fitting parameters and results.
#'    \itemize{
#'      \item Condition ID: (list length number of genes)
#'      List of fits for each gene to the samples of this condition.
#'      One entry of this format for all conditions fit.
#'      \itemize{
#'        \item Gene ID: (list length 2)
#'        Impulse, constant and sigmoidal model fit to gene observations.
#'        One entry of this format for all gene IDs.
#'        \itemize{
#'          \item lsImpulseFit: (list) List of impulse fit parameters and results.
#'          \itemize{
#'            \item vecImpulseParam: (numeric vector length 6)
#'            {beta, h0, h1, h2, t1, t2}
#'            Maximum likelihood estimators of impulse model parameters.
#'            \item vecImpulseValue: (numeric vector length number of time points)
#'            Values of impulse model fit at time points used for fit.
#'            \item lsvecBatchFactors: (list length number of confounders)
#'            List of vectors of scalar batch correction factors for each sample.
#'            These are also maximum likelihood estimators.
#'            NULL if no confounders given.
#'            \item scaDispParam: (scalar) Dispersion parameter estimate
#'            used in fitting (hyper-parameter).
#'            \item scaLL: (scalar) Loglikelihood of data under maximum likelihood
#'            estimator model.
#'            \item scaConvergence: (scalar) 
#'            Convergence status of optim on impulse model.
#'          }
#'          \item lsConstFit: (list) List of constant fit parameters and results.
#'          \itemize{
#'            \item scaMu: (scalar) Maximum likelihood estimator of
#'            negative binomial mean parameter.
#'            \item lsvecBatchFactors: (list length number of confounders)
#'            List of vectors of scalar batch correction factors for each sample.
#'            These are also maximum likelihood estimators.
#'            NULL if no confounders given.
#'            \item scaDispParam: (scalar) Dispersion parameter estimate
#'            used in fitting (hyper-parameter).
#'            \item scaLL: (scalar) Loglikelihood of data under maximum likelihood
#'            estimator model.
#'            \item scaConvergence: (scalar) 
#'            Convergence status of optim on constant model.
#'          }
#'          \item ls SigmoidFit: (list) List of sigmoidal fit parameters and results.
#'          \itemize{
#'            \item vecSigmoidParam: (numeric vector length 4)
#'            {beta, h0, h1, t}
#'            Maximum likelihood estimators of sigmoidal model parameters.
#'            \item vecSigmoidValue: (numeric vector length number of time points)
#'            Values of sigmoid model fit at time points used for fit.
#'            \item lsvecBatchFactors: (list length number of confounders)
#'            List of vectors of scalar batch correction factors for each sample.
#'            These are also maximum likelihood estimators.
#'            NULL if no confounders given.
#'            \item scaDispParam: (scalar) Dispersion parameter estimate
#'            used in fitting (hyper-parameter).
#'            \item scaLL: (scalar) Loglikelihood of data under maximum likelihood
#'            estimator model.
#'            \item scaConvergence: (scalar) 
#'            Convergence status of optim on sigmoidal model.
#'        }
#'      }
#'    }
#' @param boolCaseCtrl: (bool) 
#' 		Whether to perform case-control analysis. Does case-only
#' 		analysis if FALSE.
#' @param vecConfounders: (vector of strings number of confounding variables)
#' 		Factors to correct for during batch correction.
#' 		Names refer to columns in dfAnnotation.
#' @param boolIdentifyTransients: (bool) [Defaul FALSE]
#'    Whether to identify transiently activated or deactivated 
#'    genes. This involves an additional fitting of sigmoidal models
#'    and hypothesis testing between constant, sigmoidal and impulse model.
#' 
#' @return dfDEAnalysis (data frame genes x reported characteristics) 
#'    Summary of fitting procedure and 
#'    differential expression results for each gene.
#'    \itemize{
#'      \item Gene: Gene ID.
#'      \item p: P-value for differential expression.
#'      \item padj: Benjamini-Hochberg false-discovery rate corrected p-value
#'      for differential expression analysis.
#'      \item loglik_full: Loglikelihood of full model.
#'      \item loglik_red: Loglikelihood of reduced model.
#'      \item df_full: Degrees of freedom of full model.
#'      \item df_red: Degrees of freedom of reduced model
#'      \item mean: Inferred mean parameter of constant model over all samples.
#'    }
#'    Entries only present in case-only DE analysis:
#'    \itemize{
#'      \item converge_impulse: Convergence status of optim for 
#'      impulse model fit (full model).
#'      \item converge_const: Convergence status of optim for 
#'      constant model fit (reduced model).
#'    }
#'    Entries only present in case-control DE analysis:
#'    \itemize{
#'      \item converge_combined: Convergence status of optim for 
#'      impulse model fit to case and control samples combined (reduced model).
#'      \item converge_case: Convergence status of optim for 
#'      impulse model fit to samples of case condition (full model 1/2).
#'      \item converge_control: Convergence status of optim for 
#'      impulse model fit to samples of control condition (full model 2/2).
#'    }
#'    Entries only present if boolIdentifyTransients is TRUE:
#'    \itemize{
#'      \item impulseTOsigmoid_p: P-value of loglikelihood ratio test
#'      impulse model fit versus sigmoidal model on samples of case condition.
#'      \item dfDEAnalysis$impulseTOsigmoid_padj: Benjamini-Hochberg 
#'      false-discovery rate corrected p-value of loglikelihood ratio test
#'      impulse model fit versus sigmoid model on samples of case condition.
#'      \item dfDEAnalysis$sigmoidTOconst_p: P-value of loglikelihood ratio test
#'      sigmoidal model fit versus constant model on samples of case condition.
#'      \item dfDEAnalysis$sigmoidTOconst_padj: Benjamini-Hochberg 
#'      false-discovery rate corrected p-value of loglikelihood ratio test
#'      sigmoidal model fit versus constant model on samples of case condition.
#'      \item dfDEAnalysis$isTransient: (bool) Whether gene is transiently
#'      activated or deactivated and differentially expressed.
#'      \item dfDEAnalysis$isMonotonous: (bool) Whether gene is not transiently
#'      activated or deactivated and differentially expressed. This scenario
#'      corresponds to a montonous expression level increase or decrease.
#'    }
#'    
#' @author David Sebastian Fischer
#' 
#' @export
runDEAnalysis <- function(matCountDataProc,
                        dfAnnotationProc, 
                        lsModelFits,
												boolCaseCtrl,
												vecConfounders,
												boolIdentifyTransients){
  
  if(!boolCaseCtrl){
    # Case-only:
    # Take Values from the only fitting performed: case
    # The full model is the impulse fit.
    vecLogLikFull <- sapply(lsModelFits$case, function(fit) fit$lsImpulseFit$scaLL )
    # The reduced model is the mean fit.
    vecLogLikRed <- sapply(lsModelFits$case, function(fit) fit$lsConstFit$scaLL )
    # Mean inferred expression:
    vecMu <- sapply(lsModelFits$case, function(fit) fit$lsConstFit$scaMu )
    if(!is.null(vecConfounders)){
    	scaNBatchFactors <- sum(sapply(vecConfounders, function(confounder){ 
    		length(unique(dfAnnotationProc[,confounder]))-1 
    	}))
    } else { scaNBatchFactors <- 0 }
    # 6 impulse model parameters, 1 dispersion estimate and 
    # 1 batch factor for each batch (except for the first one) for each confounder.
    scaDegFreedomFull <- 6 + 1 + scaNBatchFactors
    # 1 constant model parameter1, 1 dispersion estimate and 
    # 1 batch factor for each batch (except for the first one) for each confounder.
    scaDegFreedomRed <- 1 + 1 + scaNBatchFactors
    
    vecConvergenceImpulse <- sapply(lsModelFits$case, function(fit) fit$lsImpulseFit$scaConvergence )
    vecConvergenceConst <- sapply(lsModelFits$case, function(fit) fit$lsConstFit$scaConvergence )
  } else {
    # Case-control:
    # Full model: Case and control model separate:
    # The log likelihood of the full model is the sum of the
    # log likelihoods of case and control fits.
    vecLogLikFull <- sapply(lsModelFits$case, function(fit) fit$lsImpulseFit$scaLL )+
      sapply(lsModelFits$control, function(fit) fit$lsImpulseFit$scaLL )
    # The reduced model is the combined data fit.
    vecLogLikRed <- sapply(lsModelFits$combined, function(fit) fit$lsImpulseFit$scaLL )
    # Mean inferred expression: On combined data
    vecMu <- sapply(lsModelFits$combined, function(fit) fit$lsConstFit$scaMu )
    if(!is.null(vecConfounders)){
    	scaNBatchFactorsFull <- sum(sapply(vecConfounders, function(confounder){ 
    		length(unique(dfAnnotationProc[dfAnnotationProc$Condition=="case",confounder]))-1+
    			length(unique(dfAnnotationProc[dfAnnotationProc$Condition=="control",confounder]))-1
    	}))
    	scaNBatchFactorsRed <- sum(sapply(vecConfounders, function(confounder){ 
    		length(unique(dfAnnotationProc[,confounder]))-1
    	}))
    } else { 
    	scaNBatchFactorsFull <- 0
    	scaNBatchFactorsRed <- 0
    }
    # 6 impulse model parameters for each case and control, 
    # 1 dispersion estimate and 
    # 1 batch factor for each batch (except for the first one) for each confounder in each condition.
    scaDegFreedomFull <- 6*2 + 1 + scaNBatchFactorsFull
    # 6 impulse model parameters, 1 dispersion estimate and 
    # 1 batch factor for each batch (except for the first one) for each confounder.
    scaDegFreedomRed <- 6 + 1 + scaNBatchFactorsRed
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
  
  scaNGenes <- dim(matCountDataProc)[1]
  if(!boolCaseCtrl){
    # Case-only
    dfDEAnalysis =   data.frame(
      Gene = row.names(matCountDataProc),
      p=vecPvalue,
      padj=vecPvalueBH,
      loglik_full=vecLogLikFull,
      loglik_red=vecLogLikRed,
      df_full=rep(scaDegFreedomFull, scaNGenes),
      df_red=rep(scaDegFreedomRed, scaNGenes),
      mean=vecMu,
      converge_impulse=vecConvergenceImpulse,
      converge_const=vecConvergenceConst,
      stringsAsFactors = FALSE)
  } else {
    # Case-control
    dfDEAnalysis =   data.frame(
      Gene = row.names(matCountDataProc),
      p=as.numeric(vecPvalue),
      padj=as.numeric(vecPvalueBH),
      loglik_full=vecLogLikFull,
      loglik_red=vecLogLikRed,
      df_full=rep(scaDegFreedomFull, scaNGenes),
      df_red=rep(scaDegFreedomRed, scaNGenes),
      mean=vecMu,
      converge_combined=vecConvergenceImpulseCombined,
      converge_case=vecConvergenceImpulseCase,
      converge_control=vecConvergenceImpulseControl,
      stringsAsFactors = FALSE)
  }
  
  # Add entries if sigmoidal fit was performed and transient tranjectories
  # are to be differentiated from monotonous ones. Do this in case condition.
  if(boolIdentifyTransients){
    # Take the sigmoid model as reference if transient are tested versus montonous model
    vecLogLikImpulse <- sapply(lsModelFits$case, function(fit) fit$lsImpulseFit$scaLL )
    vecLogLikSigmoid <- sapply(lsModelFits$case, function(fit) fit$lsSigmoidFit$scaLL )
    vecLogLikConst <- sapply(lsModelFits$case, function(fit) fit$lsConstFit$scaLL )
    scaDegFreedomImpulse <- 6 + 1 + scaNBatchFactors
    scaDegFreedomSigmoid <- 4 + 1 + scaNBatchFactors
    scaDegFreedomConst <- 1 + 1 + scaNBatchFactors
    
    # Compare impulse to sigmoid: Trajectory is transient if this is significant
    vecDevianceImpulseSigmoid <- 2*(vecLogLikImpulse - vecLogLikSigmoid)
    vecPvalueImpulseSigmoid <- pchisq(vecDevianceImpulseSigmoid,
                        df=scaDegFreedomImpulse-scaDegFreedomSigmoid,
                        lower.tail=FALSE)
    vecPvalueImpulseSigmoidBH = p.adjust(vecPvalueImpulseSigmoid, method = "BH")
    
    # Compare sigmoid to constant: Trajectory is constant if this is not significant
    # Trajectory is montonous if this is significant but vecPvalueImpulseSigmoidBH is not
    vecDevianceSigmoidConst <- 2*(vecLogLikSigmoid - vecLogLikConst)
    vecPvalueSigmoidConst <- pchisq(vecDevianceSigmoidConst,
                                      df=scaDegFreedomSigmoid-scaDegFreedomConst,
                                      lower.tail=FALSE)
    vecPvalueSigmoidConstBH = p.adjust(vecPvalueSigmoidConst, method = "BH")
    
    # Add entries into DE table
    dfDEAnalysis$impulseTOsigmoid_p <- as.numeric(vecPvalueImpulseSigmoid)
    dfDEAnalysis$impulseTOsigmoid_padj <- as.numeric(vecPvalueImpulseSigmoidBH)
    dfDEAnalysis$sigmoidTOconst_p <- as.numeric(vecPvalueSigmoidConst)
    dfDEAnalysis$sigmoidTOconst_padj <- as.numeric(vecPvalueSigmoidConstBH)
    dfDEAnalysis$isTransient <- vecPvalueImpulseSigmoidBH <= 0.01
    dfDEAnalysis$isMonotonous <- (vecPvalueSigmoidConst <= 0.01) &
      (vecPvalueImpulseSigmoidBH > 0.01)
  }
  
  return(dfDEAnalysis)
}