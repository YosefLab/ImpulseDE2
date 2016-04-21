#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++     DE analysis   ++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Compute p-value of differential expression and report summary of statistics.

# INPUT:
#   arr3DCountData: (Numeric 3D array genes x samples x replicates)
#       Contains expression values or similar locus-specific read-outs.
#   dfAnnotationRed: (Table samples x 2[time and condition]) 
#       Co-variables for the samples including condition and time points.
#       Time points vecMust be numeric.
#   imp_fit_results: (list runs*2 ["impulse_parameters","impulse_fits"])
#       impulse_parameters: (matrix genes x (NPARAM+1)*runs) +1 is value of 
#         objective of optimisation (i.e. sum of squares or weighted SS)
#       impulse_fits: (matrix genes x timepoints) model values for gene data
#   control_name: (str) name of the control condition in annotation_table.
# OUTPUT:
#   dfDEAnalysis: (data frame genes x fitting characteristics) Summary of fitting
#       procedure for each gene.

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

computePval <- function(arr3DCountData,dfAnnotationRed,vecDispersions,
  lsImpulseFits,
  strCaseName=NULL, strControlName=NULL, NPARAM=6){
  
  if(is.null(strControlName)){
    # Without control data:
    # Take Values from the only fitting performed: case
    # The full model is the impulse fit.
    vecLogLikFull <- lsImpulseFits$parameters_case[,"logL_H1"]
    # The reduced model is the mean fit.
    vecLogLikRed <- lsImpulseFits$parameters_case[,"logL_H0"]
    # Mean inferred expression:
    vecMu <- lsImpulseFits$parameters_case[,"vecMu"]
    # Parameters and 1 dispersion estimate
    scaDegFreedomFull <- NPARAM + 1
    # 1 mean and 1 dispersion estimate
    scaDegFreedomRed <- 1 + 1
  } else {
    # With control data:
    # Full model: Case and control model separate:
    # The log likelihood of the full model is the sum of the
    # log likelihoods of case and control fits.
    vecLogLikFull <- lsImpulseFits$parameters_case[,"logL_H1"] + lsImpulseFits$parameters_control[,"logL_H1"]
    # The reduced model is the combined data fit.
    vecLogLikRed <- lsImpulseFits$parameters_combined[,"logL_H1"]
    # Mean inferred expression: On combined data
    vecMu <- lsImpulseFits$parameters_combined[,"vecMu"]
    # Parameters of both models (case and control) and 1 dispersion estimate
    scaDegFreedomFull <- NPARAM*2 + 1
    # Parameters of one model (combined) and 1 dispersion estimate
    scaDegFreedomRed <- NPARAM + 1
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
      "loglik_full"=round(vecLogLikFull),
      "loglik_red"=round(vecLogLikRed),
      "deviance"=round(vecDeviance),
      "vecMu"=round(vecMu),
      "size"=round(vecDispersions),
      "converge_impulse"=impulse_fit_dfDEAnalysis$parameters_case[,"converge_H1"],
      "converge_mean"=impulse_fit_dfDEAnalysis$parameters_case[,"converge_H0"],
      stringsAsFactors = FALSE))
  } else {
    # With control data:
    dfDEAnalysis =   as.data.frame(cbind(
      "Gene" = row.names(arr3DCountData),
      "p"=as.numeric(vecPvalue),
      "adj.p"=as.numeric(vecPvalueBH),
      "loglik_full"=round(vecLogLikFull),
      "loglik_red"=round(vecLogLikRed),
      "deviance"=round(vecDeviance),
      "vecMu"=round(vecMu),
      "size"=round(vecDispersions),
      "converge_combined"=impulse_fit_dfDEAnalysis$impulse_parameters_combined[,"converge_H1"],
      "converge_case"=impulse_fit_dfDEAnalysis$impulse_parameters_case[,"converge_H1"],
      "converge_control"=impulse_fit_dfDEAnalysis$impulse_parameters_control[,"converge_H1"],
      "converge_mean"=impulse_fit_dfDEAnalysis$impulse_parameters_combined[,"converge_H0"],
      stringsAsFactors = FALSE))
  }
  
  # Order data frame by adjusted p-value
  dfDEAnalysis$adj.p <- as.numeric(as.character(dfDEAnalysis$adj.p))
  dfDEAnalysis = dfDEAnalysis[order(dfDEAnalysis$adj.p),]
  
  return(dfDEAnalysis)
}