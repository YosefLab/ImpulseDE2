################################################################################
########################     ImpulseDE2 package     ############################
################################################################################

### Version 1.0
### Author David Sebastian Fischer

################################################################################
### Libraries and source code
################################################################################

library(compiler)
library(DESeq2)
library(BiocParallel)
library(Matrix)
library(ComplexHeatmap)
library(circlize)

# Source functions in .R files from same directory as this function.
#setwd("/Users/davidsebastianfischer/gitDevelopment/ImpulseDE2/R")
#setwd("/data/yosef2/users/fischerd/code/ImpulseDE2/R")
setwd("/home/david/gitDevelopment/ImpulseDE2/R")

source("srcImpulseDE2_evalImpulse.R")
source("srcImpulseDE2_evalSigmoid.R")
source("srcImpulseDE2_compareDEMethods.R")
source("srcImpulseDE2_computeNormConst.R")
source("srcImpulseDE2_CostFunctionsFit.R")
source("srcImpulseDE2_fitImpulse.R")
source("srcImpulseDE2_fitSigmoid.R")
source("srcImpulseDE2_plotGenes.R")
source("srcImpulseDE2_plotHeatmap.R")
source("srcImpulseDE2_processData.R")
source("srcImpulseDE2_runDEAnalysis.R")
source("srcImpulseDE2_runDESeq2.R")
source("srcImpulseDE2_simulateDataSet.R")

# Compile functions
evalImpulse_comp <- cmpfun(evalImpulse)
evalSigmoid_comp <- cmpfun(evalSigmoid)
evalLogLikMu_comp <- cmpfun(evalLogLikMu)
evalLogLikImpulse_comp <- cmpfun(evalLogLikImpulse)
evalLogLikSigmoid_comp <- cmpfun(evalLogLikSigmoid)

################################################################################
### Main function / Wrapper function
################################################################################

#' ImpulseDE2 wrapper
#'
#' Wrapper to run ImpulseDE2 on bulk omics count data.
#' This wrapper can perform the entire analysis pipeline of 
#' ImpulseDE2 on its own if the right parameters are supplied.
#' Core results are returned, additional objects may be saved to
#' a directory supplied by the user.
#' To run ImpulseDE2 on bulk omics count data, use the minimal
#' parameter set:
#' \itemize{
#'    \item matCountData
#'    \item dfAnnotation
#'    \item boolCaseCtrl
#'    \itme vecConfounders
#' }
#' Additionally, you can provide:
#' \itemize{
#'    \item scaNProc to set the number of processes for parallelisation.
#'    \item scaQThres to set the cut off for your DE gene list. 
#'    \item boolPlotting and boolSimplePlot to control plotting behaviour.
#'    \item vecDispersionsExternal to supply external dispersion parameters
#'    which may be necessary depending on your confounding factors (runImpulseDE2
#'    will tell you if it is necessary).
#'    \item vecSizeFactorsExternal to supply external size factors.
#'    \item dirTemp to save temporary data to the given directory 
#'    (recommended for trouble shooting and reproducibility of results).
#' }
#' 
#' @details ImpulseDE2 is based on the impulse model proposed by
#' Chechik and Koller (Chechik and Koller, 2009).
#' The computational complexity of ImpulseDE2 is linear in the
#' number of genes and linear in the number of samples.
#' 
#' @aliases ImpulseDE2 wrapper
#' 
#' @seealso Calls the following functions:
#' \code{\link{processData}}, 
#' \code{\link{runDESeq2}},
#' \code{\link{computeNormConst}},
#' \code{\link{fitModels},
#' \code{\link{fitSigmoidModels},
#' \code{\link{runDEAnalysis}}, 
#' \code{\link{plotGenes}}.
#' The following functions are additionally available to the user:
#' \code{\link{compareDEMethods}},
#' \code{\link{fitSigmoidModels}},
#' \code{\link{plotGenes}},  
#' \code{\link{plotHeatmap}},
#' \code{\link{runDEAnalysis}},   
#' \code{\link{simulateDataSet}}.
#' 
#' @param matCountData: (matrix genes x samples) [Default NULL] 
#'    Read count data, unobserved entries are NA.
#' @param dfAnnotation: (data frame samples x covariates) 
#'    {Sample, Condition, Time (numeric), TimeCateg (str)
#'    (and confounding variables if given).}
#'    Annotation table with covariates for each sample.
#' @param boolCaseCtrl: (bool) 
#' 		Whether to perform case-control analysis. Does case-only
#' 		analysis if FALSE.
#' @param vecConfounders: (vector of strings number of confounding variables)
#' 		Factors to correct for during batch correction. Have to 
#' 		supply dispersion factors if more than one is supplied.
#' 		Names refer to columns in dfAnnotation.
#' @param scaNProc: (scalar) [Default 1] Number of processes for 
#'    parallelisation.
#' @param scaQThres: (scalar) [Default NULL] 
#'    FDR-corrected p-value cutoff for significance.
#' @param vecDispersionsExternal: (vector length number of
#'    genes in matCountData) [Default NULL]
#'    Externally generated list of gene-wise dispersion factors
#'    which overides DESeq2 generated dispersion factors.
#' @param vecSizeFactorsExternal: (vector length number of
#'    cells in matCountData) [Default NULL]
#'    Externally generated list of size factors which override
#'    size factor computation in ImpulseDE2.
#' @param boolIdentifyTransients: (bool) [Defaul FALSE]
#'    Whether to identify transiently activated or deactivated 
#'    genes. This involves an additional fitting of sigmoidal models
#'    and hypothesis testing between constant, sigmoidal and impulse model.
#' @param boolPlotting: (bool) [Default FALSE] 
#'    Whether to plot significant DE genes into output pdf.
#'    Consider setting FALSE for large data sets with many hits.
#' @param boolSimplePlot: (bool) [Default FALSE]
#'    Whether to reduce plot to data points and impulse trace without
#'    batch information.
#' @param dirTemp: (dir) Directory to which temporary results are saved.
#' 
#' @return (list length 4)
#' \itemize{
#'    \item vecDEGenes: (list number of genes) Genes IDs identified
#'    as differentially expressed by ImpulseDE2 at threshold \code{scaQThres}.
#'    \item dfDEAnalysis (data frame genes x reported characteristics) 
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
#'    \item lsModelFits: (list length number of conditions fit (1 or 3))
#'    {"case"} or {"case", "control", "combined"}
#'    One model fitting object for each condition:
#'    In case-only DE analysis, only the condition {"case"} is fit.
#'    In case-control DE analysis, the conditions 
#'    {"case", "control","combined} are fit.
#'    Each condition entry is a list of model fits for each gene.
#'    Each gene entry is a list of model fits to the individual models:
#'    Impulse model and constant model (if boolFitConst is TRUE).
#'    At this level, the sigmoid model fit can be added later.
#'    Each model fit per gene is a list of fitting parameters and results.
#'    \itemize{
#'      \item Condition ID: (list length number of genes)
#'      List of fits for each gene to the samples of this condition.
#'      One entry of this format for all conditions fit.
#'      \itemize{
#'        \item Gene ID: (list length 2)
#'        Impulse and constant model fit to gene observations.
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
#'          NULL if boolIdentifyTransients is FALSE.
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
#'          }
#'        }
#'      }
#'    }
#' }
#' Additionally, \code{ImpulseDE2} saves the following objects and tables into
#' the working directory:
#' \itemize{
#'    \item \code{ImpulseDE2_matCountDataProc.RData}
#'    matCountDataProc: (matrix genes x samples) [Default NULL] 
#'    Processed read count data, unobserved entries are NA.
#'    \item \code{ImpulseDE2_dfAnnotationProc.RData} 
#'    dfAnnotationProc: (data frame samples x covariates) 
#'    {Sample, Condition, Time (numeric), TimeCateg (str)
#'    (and confounding variables if given).}
#'    Annotation table with covariates for each sample.
#'    \item \code{ImpulseDE2_vecSizeFactors.RData}
#'    vecSizeFactors: (numeric vector number of samples) 
#'    Model scaling factors for each sample which take
#'    sequencing depth into account (size factors).
#'    \item \code{ImpulseDE2_vecDispersions.RData}
#'    vecDispersions: (vector number of genes) Gene-wise 
#'    negative binomial dispersion hyper-parameter.
#'    \item \code{ImpulseDE2_lsModelFits.RData} (list length number of conditions fit (1 or 3))
#'    {"case"} or {"case", "control", "combined"}
#'    One model fitting object for each condition:
#'    In case-only DE analysis, only the condition {"case"} is fit.
#'    In case-control DE analysis, the conditions 
#'    {"case", "control","combined} are fit.
#'    Each condition entry is a list of model fits for each gene.
#'    Each gene entry is a list of model fits to the individual models:
#'    Impulse model and constant model (if boolFitConst is TRUE).
#'    At this level, the sigmoid model fit can be added later.
#'    Each model fit per gene is a list of fitting parameters and results.
#'    Read the return-value section of runImpulseDE2 for details.
#'    \item \code{ImpulseDE2_dfImpulseDE2Results.RData} dfDEAnalysis (data frame genes x reported characteristics) 
#'    Summary of fitting procedure and 
#'    differential expression results for each gene.
#'    Read the return-value section of runImpulseDE2 for details.
#'    \item \code{ImpulseDE2_vecDEGenes.RData} vecDEGenes: (list number of genes) Genes IDs identified
#'    as differentially expressed by ImpulseDE2 at threshold \code{scaQThres}.
#' }
#' 
#' @author David Sebastian Fischer
#' 
#' @export
runImpulseDE2 <- function(matCountData=NULL, 
  dfAnnotation=NULL,
  boolCaseCtrl=NULL,
  vecConfounders=NULL,
  scaNProc=1, 
  scaQThres=NULL,
  vecDispersionsExternal=NULL,
  vecSizeFactorsExternal=NULL,
  boolIdentifyTransients=FALSE,
  boolPlotting=FALSE,
  boolSimplePlot=FALSE,
  dirTemp=NULL ){
  
  print("ImpulseDE2 v1.0 for count data")
  
  tm_runImpulseDE2 <- system.time({    
    # 1. Process input data 
    print("# Process input")
    lsProcessedData <- processData(
      dfAnnotation=dfAnnotation,
      matCountData=matCountData,
      boolCaseCtrl=boolCaseCtrl, 
      vecConfounders=vecConfounders,
      vecDispersionsExternal=vecDispersionsExternal,
      vecSizeFactorsExternal=vecSizeFactorsExternal)
    
    matCountDataProc <- lsProcessedData$matCountDataProc
    dfAnnotationProc <- lsProcessedData$dfAnnotationProc
    vecSizeFactorsExternalProc <- lsProcessedData$vecSizeFactorsExternalProc
    vecDispersionsExternalProc <- lsProcessedData$vecDispersionsExternalProc
    
    # Set parallelisation
    if(scaNProc > 1){ register(MulticoreParam(workers=scaNProc)) 
    } else { register(SerialParam()) }
    
    if(!is.null(dirTemp)){
      save(matCountDataProc,file=file.path(dirTemp,"ImpulseDE2_matCountDataProc.RData"))
      save(dfAnnotationProc,file=file.path(dirTemp,"ImpulseDE2_dfAnnotationProc.RData"))
    }
    
    # 2. Run DESeq2
    # Use dispersion factors from DESeq2 if
    if(is.null(vecDispersionsExternal)){
      print("# Run DESeq2: Using dispersion factors computed by DESeq2.")
      tm_runDESeq2 <- system.time({
        vecDispersions <- runDESeq2(
          dfAnnotationProc=dfAnnotationProc,
          matCountDataProc=matCountDataProc,
          boolCaseCtrl=boolCaseCtrl,
          vecConfounders=vecConfounders)
      })
      print(paste0("Consumed time: ",round(tm_runDESeq2["elapsed"]/60,2), " min."))
    } else {
      # Use externally provided dispersions and reorder
      print("# Using externally supplied dispersion factors.")
      vecDispersions <- vecDispersionsExternalProc
    }
    if(!is.null(dirTemp)){
      save(vecDispersions,file=file.path(dirTemp,"ImpulseDE2_vecDispersions.RData"))
    }
    
    # 3. Compute size factors
    print("# Compute size factors")
    vecSizeFactors <- computeNormConst(
    	matCountDataProc=matCountDataProc,
      vecSizeFactorsExternal=vecSizeFactorsExternalProc )
    if(!is.null(dirTemp)){
      save(vecSizeFactors,file=file.path(dirTemp,"ImpulseDE2_vecSizeFactors.RData"))
    }
    
    ###  4. Fit null and alternative model to each gene
    print("# Fitting null and alternative model to the genes")
    tm_fitImpulse <- system.time({
      lsModelFits <- fitModels(
        matCountDataProc=matCountDataProc, 
        vecDispersions=vecDispersions,
        vecSizeFactors=vecSizeFactors,
        dfAnnotationProc=dfAnnotationProc,
        vecConfounders=vecConfounders,
        boolCaseCtrl=boolCaseCtrl)
    })
    if(!is.null(dirTemp)){
      save(lsModelFits,file=file.path(dirTemp,"ImpulseDE2_lsModelFits.RData"))
    }
    print(paste0("Consumed time: ",round(tm_fitImpulse["elapsed"]/60,2)," min."))
    
    ###  5. Fit sigmoid model to case condition if desired
    if(boolIdentifyTransients){
      print("# Fitting sigmoid model to case condition")
      tm_fitSigmoid <- system.time({
        lsModelFits <- fitSigmoidModels(
          matCountDataProc=matCountDataProc,
          lsModelFits=lsModelFits,
          vecDispersions=vecDispersions,
          vecSizeFactors=vecSizeFactors,
          dfAnnotationProc=dfAnnotationProc,
          vecConfounders=vecConfounders,
          strCondition="case")
      })
      if(!is.null(dirTemp)){
        save(lsModelFits,file=file.path(dirTemp,"ImpulseDE2_lsModelFits.RData"))
      }
      print(paste0("Consumed time: ",round(tm_fitSigmoid["elapsed"]/60,2)," min."))
    }
    
    ### 6. Differentially expression analysis based on model fits
    print("# Differentially expression analysis based on model fits")
    dfImpulseDE2Results <- runDEAnalysis(
      matCountDataProc=matCountDataProc,
      dfAnnotationProc=dfAnnotationProc,
      lsModelFits=lsModelFits,
      boolCaseCtrl=boolCaseCtrl,
      vecConfounders=vecConfounders,
      boolIdentifyTransients=boolIdentifyTransients)
    
    if(!is.null(scaQThres)){
      vecDEGenes <- as.vector( 
        dfImpulseDE2Results[as.numeric(dfImpulseDE2Results$padj) <= scaQThres,"Gene"] )
      print(paste0("Found ", length(vecDEGenes)," DE genes",
        " at a FDR corrected p-value cut off of ", scaQThres, "."))
    } else {
      vecDEGenes <- NULL
    }
    if(!is.null(dirTemp)){
      save(dfImpulseDE2Results,file=file.path(dirTemp,"ImpulseDE2_dfImpulseDE2Results.RData"))
      if(!is.null(scaQThres)){
        save(vecDEGenes,file=file.path(dirTemp,"ImpulseDE2_vecDEGenes.RData"))
      }
    }
    
    ### 7. Plot differentially expressed genes
    if(boolPlotting){
      print("# Plot differentially expressed genes")
      if(length(vecDEGenes) > 500){vecDEGenesPlot <- vecDEGenes[1:500]
      } else {vecDEGenesPlot <- vecDEGenes}
      tm_plotGenes <- system.time({
        plotGenes(
          vecGeneIDs=vecDEGenesPlot,
          matCountDataProc=matCountDataProc,
          matBatchFactors=matBatchFactors,
          matSizeFactors=matSizeFactors,
          dfAnnotationProc=dfAnnotationProc, 
          lsModelFits=lsModelFits,
          dfImpulseDE2Results=dfImpulseDE2Results,
          vecRefPval=vecRefPval,
          strMode=strMode,
          strCaseName=strCaseName, 
          strControlName=strControlName, 
          strFileNameSuffix="DEgenes", 
          strPlotTitleSuffix="", 
          strPlotSubtitle="",
          strNameMethod2=strRefMethod,
          boolSimplePlot=boolSimplePlot,
          boolLogPlot=boolLogPlot)
      })
      print(paste0("Consumed time: ",round(tm_plotGenes["elapsed"]/60,2)," min."))
    }
  })
  print("Finished running ImpulseDE2.")
  print(paste0("TOTAL consumed time: ",round(tm_runImpulseDE2["elapsed"]/60,2)," min."))
  
  return(list(
    dfImpulseDE2Results=dfImpulseDE2Results,
    vecDEGenes=vecDEGenes,
    lsModelFits=lsModelFits
  ))
}
