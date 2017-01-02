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

source("srcImpulseDE2_classUnions.R")
source("srcImpulseDE2_classImpulseDE2Object.R")
source("srcImpulseDE2_evalImpulse.R")
source("srcImpulseDE2_evalSigmoid.R")
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
#' @param boolVerbose: (bool) [Default TRUE] Whether to print
#'    progress to stdout.
#' @param dirTemp: (dir) Directory to which temporary results are saved.
#' 
#' @return (object of class ImpulseDE2Object)
#' This object can be treated as a list with 2 elements:
#' (list length 2)
#' \itemize{
#'    \item vecDEGenes: (list number of genes) Genes IDs identified
#'    as differentially expressed by ImpulseDE2 at threshold \code{scaQThres}.
#'    \item dfDEAnalysis (data frame samples x reported characteristics) 
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
#'      \item allZero: (bool) Whether there were no observed non-zero observations of this gene.
#'      If TRUE, fitting and DE analsysis were skipped and entry is NA.
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
#'      \item converge_sigmoid: Convergence status of optim for 
#'      sigmoid model fit to samples of case condition.
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
  boolVerbose=TRUE,
  dirTemp=NULL ){
  
  strMessage <- "ImpulseDE2 v1.0 for count data"
  if(boolVerbose) print(strMessage)
  strReport <- strMessage
  
  tm_runImpulseDE2 <- system.time({    
    # 1. Process input data 
    strMessage <- "# Process input"
    if(boolVerbose) print(strMessage)
    strReport <- paste0(strReport, "\n", strMessage)
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
    if(boolVerbose) write(lsProcessedData$strReportProcessing, file="", ncolumns=1)
    strReport <- paste0(strReport, lsProcessedData$strReportProcessing)
    
    # Attempt accessing temp directory
    if(!is.null(dirTemp)){
      dirCurrent <- getwd()
      tryCatch({
        setwd(dirTemp)
      }, error=function(strErrorMsg){
        stop(paste0("ERROR: dirTemp not available: ", strErrorMsg))
      })
      setwd(dirCurrent)
    }
    
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
      strMessage <- "# Run DESeq2: Using dispersion factors computed by DESeq2."
      if(boolVerbose) print(strMessage)
      strReport <- paste0(strReport, "\n", strMessage)
      tm_runDESeq2 <- system.time({
        vecDispersions <- runDESeq2(
          dfAnnotationProc=dfAnnotationProc,
          matCountDataProc=matCountDataProc,
          boolCaseCtrl=boolCaseCtrl,
          vecConfounders=vecConfounders)
      })
      strMessage <- paste0("Consumed time: ",round(tm_runDESeq2["elapsed"]/60,2), " min.")
      if(boolVerbose) print(strMessage)
      strReport <- paste0(strReport, "\n", strMessage)
    } else {
      # Use externally provided dispersions and reorder
      strMessage <- "# Using externally supplied dispersion factors."
      if(boolVerbose) print(strMessage)
      strReport <- paste0(strReport, "\n", strMessage)
      vecDispersions <- vecDispersionsExternalProc
    }
    if(!is.null(dirTemp)){
      save(vecDispersions,file=file.path(dirTemp,"ImpulseDE2_vecDispersions.RData"))
    }
    
    # 3. Compute size factors
    strMessage <- "# Compute size factors"
    if(boolVerbose) print(strMessage)
    strReport <- paste0(strReport, "\n", strMessage)
    vecSizeFactors <- computeNormConst(
    	matCountDataProc=matCountDataProc,
      vecSizeFactorsExternal=vecSizeFactorsExternalProc )
    if(!is.null(dirTemp)){
      save(vecSizeFactors,file=file.path(dirTemp,"ImpulseDE2_vecSizeFactors.RData"))
    }
    
    # 4. Create instance of ImpulseDE2Object
    # Create ImpulseDE2 object
    objectImpulseDE2 <- new('ImpulseDE2Object',
                            dfImpulseDE2Results = NULL,
                            vecDEGenes          = NULL,
                            lsModelFits         = NULL,
                            matCountDataProc    = matCountDataProc,
                            dfAnnotationProc    = dfAnnotationProc,
                            vecSizeFactors      = vecSizeFactors,
                            vecDispersions      = vecDispersions,
                            boolCaseCtrl        = boolCaseCtrl,
                            vecConfounders      = vecConfounders,
                            scaNProc            = scaNProc, 
                            scaQThres           = scaQThres,
                            strReport           = strReport )
    
    #  5. Fit null and alternative model to each gene
    strMessage <- "# Fitting null and alternative model to the genes"
    if(boolVerbose) print(strMessage)
    objectImpulseDE2@strReport <- paste0(objectImpulseDE2@strReport, "\n", strMessage)
    tm_fitImpulse <- system.time({
      objectImpulseDE2 <- fitModels(
        objectImpulseDE2=objectImpulseDE2,
        vecConfounders=vecConfounders,
        boolCaseCtrl=boolCaseCtrl)
    })
    if(!is.null(dirTemp)){
      save(lsModelFits,file=file.path(dirTemp,"ImpulseDE2_lsModelFits.RData"))
    }
    strMessage <- paste0("Consumed time: ",round(tm_fitImpulse["elapsed"]/60,2)," min.")
    if(boolVerbose) print(strMessage)
    objectImpulseDE2@strReport <- paste0(objectImpulseDE2@strReport, "\n", strMessage)
    
    # 6. Fit sigmoid model to case condition if desired
    if(boolIdentifyTransients){
      strMessage <- "# Fitting sigmoid model to case condition"
      if(boolVerbose) print(strMessage)
      objectImpulseDE2@strReport <- paste0(objectImpulseDE2@strReport, "\n", strMessage)
      tm_fitSigmoid <- system.time({
        objectImpulseDE2 <- fitSigmoidModels(
          objectImpulseDE2=objectImpulseDE2,
          vecConfounders=vecConfounders,
          strCondition="case")
      })
      if(!is.null(dirTemp)){
        save(objectImpulseDE2@lsModelFits,file=file.path(dirTemp,"ImpulseDE2_lsModelFits.RData"))
      }
      strMessage <- paste0("Consumed time: ",round(tm_fitSigmoid["elapsed"]/60,2)," min.")
      if(boolVerbose) print(strMessage)
      objectImpulseDE2@strReport <- paste0(objectImpulseDE2@strReport, "\n", strMessage)
    }
    
    # 7. Differentially expression analysis based on model fits
    strMessage <- "# Differentially expression analysis based on model fits"
    if(boolVerbose) print(strMessage)
    objectImpulseDE2@strReport <- paste0(objectImpulseDE2@strReport, "\n", strMessage)
    dfImpulseDE2Results <- runDEAnalysis(
      objectImpulseDE2=objectImpulseDE2,
      vecAllIDs=rownames(matCountData),
      boolCaseCtrl=boolCaseCtrl,
      vecConfounders=vecConfounders,
      boolIdentifyTransients=boolIdentifyTransients)
    
    if(!is.null(scaQThres)){
      vecDEGenes <- as.vector( 
        dfImpulseDE2Results[as.numeric(dfImpulseDE2Results$padj) <= scaQThres,"Gene"] )
      strMessage <- paste0("Found ", length(vecDEGenes)," DE genes",
        " at a FDR corrected p-value cut off of ", scaQThres, ".")
      if(boolVerbose) print(strMessage)
      objectImpulseDE2@strReport <- paste0(objectImpulseDE2@strReport, "\n", strMessage)
    } else {
      vecDEGenes <- NULL
    }
    if(!is.null(dirTemp)){
      save(dfImpulseDE2Results,file=file.path(dirTemp,"ImpulseDE2_dfImpulseDE2Results.RData"))
      if(!is.null(scaQThres)){
        save(vecDEGenes,file=file.path(dirTemp,"ImpulseDE2_vecDEGenes.RData"))
      }
    }
    objectImpulseDE2@dfImpulseDE2Results <- dfImpulseDE2Results
    objectImpulseDE2@vecDEGenes <- vecDEGenes
  })
  strMessage <- "Finished running ImpulseDE2."
  if(boolVerbose) print(strMessage)
  objectImpulseDE2@strReport <- paste0(objectImpulseDE2@strReport, "\n", strMessage)
  strMessage <- paste0("TOTAL consumed time: ",round(tm_runImpulseDE2["elapsed"]/60,2)," min.")
  if(boolVerbose) print(strMessage)
  objectImpulseDE2@strReport <- paste0(objectImpulseDE2@strReport, "\n", strMessage)
  
  return(objectImpulseDE2)
}
