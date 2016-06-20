################################################################################
########################     ImpulseDE2 package     ############################
################################################################################

### Version 1.0
### Author David Sebastian Fischer

################################################################################
### Libraries and source code
################################################################################

library(compiler)
library(parallel)
library(DESeq2)
library(BiocParallel)

# Source functions in .R files from same directory as this function.
setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/building/code_files")
source("srcImpulseDE2_processData.R")
source("srcImpulseDE2_computeNormConst.R")
source("srcImpulseDE2_runDESeq2.R")
source("srcImpulseDE2_calcImpulse.R")
# Compile function
calcImpulse_comp <- cmpfun(calcImpulse)
source("srcImpulseDE2_CostFunctionsFit.R")
# Compile functions
evalLogLikImpulseBatch_comp <- cmpfun(evalLogLikImpulseBatch)
evalLogLikZINB_comp <- cmpfun(evalLogLikZINB)
evalLogLikImpulseSC_comp <- cmpfun(evalLogLikImpulseSC)
source("srcImpulseDE2_fitImpulse.R")
source("srcImpulseDE2_computePval.R")
source("srcImpulseDE2_plotDEGenes.R")

################################################################################
### Main function
################################################################################

#' Differential expression analysis using impulse models
#'
#' Fits an impulse model to time course data and uses this model as a basis
#' to detect differentially expressed genes. Differential expression is
#' either differential expression of a gene over time within one condition
#' or differential expression of a gene over time between two conditions (
#' case and control). With a single condition, the alternative model is the
#' impulse fit to the time course data and the null model is the mean fit.
#' The mean fit models no differential behaviour over time. With two 
#' conditions, the alternative model is separate impulse fits to case and
#' control data and the null model is an impulse fit to the combined data.
#' Here, the impulse fit to the combined data models no differential expression
#' between the conditions.
#' 
#' @details \code{ImpulseDE2} is based on the impulse model proposed by
#' Chechik and Koller (Chechik and Koller, 2009). The impulse model
#' models the response of gene activity read outs (such as RNAseq counts)
#' to environmental or developmental stimuli as the product
#' of two simgoids. This model can capture simple time course
#' patterns, such as plateaus, increase and decrease. ImpulseDE2 uses
#' the impulse model to identify differential activity over time on any
#' type of count data which follows the negative binomial distribution
#' (as frequently encountered in sequencing data). ImpulseDE2 performs
#' fitting of the impulse model and a mean model to data and evaluates
#' the fit. The computational complexity of ImpulseDE2 is linear in the
#' number of genes and linear in the number of samples.
#' \enumerate{
#'  \item Impulse fitting: The impulse model is fitted based on the assumption
#'    that the input count data follow a negative binomial distribution with
#'    dispersion as identified by DESeq2. The impulse model does not have
#'    a closed form maximum likelihood parameter estimate and must therefore
#'    be inferred from numerical optimisation.
#'  \enumerate{
#'    \item Initialisation: Initialisation is performed twice for each gene,
#'      based on a peak and a valley model. The parameters representing these 
#'      models reflect the form that these two models would have given the count
#'      data of each gene and are specific to each gene.
#'    \item Optimisation: The cost function for the fit is the log likelihood
#'      of the data which is evaluated based on negative binomial likelihoods
#'      at each observed time point, with the value of the impulse model as the
#'      mean and the gene dispersion inferred using DESeq2 as the 
#'      dispersion. Numerical optimisation is performed using the
#'      BFGS algorithm. In analogy to generalised linear models for count data,
#'      the fitting of the parameters representing limit behaviours of the two
#'      sigmoids (which are counts) are fitted in log space so that they cannot
#'      adopt negative values.
#'    \item Fit selection: The fit with the higher log likelihood of the two
#'      initialisation is selected and kept as a maximum likelihood estimate.
#'  }
#'  \item Mean fitting: The mean model is a single negative binomial and
#'    serves as the null model in the case of differential expression over
#'    time within a single condition. The maximum likelihood estimate of
#'    the mean of a negative binomial distribution given the data  and 
#'    the dispersion is the sample average. This closed form solution 
#'    is used here.
#'  \item Fit evaluation: The model comparison statistic is the deviance
#'    2 * (loglikelihood(H1) - loglikelihood(H0)). The deviance is chi-squared
#'    distributed with the difference in degrees of freedom of both models
#'    as degrees of freedom if the null model is contained in the alternative
#'    model, which is given in both modes of differential expression analysis
#'    with ImpulseDE (with and without control). Therefore p-values for 
#'    differential expression are computed based on the chi-squared distribution.
#'    The p-values are the FDR corrected (Benjamini and Hochberg, 1995).
#' }
#' 
#' @aliases ImpulseDE2
#' 
#' @param matCountData: (matrix genes x replicates) [Default NULL] 
#'    Count data of all conditions, unobserved entries are NA.
#' @param dfAnnotation: (Table) [Default NULL] 
#'    Lists co-variables of samples: 
#'    Sample, Condition, Time (numeric), (and Timecourse).
#' @param strCaseName: (str) [Default NULL] 
#'    Name of the case condition in \code{dfAnnotation}.
#' @param strControlName: (str) [Default NULL] 
#'    Name of the control condition in \code{dfAnnotation}.
#' @param strMode: (str) [Default "batch"] 
#'    {"batch","longitudinal","singlecell"}
#'    Mode of model fitting.
#' @param nProc: (scalar) [Default 1] Number of processes for 
#'    parallelisation.
#' @param Q_value: (scalar) [Default 0.01] 
#'    FDR-corrected p-value cutoff for significance.
#' @param scaSmallRun: (integer) [Default NULL] Number of rows
#'    on which ImpulseDE2 is supposed to be run, the full
#'    data set is only used for size factor estimation.
#' @param boolPlotting: (bool) [TRUE] 
#'    Whether to plot significant DE genes into output pdf.
#'    Consider setting FALSE for large data sets with many hits.
#' @param lsPseudoDE: (list) [Defaul NULL] 
#' @param vecDispersionsExternal: (vector length number of
#'    genes in matCountData) [Default NULL]
#'    Externally generated list of gene-wise dispersion factors
#'    which overides DESeq2 generated dispersion factors.
#' @param vecSizeFactorsExternal: (vector length number of
#'    cells in matCountData) [Default NULL]
#'    Externally generated list of size factors which override
#'    size factor computation in ImpulseDE2.
#' @param boolRunDESeq2: (bool) [Default TRUE]
#'    Whether to run DESeq2.
#' @param boolSimplePlot: (bool) [Default FALSE]
#'    Whether to reduce plot to data points and impulse trace.
#' @param boolLogPlot: (bool) [Default FALSE]
#'    Whether to plot in counts in log space.
#' 
#' @return (list length 4)
#' \itemize{
#'    \item vecDEGenes: (list number of genes) Genes IDs identified
#'        as differentially expressed by ImpulseDE2 at threshold \code{Q_value}.
#'    \item dfImpulseResults: (data frame) ImpulseDE2 results.
#'    \item lsImpulseFits: (list) List of matrices which
#'        contain parameter fits and model values for given time course for the
#'        case condition (and control and combined if control is present).
#'        Each parameter matrix is called parameter_'condition' and has the form
#'        (genes x [beta, h0, h1, h2, t1, t2, logL_H1, converge_H1, mu, logL_H0, 
#'        converge_H0]) where beta to t2 are parameters of the impulse
#'        model, mu is the single parameter of the mean model, logL are
#'        log likelihoods of full (H1) and reduced model (H0) respectively, converge
#'        is convergence status of numerical optimisation of model fitting by
#'        \code{optim} from \code{stats} of either model. Each value matrix is called
#'        value_'condition' and has the form (genes x time points) and contains the
#'        counts predicted by the impulse model at the observed time points.
#'    \item dfDESeq2Results: (data frame) DESeq2 results. NULL if DESeq2 is not run.
#' }
#' Additionally, \code{ImpulseDE2} saves the following objects and tables into
#' the working directory:
#' \itemize{
#'    \item \code{ImpulseDE2_matCountDataProc.RData} (2D array genes x replicates) 
#'        Count data: Reduced version of \code{matCountData}. For internal use.
#'    \item \code{ImpulseDE2_dfAnnotationProc.RData} (data frame) Annotation table.
#'    \item \code{ImpulseDE2_matSizeFactors.RData} (numeric matrix genes x samples) 
#'        Model scaling factors for each observation which take
#'        sequencing depth into account (size factors). One size
#'        factor per sample - rows of this matrix are equal.
#'    \item \code{ImpulseDE2_lsMatTranslationFactors.RData} 
#'      (list {case} or {case, ctrl, combined}) List of
#'      translation factor matrices. NULL if strMode not
#'      "longitudinal".
#'      \itemize{
#'        \item matTranslationFactorsAll: (numeric matrix genes x samples) 
#'        Model scaling factors for each observation which take
#'        longitudinal time series mean within a gene into account 
#'        (translation factors). Computed based based on all samples.
#'        \item matTranslationFactorsCase: As matTranslationFactorsAll
#'        for case samples only, non-case samples set NA.
#'        \item matTranslationFactorsCtrl: As matTranslationFactorsAll
#'        for control samples only, non-control samples set NA.
#'      }
#'    \item \code{ImpulseDE2_vecDispersions.RData} (vector number of genes) Inverse 
#'        of gene-wise negative binomial dispersion coefficients computed by DESeq2.
#'    \item \code{ImpulseDE2_dfDESeq2Results.RData} (data frame) DESeq2 results.
#'    \item \code{ImpulseDE2_lsImpulseFits.RData} (list) List of matrices which
#'        contain parameter fits and model values for given time course for the
#'        case condition (and control and combined if control is present).
#'        Each parameter matrix is called parameter_'condition' and has the form
#'        (genes x [beta, h0, h1, h2, t1, t2, logL_H1, converge_H1, mu, logL_H0, 
#'        converge_H0]) where beta to t2 are parameters of the impulse
#'        model, mu is the single parameter of the mean model, logL are
#'        log likelihoods of full (H1) and reduced model (H0) respectively, converge
#'        is convergence status of numerical optimisation of model fitting by
#'        \code{optim} from \code{stats} of either model. Each value matrix is called
#'        value_'condition' and has the form (genes x time points) and contains the
#'        counts predicted by the impulse model at the observed time points.
#'    \item \code{ImpulseDE2_dfImpulseResults.RData} (data frame) ImpulseDE2 results.
#'    \item \code{ImpulseDE2_vecDEGenes.RData} (list number of genes) Genes IDs identified
#'        as differentially expressed by ImpulseDE2 at threshold \code{Q_value}.
#'    \item \code{ImpulseDE2_ClusterOut.txt} Text-file with stdout and stderr from
#'        cluster created in \code{fitImpulse}.
#' }
#' 
#' @seealso Calls the following functions:
#' \code{\link{processData}}, \code{\link{runDESeq2}},
#' \code{\link{fitImpulse}},
#' \code{\link{computePval}}, \code{\link{plotDEGenes}}.
#' 
#' @author David Sebastian Fischer
#' 
#' @export

runImpulseDE2 <- function(matCountData=NULL, 
  dfAnnotation=NULL,
  strCaseName = NULL, 
  strControlName=NULL, 
  strMode="batch",
  strSCMode="clustered",
  nProc=1, 
  Q_value=0.01,
  scaSmallRun=NULL,
  boolPlotting=TRUE,
  lsPseudoDE=NULL, 
  vecDispersionsExternal=NULL,
  vecSizeFactorsExternal=NULL,
  boolRunDESeq2=TRUE,
  boolSimplePlot=FALSE, 
  boolLogPlot=FALSE ){
  
  NPARAM=6
  
  print("ImpulseDE2 v1.0 for count data")
  
  tm_runImpulseDE2 <- system.time({    
    # 1. Process input data 
    print("1. Prepare data")
    lsProcessedData <- processData(
      dfAnnotation=dfAnnotation,
      matCountData=matCountData,
      scaSmallRun=scaSmallRun,
      strControlName=strControlName, 
      strCaseName=strCaseName,
      strMode=strMode,
      lsPseudoDE=lsPseudoDE,
      vecDispersionsExternal=vecDispersionsExternal,
      vecSizeFactorsExternal=vecSizeFactorsExternal,
      boolRunDESeq2=boolRunDESeq2 )
    
    matCountDataProc <- lsProcessedData$matCountDataProc
    matCountDataProcFull <- lsProcessedData$matCountDataProcFull
    dfAnnotationProc <- lsProcessedData$dfAnnotationProc
    matProbNB <- lsProcessedData$matProbNB
    matDropoutRate <- lsProcessedData$matDropout
    matMuCluster <- lsProcessedData$matMuCluster
    vecClusterAssignments <- lsProcessedData$vecClusterAssignments
    vecCentroids <- lsProcessedData$vecCentroids
    
    save(matCountDataProc,file=file.path(getwd(),"ImpulseDE2_matCountDataProc.RData"))
    save(matCountDataProcFull,file=file.path(getwd(),"ImpulseDE2_matCountDataProcFull.RData"))
    save(dfAnnotationProc,file=file.path(getwd(),"ImpulseDE2_dfAnnotationProc.RData"))
    save(matProbNB,file=file.path(getwd(),"ImpulseDE2_matProbNB.RData"))
    save(matDropoutRate,file=file.path(getwd(),"ImpulseDE2_matDropoutRate.RData"))
    save(vecClusterAssignments,file=file.path(getwd(),"ImpulseDE2_vecClusterAssignments.RData"))
    save(vecCentroids,file=file.path(getwd(),"ImpulseDE2_vecCentroids.RData"))
    
    # 2. Compute normalisation constants
    print("2. Compute Normalisation constants")
    lsNormConst <- computeNormConst(
      matCountDataProcFull=matCountDataProcFull,
      matCountDataProc=matCountDataProc,
      scaSmallRun=scaSmallRun,
      dfAnnotationProc=dfAnnotationProc,
      strCaseName=strCaseName,
      strControlName=strControlName, 
      strMode=strMode,
      vecSizeFactorsExternal=vecSizeFactorsExternal)
    lsMatTranslationFactors <- lsNormConst$lsMatTranslationFactors
    matSizeFactors <- lsNormConst$matSizeFactors
    save(lsMatTranslationFactors,file=file.path(getwd(),"ImpulseDE2_lsMatTranslationFactors.RData"))
    save(matSizeFactors,file=file.path(getwd(),"ImpulseDE2_matSizeFactors.RData"))
    
    # 3. Run DESeq2
    if(boolRunDESeq2){
      print("3. Run DESeq2")
      tm_runDESeq2 <- system.time({
        lsDESeq2Results <- runDESeq2(
          dfAnnotationProc=dfAnnotationProc,
          matCountDataProc=matCountDataProc,
          nProc=nProc,
          strControlName=strControlName,
          strMode=strMode)
      })
      print(paste("Consumed time: ",round(tm_runDESeq2["elapsed"]/60,2),
        " min",sep=""))
      dfDESeq2Results <- lsDESeq2Results[[2]]
      vecRefPval <- dfDESeq2Results$padj
      names(vecRefPval) <- rownames(dfDESeq2Results)
      strRefMethod <- "DESeq2"
    } else { 
      print("3. Not running DESeq2")
      dfDESeq2Results <- NULL
      vecRefPval <- NULL
      strRefMethod <- NULL
    }
    if(is.null(vecDispersionsExternal) & boolRunDESeq2){
      # Use dispersions generated by DESeq2
      print("Using dispersion factors computed by DESeq2.")
      vecDispersions <- lsDESeq2Results[[1]]
      
      # is.null(vecDispersionsExternal) & !boolRunDESeq2 exception 
      # is caught in processData()
    } else {
      # Use externally provided dispersions and reorder
      print("Using externally supplied dispersion factors.")
      vecDispersions <- vecDispersionsExternal[rownames(matCountDataProc)]
    }
    save(vecDispersions,file=file.path(getwd(),"ImpulseDE2_vecDispersions.RData"))
    save(dfDESeq2Results,file=file.path(getwd(),"ImpulseDE2_dfDESeq2Results.RData"))
    
    ###  4. Fit Impule model to each gene 
    print("4. Fitting Impulse model to the genes")
    tm_fitImpulse <- system.time({
      lsImpulseFits <- fitImpulse(
        matCountDataProc=matCountDataProc, 
        vecDispersions=vecDispersions,
        matDropoutRate=matDropoutRate,
        matProbNB=matProbNB,
        vecClusterAssignments=vecClusterAssignments,
        matSizeFactors=matSizeFactors,
        lsMatTranslationFactors=lsMatTranslationFactors,
        dfAnnotationProc=dfAnnotationProc, 
        strCaseName=strCaseName, 
        strControlName=strControlName,
        strMode=strMode,
        strSCMode=strSCMode,
        nProc=nProc, 
        NPARAM=NPARAM)
    })
    save(lsImpulseFits,file=file.path(getwd(),"ImpulseDE2_lsImpulseFits.RData"))
    print(paste("Consumed time: ",round(tm_fitImpulse["elapsed"]/60,2),
      " min",sep=""))
    
    ### 5. Detect differentially expressed genes
    print("5. DE analysis")
    dfImpulseResults <- computePval(
      matCountDataProc=matCountDataProc,
      vecDispersions=vecDispersions,
      dfAnnotationProc=dfAnnotationProc,
      lsImpulseFits=lsImpulseFits,
      strCaseName=strCaseName, 
      strControlName=strControlName,
      strMode=strMode, 
      NPARAM=NPARAM)
    
    vecDEGenes <- as.character(as.vector( 
      dfImpulseResults[as.numeric(dfImpulseResults$adj.p) <= Q_value,"Gene"] ))
    save(dfImpulseResults,file=file.path(getwd(),"ImpulseDE2_dfImpulseResults.RData"))
    save(vecDEGenes,file=file.path(getwd(),"ImpulseDE2_vecDEGenes.RData"))
    print(paste("Found ", length(vecDEGenes)," DE genes",sep=""))
    
    ### 6. Plot differentially expressed genes
    if(boolPlotting){
      print("6. Plot differentially expressed genes")
      tm_plotDEGenes <- system.time({
        plotDEGenes(
          vecGeneIDs=vecDEGenes,
          matCountDataProc=matCountDataProc,
          lsMatTranslationFactors=lsMatTranslationFactors,
          matSizeFactors=matSizeFactors,
          dfAnnotationProc=dfAnnotationProc, 
          lsImpulseFits=lsImpulseFits,
          matMuCluster=matMuCluster,
          vecCentroids=vecCentroids,
          vecClusterAssignments=vecClusterAssignments,
          dfImpulseResults=dfImpulseResults,
          vecRefPval=vecRefPval,
          strMode=strMode,
          strSCMode=strSCMode,
          strCaseName=strCaseName, 
          strControlName=strControlName, 
          strFileNameSuffix="DEgenes", 
          strPlotTitleSuffix="", 
          strPlotSubtitle="",
          strNameMethod2=strRefMethod,
          boolSimplePlot=boolSimplePlot,
          boolLogPlot=boolLogPlot,
          NPARAM=NPARAM)
      })
      print(paste("Consumed time: ",round(tm_plotDEGenes["elapsed"]/60,2),
        " min",sep=""))
    }
  })
  print("Finished ImpulseDE2.")
  print(paste("TOTAL consumed time: ",round(tm_runImpulseDE2["elapsed"]/60,2),
    " min",sep=""))
  
  return(list(
    "vecDEGenes"=vecDEGenes,
    "dfImpulseResults"=dfImpulseResults,
    "lsImpulseFits"=lsImpulseFits,
    "dfDESeq2Results"=dfDESeq2Results))
}
