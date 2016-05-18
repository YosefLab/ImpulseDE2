################################################################################
########################     ImpulseDE2 package     ############################
################################################################################

### Date laste change:     May 15th 2016
### Author David Sebastian Fischer

################################################################################
### Libraries and source code
################################################################################

library(compiler)
library(parallel)
library(DESeq2)
library(BiocParallel)

# Source functions in .R files from same directory as this function.
setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files")
#setwd("/data/yosef2/users/fischerd/code/ImpulseDE2")
# Prepares the data and annotation for internal use
source("srcImpulseDE2_processData.R")
# Compute normalisation constants
source("srcImpulseDE2_computeNormConst.R")
# Wrapper fur running DESeq2
source("srcImpulseDE2_runDESeq2.R")
# Compute value of impulse function given parameters
source("srcImpulseDE2_calcImpulse.R")

# Compile function
calcImpulse_comp <- cmpfun(calcImpulse)

# Cost functions for fitting
source("srcImpulseDE2_CostFunctionsFit.R")

# Compile functions
evalLogLikImpulseBatch_comp <- cmpfun(evalLogLikImpulseBatch)
evalLogLikImpulseByTC_comp <- cmpfun(evalLogLikImpulseByTC)
evalLogLikImpulseSC_comp <- cmpfun(evalLogLikImpulseSC)

# Fit impulse model to a timecourse dataset
source("srcImpulseDE2_fitImpulse.R")
# Detect differentially expressed genes over time
source("srcImpulseDE2_computePval.R")
# Plot the impulse fits and input data
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
#' the fit. The computational complexity of ImpulseDE2 is O(N), where
#' N is the number of genes or regions observed.
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
#' @param matCountData (matrix genes x replicates) [Default NULL] Count data of all conditions, 
#'    unobserved entries are NA. Column labels are replicate names, row labels
#'    gene names.
#' @param dfAnnotationFull (Table) [Default NULL] Lists co-variables of individual replicates: 
#'    Replicate, Sample, Condition, Time. Time must be numeric.
#' @param strCaseName (str) [Default NULL] Name of the case condition in \code{dfAnnotationFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#'    \code{dfAnnotationFull}.
#' @param nProc (scalar) [Default 3] Number of processes for parallelisation. The
#'    specified value is internally changed to \code{min(detectCores() - 1, nProc)} 
#'    using the \code{detectCores} function from the package \code{parallel} to avoid overload.
#' @param Q_value (scalar) [Default 0.01] FDR-corrected p-value cutoff for significance.
#' @param boolPlotting (bool) [TRUE] Whether to plot significant DE genes into output pdf.
#'    Consider setting FALSE for large data sets with many hits.
#' 
#' @return (list length 4) with the following elements:
#' \itemize{
#'    \item \code{lsDEGenes} (list number of genes) Genes IDs identified
#'        as differentially expressed by ImpulseDE2 at threshold \code{Q_value}.
#'    \item \code{dfImpulseResults} (data frame) ImpulseDE2 results.
#'    \item \code{lsImpulseFits} (list) List of matrices which
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
#'    \item \code{dfDESeq2Results} (data frame) DESeq2 results.
#' }
#' Additionally, \code{ImpulseDE2} saves the following objects and tables into
#' the working directory:
#' \itemize{
#'    \item \code{ImpulseDE2_arr2DCountData.RData} (2D array genes x replicates) 
#'        Count data: Reduced version of \code{matCountData}. For internal use.
#'    \item \code{ImpulseDE2_dfAnnotationFull.RData} (data frame) Annotation table.
#'    \item \code{ImpulseDE2_vecNormConst.RData} (matrix samples x replicates) Normalisation
#'    constants for each replicate. Missing samples are set NA.
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
#'    \item \code{ImpulseDE2_lsDEGenes.RData} (list number of genes) Genes IDs identified
#'        as differentially expressed by ImpulseDE2 at threshold \code{Q_value}.
#'    \item \code{ImpulseDE2_ClusterOut.txt} Text-file with stdout and stderr from
#'        cluster created in \code{fitImpulse}.
#' }
#' 
#' @seealso Coordinates the following source functions:
#' \code{\link{processData}}, \code{\link{runDESeq2}},
#' \code{\link{fitImpulse}}, \code{\link{evalLogLikImpulse}}, 
#' \code{\link{evalLogLikMean}}, \code{\link{calcImpulse}},
#' \code{\link{computePval}}, \code{\link{plotDEGenes}}.
#' Calls directly: \code{\link{processData}}, \code{\link{runDESeq2}},
#' \code{\link{fitImpulse}}, \code{\link{computePval}},
#' \code{\link{plotDEGenes}}.
#' 
#' @author David Sebastian Fischer
#' 
#' @references Benjamini, Y. and Hochberg, Y. (1995) Controlling the false
#' discovery rate: a practical and powerful approach to multiple testing.
#' J. R. Stat. Soc. Series B Stat. Methodol., 57, 289-300.
#' @references Storey, J.D. et al. (2005) Significance analysis of time course
#' microarray experiments. Proc. Natl. Acad. Sci. USA, 102, 12837-12841.
#' @references Rangel, C., Angus, J., Ghahramani, Z., Lioumi, M., Sotheran, E.,
#' Gaiba, A., Wild, D.L., Falciani, F. (2004) Modeling T-cell activation using
#' gene expression profiling and state-space models. Bioinformatics, 20(9),
#' 1361-72.
#' @references Chechik, G. and Koller, D. (2009) Timing of Gene Expression
#' Responses to Envi-ronmental Changes. J. Comput. Biol., 16, 279-290.
#' @references Yosef, N. et al. (2013) Dynamic regulatory network controlling
#' TH17 cell differentiation. Nature, 496, 461-468.
#' @export

runImpulseDE2 <- function(matCountData=NULL, dfAnnotationFull=NULL,
  strCaseName = NULL, strControlName=NULL, strMode="batch",
  nProc=3, Q_value=0.01, boolPlotting=TRUE,
  lsPseudoDE=NULL, 
  vecDispersionsExternal=NULL, boolRunDESeq2=TRUE ){
  
  NPARAM=6
  
  print("Impulse v1.3 for count data")
  
  tm_runImpulseDE2 <- system.time({    
    # 1. Process input data 
    print("1. Prepare data")
    lsProcessedData <- processData(
      dfAnnotationFull=dfAnnotationFull,
      matCountData=matCountData,
      strControlName=strControlName, 
      strCaseName=strCaseName,
      strMode=strMode,
      lsPseudoDE=lsPseudoDE,
      vecDispersionsExternal=vecDispersionsExternal,
      boolRunDESeq2=boolRunDESeq2 )
    arr2DCountData <- lsProcessedData$arr2DCountData
    matProbNB <- lsProcessedData$matProbNB
    matDropoutRate <- lsProcessedData$matDropout
    
    save(arr2DCountData,file=file.path(getwd(),"ImpulseDE2_arr2DCountData.RData"))
    save(matProbNB,file=file.path(getwd(),"ImpulseDE2_matProbNB.RData"))
    save(matDropoutRate,file=file.path(getwd(),"ImpulseDE2_matDropoutRate.RData"))
    save(dfAnnotationFull,file=file.path(getwd(),"ImpulseDE2_dfAnnotationFull.RData"))
    
    # 2. Compute normalisation constants
    print("2. Compute Normalisation constants")
    vecNormConst <- computeNormConst(
      arr2DCountData=arr2DCountData,
      dfAnnotationFull=dfAnnotationFull,
      matProbNB=matProbNB,
      strMode=strMode)
    save(vecNormConst,file=file.path(getwd(),"ImpulseDE2_vecNormConst.RData"))
    
    # 3. Run DESeq2
    if(boolRunDESeq2){
      print("3. Run DESeq2")
      tm_runDESeq2 <- system.time({
        lsDESeq2Results <- runDESeq2(
          dfAnnotationFull=dfAnnotationFull,
          arr2DCountData=arr2DCountData,
          nProcessesAssigned=nProc,
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
      vecRefPval <- NULL
      strRefMethod <- NULL
    }
    if(is.null(vecDispersionsExternal) & boolRunDESeq2){
      # Use dispersions generated by DESeq2
      vecDispersions <- lsDESeq2Results[[1]]
      
      # is.null(vecDispersionsExternal) & !boolRunDESeq2 exception 
      # is caught in processData()
    } else {
      # Use externally provided dispersions and reorder
      vecDispersions <- vecDispersionsExternal[rownames(arr2DCountData)]
      dfDESeq2Results <- NULL
    }
    save(vecDispersions,file=file.path(getwd(),"ImpulseDE2_vecDispersions.RData"))
    save(dfDESeq2Results,file=file.path(getwd(),"ImpulseDE2_dfDESeq2Results.RData"))
    
    ###  4. Fit Impule model to each gene 
    print("4. Fitting Impulse model to the genes")
    tm_fitImpulse <- system.time({
      lsImpulseFits <- fitImpulse(
        arr2DCountData=arr2DCountData, 
        vecDispersions=vecDispersions,
        matDropoutRate=matDropoutRate,
        matProbNB=matProbNB,
        vecNormConst=vecNormConst,
        dfAnnotationFull=dfAnnotationFull, 
        strCaseName=strCaseName, 
        strControlName=strControlName,
        strMode=strMode,
        nProcessesAssigned=nProc, 
        NPARAM=NPARAM)
    })
    save(lsImpulseFits,file=file.path(getwd(),"ImpulseDE2_lsImpulseFits.RData"))
    print(paste("Consumed time: ",round(tm_fitImpulse["elapsed"]/60,2),
      " min",sep=""))
    
    ### 5. Detect differentially expressed genes
    print("5. DE analysis")
    dfImpulseResults <- computePval(
      arr2DCountData=arr2DCountData,
      vecDispersions=vecDispersions,
      dfAnnotationFull=dfAnnotationFull,
      lsImpulseFits=lsImpulseFits,
      strCaseName=strCaseName, 
      strControlName=strControlName,
      strMode=strMode, 
      NPARAM=NPARAM)
    
    lsDEGenes <- as.character(as.vector( 
      dfImpulseResults[as.numeric(dfImpulseResults$adj.p) <= Q_value,"Gene"] ))
    save(dfImpulseResults,file=file.path(getwd(),"ImpulseDE2_dfImpulseResults.RData"))
    save(lsDEGenes,file=file.path(getwd(),"ImpulseDE2_lsDEGenes.RData"))
    print(paste("Found ", length(lsDEGenes)," DE genes",sep=""))
    
    ### 6. Plot the top DE genes
    if(boolPlotting){
      print("6. Plot top DE genes")
      tm_plotDEGenes <- system.time({
        plotDEGenes(
          lsGeneIDs=lsDEGenes,
          arr2DCountData=arr2DCountData,
          vecNormConst=vecNormConst,
          dfAnnotationFull=dfAnnotationFull, 
          lsImpulseFits=lsImpulseFits,
          strCaseName=strCaseName, 
          strControlName=strControlName, 
          strFileNameSuffix="DE", 
          strPlotTitleSuffix="", 
          strPlotSubtitle="",
          dfImpulseResults=dfImpulseResults,
          vecRefPval=vecRefPval,
          strMode=strMode,
          strNameMethod2=strRefMethod,
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
    "lsDEGenes"=lsDEGenes,
    "dfImpulseResults"=dfImpulseResults,
    "lsImpulseFits"=lsImpulseFits,
    "dfDESeq2Results"=dfDESeq2Results))
}
