################################################################################
########################     Impulse DE package     ############################
################################################################################

### Version:  1.3
### Date:     2016
### Author v1.0:  Jil Sander
### Author v1.1:  David Fischer --Divide file into src files, annotation
### Author v1.2:  David Fischer --Support replicate samples, WLS fitting and residual scaling
### Author v1.3.1:  David Fischer --WLS fitting with CV weights
### Author v1.3.2:  David Fischer --NB fitting

################################################################################
### Libraries and source code
################################################################################

library(compiler)
library(parallel)
library(DESeq2)

setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building")
#Prepares the data and annotation for internal use
source("srcImpulseDE2_processData.R")
# Wrapper fur running DESeq2
source("srcImpulseDE2_runDESeq2.R")
# Compute value of impulse function given parameters
source("srcImpulseDE2_calcImpulse.R")
# Cost functions for fitting
source("srcImpulseDE2_CostFunctionsFit.R")
# Fit impulse model to a timecourse dataset
source("srcImpulseDE2_fitImpulse.R")
# Detect differentially expressed genes over time
source("srcImpulseDE2_computePval.R")
# Plot the impulse fits and input data
source("srcImpulseDE2_plotDEGenes.R")

################################################################################
### Compiling of frequently used functions
################################################################################

# Compile fitting simple functions to make them quicker
calcImpulse_comp <- cmpfun(calcImpulse)
evalLogLikImpulse_comp <- cmpfun(evalLogLikImpulse)
evalLogLikMean_comp <- cmpfun(evalLogLikMean)

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
#' Here, the impuls fit to the combined data models no differential expression
#' between the conditions.
#' @aliases ImpulseDE2
#' @param matCountData (matrix genes x replicates) Count data of all conditions, 
#' unobserved entries are NA. Column labels are replicate names, row labels
#' gene names.
#' @param dfAnnotationFull (Table) Lists co-variables of individual replicates: 
#' Sample, Condition, Time. Time must be numeric.
#' @param strCaseName (str) Name of the case condition in \code{dfAnnotationFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#' \code{dfAnnotationFull}.
#' @param nProc (scalar) [Default 3] Number of processes for parallelisation. The
#' specified value is internally changed to \code{min(detectCores() - 1, nProc)} 
#' using the \code{detectCores} function from the package \code{parallel} to avoid overload.
#' @param Q_value (scalar) [Default 0.01] FDR-corrected p-value cutoff for significance.
#' @return List containing the following elements:
#' \itemize{
#' \item \code{lsDEGenes} (list number of genes) Genes IDs identified
#' as differentially expressed by ImpulseDE2 at threshold \code{Q_value}.
#' \item \code{dfImpulseResults} (data frame) ImpulseDE2 results.
#' \item \code{lsImpulseFits} 
#' \item \code{dfDESeq2Results} (data frame) DESeq2 results.
#' }
#' Additionally, \code{ImpulseDE2} saves the following objects and tables into
#' the working directory:
#' \itemize{
#' \item \code{ImpulseDE2_arr2DCountData.RData} (2D array genes x replicates) Reduced 
#' version of \code{matCountData}. For internal use.
#' \item \code{ImpulseDE2_arr3DCountData.RData} (3D array genes x samples x replicates)
#' \code{matCountData} reshaped into a 3D array. For internal use.
#' \item \code{ImpulseDE2_dfAnnotationRed.RData} (data frame) Reduced version of 
#' \code{annotation_table}. For internal use.
#' \item \code{ImpulseDE2_vecDESeq2Dispersions.RData} (vector number of genes) Inverse 
#' of gene-wise negative binomial dispersion coefficients computed by DESeq2.
#' \item \code{ImpulseDE2_dfDESeq2Results.RData} (data frame) DESeq2 results.
#' \item \code{ImpulseDE2_lsImpulseFits.RData} (list) 
#' \item \code{ImpulseDE2_dfImpulseResults.RData} (data frame) ImpulseDE2 results.
#' \item \code{ImpulseDE2_lsDEGenes.RData} (list number of genes) Genes IDs identified
#' as differentially expressed by ImpulseDE2 at threshold \code{Q_value}.
#' \item \code{ImpulseDE2_ClusterOut.txt} Text-file with stdout and stderr from
#' cluster created in \code{fitImpulse}.
#' }
#' @details \code{ImpulseDE2} is based on the impulse model proposed by
#' Chechik and Koller, which reflects a two-step behavior of genes within a cell
#' responding to environmental changes (Chechik and Koller, 2009). To detect
#' differentially expressed genes, a five-step workflow is followed:
#' \enumerate{
#' \item The genes are clustered into a limited number of groups using k-means
#' clustering. If \code{plot_clusters} = \code{TRUE}, PDF documents are
#' generated, which contain plots of each cluster. Additionally, a text-file is
#' produced containing each gene together with its corresponding cluster number.
#' \item The impulse model is fitted to the mean expression profiles of the
#' clusters. The best parameter sets are then used for the next step.
#' \item The impulse model is fitted to each gene separately using the parameter
#' sets from step 2 as optimal start point guesses.
#' \item The impulse model is fitted to a randomized dataset (bootstrap), which
#' is essential to detect significantly differentially expressed genes
#' (Storey et al., 2005).
#' \item Detection of differentially expressed genes utilizing the fits to the
#' real and randomized data sets. FDR-correction is performed to obtain adjusted
#' p-values (Benjamini and Hochberg, 1995).
#' }
#' @examples
#' \dontrun{
#' #' Install package longitudinal and load it
#' library(longitudinal)
#' #' Attach datasets
#' data(tcell)
#' #' check dimension of data matrix of interest
#' dim(tcell.10)
#' #' generate a proper annotation table
#' annot <- as.data.frame(cbind("Time" =
#'    sort(rep(get.time.repeats(tcell.10)$time,10)),
#'    "Condition" = "activated"), stringsAsFactors = FALSE)
#' #' Time columns must be numeric
#' annot$Time <- as.numeric(annot$Time)
#' #' rownames of annotation table must appear in data table
#' rownames(annot) = rownames(tcell.10)
#' #' apply ImpulseDE2 in single time course mode
#' #' since genes must be in rows, transpose data matrix using t()
#' #' For the example, reduce random iterations to 100 and number of
#' #' used processors to 1
#' impulse_results <- impulse_DE(t(tcell.10), annot, "Time", "Condition",
#'    n_randoms = 50, nProc = 1)
#' }
#' @seealso \code{\link{processData}}, \code{\link{runDESeq2}},
#' \code{\link{fitImpulse}}, \code{\link{evalLogLikImpulse}}, 
#' \code{\link{evalLogLikMean}}, \code{\link{calcImpulse}},
#' \code{\link{computePval}}, \code{\link{plotDEGenes}}.
#' @author David Fischer
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

# INPUT
#   matCountData: (matrix genes x replicates) Count data of all conditions, 
#       unobserved entries are NA. Column labels are replicate names, row labels
#       gene names.
#   dfAnnotationFull: (Table) Lists co-variables of individual replicates:
#       Sample, Condition, Time. Time must be numeric.
#   strCaseName: (str) Name of the case condition in dfAnnotationFull
#   strControlName: (str) Name of the control condition in dfAnnotationFull
#   nProc: (scalar) [Default 3] Number of processes for parallelisation.
#   Q_value: (scalar) [default 0.01] FDR-corrected p-value cutoff for significance.
# OUTPUT
#   lsImpulseFits: (list) List containing fitted values and model parameters
#   lsDEGenes: names of genes being called as differentially expressed
#   ...more output saved to working directory, see documentation

runImpulseDE2 <- function(matCountData=NULL, dfAnnotationFull=NULL,
  strCaseName = NULL, strControlName=NULL, nProc=3, Q_value=0.01){
  
  print("###################################################################")
  print("Impulse v1.3 for count data")
  print("###################################################################")
  
  NPARAM=6
  
  tm_runImpulseDE2 <- system.time({
    
    # 1. Process input data 
    print("1. Prepare data:")
    tm_processData <- system.time({
      lsProcessedData <- processData(
        dfAnnotationFull=dfAnnotationFull,matCountData=matCountData,
        strControlName=strControlName, strCaseName=strCaseName)
    })
    arr2DCountData <- lsProcessedData[[1]]
    arr3DCountData <- lsProcessedData[[2]]
    dfAnnotationRed <- lsProcessedData[[3]]
    save(arr2DCountData,file=file.path(getwd(),"ImpulseDE2_arr2DCountData.RData"))
    save(arr3DCountData,file=file.path(getwd(),"ImpulseDE2_arr3DCountData.RData"))
    save(dfAnnotationRed,file=file.path(getwd(),"ImpulseDE2_dfAnnotationRed.RData"))
    print("DONE")
    print(paste("Consumed time: ",round(tm_processData["elapsed"]/60,2),
      " min",sep=""))
    print("###################################################################")
    
    # 2. Run DESeq2
    print("2. Run DESeq2:")
    tm_runDESeq2 <- system.time({
      lsDESeq2Results <- runDESeq2(dfAnnotationFull=dfAnnotationFull,arr2DCountData=arr2DCountData)
    })
    vecDESeq2Dispersions <- lsDESeq2Results[[1]]
    dfDESeq2Results <- lsDESeq2Results[[2]]
    save(vecDESeq2Dispersions,file=file.path(getwd(),"ImpulseDE2_vecDESeq2Dispersions.RData"))
    save(dfDESeq2Results,file=file.path(getwd(),"ImpulseDE2_dfDESeq2Results.RData"))
    print("DONE")
    print(paste("Consumed time: ",round(tm_runDESeq2["elapsed"]/60,2),
      " min",sep=""))
    print("###################################################################")
    
    ###  3. Fit Impule model to each gene 
    print("3. Fitting Impulse model to the genes")
    tm_fitImpulse <- system.time({
      lsImpulseFits <- fitImpulse(arrCountData=arr3DCountData, 
        vecDispersions=vecDESeq2Dispersions, dfAnnotation=dfAnnotationRed, 
        strCaseName=strCaseName, strControlName=strControlName, 
        nProcessesAssigned=nProc, NPARAM=NPARAM)
    })
    save(lsImpulseFits,file=file.path(getwd(),"ImpulseDE2_lsImpulseFits.RData"))
    print("DONE")
    print(paste("Consumed time: ",round(tm_fitImpulse["elapsed"]/60,2),
      " min",sep=""))
    print("###################################################################")
    
    ### 4. Detect differentially expressed genes
    print("4. DE analysis")
    tm_DE <- system.time({
      dfImpulseResults <- computePval(
        arr3DCountData=arr3DCountData,vecDispersions=vecDESeq2Dispersions,
        dfAnnotationRed=dfAnnotationRed,lsImpulseFits=lsImpulseFits,
        strCaseName=strCaseName, strControlName=strControlName,
        NPARAM=NPARAM)
      
      lsDEGenes <- as.character(as.vector( 
        dfImpulseResults[as.numeric(dfImpulseResults$adj.p) <= Q_value,"Gene"] ))
    })
    save(dfImpulseResults,file=file.path(getwd(),"ImpulseDE2_dfImpulseResults.RData"))
    save(lsDEGenes,file=file.path(getwd(),"ImpulseDE2_lsDEGenes.RData"))
    print(paste("Found ", length(lsDEGenes)," DE genes",sep=""))
    print("DONE")
    print(paste("Consumed time: ",round(tm_DE["elapsed"]/60,2),
      " min",sep=""))
    print("###################################################################")
    
    ### 5. Plot the top DE genes
    print("5. Plot top DE genes")
    tm_plotDEGenes <- system.time({
      plotDEGenes(lsGeneIDs=lsDEGenes,
        arr3DCountData=arr3DCountData, dfAnnotationRed=dfAnnotationRed, 
        lsImpulseFits=lsImpulseFits,
        strCaseName=strCaseName, strControlName=strControlName, 
        strFileNameSuffix="DE", strPlotTitleSuffix="", strPlotSubtitle="",
        dfImpulseResults=dfImpulseResults,dfDESeq2Results=dfDESeq2Results,
        NPARAM=NPARAM)
    })
    print("DONE")
    print(paste("Consumed time: ",round(tm_plotDEGenes["elapsed"]/60,2),
      " min",sep=""))
    print("##################################################################")
  })
  print("Finished ImpulseDE2.")
  print(paste("TOTAL consumed time: ",round(tm_runImpulseDE2["elapsed"]/60,2),
    " min",sep=""))
  print("##################################################################")
  
  return(list(
    "lsDEGenes"=lsDEGenes,
    "dfImpulseResults"=dfImpulseResults,
    "lsImpulseFits"=lsImpulseFits,
    "dfDESeq2Results"=dfDESeq2Results))
}
