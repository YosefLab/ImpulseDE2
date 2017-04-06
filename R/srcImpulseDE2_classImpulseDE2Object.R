### ImpulseDE2 output container class

### 1. Define output container class

# Define class unions for slots DO NOT USE SPACE BETWEEN members=c(..)
# THAT CAUSES ERROR DURING R CMD BUILD

#' @importFrom methods setClassUnion
setClassUnion("numericORNULL", members = c("numeric", "NULL"))
setClassUnion("characterORNULL", members = c("character", "NULL"))
setClassUnion("listORNULL", members = c("list", "NULL"))
setClassUnion("data.frameORNULL", members = c("data.frame", "NULL"))

#' Container class for ImpulseDE2 output
#' 
#' ImpulseDE2 output and intermediate results such as model fits.
#' 
#' @slot dfDEAnalysis (data frame samples x reported characteristics) 
#' Summary of fitting procedure and 
#' differential expression results for each gene.
#' \itemize{
#' \item Gene: Gene ID.
#' \item p: P-value for differential expression.
#' \item padj: Benjamini-Hochberg false-discovery rate corrected p-value
#' for differential expression analysis.
#' \item loglik_full: Loglikelihood of full model.
#' \item loglik_red: Loglikelihood of reduced model.
#' \item df_full: Degrees of freedom of full model.
#' \item df_red: Degrees of freedom of reduced model
#' \item mean: Inferred mean parameter of constant model of first batch.
#' From combined samples in case-ctrl. 
#' \item allZero (bool) Whether there were no observed 
#' non-zero observations of this gene.
#' If TRUE, fitting and DE analsysis were skipped and entry is NA.
#' }
#' Entries only present in case-only DE analysis:
#' \itemize{
#' \item converge_impulse: Convergence status of optim for 
#' impulse model fit (full model).
#' \item converge_const: Convergence status of optim for 
#' constant model fit (reduced model).
#' }
#' Entries only present in case-control DE analysis:
#' \itemize{
#' \item converge_combined: Convergence status of optim for 
#' impulse model fit to case and control samples combined (reduced model).
#' \item converge_case: Convergence status of optim for 
#' impulse model fit to samples of case condition (full model 1/2).
#' \item converge_control: Convergence status of optim for 
#' impulse model fit to samples of control condition (full model 2/2).
#' }
#' Entries only present if boolIdentifyTransients is TRUE:
#' \itemize{
#' \item converge_sigmoid: Convergence status of optim for 
#' sigmoid model fit to samples of case condition.
#' \item impulseTOsigmoid_p: P-value of loglikelihood ratio test
#' impulse model fit versus sigmoidal model on samples of case condition.
#' \item impulseTOsigmoid_padj: Benjamini-Hochberg 
#' false-discovery rate corrected p-value of loglikelihood ratio test
#' impulse model fit versus sigmoid model on samples of case condition.
#' \item sigmoidTOconst_p: P-value of loglikelihood ratio test
#' sigmoidal model fit versus constant model on samples of case condition.
#' \item sigmoidTOconst_padj: Benjamini-Hochberg 
#' false-discovery rate corrected p-value of loglikelihood ratio test
#' sigmoidal model fit versus constant model on samples of case condition.
#' \item isTransient (bool) Whether gene is transiently
#' activated or deactivated and differentially expressed.
#' \item isMonotonous (bool) Whether gene is not transiently
#' activated or deactivated and differentially expressed. This scenario
#' corresponds to a montonous expression level increase or decrease.
#' }
#' @slot vecDEGenes (list number of genes) Genes IDs identified
#' as differentially expressed by ImpulseDE2 at threshold \code{scaQThres}.
#' @slot lsModelFits (list length number of conditions fit (1 or 3))
#' {'case'} or {'case', 'control', 'combined'}
#' One model fitting object for each condition:
#' In case-only DE analysis, only the condition {'case'} is fit.
#' In case-control DE analysis, the conditions 
#' {'case', 'control','combined} are fit.
#' Each condition entry is a list of model fits for each gene.
#' Each gene entry is a list of model fits to the individual models:
#' Impulse model and constant model (if boolFitConst is TRUE).
#' At this level, the sigmoid model fit can be added later.
#' Each model fit per gene is a list of fitting parameters and results.
#' \itemize{
#' \item IdxGroups (list length number of conditions)
#' Samples grouped by time points and by batches and time point vectors. 
#' Sample groups are stored in the form of index vectors in which
#' samples of the same time point or batch have the same index.
#' \itemize{
#' \item Condition ID (list length 3)
#' List of index vectors and time points.
#' One entry of this format for each condition.
#' \itemize{
#' \item vecTimepointsUnique (numeric vector length number of unique
#' timepoints) Vector of unique time coordinates observed in this condition.
#' \item vecidxTimepoint (idx vector length number of samples)
#' Index of the time coordinates of each sample (reference is
#' vecTimepointsUnique).
#' \item lsvecBatchUnique (list number of confounders)
#' List of string vectors. One vector per confounder: vector of unique batches
#' in this confounder.
#' \item lsvecidxBatches (idx list length number of confounding variables)
#' List of index vectors. 
#' One vector per confounding variable.
#' Each vector has one entry per sample with the index of the batch ID
#' within the given confounding variable of the given sample. Reference
#' is the list of unique batch ids for each confounding variable.
#' }
#' }
#' \item Condition ID (list length number of genes)
#' List of fits for each gene to the samples of this condition.
#' One entry of this format for all conditions fit.
#' \itemize{
#' \item Gene ID (list length 2)
#' Impulse and constant model fit to gene observations.
#' One entry of this format for all gene IDs.
#' \itemize{
#' \item lsImpulseFit (list) List of impulse fit parameters and results.
#' \itemize{
#' \item vecImpulseParam (numeric vector length 6)
#' {beta, h0, h1, h2, t1, t2}
#' Maximum likelihood estimators of impulse model parameters.
#' \item vecImpulseValue (numeric vector length number of time points)
#' Values of impulse model fit at time points used for fit.
#' \item lsvecBatchFactors (list length number of confounders)
#' List of vectors of scalar batch correction factors for each sample.
#' These are also maximum likelihood estimators.
#' NULL if no confounders given.
#' \item scaDispParam (scalar) Dispersion parameter estimate
#' used in fitting (hyper-parameter).
#' \item scaLL (scalar) Loglikelihood of data under maximum likelihood
#' estimator model.
#' \item scaConvergence (scalar) 
#' Convergence status of optim on impulse model.
#' }
#' \item lsConstFit (list) List of constant fit parameters and results.
#' \itemize{
#' \item scaMu (scalar) Maximum likelihood estimator of
#' negative binomial mean parameter.
#' \item lsvecBatchFactors (list length number of confounders)
#' List of vectors of scalar batch correction factors for each sample.
#' These are also maximum likelihood estimators.
#' NULL if no confounders given.
#' \item scaDispParam (scalar) Dispersion parameter estimate
#' used in fitting (hyper-parameter).
#' \item scaLL (scalar) Loglikelihood of data under maximum likelihood
#' estimator model.
#' \item scaConvergence (scalar) 
#' Convergence status of optim on constant model.
#' }
#' \item ls SigmoidFit (list) List of sigmoidal fit parameters and results.
#' NULL if boolIdentifyTransients is FALSE.
#' \itemize{
#' \item vecSigmoidParam (numeric vector length 4)
#' {beta, h0, h1, t}
#' Maximum likelihood estimators of sigmoidal model parameters.
#' \item vecSigmoidValue (numeric vector length number of time points)
#' Values of sigmoid model fit at time points used for fit.
#' \item lsvecBatchFactors (list length number of confounders)
#' List of vectors of scalar batch correction factors for each sample.
#' These are also maximum likelihood estimators.
#' NULL if no confounders given.
#' \item scaDispParam (scalar) Dispersion parameter estimate
#' used in fitting (hyper-parameter).
#' \item scaLL (scalar) Loglikelihood of data under maximum likelihood
#' estimator model.
#' \item scaConvergence (scalar) 
#' Convergence status of optim on sigmoidal model.
#' }
#' }
#' }
#' }
#' @slot matCountDataProc (matrix genes x samples) [Default NULL] 
#' Read count data, unobserved entries are NA. Processed matrix.
#' @slot dfAnnotationProc (data frame samples x covariates) 
#' {Sample, Condition, Time (numeric), TimeCateg (str)
#' (and confounding variables if given).}
#' Annotation table with covariates for each sample.
#' Processed table.
#' @slot vecDispersions (numeric vector number of samples) 
#' Gene-wise negative binomial dispersion hyper-parameters.
#' @slot vecSizeFactors (numeric vector number of samples) 
#' Model scaling factors for each sample which take
#' sequencing depth into account (size factors).
#' @slot boolCaseCtrl (bool) 
#' Whether to perform case-control analysis. Does case-only
#' analysis if FALSE.
#' @slot vecConfounders (vector of strings number of confounding variables)
#' Factors to correct for during batch correction. Have to 
#' supply dispersion factors if more than one is supplied.
#' Names refer to columns in dfAnnotation.
#' @slot scaNProc (scalar) Number of processes for 
#' parallelisation.
#' @slot scaQThres (scalar)
#' FDR-corrected p-value cutoff for significance.
#' @slot strReport (str)
#' ImpulseDE2 stdout report.
#' 
#' @name ImpulseDE2Object-class
#' 
#' @author David Sebastian Fischer
setClass("ImpulseDE2Object", slots = c(
    dfImpulseDE2Results = "data.frameORNULL", 
    vecDEGenes = "characterORNULL", lsModelFits = "listORNULL", matCountDataProc = "matrix", 
    vecAllIDs = "characterORNULL", dfAnnotationProc = "data.frame", vecSizeFactors = "numeric", 
    vecDispersions = "numeric", boolCaseCtrl = "logical", vecConfounders = "characterORNULL", 
    scaNProc = "numeric", scaQThres = "numericORNULL", strReport = "characterORNULL"))

### 2. Enable accession of private elements via functions

#' ImpulseDE2Object "get" accession functions
#' 
#' Get internal data of ImpulseDE2 output object.
#' 
#' @param obj (ImpulseDE2Object)  A ImpulseDE2 output object.
#' 
#' @return The internal data object specified by the function.
#' 
#' @name get_accessors
#' @rdname get_accessors
#' @aliases 
#' get_matCountDataProc 
#' get_dfAnnotationProc 
#' get_vecSizeFactors
#' get_vecDispersions 
#' get_boolCaseCtrl
#' get_vecConfounders 
#' get_scaNProc
#' get_scaQThres
#' get_strReport
#'
#' @examples    
#' lsSimulatedData <- simulateDataSetImpulseDE2(
#' vecTimePointsA   = rep(seq(1,8),3),
#' vecTimePointsB   = NULL,
#' vecBatchesA      = NULL,
#' vecBatchesB      = NULL,
#' scaNConst        = 30,
#' scaNImp          = 10,
#' scaNLin          = 10,
#' scaNSig          = 10)
#' objectImpulseDE2 <- runImpulseDE2(
#' matCountData    = lsSimulatedData$matObservedCounts, 
#' dfAnnotation    = lsSimulatedData$dfAnnotation,
#' boolCaseCtrl    = FALSE,
#' vecConfounders  = NULL,
#' scaNProc        = 1 )
#' # Extract hidden auxillary result and processed input objects.
#' lsModelFits <- get_lsModelFits(objectImpulseDE2)
#' matCountDataProc <- get_matCountDataProc(objectImpulseDE2)
#' dfAnnotationProc <- get_dfAnnotationProc(objectImpulseDE2)
#' vecAllIDs <- get_vecAllIDs(objectImpulseDE2)
#' vecSizeFactors <- get_vecSizeFactors(objectImpulseDE2)
#' vecDispersions <- get_vecDispersions(objectImpulseDE2)
#' boolCaseCtrl <- get_boolCaseCtrl(objectImpulseDE2)
#' vecConfounders <- get_vecConfounders(objectImpulseDE2)
#' scaNProc <- get_scaNProc(objectImpulseDE2)
#' scaQThres <- get_scaQThres(objectImpulseDE2)
#' strReport <- get_strReport(objectImpulseDE2)
#' 
#' @author David Sebastian Fischer
NULL

#' @rdname get_accessors
#' @export
get_lsModelFits <- function(obj) 
    return(obj@lsModelFits)

#' @rdname get_accessors
#' @export
get_matCountDataProc <- function(obj) 
    return(obj@matCountDataProc)

#' @rdname get_accessors
#' @export
get_dfAnnotationProc <- function(obj) 
    return(obj@dfAnnotationProc)

#' @rdname get_accessors
#' @export
get_vecAllIDs <- function(obj) 
    return(obj@vecAllIDs)

#' @rdname get_accessors
#' @export
get_vecSizeFactors <- function(obj) 
    return(obj@vecSizeFactors)

#' @rdname get_accessors
#' @export
get_vecDispersions <- function(obj) 
    return(obj@vecDispersions)

#' @rdname get_accessors
#' @export
get_boolCaseCtrl <- function(obj) 
    return(obj@boolCaseCtrl)

#' @rdname get_accessors
#' @export
get_vecConfounders <- function(obj) 
    return(obj@vecConfounders)

#' @rdname get_accessors
#' @export
get_scaNProc <- function(obj) 
    return(obj@scaNProc)

#' @rdname get_accessors
#' @export
get_scaQThres <- function(obj) 
    return(obj@scaQThres)

#' @rdname get_accessors
#' @export
get_strReport <- function(obj) 
    return(obj@strReport)

#' ImpulseDE2Object "set" accession functions
#' 
#' Set internal data of ImpulseDE2 output object.
#' 
#' @param obj (ImpulseDE2Object)  A ImpulseDE2 output object.
#' @param element (type depends on element)  
#' An element of the ImpulseDE2 output object which is to be substituted.
#' 
#' @return (ImpulseDE2Object)
#' The ImpulseDE2 object with the target element substituted.
#' 
#' @aliases 
#' set_boolCaseCtrl
#' set_dfAnnotationProc
#' set_dfImpulseDE2Results
#' set_lsModelFits
#' set_matCountDataProc
#' set_vecAllIDs
#' set_vecConfounders
#' set_vecDEGenes
#' set_vecDispersions
#' set_vecSizeFactors 
#' set_scaNProc
#' set_scaQThres
#' set_strReport
#' 
#' @name set_accessors
#' 
#' @author David Sebastian Fischer
NULL


#' @name set_accessors
set_boolCaseCtrl <- function(obj,element) {
    obj@boolCaseCtrl <- element
    return(obj)
}

#' @name set_accessors
set_dfAnnotationProc <- function(obj,element) {
    obj@dfAnnotationProc <- element
    return(obj)
}

#' @name set_accessors
set_dfImpulseDE2Results <- function(obj,element) {
    obj@dfImpulseDE2Results <- element
    return(obj)
}

#' @name set_accessors
set_lsModelFits <- function(obj,element) {
    obj@lsModelFits <- element
    return(obj)
}

#' @name set_accessors
set_matCountDataProc <- function(obj,element) {
    obj@matCountDataProc <- element
    return(obj)
}

#' @name set_accessors
set_scaNProc <- function(obj,element) {
    obj@scaNProc <- element
    return(obj)
}

#' @name set_accessors
set_scaQThres <- function(obj,element) {
    obj@scaQThres <- element
    return(obj)
}

#' @name set_accessors
set_strReport <- function(obj,element) {
    obj@strReport <- element
    return(obj)
}

#' @name set_accessors
set_vecAllIDs <- function(obj,element) {
    obj@vecAllIDs <- element
    return(obj)
}

#' @name set_accessors
set_vecConfounders <- function(obj,element) {
    obj@vecConfounders <- element
    return(obj)
}

#' @name set_accessors
set_vecDEGenes <- function(obj,element) {
    obj@vecDEGenes <- element
    return(obj)
}

#' @name set_accessors
set_vecDispersions <- function(obj,element) {
    obj@vecDispersions <- element
    return(obj)
}

#' @name set_accessors
set_vecSizeFactors <- function(obj,element) {
    obj@vecSizeFactors <- element
    return(obj)
}

### 3. Enable accession of public elements via list-like properties of
### ImpulseDE2Object.

#' List-like accessor methods for ImpulseDE2Object
#' 
#' Allow usage of ImpulseDE2 ouput object like a list with
#' respect to the core output:
#' dfImpulseDE2Results and vecDEGenes.
#' 
#' @param x (ImpulseDE2Object) ImpulseDE2 output object.
#' @param i,name (idx or str) Name or index of core output 
#' element of ImpulseDE2Object.
#' @param j       Not used, only vectors.
#' @param ...     Not used.
#' 
#' @examples
#' lsSimulatedData <- simulateDataSetImpulseDE2(
#' vecTimePointsA   = rep(seq(1,8),3),
#' vecTimePointsB   = NULL,
#' vecBatchesA      = NULL,
#' vecBatchesB      = NULL,
#' scaNConst        = 30,
#' scaNImp          = 10,
#' scaNLin          = 10,
#' scaNSig          = 10)
#' objectImpulseDE2 <- runImpulseDE2(
#' matCountData    = lsSimulatedData$matObservedCounts, 
#' dfAnnotation    = lsSimulatedData$dfAnnotation,
#' boolCaseCtrl    = FALSE,
#' vecConfounders  = NULL,
#' scaNProc        = 1 )
#' names(objectImpulseDE2) # Display core output
#' # With respect to this core output, objectImpulseDE2
#' # can be treated like a list.
#' head(objectImpulseDE2[['dfImpulseDE2Results']])
#' head(objectImpulseDE2$dfImpulseDE2Results)
#' head(objectImpulseDE2[['vecDEGenes']])
#' head(objectImpulseDE2$vecDEGenes)
#' 
#' @name list_accession
#' @aliases names.ImpulseDE2Object 
#' names,ImpulseDE2Object-method          
#' $,ImpulseDE2Object-method          
#' [[,ImpulseDE2Object,character,missing-method
#' 
#' @author David Sebastian Fischer
NULL

#' @return Names of core output in ImpulseDE2Object.
#' @name list_accession
#' @export
setMethod("names", "ImpulseDE2Object", function(x) {
    return(c("dfImpulseDE2Results", "vecDEGenes"))
})

#' @return Target element from ImpulseDE2Object.
#' @name list_accession
#' @export
setMethod("[[", c("ImpulseDE2Object", "character", "missing"), 
          function(x, i, j, ...) {
              if (identical(i, "dfImpulseDE2Results")) {
                  return(x@dfImpulseDE2Results)
              } else if (identical(i, "vecDEGenes")) {
                  return(x@vecDEGenes)
              } else {
                  return(NULL)
              }
          })

#' @return Target element from ImpulseDE2Object.
#' @name list_accession
#' @export
setMethod("$", "ImpulseDE2Object", function(x, name) x[[name]])

### 4. Functions on ImpulseDE2Object that perform specific tasks

#' Append string to strReport in ImpulseDE2Object
#' 
#' Append string to strReport in ImpulseDE2Object.
#' 
#' @param obj (ImpulseDE2Object)  A ImpulseDE2 output object.
#' @param s (str) String to be appended.
#' 
#' @return (obj)
#' The ImpulseDE2 object with the str appended in strReport.
#' 
#' @author David Sebastian Fischer
append_strReport <- function(obj,s) {
    obj@strReport <- paste0(obj@strReport,"\n", s)
    return(obj)
}

#' Print ImpulseDE2 report to .txt file
#' 
#' Print ImpulseDE2 report to .txt file.
#'
#' @param object (ImpulseDE2Object) Output object of ImpulseDE2.
#' @param fileReport (file) File to print report to.
#' 
#' @return No return.
#' 
#' @examples
#' dirPWD <- getwd() # Will save into current working directory.
#' lsSimulatedData <- simulateDataSetImpulseDE2(
#' vecTimePointsA   = rep(seq(1,8),3),
#' vecTimePointsB   = NULL,
#' vecBatchesA      = NULL,
#' vecBatchesB      = NULL,
#' scaNConst        = 30,
#' scaNImp          = 10,
#' scaNLin          = 10,
#' scaNSig          = 10)
#' objectImpulseDE2 <- runImpulseDE2(
#' matCountData    = lsSimulatedData$matObservedCounts, 
#' dfAnnotation    = lsSimulatedData$dfAnnotation,
#' boolCaseCtrl    = FALSE,
#' vecConfounders  = NULL,
#' scaNProc        = 1 )
#' # Uncomment to run:
#' #writeReportToFile(
#' #object=objectImpulseDE2,
#' #file=paste0(dirPWD, 'ImpulseDE2Report.txt')
#' #)
#' 
#' @author David Sebastian Fischer
#' @export
writeReportToFile <- function(object, fileReport) {
    write(object@strReport, file = fileReport, ncolumns = 1)
}
