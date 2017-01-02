################################################################################
#################     ImpulseDE2 output container class     ####################
################################################################################

### 1. Define output container class

#' Container class for ImpulseDE2 output
#' 
#' ImpulseDE2 output and intermediate results such as model fits.
#' 
#' @slot dfDEAnalysis (data frame samples x reported characteristics) 
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
#' @slot vecDEGenes: (list number of genes) Genes IDs identified
#'    as differentially expressed by ImpulseDE2 at threshold \code{scaQThres}.
#' @slot lsModelFits: (list length number of conditions fit (1 or 3))
#'    {"case"} or {"case", "control", "combined"}
#'    One model fitting object for each condition:
#'    In case-only DE analysis, only the condition {"case"} is fit.
#'    In case-control DE analysis, the conditions 
#'    {"case", "control","combined} are fit.
#'    Each condition entry is a list of model fits for each gene.
#'       Each gene entry is a list of model fits to the individual models:
#'     Impulse model and constant model (if boolFitConst is TRUE).
#'    At this level, the sigmoid model fit can be added later.
#'    Each model fit per gene is a list of fitting parameters and results.
#'    \itemize{
#'       \item IdxGroups: (list length number of conditions)
#'       Samples grouped by time points and by batches and time point vectors. 
#'       Sample groups are stored in the form of index vectors in which
#'       samples of the same time point or batch have the same index.
#'       \itemize{
#'         \item Condition ID: (list length 3)
#'         List of index vectors and time points.
#'         One entry of this format for each condition.
#'         \itemize{
#'           \item vecTimepointsUnique: (numeric vector length number of unique
#'         timepoints) Vector of unique time coordinates observed in this condition.
#'           \item vecidxTimepoint: (idx vector length number of samples)
#'         Index of the time coordinates of each sample (reference is
#'         vecTimepointsUnique).
#'           \item lsvecBatchUnique: (list number of confounders)
#'         List of string vectors. One vector per confounder: vector of unique batches
#'         in this confounder.
#'           \item lsvecidxBatches: (idx list length number of confounding variables)
#'    		   List of index vectors. 
#'    		   One vector per confounding variable.
#'    		   Each vector has one entry per sample with the index of the batch ID
#'    		   within the given confounding variable of the given sample. Reference
#'    		   is the list of unique batch ids for each confounding variable.
#'         }
#'       }
#'       \item Condition ID: (list length number of genes)
#'       List of fits for each gene to the samples of this condition.
#'       One entry of this format for all conditions fit.
#'       \itemize{
#'         \item Gene ID: (list length 2)
#'         Impulse and constant model fit to gene observations.
#'         One entry of this format for all gene IDs.
#'         \itemize{
#'           \item lsImpulseFit: (list) List of impulse fit parameters and results.
#'           \itemize{
#'             \item vecImpulseParam: (numeric vector length 6)
#'           {beta, h0, h1, h2, t1, t2}
#'           Maximum likelihood estimators of impulse model parameters.
#'             \item vecImpulseValue: (numeric vector length number of time points)
#'           Values of impulse model fit at time points used for fit.
#'             \item lsvecBatchFactors: (list length number of confounders)
#'           List of vectors of scalar batch correction factors for each sample.
#'           These are also maximum likelihood estimators.
#'           NULL if no confounders given.
#'             \item scaDispParam: (scalar) Dispersion parameter estimate
#'           used in fitting (hyper-parameter).
#'             \item scaLL: (scalar) Loglikelihood of data under maximum likelihood
#'           estimator model.
#'             \item scaConvergence: (scalar) 
#'           Convergence status of optim on impulse model.
#'           }
#'           \item lsConstFit: (list) List of constant fit parameters and results.
#'           \itemize{
#'             \item scaMu: (scalar) Maximum likelihood estimator of
#'           negative binomial mean parameter.
#'             \item lsvecBatchFactors: (list length number of confounders)
#'           List of vectors of scalar batch correction factors for each sample.
#'           These are also maximum likelihood estimators.
#'           NULL if no confounders given.
#'             \item scaDispParam: (scalar) Dispersion parameter estimate
#'           used in fitting (hyper-parameter).
#'             \item scaLL: (scalar) Loglikelihood of data under maximum likelihood
#'           estimator model.
#'             \item scaConvergence: (scalar) 
#'           Convergence status of optim on constant model.
#'           }
#'           \item ls SigmoidFit: (list) List of sigmoidal fit parameters and results.
#'           NULL if boolIdentifyTransients is FALSE.
#'           \itemize{
#'             \item vecSigmoidParam: (numeric vector length 4)
#'           {beta, h0, h1, t}
#'           Maximum likelihood estimators of sigmoidal model parameters.
#'             \item vecSigmoidValue: (numeric vector length number of time points)
#'           Values of sigmoid model fit at time points used for fit.
#'             \item lsvecBatchFactors: (list length number of confounders)
#'           List of vectors of scalar batch correction factors for each sample.
#'           These are also maximum likelihood estimators.
#'           NULL if no confounders given.
#'             \item scaDispParam: (scalar) Dispersion parameter estimate
#'           used in fitting (hyper-parameter).
#'             \item scaLL: (scalar) Loglikelihood of data under maximum likelihood
#'           estimator model.
#'             \item scaConvergence: (scalar) 
#'           Convergence status of optim on sigmoidal model.
#'           }
#'         }
#'       }
#'     }
#' @slot matCountDataProc: (matrix genes x samples) [Default NULL] 
#'    Read count data, unobserved entries are NA. Processed matrix.
#' @slot dfAnnotationProc: (data frame samples x covariates) 
#'    {Sample, Condition, Time (numeric), TimeCateg (str)
#'    (and confounding variables if given).}
#'    Annotation table with covariates for each sample.
#'    Processed table.
#' @slot vecDispersions: (numeric vector number of samples) 
#'    Gene-wise negative binomial dispersion hyper-parameters.
#' @slot vecSizeFactors: (numeric vector number of samples) 
#'    Model scaling factors for each sample which take
#'    sequencing depth into account (size factors).
#' @slot boolCaseCtrl: (bool) 
#' 		Whether to perform case-control analysis. Does case-only
#' 		analysis if FALSE.
#' @slot vecConfounders: (vector of strings number of confounding variables)
#' 		Factors to correct for during batch correction. Have to 
#' 		supply dispersion factors if more than one is supplied.
#' 		Names refer to columns in dfAnnotation.
#' @slot scaNProc: (scalar) Number of processes for 
#'    parallelisation.
#' @slot scaQThres: (scalar)
#'    FDR-corrected p-value cutoff for significance.
#' @slot strReport: (str)
#'    ImpulseDE2 stdout report.
#'      
#' @author David Sebastian Fischer
setClass(
  'ImpulseDE2Object',
  slots = c(
    dfImpulseDE2Results = "data.frameORNULL",
    vecDEGenes          = "characterORNULL",
    lsModelFits         = "listORNULL",
    matCountDataProc    = "matrix",
    dfAnnotationProc    = "data.frame",
    vecSizeFactors      = "numeric",
    vecDispersions      = "numeric",
    boolCaseCtrl        = "logical",
    vecConfounders      = "characterORNULL",
    scaNProc            = "numeric", 
    scaQThres           = "numericORNULL",
    strReport           = "characterORNULL"
  )
)

### 2. Enable accession of private elements via functions
### which carry the same name as the element.
### I. Set generic function which defines string as a function:
### setGeneric('funName', function(object) standardGeneric('funName'), valueClass = 'funOutputClass')
### II. Define function on ImpulseDE2Object:
### setMethod('funName', 'ImpulseDE2Object', function(object) object@funName)

#' @return \code{lsModelFits} retrieves lsModelFits
#' @name ImpulseDEObject generics
#' @export
setGeneric('lsModelFits', function(object) standardGeneric('lsModelFits'), valueClass = 'listORNULL')
#' @name ImpulseDEObject accessors
#' @export
setMethod('lsModelFits', 'ImpulseDE2Object', function(object) object@lsModelFits)

#' @return \code{matCountDataProc} retrieves matCountDataProc
#' @name ImpulseDEObject generics
#' @export
setGeneric('matCountDataProc', function(object) standardGeneric('matCountDataProc'), valueClass = 'matrix')
#' @name ImpulseDEObject accessors
#' @export
setMethod('matCountDataProc', 'ImpulseDE2Object', function(object) object@matCountDataProc)

#' @return \code{dfAnnotationProc} retrieves dfAnnotationProc
#' @name ImpulseDEObject generics
#' @export
setGeneric('dfAnnotationProc', function(object) standardGeneric('dfAnnotationProc'), valueClass = 'data.frame')
#' @name ImpulseDEObject accessors
#' @export
setMethod('dfAnnotationProc', 'ImpulseDE2Object', function(object) object@dfAnnotationProc)

#' @return \code{vecSizeFactors} retrieves vecSizeFactors
#' @name ImpulseDEObject generics
#' @export
setGeneric('vecSizeFactors', function(object) standardGeneric('vecSizeFactors'), valueClass = 'numeric')
#' @name ImpulseDEObject accessors
#' @export
setMethod('vecSizeFactors', 'ImpulseDE2Object', function(object) object@vecSizeFactors)

#' @return \code{vecDispersions} retrieves vecDispersions
#' @name ImpulseDEObject generics
#' @export
setGeneric('vecDispersions', function(object) standardGeneric('vecDispersions'), valueClass = 'numeric')
#' @name ImpulseDEObject accessors
#' @export
setMethod('vecDispersions', 'ImpulseDE2Object', function(object) object@vecDispersions)

#' @return \code{boolCaseCtrl} retrieves boolCaseCtrl
#' @name ImpulseDEObject generics
#' @export
setGeneric('boolCaseCtrl', function(object) standardGeneric('boolCaseCtrl'), valueClass = 'logical')
#' @name ImpulseDEObject accessors
#' @export
setMethod('boolCaseCtrl', 'ImpulseDE2Object', function(object) object@boolCaseCtrl)

#' @return \code{vecConfounders} retrieves vecConfounders
#' @name ImpulseDEObject generics
#' @export
setGeneric('vecConfounders', function(object) standardGeneric('vecConfounders'), valueClass = 'characterORNULL')
#' @name ImpulseDEObject accessors
#' @export
setMethod('vecConfounders', 'ImpulseDE2Object', function(object) object@vecConfounders)

#' @return \code{scaNProc} retrieves scaNProc
#' @name ImpulseDEObject generics
#' @export
setGeneric('scaNProc', function(object) standardGeneric('scaNProc'), valueClass = 'numeric')
#' @name ImpulseDEObject accessors
#' @export
setMethod('scaNProc', 'ImpulseDE2Object', function(object) object@scaNProc)

#' @return \code{scaQThres} retrieves scaQThres
#' @name ImpulseDEObject generics
#' @export
setGeneric('scaQThres', function(object) standardGeneric('scaQThres'), valueClass = 'numericOrNULL')
#' @name ImpulseDEObject accessors
#' @export
setMethod('scaQThres', 'ImpulseDE2Object', function(object) object@scaQThres)

#' @return \code{strReport} retrieves strReport
#' @name ImpulseDEObject generics
#' @export
setGeneric('strReport', function(object) standardGeneric('strReport'), valueClass = 'characterORNULL')
#' @name ImpulseDEObject accessors
#' @export
setMethod('strReport', 'ImpulseDE2Object', function(object) object@strReport)

### 2. Enable accession of public elements via list-like
### properties of ImpulseDE2Object.

# a) Enable names()
#' @name extractions
#' @export
setMethod('names', 'ImpulseDE2Object', function(x) {
  return( c("dfImpulseDE2Results", "vecDEGenes") )
})

# b) Enable object[[ element ]] operator
#' @name extractions
#' @export
setMethod('[[', c('ImpulseDE2Object', 'character', 'missing'), function(x, i, j, ...) {
  if(identical(i, "dfImpulseDE2Results")){ return(x@dfImpulseDE2Results)
  } else if(identical(i, "vecDEGenes")){ return(x@vecDEGenes)
  } else { return(NULL) }
})

# c) Enable object$element operator, which relies on [[ ]]
#' @name extractions
#' @export
setMethod('$', 'ImpulseDE2Object', function(x, name) x[[name]] )

### 3. Enable printing of report to .txt file

#' @return \code{writeReportToFile} writes strReport to file
#' @name ImpulseDEObject generics
#' @export
setGeneric('writeReportToFile', function(object, fileReport) standardGeneric('writeReportToFile'), valueClass = 'NULL')
#' @export
setMethod('writeReportToFile', signature(object='ImpulseDE2Object', fileReport='character'), 
          function(object, fileReport) write(object@strReport, file=fileReport, ncolumns=1) 
)
