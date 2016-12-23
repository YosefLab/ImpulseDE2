################################################################################
#################     ImpulseDE2 output container class     ####################
################################################################################

# Define class unions for slots
setClassUnion('numericORNULL', members = c('numeric', 'NULL'))

### 1. Define output container class

setClass(
  'ImpulseDE2Object',
  slots = c(
    dfImpulseDE2Results = "data.frame",
    vecDEGenes          = "characterORNULL",
    lsModelFits         = "list",
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
setGeneric('lsModelFits', function(object) standardGeneric('lsModelFits'), valueClass = 'list')
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
setGeneric('vecConfounders', function(object) standardGeneric('vecConfounders'), valueClass = 'characterOrnumericOrNULL')
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
setGeneric('scaQThres', function(object) standardGeneric('scaQThres'), valueClass = 'characterOrnumericOrNULL')
#' @name ImpulseDEObject accessors
#' @export
setMethod('scaQThres', 'ImpulseDE2Object', function(object) object@scaQThres)

#' @return \code{strReport} retrieves strReport
#' @name ImpulseDEObject generics
#' @export
setGeneric('strReport', function(object) standardGeneric('strReport'), valueClass = 'characterOrnumericOrNULL')
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
