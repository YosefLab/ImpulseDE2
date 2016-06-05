#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++     Process single cell data    +++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Prepare single cell data for analysis
#' 
#' Check that input is correctly supplied and formatted. Then process input
#' data for analysis.
#' (I) Helper functions:
#'    checkNull() Check whether object was supplied (is not NULL).
#'    checkCounts() Checks whether elements are count data.
#'    checkNumeric() Checks whether elements are numeric.
#'    checkCounts() Checks whether elements are count data.
#' 
#' @seealso Called by \code{runPseudoDE}.
#' 
#' @param matCounts: (matrix genes x samples)
#'    Count data of all cells, unobserved entries are NA.
#' @param vecPseudotime: (numerical vector length number of cells)
#'    Pseudotime coordinates (1D) of cells: One scalar per cell.
#'    Has to be named: Names of elements are cell names.
#' 
#' @return matCountsProc: (matrix genes x samples)
#'    Processed count data of all cells, unobserved entries are NA.
#' @return vecPseudotimeProc: (numerical vector length number of cells)
#'    Processed pseudotime coordinates (1D) of cells: One scalar per cell.
#'    Names of elements are cell names.
#' @export

processSCData <- function(matCounts,
  vecPseudotime,
  boolDEAnalysisImpulseModel,
  boolDEAnalysisModelFree ){
  
  # Check whether object was supplied (is not NULL).
  checkNull <- function(objectInput){
    if(is.null(objectInput)){
      stop(paste0( "ERROR: ", objectInput," was not given as input." ))
    }
  }
  
  # Checks whether elements are numeric
  checkNumeric <- function(matInput){
    if(any(!is.numeric(matInput))){
      stop(paste0( "ERROR: ", matInput, " contains non-numeric elements. Requires count data." ))
    }
  }
  # Checks whether elements are count data: non-negative integer finite numeric elements.
  # Note that NA are allowed.
  checkCounts <- function(matInput){
    checkNumeric(matInput)
    if(any(matInput %% 1 != 0)){
      stop(paste0( "ERROR: ", matInput, " contains non-integer elements. Requires count data." ))
    }
    if(any(!is.finite(matInput))){
      stop(paste0( "ERROR: ", matInput, " contains infinite elements. Requires count data." ))
    }
    if(any(matInput<0)){
      stop(paste0( "ERROR: ", matInput, " contains negative elements. Requires count data." ))
    }
  }
  
  # (I) Check input
  # 1. matCounts
  checkNull(matCounts)
  checkCounts(matCounts)
  
  # 2. vecPseudotime
  checkNull(vecPseudotime)
  checkNumeric(vecPseudotime)
  if(!all(names(vecPseudotime) %in% colnames(matCounts))){
    stop("ERROR: Not all cells in vecPseudotime are given matCounts.")
  }
  
  # 3. boolDEAnalysisImpulseModel and boolDEAnalysisModelFree
  if(!boolDEAnalysisImpulseModel & !boolDEAnalysisModelFree){
    stop(paste0("ERROR: No differential analysis mode selected.",
      "Set either boolDEAnalysisImpulseModel or boolDEAnalysisModelFree as TRUE."))
  }
  
  # (II) Process data
  # Convert from data frame to matrix for zinb()
  if(is.list(matCounts)){
    matCounts <- data.matrix(matCounts)
  }
  # Take out NA cells
  vecPseudotimeProc <- vecPseudotime[!is.na(vecPseudotime)]
  # Adjust ording of cells in objects
  matCountsProc <- matCounts[,names(vecPseudotime)]
  # Remove all zero genes: Mu is initialised
  # as log(sum counts)
  matCountsProc <- matCountsProc[apply(matCountsProc,1,function(gene){any(gene!=0)}),]
  
  return(list(matCountsProc=matCountsProc,
    vecPseudotimeProc=vecPseudotimeProc))
}