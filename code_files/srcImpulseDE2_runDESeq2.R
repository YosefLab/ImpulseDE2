#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++     Run DESeq2    +++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Wrapper function for running DESeq2
#' 
#' Run DESeq2 and extract differential expression analysis results and 
#' overdispersion coefficients.
#' 
#' @seealso Called by \code{runImpulseDE2}.
#' 
#' @param dfAnnotationFull (Table) Lists co-variables of individual replicates:
#'    Replicate, Sample, Condition, Time. Time must be numeric.
#' @param arr2DCountData (2D array genes x replicates) Count data: Reduced 
#'    version of \code{matCountData}. For internal use.
#' @param nProcessesAssigned: (scalar) [Default 1] Number of processes for parallelisation. The
#'    specified value is internally changed to \code{min(detectCores() - 1, nProc)} 
#'    using the \code{detectCores} function from the package \code{parallel} to 
#'    avoid overload.
#' @param strMode: (str) [Default "batch"] {"batch","timecourses","singlecell"}
#'    Mode of model fitting.
#'    
#' @return (list length 2) with the following elements:
#' \itemize{
#'  \item \code{dds_dispersions} (vector number of genes) Inverse of gene-wise 
#'      negative binomial dispersion coefficients computed by DESeq2.
#'  \item \code{dds_resultsTable} (data frame) DESeq2 results.
#' }
#' @export

runDESeq2 <- function(dfAnnotationFull, arr2DCountData,
  nProcessesAssigned=1, strControlName=NULL, strMode="batch"){
  
  # Set number of processes to number of cores assigned if available
  nProcesses <- min(detectCores() - 1, nProcessesAssigned)
  register(MulticoreParam(nProcesses))
  
  if(is.null(strControlName)){
    # Without control data
    
    # The covariate Timecourses, indicating the time series
    # experiment a sample belongs to, is ignored in batch mode,
    # in which replicates of one sample from different time series
    # experimentes are assumed to be i.i.d. In batch mode,
    # all samples should originate form the same series of 
    # independent samples.
    if(strMode=="batch" | strMode=="singlecell"){
      # Create DESeq2 data object
      dds <- DESeqDataSetFromMatrix(countData = arr2DCountData,
        colData = dfAnnotationFull,
        design = ~ Sample)
      # Run DESeq2
      ddsDESeqObject <- DESeq(dds, test = "LRT", 
        full = ~ Sample, reduced = ~ 1,
        parallel=TRUE)
    } else if(strMode=="timecourses"){
      # Create DESeq2 data object
      dds <- DESeqDataSetFromMatrix(countData = arr2DCountData,
        colData = dfAnnotationFull,
        design = ~ Sample + Timecourse)
      # Run DESeq2
      ddsDESeqObject <- DESeq(dds, test = "LRT", 
        full = ~ Sample + Timecourse, reduced = ~ Timecourse,
        parallel=TRUE)
    } else {
      stop(paste0("ERROR: Unrecognised strMode in runDESeq2(): ",strMode))
    }
    # Get gene-wise dispersion estimates
    # var = mean + alpha * mean^2, alpha is dispersion
    # DESeq2 dispersion is 1/size used dnbinom (used in cost function
    # for evaluation of likelihood)
    dds_dispersions <- 1/dispersions(ddsDESeqObject) 
    # DESeq results for comparison
    dds_resultsTable <- results(ddsDESeqObject)
  } else {
    # With control data: 
    # 1. Fit dispersion with full information, i.e. sample 
    # names.  
    # 2. Evaluate p-value for differential expression based 
    # on comparison condition versus no-contition information. 
    
    # Create DESeq2 data object
    dds <- DESeqDataSetFromMatrix(countData = arr2DCountData,
      colData = dfAnnotationFull,
      design = ~ Sample)
    # Run DESeq2
    ddsDESeqObjectSample <- DESeq(dds, test = "LRT", 
      full = ~ Sample, reduced = ~ 1,
      parallel=TRUE)
    # Get gene-wise dispersion estimates
    # var = mean + alpha * mean^2, alpha is dispersion
    # DESeq2 dispersion is 1/size used dnbinom (used in cost function
    # for evaluation of likelihood)
    dds_dispersions <- 1/dispersions(ddsDESeqObjectSample)
    
    dds2 <- DESeqDataSetFromMatrix(countData = arr2DCountData,
      colData = dfAnnotationFull,
      design = ~ Condition)
    # Run DESeq2
    ddsDESeqObjectCond <- DESeq(dds2, test = "LRT", 
      full = ~ Condition, reduced = ~ 1,
      parallel=TRUE)
    # DESeq results for comparison
    dds_resultsTable <- results(ddsDESeqObjectCond)
  }      
  return(list(dds_dispersions,dds_resultsTable))
}