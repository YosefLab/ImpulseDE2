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
#' @return (list length 2) with the following elements:
#' \itemize{
#'  \item \code{dds_dispersions} (vector number of genes) Inverse of gene-wise 
#'      negative binomial dispersion coefficients computed by DESeq2.
#'  \item \code{dds_resultsTable} (data frame) DESeq2 results.
#' }
#' @export

runDESeq2 <- function(dfAnnotationFull, arr2DCountData, 
  nProcessesAssigned=1, strMode="batch"){
  
  # Set number of processes to number of cores assigned if available
  nProcesses <- min(detectCores() - 1, nProcessesAssigned)
  register(MulticoreParam(nProcesses))
  
  dfCountData <- arr2DCountData[,colnames(arr2DCountData) %in% dfAnnotationFull$Replicate]
  colnames(dfCountData) <- dfAnnotationFull$Sample[
    match(colnames(dfCountData),dfAnnotationFull$Replicate)]
  # The covariate Timecourses, indicating the time series
  # experiment a sample belongs to, is ignored in batch mode,
  # in which replicates of one sample from different time series
  # experimentes are assumed to be i.i.d. In batch mode,
  # all samples should originate form the same series of 
  # independent samples.
  if(strMode=="batch" | strMode=="singlecell"){
    # Create DESeq2 data object
    dds <- DESeqDataSetFromMatrix(countData = dfCountData,
      colData = dfAnnotationFull,
      design = ~ Sample)
    # Run DESeq2
    ddsDESeqObject <- DESeq(dds, test = "LRT", 
      full = ~ Sample, reduced = ~ 1,
      parallel=TRUE)
  } else if(strMode=="timecourses"){
    # Create DESeq2 data object
    dds <- DESeqDataSetFromMatrix(countData = dfCountData,
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
  
  return(list(dds_dispersions,dds_resultsTable))
}