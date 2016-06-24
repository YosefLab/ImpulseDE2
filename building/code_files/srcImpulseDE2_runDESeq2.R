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
#' @param matCountDataProc: (matrix genes x samples)
#'    Count data: Reduced version of \code{matCountData}. 
#'    For internal use.
#' @param dfAnnotationProc: (Table) Processed annotation table. 
#'    Lists co-variables of samples: 
#'    Sample, Condition, Time (numeric), TimeCateg (str)
#'    (and Timecourse). For internal use.
#' @param nProc: (scalar) [Default 1] 
#'    Number of processes for parallelisation.
#' @param strMode: (str) [Default "batch"] 
#'    {"batch","longitudinal","singlecell"}
#'    Mode of model fitting.
#'    
#' @return (list length 2) with the following elements:
#' \itemize{
#'  \item \code{dds_dispersions} (vector number of genes) Inverse of gene-wise 
#'      negative binomial dispersion coefficients computed by DESeq2.
#'  \item \code{dds_resultsTable} (data frame) DESeq2 results.
#' }
#' @export

runDESeq2 <- function(dfAnnotationProc, 
  matCountDataProc,
  nProc=1, 
  strControlName=NULL, 
  strMode="batch"){
  
  # Set number of processes for parallelisation
  register(MulticoreParam(nProc))
  
  if(is.null(strControlName)){
    # Without control data:
    # The covariate LongitudinalSeries, indicating the time series
    # experiment a sample belongs to, is ignored in batch mode,
    # in which replicates of one sample from different time series
    # experimentes are assumed to be i.i.d. In batch mode,
    # all samples should originate form the same series of 
    # independent samples.
    
    if(strMode=="batch" | strMode=="singlecell"){
      # Create DESeq2 data object
      dds <- suppressWarnings( DESeqDataSetFromMatrix(countData = matCountDataProc,
        colData = dfAnnotationProc,
        design = ~ TimeCateg) )
      # Run DESeq2
      ddsDESeqObject <- DESeq(dds, test = "LRT", 
        full = ~ TimeCateg, reduced = ~ 1,
        parallel=TRUE)
      
    } else if(strMode=="longitudinal"){
      # Create DESeq2 data object
      dds <- suppressWarnings( DESeqDataSetFromMatrix(countData = matCountDataProc,
        colData = dfAnnotationProc,
        design = ~ TimeCateg + LongitudinalSeries) )
      # Run DESeq2
      ddsDESeqObject <- DESeq(dds, test = "LRT", 
        full = ~ TimeCateg + LongitudinalSeries, reduced = ~ LongitudinalSeries,
        parallel=TRUE)
      
    } else {
      stop(paste0("ERROR: Unrecognised strMode in runDESeq2(): ",strMode))
    }
    # Get gene-wise dispersion estimates
    # var = mean + alpha * mean^2, alpha is dispersion
    # DESeq2 dispersion is 1/size used dnbinom (used in cost function
    # for evaluation of likelihood)
    dds_dispersions <- 1/dispersions(ddsDESeqObject) 
    names(dds_dispersions) <- rownames(ddsDESeqObject)
    # DESeq results for comparison
    dds_resultsTable <- results(ddsDESeqObject)
  } else {
    # With control data:
    
    # Hypothesis testing operates under the same model
    # for batch, singlecell and longitudinal: Condition
    # is a predictor included in LongitudinalSeries and therefore
    # not allowed by DESeq2. Note that that the dispersions
    # are fit on the more exact full model for LongitudinalSeries
    
    if(strMode=="batch" | strMode=="singlecell"){
      # Create DESeq2 data object
      dds <- suppressWarnings( DESeqDataSetFromMatrix(countData = matCountDataProc,
        colData = dfAnnotationProc,
        design = ~ TimeCateg + Condition) )
      # Run DESeq2
      ddsDESeqObject <- DESeq(dds, test = "LRT", 
        full = ~ TimeCateg + Condition, reduced = ~ TimeCateg,
        parallel=TRUE)
      
      # Get gene-wise dispersion estimates
      # var = mean + alpha * mean^2, alpha is dispersion
      # DESeq2 dispersion is 1/size used dnbinom (used in cost function
      # for evaluation of likelihood)
      dds_dispersions <- 1/dispersions(ddsDESeqObject)
      names(dds_dispersions) <- rownames(ddsDESeqObject)
      # DESeq results for comparison
      dds_resultsTable <- results(ddsDESeqObject)
      
    } else if(strMode=="longitudinal"){
      # Create DESeq2 data object
      # Define linear model suited to hypothesis testing
      dds <- suppressWarnings( DESeqDataSetFromMatrix(countData = matCountDataProc,
        colData = dfAnnotationProc,
        design = ~ TimeCateg + Condition) )
      # Run DESeq2
      ddsDESeqObjectTest <- DESeq(dds, test = "LRT", 
        full = ~ TimeCateg + Condition, reduced = ~ TimeCateg,
        parallel=TRUE)
      # DESeq results for comparison
      dds_resultsTable <- results(ddsDESeqObjectTest)
      
      # Define linear model suited to overdispersion fitting
      # Only need dispersions now, but the run time of the rest
      # of DESeq2 is negligible.
      dds <- suppressWarnings( DESeqDataSetFromMatrix(countData = matCountDataProc,
        colData = dfAnnotationProc,
        design = ~ TimeCateg + LongitudinalSeries) )
      # Run DESeq2
      ddsDESeqObjectFit <- DESeq(dds, test = "LRT", 
        full = ~ TimeCateg + LongitudinalSeries, reduced = ~ TimeCateg,
        parallel=TRUE)
      # Get gene-wise dispersion estimates
      # var = mean + alpha * mean^2, alpha is dispersion
      # DESeq2 dispersion is 1/size used dnbinom (used in cost function
      # for evaluation of likelihood)
      dds_dispersions <- 1/dispersions(ddsDESeqObjectFit)
      names(dds_dispersions) <- rownames(ddsDESeqObjectFit)
      
    } else {
      stop(paste0("ERROR: Unrecognised strMode in runDESeq2(): ",strMode))
    }
  }
  
  return(list(dds_dispersions,dds_resultsTable))
}