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
#' @param strCaseName: (str) [Default NULL] 
#'    Name of the case condition in \code{dfAnnotation}.
#' @param strControlName: (str) [Default NULL] 
#'    Name of the control condition in \code{dfAnnotation}.
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
  strCaseName=NULL,
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
    
    if(strMode=="batch"){
      # Create DESeq2 data object
      dds <- suppressWarnings( DESeqDataSetFromMatrix(countData = matCountDataProc,
        colData = dfAnnotationProc,
        design = ~Condition + Condition:TimeCateg) )
      # Run DESeq2
      ddsDESeqObject <- DESeq(dds, test = "LRT", 
        full = ~Condition + Condition:TimeCateg,
        reduced = ~ TimeCateg,
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
      # Note: Have to distinguish model formula between
      # 1. One condition only has one series (i.e. that condition is 
      # a linear combination of that series)
      # 2. Both conditions have multiple series and are not linear
      # combinations of any one series.
      # Note that Condition is in both null models as it complements
      # the longitudinal series information in LSnested.
      boolSingleSeriesCond <- FALSE
      if(length(unique( dfAnnotationProc[dfAnnotationProc$Condition==strCaseName,]$LongSerNested ))==1){
        boolSingleSeriesCond <- TRUE
      }
      if(length(unique( dfAnnotationProc[dfAnnotationProc$Condition==strControlName,]$LongSerNested ))==1){
        boolSingleSeriesCond <- TRUE
      }
      if(boolSingleSeriesCond){
        # Create DESeq2 data object: Don't need interaction of
        # Condition and LongSerNested as longitudinal label is irrelevant
        # for one condition as it covers the same sample range as this
        # condition.
        dds <- suppressWarnings( DESeqDataSetFromMatrix(countData = matCountDataProc,
          colData = dfAnnotationProc,
          design = ~Condition + LongSerNested + Condition:TimeCateg) )
        # Run DESeq2
        ddsDESeqObject <- DESeq(dds, test = "LRT", 
          full = ~Condition + LongSerNested + Condition:TimeCateg,
          reduced = ~Condition + LongSerNested + TimeCateg,
          parallel=TRUE)
      } else {
        # Create DESeq2 data object: Need interaction between 
        # Condition and LongSerNested if each condition has multiple series.
        dds <- suppressWarnings( DESeqDataSetFromMatrix(countData = matCountDataProc,
          colData = dfAnnotationProc,
          design = ~Condition + Condition:LongSerNested + Condition:TimeCateg) )
        # Run DESeq2
        ddsDESeqObject <- DESeq(dds, test = "LRT", 
          full = ~Condition + Condition:LongSerNested + Condition:TimeCateg,
          reduced = ~Condition + Condition:LongSerNested + TimeCateg,
          parallel=TRUE)
      }
      # Get gene-wise dispersion estimates
      # var = mean + alpha * mean^2, alpha is dispersion
      # DESeq2 dispersion is 1/size used dnbinom (used in cost function
      # for evaluation of likelihood)
      dds_dispersions <- 1/dispersions(ddsDESeqObject)
      names(dds_dispersions) <- rownames(ddsDESeqObject)
      
    } else {
      stop(paste0("ERROR: Unrecognised strMode in runDESeq2(): ",strMode))
    }
  }
  
  return(list(dds_dispersions,dds_resultsTable))
}