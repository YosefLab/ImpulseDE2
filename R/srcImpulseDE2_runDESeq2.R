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
#'    
#' @return (list length 2) with the following elements:
#' \itemize{
#'  \item \code{vecDispersions} (vector number of genes) Inverse of gene-wise 
#'      negative binomial dispersion coefficients computed by DESeq2.
#'  \item \code{ddsResults} (data frame) DESeq2 results.
#' }
#' @export

runDESeq2 <- function(dfAnnotationProc, 
  matCountDataProc,
  boolCaseCtrl,
  vecConfounders){
  
  # Get gene-wise dispersion estimates
  # var = mean + alpha * mean^2, alpha is dispersion
  # DESeq2 dispersion is 1/size used dnbinom (used in cost function
  # for evaluation of likelihood)
  
  # Conditions to be fitted separately
  if(boolCaseCtrl){
    vecLabels <- c("combined","case","control")
  } else {
    vecLabels <- c("case")
  }
  
  if(boolCaseCtrl){    
    lsSamplesByCond <- list(
      combined=dfAnnotationProc$Sample,
      case=dfAnnotationProc[dfAnnotationProc$Condition=="case",]$Sample,
      control=dfAnnotationProc[dfAnnotationProc$Condition=="control",]$Sample )
  } else {
    lsSamplesByCond <- list(
      case=dfAnnotationProc[dfAnnotationProc$Condition=="case",]$Sample )
  }
  
  lsvecDispersions <- lapply(vecLabels, function(label){
    # Run DESeq2 on samples from ONE CONDITION
    # Estimate one dispersion parameter PER CONDITION.
    print(paste0("Run DESeq2 on condition ", match(label, vecLabels), " (", label,
                 ") out of ", length(vecLabels), " condition(s)."))
    vecSamples <- lsSamplesByCond[[label]]
    dds <- NULL
    ddsDESeqObject <- NULL
    
    if(is.null(vecConfounders)){
      # No batch correction
      dds <- suppressWarnings( DESeqDataSetFromMatrix(
        countData = matCountDataProc[,vecSamples],
        colData = dfAnnotationProc[vecSamples,],
        design = ~ TimeCateg) )
      ddsDESeqObject <- DESeq(dds, test = "LRT", 
                              quiet=TRUE, parallel=TRUE,           
                              full = ~ TimeCateg, 
                              reduced = ~ 1)
      
    } else {
      # With batch correction
      # Catch case in which at least one batch has unique time points.
      tryCatch({
        dds <- suppressWarnings( DESeqDataSetFromMatrix(
          countData = matCountDataProc[,vecSamples],
          colData = dfAnnotationProc[vecSamples,],
          design = ~TimeCateg + Batch  ) )
        ddsDESeqObject <- DESeq(dds, test = "LRT",
                                quiet=TRUE, parallel=TRUE,           
                                full = ~TimeCateg + Batch, 
                                reduced = ~Batch)
      }, error=function(strErrorMsg){
        print(strErrorMsg)
        print(paste0("WARNING: DESeq2 failed on full model - dispersions may be inaccurate.",
                     "Estimating dispersions on reduced model formulation [full = ~Batch",
                     " reduced = ~1]. Supply externally generated dispersion parameters via ",
                     "lsvecDispersionsExternal if there is a more accurate model for your data set."))
        warning("Warning generated in dispersion factor estimation, read stdout.")
      }, finally={
        if(is.null(dds)){
          dds <- suppressWarnings( DESeqDataSetFromMatrix(
            countData = matCountDataProc[,vecSamples],
            colData = dfAnnotationProc[vecSamples,],
            design = ~Batch) )
          ddsDESeqObject <- DESeq(dds, test = "LRT",
                                  quiet=TRUE, parallel=TRUE,                      
                                  full = ~Batch, 
                                  reduced = ~1)
        }
      })
      
    }
    
    vecDispersions <- 1/dispersions(ddsDESeqObject)
    names(vecDispersions) <- rownames(ddsDESeqObject)
    return(vecDispersions)
  })
  names(lsvecDispersions) <- vecLabels

  return(lsvecDispersions)
}