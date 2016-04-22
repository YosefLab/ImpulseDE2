#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++     Run DESeq2    +++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Run DESeq2 and extract differential expression analysis results and 
### overdispersion coefficients.

# INPUT:
#   dfAnnotationFull: (Table) Lists co-variables of individual replicates:
#       Sample, Condition, Time. Time must be numeric.
#   arr2DCountData: (2D array genes x replicates) Count data: Reduced 
#       version of \code{matCountData}. For internal use.
# OUTPUT:
#   List of:
#   dds_dispersions: (vector number of genes) Inverse of gene-wise 
#       negative binomial dispersion coefficients computed by DESeq2.
#   dds_resultsTable: (data frame) DESeq2 results.

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

runDESeq2 <- function(dfAnnotationFull, arr2DCountData){
  
  dfCountData <- arr2DCountData[,colnames(arr2DCountData) %in% dfAnnotationFull$Replicate]
  colnames(dfCountData) <- dfAnnotationFull$Sample[match(colnames(dfCountData),dfAnnotationFull$Replicate)]
  # Create DESeq2 data object
  dds <- DESeqDataSetFromMatrix(countData = dfCountData,
    colData = dfAnnotationFull,
    design = ~ Sample)
  # Run DESeq2
  ddsDESeqObject <- DESeq(dds, test = "LRT", full = ~ Sample, reduced = ~1)
  # Get gene-wise dispersion estimates
  # var = mean + alpha * mean^2, alpha is dispersion
  # DESeq2 dispersion is 1/size used dnbinom (used in cost function
  # for evaluation of likelihood)
  dds_dispersions <- 1/dispersions(ddsDESeqObject) 
  # DESeq results for comparison
  dds_resultsTable <- results(ddsDESeqObject)
  
  return(list(dds_dispersions,dds_resultsTable))
}