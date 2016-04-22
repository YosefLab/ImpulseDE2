#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Annotation preparation    +++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# Prepare annotation table for internal use

# INPUT:
#   dfAnnotationFull: (Table samples x 2 [time and condition]) providing 
#       co-variables for the samples including condition and time points.
#       Time points must be numeric numbers.
#   arr2DCountDatas: (Numeric 3D array genes x samples x replicates).
#       Contains expression values or similar locus-specific read-outs.
# OUTPUT:
#   annot: (Table samples x 2[time and condition]) dfAnnotationFull reduced to 
#       target samples 

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