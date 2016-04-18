#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Annotation preparation    +++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# Prepare annotation table for internal use

# INPUT:
#   annotation_table: (Table samples x 2 [time and condition]) providing 
#       co-variables for the samples including condition and time points.
#       Time points must be numeric numbers.
#   data_tables: (Numeric 3D array genes x samples x replicates).
#       Contains expression values or similar locus-specific read-outs.
# OUTPUT:
#   annot: (Table samples x 2[time and condition]) annotation_table reduced to 
#       target samples 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

runDESeq2 <- function(annotation_table, data_table){
  
  dfCountData <- data_table[,colnames(data_table) %in% annotation_table$Replicate_name]
  colnames(dfCountData) <- annotation_table$Sample[match(colnames(dfCountData),annotation_table$Replicate_name)]
  # Create DESeq2 data object
  dds <- DESeqDataSetFromMatrix(countData = dfCountData,
    colData = annotation_table,
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