# Local
rm(list = ls())
boolCluster <- FALSE

print("Process data RNAseq Jankovic")
# Load Data set RNAseq
# 1. Counts
dfRNA <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/rsem_readCountsTable.txt",
  sep="\t",header=F,colClasses=c(
    "numeric","numeric","numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric","numeric","numeric"))
dfGeneIDs <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/gene_list.txt",sep="\t",header=F)
dfCellIDs <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/cell_list.txt",sep="\t",header=F)
vecSamples <- apply(dfCellIDs, 1, function(name){ unlist(strsplit(unlist(strsplit(name,"_1"))[1],"/"))[2] })

matDataA <- round(t(apply(dfRNA,1,as.numeric)))
vecGeneIDs <- as.vector(dfGeneIDs[,2])
vecboolidxDupIDs <- duplicated(vecGeneIDs)
if(sum(vecboolidxDupIDs)>0){
  vecGeneIDs[vecboolidxDupIDs] <- paste0(rep("DuplicatedGene_",sum(vecboolidxDupIDs)),seq(1,sum(vecboolidxDupIDs)))
}
rownames(matDataA) <- vecGeneIDs
colnames(matDataA) <- vecSamples
matDataA <- matDataA

# 2. Annotation
dfAnnotationRNA <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/AnnotationTable_RNAseqJankovic.tab",header=T)
dfAnnotationRNA$TimeCateg <- paste0(rep("_",length(dfAnnotationRNA$Time)),dfAnnotationRNA$Time)
dfAnnotationA <- dfAnnotationRNA

matDataA <- matDataA[,as.vector(dfAnnotationA$Sample)]
colnames(matDataA) <- dfAnnotationA$Sample

# All zero rows
matDataA[!is.finite(matDataA)] <- NA
vecboolNonzeroA <- apply(matDataA,1,function(gene){any(gene>0 & !is.na(gene) & is.finite(gene))})
matDataA <- matDataA[vecboolNonzeroA,]

# Summary matrices
matQval <- matrix(c("NA",rep(NA,4)),nrow=dim(matDataA)[1],ncol=4+1,byrow=TRUE)
matPval <- matrix(c("NA",rep(NA,4)),nrow=dim(matDataA)[1],ncol=4+1,byrow=TRUE)
matRunTime_RNAseqData <- array(NA,c(4))

rownames(matQval) <- rownames(matDataA)
colnames(matQval) <- c("Gene",paste0(
  rep("A_",4), 
  c("ImpulseDE2", "ImpulseDE", "DESeq2", "edge")
))
matQval[,"Gene"] <- rownames(matDataA)

rownames(matPval) <- rownames(matDataA)
colnames(matPval) <- c("Gene",paste0(
  rep("A_",4), 
  c("ImpulseDE2", "ImpulseDE", "DESeq2", "edge")
))
matPval[,"Gene"] <- rownames(matDataA)

names(matRunTime_RNAseqData) <- c("ImpulseDE2", "ImpulseDE", "DESeq2", "edge")


################################################################################
# ImpulseDE2
print("Run ImpulseDE2")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/building/code_files/ImpulseDE2_main.R")

# Create input data set
# Only retain non zero
matDataA_ImpulseDE2 <- matDataA

tm_ImpulseDE2A <- system.time({
  setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/ctrl/ImpulseDE2")
  strControlName = "ctrl"
  strCaseName = "case"
  lsImpulseDE_resultsA <- runImpulseDE2(matCountData=matDataA_ImpulseDE2, 
    dfAnnotation=dfAnnotationA,
    strCaseName = strCaseName, 
    strControlName=strControlName, 
    strMode="longitudinal",
    nProc=2, 
    Q_value=10^(-3),
    boolPlotting=FALSE)
  dfImpulseResultsA <- lsImpulseDE_resultsA$dfImpulseResults
  qvals_A <- dfImpulseResultsA$adj.p
  pvals_A <- dfImpulseResultsA$p
  ids_A <- as.vector(dfImpulseResultsA$Gene)
})

# Extract Results
matQval[ids_A,"A_ImpulseDE2"] <- qvals_A
matPval[ids_A,"A_ImpulseDE2"] <- pvals_A

matRunTime_RNAseqData["ImpulseDE2"] <- round(tm_ImpulseDE2A["elapsed"])

lsResDEcomparison_RNAseqdata <- list("matQval"=matQval, 
  "matPval"=matPval, 
  "matRunTime_RNAseqData"=matRunTime_RNAseqData)

save(lsResDEcomparison_RNAseqdata,file=file.path("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/lsResDEcomparison_RNAseqdata.RData"))
print("Finished ImpulseDE2")

################################################################################
# DESeq2
print("Run DESeq2")
library(DESeq2)
library(BiocParallel)
# Parallelisation
nProcessesAssigned <- 3
nProcesses <- min(detectCores() - 1, nProcessesAssigned)
register(MulticoreParam(nProcesses))

# Create input data set
# Only retain non zero
matDataA_DESeq2 <- matDataA

tm_DESeq2A <- system.time({
  # Create DESeq2 data object
  dds <- DESeqDataSetFromMatrix(countData = matDataA_DESeq2,
    colData = dfAnnotationA,
    design = ~TimeCateg + Condition + LongitudinalSeries)
  # Run DESeq2
  ddsDESeqObjectA <- DESeq(dds, test = "LRT", 
    full = ~TimeCateg + Condition + LongitudinalSeries, reduced = ~ TimeCateg + LongitudinalSeries,
    parallel=TRUE)
  
  # DESeq results for comparison
  dds_resultsTableA <- results(ddsDESeqObjectA)
  qvals_A <- dds_resultsTableA$padj
  pvals_A <- dds_resultsTableA$pvalue
  ids_A <- rownames(dds_resultsTableA)
})

# Extract Results
matQval[ids_A,"A_DESeq2"] <- qvals_A
matPval[ids_A,"A_DESeq2"] <- pvals_A
matRunTime_RNAseqData["DESeq2"] <- round(tm_DESeq2A["elapsed"])

lsResDEcomparison_RNAseqdata <- list("matQval"=matQval, 
  "matPval"=matPval, 
  "matRunTime_RNAseqData"=matRunTime_RNAseqData)

save(lsResDEcomparison_RNAseqdata,file=file.path("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/lsResDEcomparison_RNAseqdata.RData"))
print("Finished DESeq2.")
if(FALSE)  {
  ################################################################################
  # ImpulseDE
  # Run on server
  print("Run ImpulseDE")
  library(compiler)
  library(grDevices)
  library(graphics)
  library(stats)
  library(stats)
  library(utils)
  library(amap)
  library(ImpulseDE)
  
  # Create input data set
  # Only retain non zero, normalise and log transform
  matDataA_ImpulseDE <- matDataA/rep(sizeFactors(ddsDESeqObjectA), each = nrow(matDataA))
  matDataA_ImpulseDE[matDataA_ImpulseDE==0] <- 1
  matDataA_ImpulseDE <- log(matDataA_ImpulseDE)/log(2)
  
  dfAnnotationA_ImpulseDE <- dfAnnotationA[dfAnnotationA$Condition=="case",c("Time","Condition")]
  dfAnnotationA_ImpulseDE$Condition <- as.vector(dfAnnotationA_ImpulseDE$Condition)
  dfAnnotationA_ImpulseDE$Time <- as.numeric(as.vector(dfAnnotationA_ImpulseDE$Time))
  rownames(dfAnnotationA_ImpulseDE) <- as.vector(dfAnnotationA$Sample)
  tm_ImpulseDEA <- system.time({
    setwd("/data/yosef2/users/fischerd/data/RNAseq_Jovanovic/comparative_analysis/output/ImpulseDE")
    strControlName = NULL
    strCaseName = "case"
    lsImpulseDE_resultsA <- impulse_DE(
      expression_table = matDataA_ImpulseDE, 
      annotation_table = dfAnnotationA_ImpulseDE,
      colname_time = "Time", 
      colname_condition = "Condition", 
      control_timecourse = FALSE,
      control_name = strControlName,
      case_name = strCaseName, 
      expr_type = "Seq",
      plot_clusters = FALSE, 
      n_iter = 100, 
      n_randoms = 50000, 
      n_process = NCORES,
      Q_value = 0.01)
    dfImpulseResultsA <- lsImpulseDE_resultsA$DE_results
    qvals_A <- dfImpulseResultsA$adj.p
    ids_A <- dfImpulseResultsA$Gene
  })
  
  # Extract Results
  matQval[ids_A,"A_ImpulseDE"] <- qvals_A
  
  matRunTime_RNAseqData["ImpulseDE"] <- round(tm_ImpulseDE2A["elapsed"])
  
  lsResDEcomparison_RNAseqdata <- list("matQval"=matQval, 
    "matPval"=matPval, 
    "matRunTime_RNAseqData"=matRunTime_RNAseqData)
  
  save(lsResDEcomparison_RNAseqdata,file=file.path("/data/yosef2/users/fischerd/data/RNAseq_Jovanovic/comparative_analysis/output/lsResDEcomparison_RNAseqdata.RData"))
  print("Finished ImpulseDE")
  
}