# Cluster
rm(list = ls())
NCORES <- 16

print("Process data RNAseq Jankovic")
# Load Data set RNAseq
# 1. Counts
dfRNA <- read.table("/data/yosef2/users/fischerd/data/RNAseq_Jovanovic/summary/rsem/rsem_readCountsTable.txt",
  sep="\t",header=F,colClasses=c(
    "numeric","numeric","numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric","numeric","numeric"))
dfGeneIDs <- read.table("/data/yosef2/users/fischerd/data/RNAseq_Jovanovic/summary/rsem/gene_list.txt",sep="\t",header=F)
dfCellIDs <- read.table("/data/yosef2/users/fischerd/data/RNAseq_Jovanovic/summary/rsem/cell_list.txt",sep="\t",header=F)
vecSamples <- apply(dfCellIDs, 1, function(name){ unlist(strsplit(unlist(strsplit(name,"_1"))[1],"/"))[2] })

matDataA <- round(t(apply(dfRNA,1,as.numeric)))
vecGeneIDs <- as.vector(dfGeneIDs[,2])
vecboolidxDupIDs <- duplicated(vecGeneIDs)
if(sum(vecboolidxDupIDs)>0){
  vecGeneIDs[vecboolidxDupIDs] <- paste0(rep("DuplicatedGene_",sum(vecboolidxDupIDs)),seq(1,sum(vecboolidxDupIDs)))
}
rownames(matDataA) <- vecGeneIDs
colnames(matDataA) <- vecSamples
#matDataA <- matDataA[1:2000,]

# 2. Annotation
dfAnnotationRNA <- read.table("/data/yosef2/users/fischerd/data/RNAseq_Jovanovic/comparative_analysis/input/AnnotationTable_RNAseqJankovic.tab",header=T)
dfAnnotationRNA$TimeCateg <- paste0(rep("_",length(dfAnnotationRNA$Time)),dfAnnotationRNA$Time)
dfAnnotationA <- dfAnnotationRNA

# Expend 0h ctrl sample to both conditions
vecExpanedSamples <- dfAnnotationA$Sample
vecExpanedSamples[vecExpanedSamples=="SRR1525500cs"] <- "SRR1525500"
vecExpanedSamples[vecExpanedSamples=="SRR1525513cs"] <- "SRR1525513"

matDataA <- matDataA[,as.vector(vecExpanedSamples)]
colnames(matDataA) <- dfAnnotationA$Sample
#only case
matDataA <- matDataA[,dfAnnotationA$Condition=="case"]
colnames(matDataA) <- dfAnnotationA[dfAnnotationA$Condition=="case",]$Sample
dfAnnotationA <- dfAnnotationA[dfAnnotationA$Condition=="case",]
rownames(dfAnnotationA) <- dfAnnotationA$Sample

# All zero rows
matDataA[!is.finite(matDataA)] <- NA
vecboolNonzeroA <- apply(matDataA,1,function(gene){any(gene>0 & !is.na(gene) & is.finite(gene))})
matDataA <- matDataA[vecboolNonzeroA,]

# Summary matrices
matQval_RNAseqData <- matrix(c("NA",rep(NA,4)),nrow=dim(matDataA)[1],ncol=4+1,byrow=TRUE)
matPval_RNAseqData <- matrix(c("NA",rep(NA,4)),nrow=dim(matDataA)[1],ncol=4+1,byrow=TRUE)
matRunTime_RNAseqData <- array(NA,c(4))

rownames(matQval_RNAseqData) <- rownames(matDataA)
colnames(matQval_RNAseqData) <- c("Gene",paste0(
  rep("A_",4), 
  c("ImpulseDE2", "ImpulseDE", "DESeq2", "edge")
))
matQval_RNAseqData[,"Gene"] <- rownames(matDataA)

rownames(matPval_RNAseqData) <- rownames(matDataA)
colnames(matPval_RNAseqData) <- c("Gene",paste0(
  rep("A_",4), 
  c("ImpulseDE2", "ImpulseDE", "DESeq2", "edge")
))
matPval_RNAseqData[,"Gene"] <- rownames(matDataA)

names(matRunTime_RNAseqData) <- c("ImpulseDE2", "ImpulseDE", "DESeq2", "edge")

################################################################################
# ImpulseDE2
print("Run ImpulseDE2")
source("/data/yosef2/users/fischerd/code/ImpulseDE2/building/code_files/ImpulseDE2_main.R")

# Create input data set
# Only retain non zero
matDataA_ImpulseDE2 <- matDataA

tm_ImpulseDE2A <- system.time({
  setwd("/data/yosef2/users/fischerd/data/RNAseq_Jovanovic/comparative_analysis/output/ImpulseDE2")
  strControlName = NULL
  strCaseName = "case"
  lsImpulseDE_resultsA <- runImpulseDE2(matCountData=matDataA_ImpulseDE2, 
    dfAnnotation=dfAnnotationA,
    strCaseName = strCaseName, 
    strControlName=strControlName, 
    strMode="batch",
    nProc=NCORES, 
    Q_value=10^(-3),
    boolPlotting=FALSE)
  dfImpulseResultsA <- lsImpulseDE_resultsA$dfImpulseResults
  qvals_A <- dfImpulseResultsA$adj.p
  pvals_A <- dfImpulseResultsA$p
  ids_A <- as.vector(dfImpulseResultsA$Gene)
})

# Extract Results
matQval_RNAseqData[ids_A,"A_ImpulseDE2"] <- qvals_A
matPval_RNAseqData[ids_A,"A_ImpulseDE2"] <- pvals_A

matRunTime_RNAseqData["ImpulseDE2"] <- round(tm_ImpulseDE2A["elapsed"])

lsResDEcomparison_RNAseqdata <- list("matQval_RNAseqData"=matQval_RNAseqData, 
  "matPval_RNAseqData"=matPval_RNAseqData, 
  "matRunTime_RNAseqData"=matRunTime_RNAseqData)

save(lsResDEcomparison_RNAseqdata,file=file.path("/data/yosef2/users/fischerd/data/RNAseq_Jovanovic/comparative_analysis/output/lsResDEcomparison_RNAseqdata.RData"))
print("Finished ImpulseDE2")

################################################################################
# DESeq2
print("Run DESeq2")
library(DESeq2)
library(BiocParallel)
# Parallelisation
nProcessesAssigned <- NCORES
register(MulticoreParam(nProcesses))

# Create input data set
# Only retain non zero
matDataA_DESeq2 <- matDataA

tm_DESeq2A <- system.time({
  # Create DESeq2 data object
  dds <- DESeqDataSetFromMatrix(countData = matDataA_DESeq2,
    colData = dfAnnotationA,
    design = ~TimeCateg)
  # Run DESeq2
  ddsDESeqObjectA <- DESeq(dds, test = "LRT", 
    full = ~ TimeCateg, reduced = ~ 1,
    parallel=TRUE)
  
  # DESeq results for comparison
  dds_resultsTableA <- results(ddsDESeqObjectA)
  qvals_A <- dds_resultsTableA$padj
  pvals_A <- dds_resultsTableA$pvalue
  ids_A <- rownames(dds_resultsTableA)
})

# Extract Results
matQval_RNAseqData[ids_A,"A_DESeq2"] <- qvals_A
matPval_RNAseqData[ids_A,"A_DESeq2"] <- pvals_A
matRunTime_RNAseqData["DESeq2"] <- round(tm_DESeq2A["elapsed"])

lsResDEcomparison_RNAseqdata <- list("matQval_RNAseqData"=matQval_RNAseqData, 
  "matPval_RNAseqData"=matPval_RNAseqData, 
  "matRunTime_RNAseqData"=matRunTime_RNAseqData)

save(lsResDEcomparison_RNAseqdata,file=file.path("/data/yosef2/users/fischerd/data/RNAseq_Jovanovic/comparative_analysis/output/lsResDEcomparison_RNAseqdata.RData"))
print("Finished DESeq2.")

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
matQval_RNAseqData[ids_A,"A_ImpulseDE"] <- qvals_A

matRunTime_RNAseqData["ImpulseDE"] <- round(tm_ImpulseDEA["elapsed"])

lsResDEcomparison_RNAseqdata <- list("matQval_RNAseqData"=matQval_RNAseqData, 
  "matPval_RNAseqData"=matPval_RNAseqData, 
  "matRunTime_RNAseqData"=matRunTime_RNAseqData)

save(lsResDEcomparison_RNAseqdata,file=file.path("/data/yosef2/users/fischerd/data/RNAseq_Jovanovic/comparative_analysis/output/lsResDEcomparison_RNAseqdata.RData"))
print("Finished ImpulseDE")