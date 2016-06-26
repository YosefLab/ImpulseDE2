rm(list = ls())
NCORES = 16

print("Process data iChIP")
# Load Data set iChIP
# 1. Counts
dfH3K4me1 <- read.table("/data/yosef2/users/fischerd/data/iChIP_Friedman/comparative_analysis/input/H3K4me1_EryLineage.tab",
  sep="\t",header=T,colClasses=c(
    "character","character","character",
    "numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric"))
vecPeaksIDs <- apply(dfH3K4me1[,1:3],1,function(coord){paste(coord,collapse="-")})
matDataA <- dfH3K4me1[,4:dim(dfH3K4me1)[2]]
matDataA <- t(apply(matDataA,1,as.numeric))

rownames(matDataA) <- vecPeaksIDs
colnames(matDataA) <- colnames(dfH3K4me1)[4:dim(dfH3K4me1)[2]]

# 2. Annotation
dfAnnotationiChIP <- read.table("/data/yosef2/users/fischerd/data/iChIP_Friedman/comparative_analysis/input/AnnotationTable_iChIPdata_H3K4me1_EryLineage.tab",header=T)
dfAnnotationiChIP$TimeCateg <- paste0(rep("_",length(dfAnnotationiChIP$Time)),dfAnnotationiChIP$Time)
rownames(dfAnnotationiChIP) <- dfAnnotationiChIP$Sample

dfAnnotationA <- dfAnnotationiChIP

# All zero rows
vecboolNonzeroA <- apply(matDataA,1,function(gene){any(gene>0)})
matDataA <- matDataA[vecboolNonzeroA,]

# Summary matrices
matQval_iChIPData <- matrix(c("NA",rep(NA,4)),nrow=dim(matDataA)[1],ncol=4+1,byrow=TRUE)
matPval_iChIPData <- matrix(c("NA",rep(NA,4)),nrow=dim(matDataA)[1],ncol=4+1,byrow=TRUE)
matRunTime_iChIPData <- array(NA,c(4))

rownames(matQval_iChIPData) <- rownames(matDataA)
colnames(matQval_iChIPData) <- c("Gene",paste0(
  rep("A_",4), 
  c("ImpulseDE2", "ImpulseDE", "DESeq2", "edge")
))
matQval_iChIPData[,"Gene"] <- rownames(matDataA)

rownames(matPval_iChIPData) <- rownames(matDataA)
colnames(matPval_iChIPData) <- c("Gene",paste0(
  rep("A_",4), 
  c("ImpulseDE2", "ImpulseDE", "DESeq2", "edge")
))
matPval_iChIPData[,"Gene"] <- rownames(matDataA)

names(matRunTime_iChIPData) <- c("ImpulseDE2", "ImpulseDE", "DESeq2", "edge")

################################################################################
# ImpulseDE2
print("Run ImpulseDE2")
source("/data/yosef2/users/fischerd/code/ImpulseDE2/building/code_files/ImpulseDE2_main.R")

# Create input data set
# Only retain non zero
matDataA_ImpulseDE2 <- matDataA

tm_ImpulseDE2A <- system.time({
  setwd("/data/yosef2/users/fischerd/data/iChIP_Friedman/comparative_analysis/output/ImpulseDE2")
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
})

# Extract Results
matQval_iChIPData[rownames(dfImpulseResultsA),"A_ImpulseDE2"] <- qvals_A
matPval_iChIPData[rownames(dfImpulseResultsA),"A_ImpulseDE2"] <- pvals_A

matRunTime_iChIPData["ImpulseDE2"] <- round(tm_ImpulseDE2A["elapsed"])

lsResDEcomparison_iChIPdataH3K4me1 <- list("matQval_iChIPData"=matQval_iChIPData, 
  "matPval_iChIPData"=matPval_iChIPData, 
  "matRunTime_iChIPData"=matRunTime_iChIPData)

save(lsResDEcomparison_iChIPdataH3K4me1,file=file.path("/data/yosef2/users/fischerd/data/iChIP_Friedman/comparative_analysis/output/lsResDEcomparison_iChIPdataH3K4me1_EryLineage.RData"))
print("Finished ImpulseDE2")

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
    design = ~TimeCateg)
  # Run DESeq2
  ddsDESeqObjectA <- DESeq(dds, test = "LRT", 
    full = ~ TimeCateg, reduced = ~ 1,
    parallel=TRUE)
  
  # DESeq results for comparison
  dds_resultsTableA <- results(ddsDESeqObjectA)
  qvals_A <- dds_resultsTableA$padj
  pvals_A <- dds_resultsTableA$pvalue
})

# Extract Results
matQval_iChIPData[rownames(dds_resultsTableA),"A_DESeq2"] <- qvals_A
matPval_iChIPData[rownames(dds_resultsTableA),"A_DESeq2"] <- pvals_A
matRunTime_iChIPData["DESeq2"] <- round(tm_DESeq2A["elapsed"])

lsResDEcomparison_iChIPdataH3K4me1 <- list("matQval_iChIPData"=matQval_iChIPData, 
  "matPval_iChIPData"=matPval_iChIPData, 
  "matRunTime_iChIPData"=matRunTime_iChIPData)

save(lsResDEcomparison_iChIPdataH3K4me1,file=file.path("/data/yosef2/users/fischerd/data/iChIP_Friedman/comparative_analysis/output/lsResDEcomparison_iChIPdataH3K4me1_EryLineage.RData"))
print("Finished DESeq2.")

# ImpulseDE
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
matDataA_ImpulseDE <- matDataA/rep(sizeFactors(ddsDESeqObjectA), each = nrow(matDataA))
matDataA_ImpulseDE <- log(matDataA_ImpulseDE+1)/log(2)

dfAnnotationA_ImpulseDE <- dfAnnotationA[,c("Time","Condition")]
dfAnnotationA_ImpulseDE$Condition <- as.vector(dfAnnotationA_ImpulseDE$Condition)
dfAnnotationA_ImpulseDE$Time <- as.numeric(as.vector(dfAnnotationA_ImpulseDE$Time))
rownames(dfAnnotationA_ImpulseDE) <- as.vector(dfAnnotationA$Sample)

tm_ImpulseDEA <- system.time({
  setwd("/data/yosef2/users/fischerd/data/iChIP_Friedman/comparative_analysis/output/ImpulseDE")
  strControlName = NULL
  strCaseName = "case"
  lsImpulseDE_resultsA <- impulse_DE(
    expression_table = matDataA_ImpulseDE, 
    annotation_table = dfAnnotationA_ImpulseDE,
    colname_time = "Time", colname_condition = "Condition", 
    control_timecourse = FALSE,
    control_name = strControlName, case_name = strCaseName, 
    expr_type = "Seq",
    plot_clusters = FALSE, 
    n_iter = 100, n_randoms = 50000, 
    n_process = 3,
    Q_value = 0.01)
  dfImpulseResultsA <- lsImpulseDE_resultsA$impulse_fit_results
  qvals_A <- dfImpulseResultsA$adj.p
  pvals_A <- dfImpulseResultsA$p
})

# Extract Results
matQval_iChIPData[rownames(dds_resultsTableA),"A_ImpulseDE"] <- qvals_A
matPval_iChIPData[rownames(dds_resultsTableA),"A_ImpulseDE"] <- pvals_A
matRunTime_iChIPData["ImpulseDE"] <- round(tm_DESeq2A["elapsed"])

lsResDEcomparison_iChIPdataH3K4me1 <- list("matQval_iChIPData"=matQval_iChIPData, 
  "matPval_iChIPData"=matPval_iChIPData, 
  "matRunTime_iChIPData"=matRunTime_iChIPData)

save(lsResDEcomparison_iChIPdataH3K4me1,file=file.path("/data/yosef2/users/fischerd/data/iChIP_Friedman/comparative_analysis/output/lsResDEcomparison_iChIPdataH3K4me1_EryLineage.RData"))
print("Finished DESeq2.")