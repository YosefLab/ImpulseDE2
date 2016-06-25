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
matDataA_ImpulseDE2 <- matDataAA

tm_ImpulseDE2A <- system.time({
  setwd("/data/yosef2/users/fischerd/data/iChIP_Friedman/comparative_analysis")
  strControlName = NULL
  strCaseName = "case"
  lsImpulseDE_resultsA <- runImpulseDE2(matCountData=matDataA_ImpulseDE2, 
    dfAnnotation=dfAnnotationA,
    strCaseName = strCaseName, 
    strControlName=strControlName, 
    strMode="batch",
    nProc=NCORES, 
    Q_value=10^(-3),
    boolPlotting=TRUE)
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

save(lsResDEcomparison_iChIPdataH3K4me1,file=file.path("/data/yosef2/users/fischerd/data/iChIP_Friedman/comparative_analysis/lsResDEcomparison_iChIPdataH3K4me1_EryLineage.RData"))
print("Finished ImpulseDE2")