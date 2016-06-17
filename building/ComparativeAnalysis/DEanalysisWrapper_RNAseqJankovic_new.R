########################################
# 1. Plot ImpulseDE2 DE genes
########################################
rm(list=ls())

source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/ImpulseDE2_main.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/srcImpulseDE2_plotDEGenes.R")
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/ImpulseDE2")
load("ImpulseDE2_matCountDataProc.RData")
load("ImpulseDE2_dfAnnotationProc.RData")
load("ImpulseDE2_dfImpulseResults.RData")
load("ImpulseDE2_vecDEGenes.RData")
load("ImpulseDE2_lsImpulseFits.RData")
load("ImpulseDE2_dfDESeq2Results.RData")
load("ImpulseDE2_lsMatTranslationFactors.RData")
load("ImpulseDE2_matSizeFactors.RData")
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/pdfs")

QplotImpulseDE2genes <- 10^(-3)
vecImpulse2_DEgenes <- as.vector(dfImpulseResults[dfImpulseResults$adj.p <= QplotImpulseDE2genes,]$Gene)
strCaseName <- "case"
strControlName <- NULL
strMode <- "batch"
if(length(vecImpulse2_DEgenes)==0){
  vecImpulse2_DEgenes <- rownames(matCountDataProc[1:32,])
}
vecRefResults <- dfDESeq2Results$padj
names(vecRefResults) <- rownames(dfDESeq2Results)
plotDEGenes(vecGeneIDs=vecImpulse2_DEgenes,
  matCountDataProc=matCountDataProc,
  lsMatTranslationFactors=lsMatTranslationFactors,
  matSizeFactors=matSizeFactors,
  dfAnnotation=dfAnnotationProc, 
  lsImpulseFits=lsImpulseFits,
  strCaseName=strCaseName, 
  strControlName=strControlName, 
  strFileNameSuffix="DEgenes", 
  strPlotTitleSuffix="", 
  strPlotSubtitle="",
  dfImpulseResults=dfImpulseResults,
  vecRefPval=vecRefResults,
  strMode=strMode, 
  NPARAM=6)

########################################
# 2. Graphical comparison
########################################
rm(list=ls())
library(gplots)
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/lsResDEcomparison_RNAseqdata.Rdata")
matQval_RNAseqData <- lsResDEcomparison_RNAseqdata$matQval_RNAseqData 
matRunTime_RNAseqData <- lsResDEcomparison_RNAseqdata$matRunTime_RNAseqData

source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/ImpulseDE2_main.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/srcImpulseDE2_plotDEGenes.R")
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/ImpulseDE2")
load("ImpulseDE2_matCountDataProc.RData")
load("ImpulseDE2_dfAnnotationProc.RData")
load("ImpulseDE2_dfImpulseResults.RData")
load("ImpulseDE2_vecDEGenes.RData")
load("ImpulseDE2_lsImpulseFits.RData")
load("ImpulseDE2_dfDESeq2Results.RData")
load("ImpulseDE2_lsMatTranslationFactors.RData")
load("ImpulseDE2_matSizeFactors.RData")

setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/pdfs")

Q <- 10^(-3)
Qdelta <- 10^(2) # difference factor required to be plotted

compareDEMethods <- function(matQval,
  strMethod2="DESeq2",
  Q = Q,
  Qdelta = Qdelta,
  matCountDataProc = matCountDataProc,
  lsMatTranslationFactors = lsMatTranslationFactors,
  matSizeFactors = matSizeFactors,
  dfAnnotationProc = dfAnnotationProc, 
  lsImpulseFits = lsImpulseFits,
  dfImpulseResults = dfImpulseResults,
  strCaseName="case", 
  strControlName = NULL, 
  strMode="batch",
  strDataDescriptionFilename="")