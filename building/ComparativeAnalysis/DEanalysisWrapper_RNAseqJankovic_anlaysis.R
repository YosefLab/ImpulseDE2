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
## Full data set
rm(list=ls())
library(gplots)
#load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/edge/lsResDEcomparison_RNAseqdata.Rdata")
#matQval_edge <- lsResDEcomparison_RNAseqdata$matQval_RNAseqData 
#matRunTime_edge <- lsResDEcomparison_RNAseqdata$matRunTime_RNAseqData

setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/building/code_files")
source("ImpulseDE2_main.R")
source("srcImpulseDE2_plotDEGenes.R")
source("srcImpulseDE2_compareDEMethods.R")
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/clusterruns/output/ImpulseDE2")
load("ImpulseDE2_matCountDataProc.RData")
load("ImpulseDE2_dfAnnotationProc.RData")
load("ImpulseDE2_dfImpulseResults.RData")
load("ImpulseDE2_vecDEGenes.RData")
load("ImpulseDE2_lsImpulseFits.RData")
load("ImpulseDE2_dfDESeq2Results.RData")
load("ImpulseDE2_lsMatTranslationFactors.RData")
load("ImpulseDE2_matSizeFactors.RData")

Q <- 10^(-3)
Qdelta <- 10^(2) # difference factor required to be plotted

# Load run data ImpulseDE2, ImpulseDE, DESeq2
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/clusterruns/output/lsResDEcomparison_RNAseqdata.RData")
matQval <- lsResDEcomparison_RNAseqdata$matQval_RNAseqData
matRunTime <- lsResDEcomparison_RNAseqdata$matRunTime_RNAseqData

colnames(matQval) <- c("Gene", sapply(colnames(matQval),function(str){unlist(strsplit(str,"A_"))[2]})[2:5])
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/pdfs/full")

compareDEMethods(matQval,
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

compareDEMethods(matQval,
  strMethod2="ImpulseDE",
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

## Full data set witout 12h
rm(list=ls())
library(gplots)

setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/building/code_files")
source("ImpulseDE2_main.R")
source("srcImpulseDE2_plotDEGenes.R")
source("srcImpulseDE2_compareDEMethods.R")
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/clusterruns/output_no12h/ImpulseDE2")
load("ImpulseDE2_matCountDataProc.RData")
load("ImpulseDE2_dfAnnotationProc.RData")
load("ImpulseDE2_dfImpulseResults.RData")
load("ImpulseDE2_vecDEGenes.RData")
load("ImpulseDE2_lsImpulseFits.RData")
load("ImpulseDE2_dfDESeq2Results.RData")
load("ImpulseDE2_lsMatTranslationFactors.RData")
load("ImpulseDE2_matSizeFactors.RData")

Q <- 10^(-3)
Qdelta <- 10^(2) # difference factor required to be plotted

# Load run data ImpulseDE2, ImpulseDE, DESeq2
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/clusterruns/output_no12h/lsResDEcomparison_RNAseqdata.RData")
matQval <- lsResDEcomparison_RNAseqdata$matQval_RNAseqData
matRunTime <- lsResDEcomparison_RNAseqdata$matRunTime_RNAseqData

colnames(matQval) <- c("Gene", sapply(colnames(matQval),function(str){unlist(strsplit(str,"A_"))[2]})[2:5])
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/pdfs/no12h")

compareDEMethods(matQval,
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

compareDEMethods(matQval,
  strMethod2="ImpulseDE",
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