### Clear
rm(list = ls())

### LOAD DATA
### Load annotation
# With control
dfAnnotationCtrl <- read.table("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/building/software_test/annotation_table_RNAseq_D50-D51-D54_D50control.tab",header=T)
rownames(dfAnnotationCtrl) <- dfAnnotationCtrl$Sample
# Without control
dfAnnotationCase <- read.table("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/building/software_test/annotation_table_RNAseq_D50-D51-D54.tab",header=T)
rownames(dfAnnotationCase) <- dfAnnotationCase$Sample

### Load count data
dfCountData <- read.table("/Users/davidsebastianfischer/MasterThesis/data/DESeq_RNAseq/raw_RNAseq/hDC_TPM_D50-D51-D54.tab",sep="\t",header=T)

### Process count data
matCountData <- dfCountData[,2:dim(dfCountData)[2]]
matCountData <- t(apply(matCountData,1,as.numeric))

rownames(matCountData) <- dfCountData[,1]
colnames(matCountData) <- colnames(dfCountData)[2:dim(dfCountData)[2]]

# Only chose high counts
#expression_table <- expression_table[apply(expression_table,1,function(gene){max(gene,na.rm=TRUE)}>500),]

source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/building/code_files/ImpulseDE2_main.R")

setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/software_test_out")
strCaseName = "case"
strControlName = NULL
if(!is.null(strControlName)){dfAnnotation=dfAnnotationCtrl
}else{dfAnnotation=dfAnnotationCase}
n_process = 3
Q_value = 10^(-2)
strMode <- "batch"
boolPlotting <- TRUE
lsImpulseDE_results <- runImpulseDE2(
  matCountData=matCountData, 
  dfAnnotation=dfAnnotation,
  strCaseName = strCaseName, 
  strControlName=strControlName, 
  strMode=strMode,
  nProc=n_process, 
  Q_value=Q_value, 
  boolPlotting=boolPlotting,
  scaSmallRun=1000)

lsDEGenes <- lsImpulseDE_results$lsDEGenes
dfImpulseResults <- lsImpulseDE_results$dfImpulseResults
lsImpulseFits <- lsImpulseDE_results$lsImpulseFits
dfDESeq2Results <- lsImpulseDE_results$dfDESeq2Results

# test plotting to change minor aspects
if(FALSE){
  # Load files from interior of ImpulseDE
  setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
  load("ImpulseDE2_matCountDataProc.RData")
  load("ImpulseDE2_dfAnnotationProc.RData")
  # Load Impulse output
  load("ImpulseDE2_dfImpulseResults.RData")
  load("ImpulseDE2_vecDEGenes.RData")
  load("ImpulseDE2_lsImpulseFits.RData")
  load("ImpulseDE2_dfDESeq2Results.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/ImpulseDE2_lsMatTranslationFactors.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/ImpulseDE2_matSizeFactors.RData")
  NPARAM <- 6
  setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files")
  source("srcImpulseDE2_plotDEGenes.R")
  setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
  plotDEGenes(
    vecGeneIDs=rownames(matCountDataProc)[1:16],
    matCountDataProc=matCountDataProc,
    lsMatTranslationFactors=lsMatTranslationFactors,
    matSizeFactors=matSizeFactors,
    dfAnnotation=dfAnnotation, 
    lsImpulseFits=lsImpulseFits,
    strCaseName=strCaseName, 
    strControlName=strControlName, 
    strFileNameSuffix="DE", 
    strPlotTitleSuffix="", 
    strPlotSubtitle="",
    dfImpulseResults=dfImpulseResults,
    vecRefPval=dfDESeq2Results$padj,
    strMode=strMode, 
    NPARAM=NPARAM)
}
