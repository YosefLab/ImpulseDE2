### Clear
rm(list = ls())

### LOAD DATA
### Load annotation
# With control
dfAnnotationCtrl <- read.table("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/software_test/annotation_table_ATACseq_E70-E72-E81-E91_E81-E91control.tab",header=T)
rownames(dfAnnotationCase) <- dfAnnotationCase$Sample
# Without control
dfAnnotationCase <- read.table("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/software_test/annotation_table_ATACseq_E70-E72-E81-E91.tab",header=T)
rownames(dfAnnotationCtrl) <- dfAnnotationCtrl$Sample

### Load count data
# All genes
setwd("/Users/davidsebastianfischer/MasterThesis/data/ATACseq_analysis")
dfCountData <- read.table("E70-E72-E81-E91_peak_counts_DESeq_E72_24h_LPS_downsampled.tab",
  sep="\t",header=T,colClasses=c(
    "character","character","character",
    "numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric"))

### Process count data
dfCountData <- dfCountData
matCountData <- dfCountData[,4:dim(dfCountData)[2]]
matCountData <- t(apply(matCountData,1,as.numeric))
lsPeaksIDs <- apply(dfCountData[,1:3],1,function(coord){paste(coord,collapse="-")})
# Rename elements for shaked's poster: hide coordinates
lsPeaksIDs <- paste0(rep("Region_",length(lsPeaksIDs)),c(1:length(lsPeaksIDs)))

rownames(matCountData) <- lsPeaksIDs
colnames(matCountData) <- colnames(dfCountData)[4:dim(dfCountData)[2]]

source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/ImpulseDE2_main.R")
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")

strCaseName = "case"
strControlName = "ctrl"
if(!is.null(strControlName)){dfAnnotation=dfAnnotationCtrl
}else{dfAnnotation=dfAnnotationCase}

n_process = 3
Q_value = 10^(-2)
strMode <- "batch"
lsImpulseDE_results <- runImpulseDE2(
  matCountData=matCountData, 
  dfAnnotation=dfAnnotation,
  strCaseName = strCaseName, 
  strControlName=strControlName, 
  strMode=strMode,
  nProc=n_process,
  Q_value=Q_value)

lsDEGenes <- lsImpulseDE_results$lsDEGenes
dfImpulseResults <- lsImpulseDE_results$dfImpulseResults
lsImpulseFits <- lsImpulseDE_results$lsImpulseFits
dfDESeq2Results <- lsImpulseDE_results$dfDESeq2Results

# test plotting to change minor aspects
if(FALSE){
  # Load files from interior of ImpulseDE
  setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
  load("ImpulseDE2_matCountData.RData")
  load("ImpulseDE2_dfAnnotation.RData")
  # Load Impulse output
  load("ImpulseDE2_dfImpulseResults.RData")
  load("ImpulseDE2_lsDEGenes.RData")
  load("ImpulseDE2_lsImpulseFits.RData")
  load("ImpulseDE2_dfDESeq2Results.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/ImpulseDE2_vecNormConst.RData")
  NPARAM <- 6
  setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files")
  source("srcImpulseDE2_plotDEGenes.R")
  setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
  plotDEGenes(
    lsGeneIDs=rownames(matCountData)[1:4],
    matCountData=matCountData,
    vecNormConst=vecNormConst,
    dfAnnotation=dfAnnotation, 
    lsImpulseFits=lsImpulseFits,
    strCaseName=strCaseName, 
    strControlName=strControlName, 
    strFileNameSuffix="DE", 
    strPlotTitleSuffix="", 
    strPlotSubtitle="",
    dfImpulseResults=dfImpulseResults,
    vecMethod2Results=dfDESeq2Results$padj,
    strMode=strMode, 
    NPARAM=NPARAM)
}