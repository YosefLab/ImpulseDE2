### Clear
rm(list = ls())

### LOAD DATA
### Load annotation
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/software_test")
dfAnnotationFull <- read.table("annotation_table_ATACseq_E70-E72-E81-E91_E81-E91control.tab",header=T)
rownames(dfAnnotationFull) <- dfAnnotationFull$Replicate

### Load count data
# All genes
setwd("/Users/davidsebastianfischer/MasterThesis/data/ATACseq_analysis")
expression_table_raw <- read.table("E70-E72-E81-E91_peak_counts_DESeq_E72_24h_LPS_downsampled.tab",
  sep="\t",header=T,colClasses=c(
    "character","character","character",
    "numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric"))

### Process count data
expression_table_raw <- expression_table_raw[232979:(232979+200),]
expression_table <- expression_table_raw[,4:dim(expression_table_raw)[2]]
expression_table <- t(apply(expression_table,1,as.numeric))
lsPeaksIDs <- apply(expression_table_raw[,1:3],1,function(coord){paste(coord,collapse="-")})
# Rename elements for shaked's poster: hide coordinates
lsPeaksIDs <- paste0(rep("Region_",length(lsPeaksIDs)),c(1:length(lsPeaksIDs)))

rownames(expression_table) <- lsPeaksIDs
colnames(expression_table) <- colnames(expression_table_raw)[4:dim(expression_table_raw)[2]]

# Only chose high counts
#expression_table <- expression_table[apply(expression_table,1,function(gene){max(gene,na.rm=TRUE)}>500),]

setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files")
source("ImpulseDE2_main.R")
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")

# mix up expression table to check indexing works
expression_table<-expression_table[,1:10,15,20,19,18,17,16,11:14]

strControlName = "ctrl"
strCaseName = "case"
n_process = 3
Q_value = 10^(-2)
strMode <- "batch"
#strMode <- "timecourses"
print("-----------------dont include too few samples here!")
lsImpulseDE_results <- runImpulseDE2(matCountData=expression_table, dfAnnotationFull=dfAnnotationFull,
  strCaseName = strCaseName, strControlName=strControlName, strMode=strMode,
  nProc=n_process, Q_value=Q_value)

lsDEGenes <- lsImpulseDE_results$lsDEGenes
dfImpulseResults <- lsImpulseDE_results$dfImpulseResults
lsImpulseFits <- lsImpulseDE_results$lsImpulseFits
dfDESeq2Results <- lsImpulseDE_results$dfDESeq2Results

# test plotting to change minor aspects
if(FALSE){
  # Load files from interior of ImpulseDE
  setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
  load("ImpulseDE2_arr2DCountData.RData")
  load("ImpulseDE2_dfAnnotationFull.RData")
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
    lsGeneIDs=rownames(expression_table)[1:10],
    arr2DCountData=arr2DCountData,
    vecNormConst=vecNormConst,
    dfAnnotationFull=dfAnnotationFull, 
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