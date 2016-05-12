### Clear
rm(list = ls())

### LOAD DATA
### Load annotation
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/software_test")
#dfAnnotationFull <- read.table("annotation_table_RNAseq_D50-D51-D54_missing.tab",header=T)
dfAnnotationFull <- read.table("annotation_table_RNAseq_D50-D51-D54.tab",header=T)
rownames(dfAnnotationFull) <- dfAnnotationFull$Replicate

### Load count data
# DE Genes (as of DESeq)
#setwd("/Users/davidsebastianfischer/MasterThesis/data/DESeq_RNAseq/DESeq_output")
#expression_table_raw <- read.table("D50-D51-D54_DEseq_RNAseq_filteredGenes_TPM_p1e-05.tab",sep="\t",header=T)
# All genes
setwd("/Users/davidsebastianfischer/MasterThesis/data/DESeq_RNAseq/raw_RNAseq")
expression_table_raw <- read.table("hDC_TPM_D50-D51-D54.tab",sep="\t",header=T)

### Process count data
expression_table <- expression_table_raw[,2:dim(expression_table_raw)[2]]
expression_table <- t(apply(expression_table,1,as.numeric))

rownames(expression_table) <- expression_table_raw[,1]
colnames(expression_table) <- colnames(expression_table_raw)[2:dim(expression_table_raw)[2]]

# Only chose high counts
#expression_table <- expression_table[apply(expression_table,1,function(gene){max(gene,na.rm=TRUE)}>500),]

# replace by simulation data
# expression_table_cut[,] <- rnbinom(dim(expression_table_cut)[1]*dim(expression_table_cut)[2],
#  size=2,mu=100)

setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files")
source("ImpulseDE2_main.R")
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")


control_timecourse = FALSE
strControlName = NULL
strCaseName = "case"
n_process = 3
Q_value = 10^(-2)
strMode <- "timecourses"
#strMode <- "batch"
boolPlotting <- TRUE
lsImpulseDE_results <- runImpulseDE2(matCountData=expression_table, dfAnnotationFull=dfAnnotationFull,
  strCaseName = strCaseName, strControlName=strControlName, strMode=strMode,
  nProc=n_process, Q_value=Q_value, boolPlotting=boolPlotting)

lsDEGenes <- lsImpulseDE_results$lsDEGenes
dfImpulseResults <- lsImpulseDE_results$dfImpulseResults
lsImpulseFits <- lsImpulseDE_results$lsImpulseFits
dfDESeq2Results <- lsImpulseDE_results$dfDESeq2Results

# test plotting to change minor aspects
if(FALSE){
  # Load files from interior of ImpulseDE
  setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
  load("ImpulseDE2_arr3DCountData.RData")
  load("ImpulseDE2_dfAnnotationRed.RData")
  # Load Impulse output
  load("ImpulseDE2_dfImpulseResults.RData")
  load("ImpulseDE2_lsDEGenes.RData")
  load("ImpulseDE2_lsImpulseFits.RData")
  load("ImpulseDE2_dfDESeq2Results.RData")
  NPARAM <- 6
  setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files")
  source("srcImpulseDE2_plotDEGenes.R")
  setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
  plotDEGenes(lsGeneIDs=rownames(expression_table)[1:100],
    arr3DCountData=arr3DCountData, dfAnnotationRed=dfAnnotationRed, 
    lsImpulseFits=lsImpulseFits,
    strCaseName=strCaseName, strControlName=strControlName, 
    strFileNameSuffix="DE", strPlotTitleSuffix="", strPlotSubtitle="",
    dfImpulseResults=dfImpulseResults,dfDESeq2Results=dfDESeq2Results,
    NPARAM=NPARAM,boolPlotting=TRUE)
}
