### Clear
rm(list = ls())

### LOAD DATA
### Load count data
# All genes
setwd("/Users/davidsebastianfischer/MasterThesis/data/Pseudotime")
expression_table_raw <- read.table("olfactoryTPMfilt_counts.txt",sep="\t",header=T)

### Process count data
expression_table <- t(apply(expression_table_raw,1,as.numeric))

rownames(expression_table) <- rownames(expression_table_raw)
colnames(expression_table) <- colnames(expression_table_raw)

matCounts <- expression_table

# Load Pseudotime data p63
load("/Users/davidsebastianfischer/MasterThesis/data/Pseudotime/EWT_ordered.Rdata")

# Investigate distribution of cells over pseudotime
vecPseudotime <- metadata$Pseudotime.1
names(vecPseudotime) <- metadata$sample
vecPseudotime <- vecPseudotime[!is.na(vecPseudotime)]

sum(colnames(matCounts) %in% names(vecPseudotime))
vecPseudotime <- vecPseudotime[names(vecPseudotime) %in% colnames(matCounts)]
matCounts <- matCounts[,names(vecPseudotime)]

library(compiler)
library(parallel)
library(DESeq2)

source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/srcImpulseDE2_runDESeq2.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/clusterCellsInPseudotime.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/formatDataClusters.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/fitHurdleModel.R")

source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/srcImpulseDE2_CostFunctionsFit.R")
evalLogLikHurdleNB_comp <- cmpfun(evalLogLikHurdleNB)
evalLogLikHurdleDrop_comp <- cmpfun(evalLogLikHurdleDrop)
evalLogLikHurdle_comp <- cmpfun(evalLogLikHurdle)

setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")

matCounts <- matCounts[apply(matCounts,1,function(gene){!all(is.na(gene) | gene==0)}),]

lsResultsClustering <- clusterCellsInPseudotime(vecPseudotime=vecPseudotime)
dfAnnotationClusters <- formatDataClusters(matCounts=matCounts,
  vecPseudotime=vecPseudotime,lsResultsClustering=lsResultsClustering)

lsHurdleParamters <- fitHurdleModel(matCounts=matCounts, 
  lsResultsClustering=lsResultsClustering, dfAnnotationClusters=dfAnnotationClusters)

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
    NPARAM=NPARAM)
}
