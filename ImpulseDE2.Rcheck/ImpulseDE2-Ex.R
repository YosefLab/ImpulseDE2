pkgname <- "ImpulseDE2"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('ImpulseDE2')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("computeNormConst")
### * computeNormConst

flush(stderr()); flush(stdout())

### Name: computeNormConst
### Title: Compute a normalisation constant for each sample
### Aliases: computeNormConst

### ** Examples

lsSimulatedData <- simulateDataSetImpulseDE2(
vecTimePointsA   = rep(seq(1,8),3),
vecTimePointsB   = NULL,
vecBatchesA      = NULL,
vecBatchesB      = NULL,
scaNConst        = 100,
scaNImp          = 200,
scaNLin          = 100,
scaNSig          = 200)
vecSizeFactors <- computeNormConst(
matCountData = lsSimulatedData$matObservedCounts)




cleanEx()
nameEx("fitSigmoidModels")
### * fitSigmoidModels

flush(stderr()); flush(stdout())

### Name: fitSigmoidModels
### Title: Fits sigmoidal models to all genes on all all samples of a
###   condition
### Aliases: fitSigmoidModels

### ** Examples

lsSimulatedData <- simulateDataSetImpulseDE2(
vecTimePointsA   = rep(seq(1,8),3),
vecTimePointsB   = NULL,
vecBatchesA      = NULL,
vecBatchesB      = NULL,
scaNConst        = 0,
scaNImp          = 20,
scaNLin          = 10,
scaNSig          = 20)
objectImpulseDE2 <- runImpulseDE2(
matCountData    = lsSimulatedData$matObservedCounts, 
dfAnnotation    = lsSimulatedData$dfAnnotation,
boolCaseCtrl    = FALSE,
vecConfounders  = NULL,
boolIdentifyTransients = FALSE,
scaNProc        = 1 )
# You could have used boolIdentifyTransients=TRUE
# to avoid the following post wrapper fitting.
objectImpulseDE2 <- fitSigmoidModels(
objectImpulseDE2 = objectImpulseDE2,
vecConfounders   = NULL,
strCondition     = "case")
objectImpulseDE2 <- updateDEAnalysis(
objectImpulseDE2=objectImpulseDE2,
scaQThresTransients=0.001)
head(objectImpulseDE2$dfImpulseDE2Results)
# dfImpulseDE2Results now contain 'transients-analysis'.
   



cleanEx()
nameEx("get_accessors")
### * get_accessors

flush(stderr()); flush(stdout())

### Name: get_accessors
### Title: ImpulseDE2Object accession methods
### Aliases: get_accessors get_lsModelFits,ImpulseDE2Object-method
###   get_matCountDataProc,ImpulseDE2Object-method
###   get_dfAnnotationProc,ImpulseDE2Object-method
###   get_vecSizeFactors,ImpulseDE2Object-method
###   get_vecDispersions,ImpulseDE2Object-method
###   get_boolCaseCtrl,ImpulseDE2Object-method
###   get_vecConfounders,ImpulseDE2Object-method
###   get_scaNProc,ImpulseDE2Object-method
###   get_scaQThres,ImpulseDE2Object-method
###   get_strReport,ImpulseDE2Object-method get_accessors get_accessors
###   get_accessors get_accessors get_accessors get_accessors get_accessors
###   get_accessors get_accessors get_accessors

### ** Examples

   
lsSimulatedData <- simulateDataSetImpulseDE2(
vecTimePointsA   = rep(seq(1,8),3),
vecTimePointsB   = NULL,
vecBatchesA      = NULL,
vecBatchesB      = NULL,
scaNConst        = 30,
scaNImp          = 10,
scaNLin          = 10,
scaNSig          = 10)
objectImpulseDE2 <- runImpulseDE2(
matCountData    = lsSimulatedData$matObservedCounts, 
dfAnnotation    = lsSimulatedData$dfAnnotation,
boolCaseCtrl    = FALSE,
vecConfounders  = NULL,
scaNProc        = 1 )
# Extract hidden auxillary result and processed input objects.
lsModelFits <- get_lsModelFits(objectImpulseDE2)
matCountDataProc <- get_matCountDataProc(objectImpulseDE2)
dfAnnotationProc <- get_dfAnnotationProc(objectImpulseDE2)
vecSizeFactors <- get_vecSizeFactors(objectImpulseDE2)
vecDispersions <- get_vecDispersions(objectImpulseDE2)
boolCaseCtrl <- get_boolCaseCtrl(objectImpulseDE2)
vecConfounders <- get_vecConfounders(objectImpulseDE2)
scaNProc <- get_scaNProc(objectImpulseDE2)
scaQThres <- get_scaQThres(objectImpulseDE2)
strReport <- get_strReport(objectImpulseDE2)




cleanEx()
nameEx("list_accession")
### * list_accession

flush(stderr()); flush(stdout())

### Name: list_accession
### Title: List-like accessor methods for ImpulseDE2Object
### Aliases: list_accession names.ImpulseDE2Object
###   names,ImpulseDE2Object-method $,ImpulseDE2Object-method
###   [[,ImpulseDE2Object,character,missing-method list_accession
###   list_accession list_accession

### ** Examples

lsSimulatedData <- simulateDataSetImpulseDE2(
vecTimePointsA   = rep(seq(1,8),3),
vecTimePointsB   = NULL,
vecBatchesA      = NULL,
vecBatchesB      = NULL,
scaNConst        = 30,
scaNImp          = 10,
scaNLin          = 10,
scaNSig          = 10)
objectImpulseDE2 <- runImpulseDE2(
matCountData    = lsSimulatedData$matObservedCounts, 
dfAnnotation    = lsSimulatedData$dfAnnotation,
boolCaseCtrl    = FALSE,
vecConfounders  = NULL,
scaNProc        = 1 )
names(objectImpulseDE2) # Display core output
# With respect to this core output, objectImpulseDE2
# can be treated like a list.
head(objectImpulseDE2[["dfImpulseDE2Results"]])
head(objectImpulseDE2$dfImpulseDE2Results)
head(objectImpulseDE2[["vecDEGenes"]])
head(objectImpulseDE2$vecDEGenes)




cleanEx()
nameEx("plotGenes")
### * plotGenes

flush(stderr()); flush(stdout())

### Name: plotGenes
### Title: Plots the impulse fits and data
### Aliases: plotGenes

### ** Examples

lsSimulatedData <- simulateDataSetImpulseDE2(
vecTimePointsA   = rep(seq(1,8),3),
vecTimePointsB   = NULL,
vecBatchesA      = NULL,
vecBatchesB      = NULL,
scaNConst        = 0,
scaNImp          = 40,
scaNLin          = 20,
scaNSig          = 40)
objectImpulseDE2 <- runImpulseDE2(
matCountData    = lsSimulatedData$matObservedCounts, 
dfAnnotation    = lsSimulatedData$dfAnnotation,
boolCaseCtrl    = FALSE,
vecConfounders  = NULL,
boolIdentifyTransients = FALSE,
scaNProc        = 1 )
lsgplotsID <- plotGenes(
scaNTopIDs=5,
objectImpulseDE2=objectImpulseDE2,
boolCaseCtrl=FALSE,
boolMultiplePlotsPerPage=TRUE,
boolSimplePlot=FALSE)
lsgplotsID[[1]]




cleanEx()
nameEx("plotHeatmap")
### * plotHeatmap

flush(stderr()); flush(stdout())

### Name: plotHeatmap
### Title: Plot structured z-value heatmaps of differentially expressed
###   genes
### Aliases: plotHeatmap

### ** Examples

library(ComplexHeatmap)
lsSimulatedData <- simulateDataSetImpulseDE2(
vecTimePointsA   = rep(seq(1,8),3),
vecTimePointsB   = NULL,
vecBatchesA      = NULL,
vecBatchesB      = NULL,
scaNConst        = 0,
scaNImp          = 50,
scaNLin          = 0,
scaNSig          = 50)
objectImpulseDE2 <- runImpulseDE2(
matCountData    = lsSimulatedData$matObservedCounts, 
dfAnnotation    = lsSimulatedData$dfAnnotation,
boolCaseCtrl    = FALSE,
vecConfounders  = NULL,
boolIdentifyTransients = TRUE,
scaNProc        = 1 )
lsHeatmaps <- plotHeatmap(
objectImpulseDE2=objectImpulseDE2,
strCondition="case",
boolIdentifyTransients=TRUE,
scaQThres=0.01)
draw(lsHeatmaps$complexHeatmapRaw)




cleanEx()
nameEx("runImpulseDE2")
### * runImpulseDE2

flush(stderr()); flush(stdout())

### Name: runImpulseDE2
### Title: ImpulseDE2 wrapper
### Aliases: runImpulseDE2 ImpulseDE2 wrapper

### ** Examples

lsSimulatedData <- simulateDataSetImpulseDE2(
vecTimePointsA   = rep(seq(1,8),3),
vecTimePointsB   = NULL,
vecBatchesA      = NULL,
vecBatchesB      = NULL,
scaNConst        = 30,
scaNImp          = 10,
scaNLin          = 10,
scaNSig          = 10)
objectImpulseDE2 <- runImpulseDE2(
matCountData    = lsSimulatedData$matObservedCounts, 
dfAnnotation    = lsSimulatedData$dfAnnotation,
boolCaseCtrl    = FALSE,
vecConfounders  = NULL,
scaNProc        = 1 )
head(objectImpulseDE2$dfImpulseDE2Results)




cleanEx()
nameEx("simulateDataSetImpulseDE2")
### * simulateDataSetImpulseDE2

flush(stderr()); flush(stdout())

### Name: simulateDataSetImpulseDE2
### Title: Simulate a data set for ImpulseDE2
### Aliases: simulateDataSetImpulseDE2

### ** Examples

lsSimulatedData <- simulateDataSetImpulseDE2(
vecTimePointsA   = rep(seq(1,8),3),
vecTimePointsB   = NULL,
vecBatchesA      = NULL,
vecBatchesB      = NULL,
scaNConst        = 30,
scaNImp          = 10,
scaNLin          = 10,
scaNSig          = 10)
   



cleanEx()
nameEx("updateDEAnalysis")
### * updateDEAnalysis

flush(stderr()); flush(stdout())

### Name: updateDEAnalysis
### Title: Update dfImpulseDE2Results after sigmoids have been fit through
###   external call
### Aliases: updateDEAnalysis

### ** Examples

lsSimulatedData <- simulateDataSetImpulseDE2(
vecTimePointsA   = rep(seq(1,8),3),
vecTimePointsB   = NULL,
vecBatchesA      = NULL,
vecBatchesB      = NULL,
scaNConst        = 0,
scaNImp          = 50,
scaNLin          = 0,
scaNSig          = 50)
objectImpulseDE2 <- runImpulseDE2(
matCountData    = lsSimulatedData$matObservedCounts, 
dfAnnotation    = lsSimulatedData$dfAnnotation,
boolCaseCtrl    = FALSE,
vecConfounders  = NULL,
boolIdentifyTransients = FALSE,
scaNProc        = 1 )
# You could have used boolIdentifyTransients=TRUE
# to avoid the following post wrapper fitting.
objectImpulseDE2 <- fitSigmoidModels(
objectImpulseDE2 = objectImpulseDE2,
vecConfounders   = NULL,
strCondition     = "case")
objectImpulseDE2 <- updateDEAnalysis(
objectImpulseDE2=objectImpulseDE2,
scaQThresTransients=0.001)
head(objectImpulseDE2$dfImpulseDE2Results)
# dfImpulseDE2Results now contain 'transients-analysis'.
   



cleanEx()
nameEx("writeReportToFile")
### * writeReportToFile

flush(stderr()); flush(stdout())

### Name: writeReportToFile
### Title: Print ImpulseDE2 report to .txt file
### Aliases: writeReportToFile

### ** Examples

dirPWD <- getwd() # Will save into current working directory.
lsSimulatedData <- simulateDataSetImpulseDE2(
vecTimePointsA   = rep(seq(1,8),3),
vecTimePointsB   = NULL,
vecBatchesA      = NULL,
vecBatchesB      = NULL,
scaNConst        = 30,
scaNImp          = 10,
scaNLin          = 10,
scaNSig          = 10)
objectImpulseDE2 <- runImpulseDE2(
matCountData    = lsSimulatedData$matObservedCounts, 
dfAnnotation    = lsSimulatedData$dfAnnotation,
boolCaseCtrl    = FALSE,
vecConfounders  = NULL,
scaNProc        = 1 )
# Uncomment to run:
#writeReportToFile(
#object=objectImpulseDE2,
#file=paste0(dirPWD, "ImpulseDE2Report.txt")
#)




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
