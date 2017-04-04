## ----options, cache=FALSE, include=FALSE, results='hide', message=FALSE, warning=FALSE----

knitr::opts_chunk$set(fig.align="center", 
                      cache=FALSE,
                      error=FALSE,
                      fig.width=6,fig.height=6,
                      autodep=TRUE,
                      out.width="600px", 
                      out.height="600px",
                      message=FALSE,
                      warning=FALSE,
                      results="hide", 
                      echo=TRUE, 
                      eval=TRUE)

options(getClass.msg=FALSE)

## ----case-only-----------------------------------------------------------
library(ImpulseDE2)
lsSimulatedData <- simulateDataSetImpulseDE2(
  vecTimePointsA   = rep(seq(1,8),3),
  vecTimePointsB   = NULL,
  vecBatchesA      = NULL,
  vecBatchesB      = NULL,
  scaNConst        = 30,
  scaNImp          = 10,
  scaNLin          = 10,
  scaNSig          = 10,
  scaMuBatchEffect = NULL,
  scaSDBatchEffect = NULL,
  dirOutSimulation = NULL)

## ----case-only-annotation, results='markdown'----------------------------
lsSimulatedData$dfAnnotation

## ----case-only2----------------------------------------------------------
objectImpulseDE2 <- runImpulseDE2(
  matCountData    = lsSimulatedData$matObservedCounts, 
  dfAnnotation    = lsSimulatedData$dfAnnotation,
  boolCaseCtrl    = FALSE,
  vecConfounders  = NULL,
  scaNProc        = 1 )

## ----case-only-results, results='markdown'-------------------------------
head(objectImpulseDE2$dfImpulseDE2Results)

## ----case-only-batch-----------------------------------------------------
library(ImpulseDE2)
lsSimulatedData <- simulateDataSetImpulseDE2(
  vecTimePointsA   = rep(seq(1,8),3),
  vecTimePointsB   = NULL,
  vecBatchesA      = c(rep("B1",8), rep("B2",8), rep("B3",8)),
  vecBatchesB      = NULL,
  scaNConst        = 30,
  scaNImp          = 10,
  scaNLin          = 10,
  scaNSig          = 10,
  scaMuBatchEffect = 1,
  scaSDBatchEffect = 0.2,
  dirOutSimulation = NULL)

## ----case-only-batch-annotation, results='markdown'----------------------
lsSimulatedData$dfAnnotation

## ----case-only-batch2----------------------------------------------------
objectImpulseDE2 <- runImpulseDE2(
  matCountData    = lsSimulatedData$matObservedCounts, 
  dfAnnotation    = lsSimulatedData$dfAnnotation,
  boolCaseCtrl    = FALSE,
  vecConfounders  = c("Batch"),
  scaNProc        = 1 )

## ----case-only-batch-results, results='markdown'-------------------------
head(objectImpulseDE2$dfImpulseDE2Results)

## ----plot-genes----------------------------------------------------------
# Continue script of "Batch effects"
library(ggplot2)
lsgplotsGenes <- plotGenes(
  vecGeneIDs       = NULL,
  scaNTopIDs       = 10,
  objectImpulseDE2 = objectImpulseDE2,
  boolCaseCtrl     = FALSE,
  dirOut           = NULL,
  strFileName      = NULL,
  vecRefPval       = NULL, 
  strNameRefMethod = NULL)
print(lsgplotsGenes[[1]])

## ----case-control-batch--------------------------------------------------
lsSimulatedData <- simulateDataSetImpulseDE2(
  vecTimePointsA   = rep(seq(1,8),3),
  vecTimePointsB   = rep(seq(1,8),3),
  vecBatchesA      = c(rep("B1",8), rep("B2",8), rep("B3",8)),
  vecBatchesB      = c(rep("C1",8), rep("C2",8), rep("C3",8)),
  scaNConst        = 30,
  scaNImp          = 10,
  scaNLin          = 10,
  scaNSig          = 10,
  scaMuBatchEffect = 1,
  scaSDBatchEffect = 0.1,
  dirOutSimulation = NULL)

## ----case-control-annotation, results='markdown'-------------------------
lsSimulatedData$dfAnnotation

## ----case-control2-------------------------------------------------------
objectImpulseDE2 <- runImpulseDE2(
  matCountData    = lsSimulatedData$matObservedCounts, 
  dfAnnotation    = lsSimulatedData$dfAnnotation,
  boolCaseCtrl    = TRUE,
  vecConfounders  = c("Batch"),
  scaNProc        = 1 )

## ----case-control-results, results='markdown'----------------------------
head(objectImpulseDE2$dfImpulseDE2Results)

## ----transients----------------------------------------------------------
library(ImpulseDE2)
lsSimulatedData <- simulateDataSetImpulseDE2(
  vecTimePointsA   = rep(seq(1,8),3),
  vecTimePointsB   = NULL,
  vecBatchesA      = c(rep("B1",8), rep("B2",8), rep("B3",8)),
  vecBatchesB      = NULL,
  scaNConst        = 0,
  scaNImp          = 100,
  scaNLin          = 0,
  scaNSig          = 0,
  scaMuBatchEffect = 1,
  scaSDBatchEffect = 0.2,
  dirOutSimulation = NULL)
objectImpulseDE2 <- runImpulseDE2(
  matCountData           = lsSimulatedData$matObservedCounts, 
  dfAnnotation           = lsSimulatedData$dfAnnotation,
  boolCaseCtrl           = FALSE,
  vecConfounders         = c("Batch"),
  boolIdentifyTransients = TRUE,
  scaNProc               = 1 )

## ----transient-results, results='markdown'-------------------------------
head(objectImpulseDE2$dfImpulseDE2Results)

## ----heatmap-------------------------------------------------------------
# Continuing script of "Transiently regulated genes"
library(ComplexHeatmap)
lsHeatmaps <- plotHeatmap(
  objectImpulseDE2       = objectImpulseDE2,
  strCondition           = "case",
  boolIdentifyTransients = TRUE,
  scaQThres              = 0.01)
draw(lsHeatmaps$complexHeatmapRaw) # Heatmap based on normalised counts

## ----session-------------------------------------------------------------
sessionInfo()

