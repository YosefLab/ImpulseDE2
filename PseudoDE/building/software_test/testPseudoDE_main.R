rm(list = ls())

# Load Pseudotime data p63
load("/Users/davidsebastianfischer/MasterThesis/data/Pseudotime/p63_ordered.Rdata")
dfCountsE4 <- read.table("/Users/davidsebastianfischer/MasterThesis/data/Pseudotime/E4_Tophatfilt_counts.txt", header=TRUE, row.names=1)
# Make sure got same cells
sum(colnames(dfCountsE4) %in% metadata$sample)
length(colnames(dfCountsE4))
length(metadata$sample)

# Investigate distribution of cells over pseudotime
lsPTpointsAll_p63PT1 <- as.vector(metadata$Pseudotime.1)
names(lsPTpointsAll_p63PT1) <- as.vector(metadata$sample)
lsPTpoints_p63PT1 <- lsPTpointsAll_p63PT1[!is.na(lsPTpointsAll_p63PT1)]
print(paste0("p63PT1: Total cells: ",length(lsPTpointsAll_p63PT1),", Non NA: ",length(lsPTpoints_p63PT1)))

matCounts <- data.matrix(dfCountsE4)
lsPTpoints_p63PT1 <- lsPTpoints_p63PT1[!is.na(lsPTpoints_p63PT1)]
matCounts <- matCounts[,colnames(matCounts) %in% names(lsPTpoints_p63PT1)]
matCounts <- matCounts[1:500,]
#matCounts <- round(counts)
nProc=2
vecPseudotime <- lsPTpoints_p63PT1[colnames(matCounts)]
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/PseudoDE_main.R")
setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
lsDEresults <- runPseudoDE(matCounts=matCounts,vecPseudotime=vecPseudotime)

if(FALSE){
  vecboolNonzeroGenes <- apply(matCounts,1,
    function(gene){mean(gene)>10})
  matCounts <- matCounts[vecboolNonzeroGenes,]
  MAXITER <- 3
  
  # Fit zinb model
  lsZINBparam <- estimate_zinb(
    Y = matCounts, 
    maxiter = MAXITER, 
    verbose = TRUE)
}

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
  strMode <- "singlecell"
  setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files")
  source("srcImpulseDE2_plotDEGenes.R")
  setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
  plotDEGenes(
    lsGeneIDs=rownames(matCountsClean),
    arr2DCountData=arr2DCountData,
    vecNormConst=vecNormConst,
    dfAnnotationFull=dfAnnotationFull, 
    lsImpulseFits=lsImpulseFits,
    strCaseName="case", 
    strControlName=NULL, 
    strFileNameSuffix="DE", 
    strPlotTitleSuffix="", 
    strPlotSubtitle="",
    dfImpulseResults=dfImpulseResults,
    vecMethod2Results=dfDESeq2Results$padj,
    strMode=strMode, 
    NPARAM=NPARAM)
}
if(FALSE){
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_dfAnnotation.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_vecDispersions.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_lsInputToImpulseDE2.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_lsResultsClustering.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_matCountsClean.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_lsZINBparam.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_matCountsClean.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_matClusterMeansFitted.RData")
  matDropout <- lsInputToImpulseDE2$matDropout
  matProbNB <- lsInputToImpulseDE2$matProbNB
  matCountsImputed <- lsInputToImpulseDE2$matCountsImputed
}
if(FALSE){
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_coefs_mu.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_fit_mu.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_Y.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_zhat.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_thetahat.RData")
  #print(coefs_mu)
  print(fit_mu)
  print(thetahat)
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_fit_pi.RData")
  print(fit_pi)
}
