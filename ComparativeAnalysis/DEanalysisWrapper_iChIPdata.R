rm(list = ls())

print("Process data iChIP")
# Load Data set iChIP
# 1. Counts
rm(list=ls())

dfH3K27ac <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/iChIP_H3K27Ac.tab.txt",sep="\t")
matCellHeaders27ac <- sapply(dfH3K27ac[1,c(14:30)],function(cell){unlist(strsplit(as.character(cell),split=" "))})
vecCellNames27ac <- as.vector(matCellHeaders27ac[1,])
vecCols27ac <- unlist(list("PeakID",(dfH3K27ac[1,c(2:13)]),vecCellNames27ac))
dfH3K27ac <- dfH3K27ac[c(2:dim(dfH3K27ac)[1]),]
colnames(dfH3K27ac) <- vecCols27ac
rownames(dfH3K27ac) <- dfH3K27ac$PeakID
vecNormFac27ac <- 10000000/as.numeric(sapply(matCellHeaders27ac[7,],function(cell){unlist(strsplit(cell,"[(]"))[2]}))
matH3K27acProc <- t(apply(dfH3K27ac[,c(14:30)],1,as.numeric))
colnames(matH3K27acProc) <- vecCellNames27ac
matH3K27acRaw <- t(apply(matH3K27acProc,1,function(gene){gene/vecNormFac27ac}))
matH3K27acCounts <- round(matH3K27acRaw)

dfH3K4me1 <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/iChIP_H3K4me1.tab.txt",sep="\t")
matCellHeaders4me1 <- sapply(dfH3K4me1[1,c(14:30)],function(cell){unlist(strsplit(as.character(cell),split=" "))})
vecCellNames4me1 <- as.vector(matCellHeaders4me1[1,])
vecCols4me1 <- unlist(list("PeakID",(dfH3K4me1[1,c(2:13)]),vecCellNames4me1))
dfH3K4me1 <- dfH3K4me1[c(2:dim(dfH3K4me1)[1]),]
colnames(dfH3K4me1) <- vecCols4me1
rownames(dfH3K4me1) <- dfH3K4me1$PeakID
vecNormFac4me1 <- 10000000/as.numeric(sapply(matCellHeaders4me1[7,],function(cell){unlist(strsplit(cell,"[(]"))[2]}))
matH3K4me1Proc <- t(apply(dfH3K4me1[,c(14:30)],1,as.numeric))
colnames(matH3K4me1Proc) <- vecCellNames4me1
matH3K4me1Raw <- t(apply(matH3K4me1Proc,1,function(gene){gene/vecNormFac4me1}))
matH3K4me1Counts <- round(matH3K4me1Raw)

dfH3K4me2 <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/iChIP_H3K4me2.tab.txt",sep="\t")
# lacks CD4 and CD8
matCellHeaders4me2 <- sapply(dfH3K4me2[1,c(14:28)],function(cell){unlist(strsplit(as.character(cell),split=" "))})
vecCellNames4me2 <- as.vector(matCellHeaders4me2[1,])
vecCols4me2 <- unlist(list("PeakID",(dfH3K4me2[1,c(2:13)]),vecCellNames4me2))
dfH3K4me2 <- dfH3K4me2[c(2:dim(dfH3K4me2)[1]),]
colnames(dfH3K4me2) <- vecCols4me2
rownames(dfH3K4me2) <- dfH3K4me2$PeakID
vecNormFac4me2 <- 10000000/as.numeric(sapply(matCellHeaders4me2[7,],function(cell){unlist(strsplit(cell,"[(]"))[2]}))
matH3K4me2Proc <- t(apply(dfH3K4me2[,c(14:28)],1,as.numeric))
colnames(matH3K4me2Proc) <- vecCellNames4me2
matH3K4me2Raw <- t(apply(matH3K4me2Proc,1,function(gene){gene/vecNormFac4me2}))
matH3K4me2Counts <- round(matH3K4me2Raw)

# Reduce to target branch
vecBranchCells <- vecCellNames27ac[c(11,17,15,5,12,7,8)]
matH3K4me1CountsBranch <- matH3K4me1Counts[,c(vecBranchCells)]
matH3K4me2CountsBranch <- matH3K4me2Counts[,c(vecBranchCells)]
matH3K27acCountsBranch <- matH3K27acCounts[,c(vecBranchCells)]

matDataA <- matH3K4me1CountsBranch[1:5000,]

# 2. Annotation
dfAnnotationiChIP <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/AnnotationTable_iChIPdata.tab",header=T)
dfAnnotationiChIP$TimeCateg <- paste0(rep("_",length(dfAnnotationiChIP$Time)),dfAnnotationiChIP$Time)
rownames(dfAnnotationiChIP) <- dfAnnotationiChIP$Sample

dfAnnotationA <- dfAnnotationiChIP

# All zero rows
vecboolNonzeroA <- apply(matDataA,1,function(gene){any(gene>0)})

# Summary matrices
matQval_iChIPData <- matrix(c("NA",rep(NA,4)),nrow=dim(matDataA)[1],ncol=4+1,byrow=TRUE)
matPval_iChIPData <- matrix(c("NA",rep(NA,4)),nrow=dim(matDataA)[1],ncol=4+1,byrow=TRUE)
matRunTime_iChIPData <- array(NA,c(4))

rownames(matQval_iChIPData) <- rownames(matDataA)
colnames(matQval_iChIPData) <- c("Gene",paste0(
  rep("A_",4), 
  c("ImpulseDE2", "ImpulseDE", "DESeq2", "edge")
))
matQval_iChIPData[,"Gene"] <- rownames(matDataA)

rownames(matPval_iChIPData) <- rownames(matDataA)
colnames(matPval_iChIPData) <- c("Gene",paste0(
  rep("A_",4), 
  c("ImpulseDE2", "ImpulseDE", "DESeq2", "edge")
))
matPval_iChIPData[,"Gene"] <- rownames(matDataA)

names(matRunTime_iChIPData) <- c("ImpulseDE2", "ImpulseDE", "DESeq2", "edge")

################################################################################
# ImpulseDE2
print("Run ImpulseDE2")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/ImpulseDE2_main.R")

# Create input data set
# Only retain non zero
matDataA_ImpulseDE2 <- matDataA[vecboolNonzeroA,]

tm_ImpulseDE2A <- system.time({
  setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/ImpulseDE2/H3K4me1")
  strControlName = NULL
  strCaseName = "case"
  lsImpulseDE_resultsA <- runImpulseDE2(matCountData=matDataA_ImpulseDE2, 
    dfAnnotation=dfAnnotationA,
    strCaseName = strCaseName, 
    strControlName=strControlName, 
    strMode="batch",
    nProc=3, 
    Q_value=10^(-3),
    boolPlotting=FALSE)
  dfImpulseResultsA <- lsImpulseDE_resultsA$dfImpulseResults
  qvals_A <- dfImpulseResultsA$adj.p
  pvals_A <- dfImpulseResultsA$p
})

# Extract Results
matQval_iChIPData[vecboolNonzeroA,"A_ImpulseDE2"] <- qvals_A
matPval_iChIPData[vecboolNonzeroA,"A_ImpulseDE2"] <- pvals_A

matRunTime_iChIPData["ImpulseDE2"] <- round(tm_ImpulseDE2A["elapsed"])

lsResDEcomparison_iChIPdataH3K4me1 <- list("matQval_iChIPData"=matQval_iChIPData, 
  "matPval_iChIPData"=matPval_iChIPData, 
  "matRunTime_iChIPData"=matRunTime_iChIPData)

save(lsResDEcomparison_iChIPdataH3K4me1,file=file.path("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/lsResDEcomparison_iChIPdataH3K4me1.RData"))
print("Finished ImpulseDE2")

################################################################################
# ImpulseDE
print("Run ImpulseDE")
source("/Users/davidsebastianfischer/MasterThesis/software/ImpulseDE/R/Impulse_DE_fin.R")

# Create input data set
# Only retain non zero
matDataA_ImpulseDE <- matDataA[vecboolNonzeroA,]

if(FALSE){
tm_ImpulseDE2A <- system.time({
  setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/ImpulseDE/H3K4me1")
  strControlName = NULL
  strCaseName = "case"
  lsImpulseDE_resultsA <- impulse_DE(
    expression_table = matDataA_ImpulseDE, 
    annotation_table = dfAnnotationA,
    colname_time = "Time", colname_condition = "Condition", 
    control_timecourse = FALSE,
    control_name = strControlName, case_name = strCaseName, 
    expr_type = "Seq",
    plot_clusters = FALSE, 
    n_iter = 100, n_randoms = 50000, 
    n_process = 3,
    Q_value = 0.01)
  dfImpulseResultsA <- lsImpulseDE_resultsA$impulse_fit_results
  qvals_A <- dfImpulseResultsA$adj.p
  pvals_A <- dfImpulseResultsA$p
})
}
################################################################################
# DESeq2
print("Run DESeq2")
library(DESeq2)
library(BiocParallel)
# Parallelisation
nProcessesAssigned <- 3
nProcesses <- min(detectCores() - 1, nProcessesAssigned)
register(MulticoreParam(nProcesses))

# Create input data set
# Only retain non zero
matDataA_DESeq2 <- matDataA[vecboolNonzeroA,]

tm_DESeq2A <- system.time({
  # Create DESeq2 data object
  dds <- DESeqDataSetFromMatrix(countData = matDataA_DESeq2,
    colData = dfAnnotationA,
    design = ~TimeCateg)
  # Run DESeq2
  ddsDESeqObjectA <- DESeq(dds, test = "LRT", 
    full = ~ TimeCateg, reduced = ~ 1,
    parallel=TRUE)
  
  # DESeq results for comparison
  dds_resultsTableA <- results(ddsDESeqObjectA)
  qvals_A <- dds_resultsTableA$padj
  pvals_A <- dds_resultsTableA$pvalue
})

# Extract Results
matQval_iChIPData[vecboolNonzeroA,"A_DESeq2"] <- qvals_A
matPval_iChIPData[vecboolNonzeroA,"A_DESeq2"] <- pvals_A
matRunTime_iChIPData["DESeq2"] <- round(tm_DESeq2A["elapsed"])

lsResDEcomparison_iChIPdataH3K4me1 <- list("matQval_iChIPData"=matQval_iChIPData, 
  "matPval_iChIPData"=matPval_iChIPData, 
  "matRunTime_iChIPData"=matRunTime_iChIPData)

save(lsResDEcomparison_iChIPdataH3K4me1,file=file.path("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/lsResDEcomparison_iChIPdataH3K4me1.RData"))
print("Finished DESeq2.")

################################################################################
# EDGE
print("Run edge")
library(edge)
library(splines)

# Create input data set
# Take out row and col names, only retain non zero
matDataA_EDGE <- matDataA[vecboolNonzeroA,]
rownames(matDataA_EDGE) <- NULL
colnames(matDataA_EDGE) <- NULL

# (I) CASE
# A)
tm_edgeA <- system.time({
  cov <- data.frame(time = dfAnnotationA$Time)
  null_model <- ~1
  full_model <- ~ns(time, df=4)
  edge_obj <- build_models(data = matDataA_EDGE, cov = cov, null.model = null_model, full.model = full_model)
  
  # Normalise
  edge_norm <- apply_snm(edge_obj, int.var=1:ncol(exprs(edge_obj)), diagnose=FALSE)
  # Adjust for unmodelled variables
  #edge_sva <- apply_sva(edge_norm)
  
  # optimal discovery procedure: Chose this one
  #edge_odp <- odp(edge_sva, bs.its = 30, verbose=FALSE)
  edge_odp <- odp(edge_norm, bs.its = 30, verbose=FALSE)
  # likelihood ratio test: throws error
  # edge_lrt <- lrt(edge_sva)
  
  # Extract results
  qval_obj_A <- qvalueObj(edge_odp)
  qvals_A <- qval_obj_A$qvalues
  pvals_A <- qval_obj_A$pvalues
  lfdr_A <- qval_obj_A$lfdr
  pi0_A <- qval_obj_A$pi0
})

# Extract Results
matQval_iChIPData[vecboolNonzeroA,"A_edge"] <- qvals_A
matPval_iChIPData[vecboolNonzeroA,"A_edge"] <- pvals_A
matRunTime_iChIPData["edge"] <- round(tm_edgeA["elapsed"])

lsResDEcomparison_iChIPdataH3K4me1 <- list("matQval_iChIPData"=matQval_iChIPData, 
  "matPval_iChIPData"=matPval_iChIPData, 
  "matRunTime_iChIPData"=matRunTime_iChIPData)

save(lsResDEcomparison_iChIPdataH3K4me1,file=file.path("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/lsResDEcomparison_iChIPdataH3K4me1.RData"))
print("Finished edge.")

#################################################
#################################################
#################################################
# Graphical analysis
#################################################
#################################################
#################################################

########################################
# 1. Plot ImpulseDE2 DE genes
########################################
rm(list=ls())

source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/ImpulseDE2_main.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/srcImpulseDE2_plotDEGenes.R")
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/ImpulseDE2/H3K4me1")
load("ImpulseDE2_matCountDataProc.RData")
load("ImpulseDE2_dfAnnotationProc.RData")
load("ImpulseDE2_dfImpulseResults.RData")
load("ImpulseDE2_vecDEGenes.RData")
load("ImpulseDE2_lsImpulseFits.RData")
load("ImpulseDE2_dfDESeq2Results.RData")
load("ImpulseDE2_lsMatTranslationFactors.RData")
load("ImpulseDE2_matSizeFactors.RData")

setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/pdfs")
QplotImpulseDE2genes <- 10^(-3)
vecImpulse2_DEgenes <- dfImpulseResults[dfImpulseResults$padj <= QplotImpulseDE2genes,]$Gene 
strCaseName <- "case"
strControlName <- NULL
strMode <- "batch"
if(length(vecImpulse2_DEgenes)==0){
  vecImpulse2_DEgenes <- rownames(matCountDataProc[1:32,])
}
plotDEGenes(vecGeneIDs=vecImpulse2_DEgenes,
    matCountDataProc=matCountDataProc,
    lsMatTranslationFactors=lsMatTranslationFactors,
    matSizeFactors=matSizeFactors,
    dfAnnotation=dfAnnotationProc, 
    lsImpulseFits=lsImpulseFits,
    strCaseName=strCaseName, 
    strControlName=strControlName, 
    strFileNameSuffix="H4K3me1_ImpulseDE2_DEgenes", 
    strPlotTitleSuffix="", 
    strPlotSubtitle="",
    dfImpulseResults=dfImpulseResults,
    vecRefPval=dfDESeq2Results$padj,
    strMode=strMode, 
    NPARAM=6)

########################################
# 2. Heatmap: Q-value thresholds
########################################
rm(list=ls())
library(gplots)
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/lsResDEcomparison_iChIPdataH3K4me1.Rdata")
matQval_iChIPData <- lsResDEcomparison_iChIPdataH3K4me1$matQval_iChIPData
matPval_iChIPData <- lsResDEcomparison_iChIPdataH3K4me1$matPval_iChIPData
matRunTime_iChIPData <- lsResDEcomparison_iChIPdataH3K4me1$matRunTime_iChIPData
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/pdfs")
########################################
# DESeq2
mat_overlap <- array(NA,c(11,11))
for(i in 0:10){
  for(j in 0:10){
    sig_Impulse <- matQval_iChIPData[,"A_ImpulseDE2"] <= 10^(-j)
    sig_Ref <- matQval_iChIPData[,"A_DESeq2"] <= 10^(-i)
    mat_overlap[i+1,j+1] <- sum(sig_Ref & sig_Impulse)/sum(sig_Ref | sig_Impulse)
  }
}
#mat_overlap[is.na(mat_overlap)] <- 0
rownames(mat_overlap) <- 0:-10
colnames(mat_overlap) <- 0:-10

graphics.off()
heatmap(mat_overlap, keep.dendro = FALSE,Rowv=NA,Colv= "Rowv",symm=FALSE,
  xlab =  paste0("ImpulseDE2 scores"), ylab = "edge scores")
breaks <- seq(0,1,by=0.01)
hm.colors <- colorpanel( length(breaks)-1, "yellow", "red" )
graphics.off()
pdf(paste('ImpulseDE2-DESeq2_JaccardCoeff_Heatmap.pdf',sep=''),width=7,height=7)
heatmap.2(mat_overlap, dendrogram="none", Rowv=FALSE,Colv=FALSE, 
  xlab =  paste0("log(p-value) ImpulseDE2"), ylab = "log(p-value) DESeq2",
  breaks=breaks,col=hm.colors, scale="none",
  trace="none",density.info="none",
  key.title = " ", key.xlab = paste0("Jaccard coefficient"), key.ylab = NULL,
  symkey=FALSE,
  cellnote=round(mat_overlap,digits=2),notecol="grey",
  lmat=rbind( c(3,4),c(2,1) ),lhei=c(1,4), lwid=c(1,4), margins=c(5,5))
dev.off()

########################################
# edge
mat_overlap <- array(NA,c(11,11))
for(i in 0:10){
  for(j in 0:10){
    sig_Impulse <- matQval_iChIPData[,"A_ImpulseDE2"] <= 10^(-j)
    sig_Ref <- matQval_iChIPData[,"A_edge"] <= 10^(-i)
    mat_overlap[i+1,j+1] <- sum(sig_Ref & sig_Impulse)/sum(sig_Ref | sig_Impulse)
  }
}
#mat_overlap[is.na(mat_overlap)] <- 0
rownames(mat_overlap) <- 0:-10
colnames(mat_overlap) <- 0:-10

graphics.off()
heatmap(mat_overlap, keep.dendro = FALSE,Rowv=NA,Colv= "Rowv",symm=FALSE,
  xlab =  paste0("ImpulseDE2 scores"), ylab = "edge scores")
breaks <- seq(0,1,by=0.01)
hm.colors <- colorpanel( length(breaks)-1, "yellow", "red" )
graphics.off()
pdf(paste('ImpulseDE2-edge_JaccardCoeff_Heatmap.pdf',sep=''),width=7,height=7)
heatmap.2(mat_overlap, dendrogram="none", Rowv=FALSE,Colv=FALSE, 
  xlab =  paste0("log(p-value) ImpulseDE2"), ylab = "log(p-value) edge",
  breaks=breaks,col=hm.colors, scale="none",
  trace="none",density.info="none",
  key.title = " ", key.xlab = paste0("Jaccard coefficient"), key.ylab = NULL,
  symkey=FALSE,
  cellnote=round(mat_overlap,digits=2),notecol="grey",
  lmat=rbind( c(3,4),c(2,1) ),lhei=c(1,4), lwid=c(1,4), margins=c(5,5))
dev.off()

########################################
# 3. Impulse fits, differentially called DE genes
########################################
rm(list=ls())
library(gplots)
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/lsResDEcomparison_iChIPdataH3K4me1.Rdata")
matQval_iChIPData <- lsResDEcomparison_iChIPdataH3K4me1$matQval_iChIPData 
matPval_iChIPData <- lsResDEcomparison_iChIPdataH3K4me1$matPval_iChIPData
matRunTime_iChIPData <- lsResDEcomparison_iChIPdataH3K4me1$matRunTime_iChIPData

source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/ImpulseDE2_main.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/srcImpulseDE2_plotDEGenes.R")
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/ImpulseDE2/H3K4me1")
load("ImpulseDE2_matCountDataProc.RData")
load("ImpulseDE2_dfAnnotationProc.RData")
load("ImpulseDE2_dfImpulseResults.RData")
load("ImpulseDE2_vecDEGenes.RData")
load("ImpulseDE2_lsImpulseFits.RData")
load("ImpulseDE2_dfDESeq2Results.RData")
load("ImpulseDE2_lsMatTranslationFactors.RData")
load("ImpulseDE2_matSizeFactors.RData")

setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/pdfs")

Q <- 10^(-3)
Qdelta <- 10^(2) # difference factor required to be plotted

########################################
# DESeq2
DEgenes_Ref_only <- matQval_iChIPData[(matQval_iChIPData[,"A_DESeq2"] < Q & matQval_iChIPData[,"A_ImpulseDE2"] >= Qdelta*Q),"Gene"]
DEgenes_Impulse_only <- matQval_iChIPData[(matQval_iChIPData[,"A_ImpulseDE2"] >= Qdelta*Q & matQval_iChIPData[,"A_DESeq2"] < Q),"Gene"]
graphics.off()
print(paste0("Number of significant DE genes at q-value of ",Q))
print(paste0("ImpulseDE only ",length(DEgenes_Impulse_only)))
print(paste0("DESeq2 only ",length(DEgenes_Ref_only,)))

vecRefResults <- matQval_iChIPData[,"A_DESeq2"]
names(vecRefResults) <- matQval_iChIPData[,"Gene"]
strControlName <- NULL
strCaseName <- "case"
strMode="batch"
NPARAM <- 6

# sort DESeq2 only genes by padj of DESeq2
dfDESeq_ImpulseDESeqOnly <- dfDESeq_Impulse[DEgenes_Ref_only,]
DEgenes_Ref_onlySorted <- dfDESeq_ImpulseDESeqOnly[order(dfDESeq_ImpulseDESeqOnly$DESeq),"Gene"]
plotDEGenes(vecGeneIDs=DEgenes_Ref_onlySorted,
    matCountDataProc=matCountDataProc,
    lsMatTranslationFactors=lsMatTranslationFactors,
    matSizeFactors=matSizeFactors,
    dfAnnotation=dfAnnotationProc, 
    lsImpulseFits=lsImpulseFits,
    strCaseName=strCaseName, 
    strControlName=strControlName, 
    strFileNameSuffix="H4K3me1_DESeq2BeatsImpulseDE2", 
    strPlotTitleSuffix="", 
    strPlotSubtitle="",
    dfImpulseResults=dfImpulseResults,
    vecRefPval=dfDESeq2Results$padj,
    strMode=strMode, 
    NPARAM=NPARAM)

plotDEGenes(vecGeneIDs=DEgenes_Impulse_only,
    matCountDataProc=matCountDataProc,
    lsMatTranslationFactors=lsMatTranslationFactors,
    matSizeFactors=matSizeFactors,
    dfAnnotation=dfAnnotationProc, 
    lsImpulseFits=lsImpulseFits,
    strCaseName=strCaseName, 
    strControlName=strControlName, 
    strFileNameSuffix="H4K3me1_ImpulseDE2BeatsDESeq2", 
    strPlotTitleSuffix="", 
    strPlotSubtitle="",
    dfImpulseResults=dfImpulseResults,
    vecRefPval=dfDESeq2Results$padj,
    strMode=strMode, 
    NPARAM=NPARAM)

########################################
# edge
DEgenes_Ref_only <- matQval_iChIPData[(matQval_iChIPData[,"A_DESeq2"] < Q & matQval_iChIPData[,"A_ImpulseDE2"] >= Qdelta*Q),"Gene"]
DEgenes_Impulse_only <- matQval_iChIPData[(matQval_iChIPData[,"A_ImpulseDE2"] >= Qdelta*Q & matQval_iChIPData[,"A_DESeq2"] < Q),"Gene"]
graphics.off()
print(paste0("Number of significant DE genes at q-value of ",Q))
print(paste0("ImpulseDE only ",length(DEgenes_Impulse_only)))
print(paste0("edge only ",length(DEgenes_edge_only.RData)))

vecRefResults <- matQval_iChIPData[,"A_edge"]
names(vecRefResults) <- matQval_iChIPData[,"Gene"]
strControlName <- NULL
strCaseName <- "case"
NPARAM <- 6

# sort DESeq2 only genes by padj of DESeq2
dfDESeq_ImpulseDESeqOnly <- dfDESeq_Impulse[DEgenes_Ref_only,]
DEgenes_Ref_onlySorted <- dfDESeq_ImpulseDESeqOnly[order(dfDESeq_ImpulseDESeqOnly$DESeq),"Gene"]
plotDEGenes(vecGeneIDs=DEgenes_Ref_onlySorted,
    matCountDataProc=matCountDataProc,
    lsMatTranslationFactors=lsMatTranslationFactors,
    matSizeFactors=matSizeFactors,
    dfAnnotation=dfAnnotationProc, 
    lsImpulseFits=lsImpulseFits,
    strCaseName=strCaseName, 
    strControlName=strControlName, 
    strFileNameSuffix="H4K3me1_EdgeBeatsImpulseDE2", 
    strPlotTitleSuffix="", 
    strPlotSubtitle="",
    dfImpulseResults=dfImpulseResults,
    vecRefPval=dfDESeq2Results$padj,
    strMode=strMode, 
    NPARAM=NPARAM)

plotDEGenes(vecGeneIDs=DEgenes_Impulse_only,
    matCountDataProc=matCountDataProc,
    lsMatTranslationFactors=lsMatTranslationFactors,
    matSizeFactors=matSizeFactors,
    dfAnnotation=dfAnnotationProc, 
    lsImpulseFits=lsImpulseFits,
    strCaseName=strCaseName, 
    strControlName=strControlName, 
    strFileNameSuffix="H4K3me1_ImpulseDE2BeatsEdge", 
    strPlotTitleSuffix="", 
    strPlotSubtitle="",
    dfImpulseResults=dfImpulseResults,
    vecRefPval=dfDESeq2Results$padj,
    strMode=strMode, 
    NPARAM=NPARAM)
