# Local
rm(list = ls())
boolCluster <- FALSE

print("Process data RNAseq Jankovic")
# Load Data set RNAseq
# 1. Counts
dfRNA <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/rsem_readCountsTable.txt",
  sep="\t",header=F,colClasses=c(
    "numeric","numeric","numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric","numeric","numeric"))
dfGeneIDs <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/gene_list.txt",sep="\t",header=F)
dfCellIDs <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/cell_list.txt",sep="\t",header=F)
vecSamples <- apply(dfCellIDs, 1, function(name){ unlist(strsplit(unlist(strsplit(name,"_1"))[1],"/"))[2] })

matDataA <- round(t(apply(dfRNA,1,as.numeric)))
vecGeneIDs <- as.vector(dfGeneIDs[,2])
vecboolidxDupIDs <- duplicated(vecGeneIDs)
if(sum(vecboolidxDupIDs)>0){
  vecGeneIDs[vecboolidxDupIDs] <- paste0(rep("DuplicatedGene_",sum(vecboolidxDupIDs)),seq(1,sum(vecboolidxDupIDs)))
}
rownames(matDataA) <- vecGeneIDs
colnames(matDataA) <- vecSamples
matDataA <- matDataA

# 2. Annotation
dfAnnotationRNA <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/AnnotationTable_RNAseqJankovic.tab",header=T)
dfAnnotationRNA$TimeCateg <- paste0(rep("_",length(dfAnnotationRNA$Time)),dfAnnotationRNA$Time)
dfAnnotationA <- dfAnnotationRNA

# Expend 0h ctrl sample to both conditions
vecExpanedSamples <- dfAnnotationA$Sample
vecExpanedSamples[vecExpanedSamples=="SRR1525500cs"] <- "SRR1525500"
vecExpanedSamples[vecExpanedSamples=="SRR1525513cs"] <- "SRR1525513"

matDataA <- matDataA[,as.vector(vecExpanedSamples)]
colnames(matDataA) <- dfAnnotationA$Sample
#only case
matDataA <- matDataA[,dfAnnotationA$Condition=="case"]
colnames(matDataA) <- dfAnnotationA[dfAnnotationA$Condition=="case",]$Sample
dfAnnotationA <- dfAnnotationA[dfAnnotationA$Condition=="case",]
rownames(dfAnnotationA) <- dfAnnotationA$Sample

# All zero rows
matDataA[!is.finite(matDataA)] <- NA
vecboolNonzeroA <- apply(matDataA,1,function(gene){any(gene>0 & !is.na(gene) & is.finite(gene))})
matDataA <- matDataA[vecboolNonzeroA,]

# Summary matrices
matQval <- matrix(c("NA",rep(NA,4)),nrow=dim(matDataA)[1],ncol=4+1,byrow=TRUE)
matPval <- matrix(c("NA",rep(NA,4)),nrow=dim(matDataA)[1],ncol=4+1,byrow=TRUE)
matRunTime_RNAseqData <- array(NA,c(4))

rownames(matQval) <- rownames(matDataA)
colnames(matQval) <- c("Gene",paste0(
  rep("A_",4), 
  c("ImpulseDE2", "ImpulseDE", "DESeq2", "edge")
))
matQval[,"Gene"] <- rownames(matDataA)

rownames(matPval) <- rownames(matDataA)
colnames(matPval) <- c("Gene",paste0(
  rep("A_",4), 
  c("ImpulseDE2", "ImpulseDE", "DESeq2", "edge")
))
matPval[,"Gene"] <- rownames(matDataA)

names(matRunTime_RNAseqData) <- c("ImpulseDE2", "ImpulseDE", "DESeq2", "edge")
if(boolCluster){
  load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/clusterruns/output/lsResDEcomparison_RNAseqdata.RData")
  matQval <- lsResDEcomparison_RNAseqdata$matQval
  matPval <- lsResDEcomparison_RNAseqdata$matPval
  matRunTime <- lsResDEcomparison_RNAseqdata$matRunTime_RNAseqData
}

################################################################################
# ImpulseDE2
if(!boolCluster){
  print("Run ImpulseDE2")
  source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/ImpulseDE2_main.R")
  
  # Create input data set
  # Only retain non zero
  matDataA_ImpulseDE2 <- matDataA
  
  tm_ImpulseDE2A <- system.time({
    setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/ImpulseDE2")
    strControlName = NULL
    strCaseName = "case"
    lsImpulseDE_resultsA <- runImpulseDE2(matCountData=matDataA_ImpulseDE2, 
      dfAnnotation=dfAnnotationA,
      strCaseName = strCaseName, 
      strControlName=strControlName, 
      strMode="longitudinal",
      nProc=3, 
      Q_value=10^(-3),
      boolPlotting=FALSE)
    dfImpulseResultsA <- lsImpulseDE_resultsA$dfImpulseResults
    qvals_A <- dfImpulseResultsA$adj.p
    pvals_A <- dfImpulseResultsA$p
    ids_A <- as.vector(dfImpulseResultsA$Gene)
  })
  
  # Extract Results
  matQval[ids_A,"A_ImpulseDE2"] <- qvals_A
  matPval[ids_A,"A_ImpulseDE2"] <- pvals_A
  
  matRunTime_RNAseqData["ImpulseDE2"] <- round(tm_ImpulseDE2A["elapsed"])
  
  lsResDEcomparison_RNAseqdata <- list("matQval"=matQval, 
    "matPval"=matPval, 
    "matRunTime_RNAseqData"=matRunTime_RNAseqData)
  
  save(lsResDEcomparison_RNAseqdata,file=file.path("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/lsResDEcomparison_RNAseqdata.RData"))
  print("Finished ImpulseDE2")
  
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
  matDataA_DESeq2 <- matDataA
  
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
    ids_A <- rownames(dds_resultsTableA)
  })
  
  # Extract Results
  matQval[ids_A,"A_DESeq2"] <- qvals_A
  matPval[ids_A,"A_DESeq2"] <- pvals_A
  matRunTime_RNAseqData["DESeq2"] <- round(tm_DESeq2A["elapsed"])
  
  lsResDEcomparison_RNAseqdata <- list("matQval"=matQval, 
    "matPval"=matPval, 
    "matRunTime_RNAseqData"=matRunTime_RNAseqData)
  
  save(lsResDEcomparison_RNAseqdata,file=file.path("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/lsResDEcomparison_RNAseqdata.RData"))
  print("Finished DESeq2.")
  
  ################################################################################
  # ImpulseDE
  # Run on server
  print("Run ImpulseDE")
  library(compiler)
  library(grDevices)
  library(graphics)
  library(stats)
  library(stats)
  library(utils)
  library(amap)
  library(ImpulseDE)
  
  # Create input data set
  # Only retain non zero, normalise and log transform
  matDataA_ImpulseDE <- matDataA/rep(sizeFactors(ddsDESeqObjectA), each = nrow(matDataA))
  matDataA_ImpulseDE[matDataA_ImpulseDE==0] <- 1
  matDataA_ImpulseDE <- log(matDataA_ImpulseDE)/log(2)
  
  dfAnnotationA_ImpulseDE <- dfAnnotationA[dfAnnotationA$Condition=="case",c("Time","Condition")]
  dfAnnotationA_ImpulseDE$Condition <- as.vector(dfAnnotationA_ImpulseDE$Condition)
  dfAnnotationA_ImpulseDE$Time <- as.numeric(as.vector(dfAnnotationA_ImpulseDE$Time))
  rownames(dfAnnotationA_ImpulseDE) <- as.vector(dfAnnotationA$Sample)
  tm_ImpulseDEA <- system.time({
    setwd("/data/yosef2/users/fischerd/data/RNAseq_Jovanovic/comparative_analysis/output/ImpulseDE")
    strControlName = NULL
    strCaseName = "case"
    lsImpulseDE_resultsA <- impulse_DE(
      expression_table = matDataA_ImpulseDE, 
      annotation_table = dfAnnotationA_ImpulseDE,
      colname_time = "Time", 
      colname_condition = "Condition", 
      control_timecourse = FALSE,
      control_name = strControlName,
      case_name = strCaseName, 
      expr_type = "Seq",
      plot_clusters = FALSE, 
      n_iter = 100, 
      n_randoms = 50000, 
      n_process = NCORES,
      Q_value = 0.01)
    dfImpulseResultsA <- lsImpulseDE_resultsA$DE_results
    qvals_A <- dfImpulseResultsA$adj.p
    ids_A <- dfImpulseResultsA$Gene
  })
  
  # Extract Results
  matQval[ids_A,"A_ImpulseDE"] <- qvals_A
  
  matRunTime_RNAseqData["ImpulseDE"] <- round(tm_ImpulseDE2A["elapsed"])
  
  lsResDEcomparison_RNAseqdata <- list("matQval"=matQval, 
    "matPval"=matPval, 
    "matRunTime_RNAseqData"=matRunTime_RNAseqData)
  
  save(lsResDEcomparison_RNAseqdata,file=file.path("/data/yosef2/users/fischerd/data/RNAseq_Jovanovic/comparative_analysis/output/lsResDEcomparison_RNAseqdata.RData"))
  print("Finished ImpulseDE")
  
}
################################################################################
# EDGE
print("Run edge")
library(edge)
library(splines)

# Create input data set
# Take out row and col names, only retain non zero
matDataA_EDGE <- matDataA
rownames(matDataA_EDGE) <- NULL
colnames(matDataA_EDGE) <- NULL

# (I) CASE
# A)
tm_edgeA <- system.time({
  cov <- data.frame(time = dfAnnotationA$Time )
  null_model <- ~1
  full_model <- ~ns(time, df=5)
  edge_obj <- build_models(data = matDataA_EDGE, cov = cov, null.model = null_model, full.model = full_model)
  
  # Normalise
  edge_norm <- apply_snm(edge_obj, int.var=1:ncol(exprs(edge_obj)), diagnose=FALSE)
  # Adjust for unmodelled variables
  #edge_sva <- apply_sva(edge_norm)
  
  # optimal discovery procedure: Chose this one
  #edge_res <- odp(edge_sva, bs.its = 30, verbose=FALSE)
  edge_res <- odp(edge_norm, bs.its = 30, verbose=FALSE)
  ## likelihood ratio test: throws error? not anymore
  #edge_res <- lrt(edge_norm)
  
  # Extract results
  qval_obj_A <- qvalueObj(edge_res)
  qvals_A <- qval_obj_A$qvalues
  pvals_A <- qval_obj_A$pvalues
  lfdr_A <- qval_obj_A$lfdr
  pi0_A <- qval_obj_A$pi0
})

# Extract Results
matQval[,"A_edge"] <- qvals_A
matPval[,"A_edge"] <- pvals_A
matRunTime_RNAseqData["edge"] <- round(tm_edgeA["elapsed"])

lsResDEcomparison_RNAseqdata <- list("matQval"=matQval, 
  "matPval"=matPval, 
  "matRunTime_RNAseqData"=matRunTime_RNAseqData)

save(lsResDEcomparison_RNAseqdata,file=file.path("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/edge/lsResDEcomparison_RNAseqdata.RData"))
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
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/clusterruns/ImpulseDE2")
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
# 2. Heatmap: Q-value thresholds
########################################
rm(list=ls())
library(gplots)
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/clusterruns/lsResDEcomparison_RNAseqdata.Rdata")
matQval <- lsResDEcomparison_RNAseqdata$matQval
matPval <- lsResDEcomparison_RNAseqdata$matPval
matRunTime_RNAseqData <- lsResDEcomparison_RNAseqdata$matRunTime_RNAseqData
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/pdfs")
########################################
# DESeq2
mat_overlap <- array(NA,c(11,11))
for(i in 0:10){
  for(j in 0:10){
    sig_Impulse <- as.numeric(matQval[,"A_ImpulseDE2"]) <= 10^(-j) & !is.na(as.numeric(matQval[,"A_ImpulseDE2"]))
    sig_Ref <- as.numeric(matQval[,"A_DESeq2"]) <= 10^(-i) & !is.na(as.numeric(matQval[,"A_DESeq2"]))
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
    sig_Impulse <- as.numeric(matQval[,"A_ImpulseDE2"]) <= 10^(-j) & !is.na(as.numeric(matQval[,"A_ImpulseDE2"]))
    sig_Ref <- as.numeric(matQval[,"A_edge"]) <= 10^(-i) & !is.na(as.numeric(matQval[,"A_edge"]))
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
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/lsResDEcomparison_RNAseqdata.Rdata")
matQval <- lsResDEcomparison_RNAseqdata$matQval 
matPval <- lsResDEcomparison_RNAseqdata$matPval
matRunTime_RNAseqData <- lsResDEcomparison_RNAseqdata$matRunTime_RNAseqData

source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/ImpulseDE2_main.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/srcImpulseDE2_plotDEGenes.R")
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/clusterruns/ImpulseDE2")
load("ImpulseDE2_matCountDataProc.RData")
load("ImpulseDE2_dfAnnotationProc.RData")
load("ImpulseDE2_dfImpulseResults.RData")
load("ImpulseDE2_vecDEGenes.RData")
load("ImpulseDE2_lsImpulseFits.RData")
load("ImpulseDE2_dfDESeq2Results.RData")
load("ImpulseDE2_lsMatTranslationFactors.RData")
load("ImpulseDE2_matSizeFactors.RData")

setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/pdfs")

Q <- 10^(-3)
Qdelta <- 10^(2) # difference factor required to be plotted

########################################
# DESeq2
# Note cannot plot genes which are NA in ImpulseDE2
vecboolImpulseFitted <- !is.na(as.numeric(matQval[,"A_DESeq2"]))
vecboolRefBeatsQ <- as.numeric(matQval[,"A_DESeq2"]) < Q & !is.na(as.numeric(matQval[,"A_DESeq2"]))
vecBoolRefBeatsImpulse <- as.numeric(matQval[,"A_ImpulseDE2"]) >= Qdelta*as.numeric(matQval[,"A_DESeq2"]) | is.na(as.numeric(matQval[,"A_DESeq2"]))
vecboolImpulseBeatsQ <- as.numeric(matQval[,"A_ImpulseDE2"]) < Q & !is.na(as.numeric(matQval[,"A_ImpulseDE2"]))
vecboolImpulseBeatsRef <- as.numeric(matQval[,"A_DESeq2"]) >= Qdelta*as.numeric(matQval[,"A_ImpulseDE2"]) | is.na(as.numeric(matQval[,"A_ImpulseDE2"]))

vecDEgenes_Ref_only <- as.vector(matQval[vecboolRefBeatsQ & vecBoolRefBeatsImpulse & vecboolImpulseFitted,"Gene"])
vecDEgenes_Impulse_only <- as.vector(matQval[vecboolImpulseBeatsQ & vecboolImpulseBeatsRef & vecboolImpulseFitted,"Gene"])
graphics.off()
print(paste0("Number of significant DE genes at q-value of ",Q))
print(paste0("ImpulseDE only ",length(vecDEgenes_Impulse_only)))
print(paste0("DESeq2 only ",length(vecDEgenes_Ref_only)))

vecRefResults <- as.numeric(matQval[,"A_DESeq2"])
names(vecRefResults) <- matQval[,"Gene"]
strControlName <- NULL
strCaseName <- "case"
strMode="batch"
NPARAM <- 6

plotDEGenes(vecGeneIDs=vecDEgenes_Ref_only,
  matCountDataProc=matCountDataProc,
  lsMatTranslationFactors=lsMatTranslationFactors,
  matSizeFactors=matSizeFactors,
  dfAnnotation=dfAnnotationProc, 
  lsImpulseFits=lsImpulseFits,
  strCaseName=strCaseName, 
  strControlName=strControlName, 
  strFileNameSuffix="RNAseq_DESeq2BeatsImpulseDE2", 
  strPlotTitleSuffix="", 
  strPlotSubtitle="",
  dfImpulseResults=dfImpulseResults,
  vecRefPval=vecRefResults,
  strNameMethod2="DESeq2",
  strMode=strMode, 
  NPARAM=NPARAM)

if(length(vecDEgenes_Impulse_only)){
  plotDEGenes(vecGeneIDs=vecDEgenes_Impulse_only,
    matCountDataProc=matCountDataProc,
    lsMatTranslationFactors=lsMatTranslationFactors,
    matSizeFactors=matSizeFactors,
    dfAnnotation=dfAnnotationProc, 
    lsImpulseFits=lsImpulseFits,
    strCaseName=strCaseName, 
    strControlName=strControlName, 
    strFileNameSuffix="RNAseq_ImpulseDE2BeatsDESeq2", 
    strPlotTitleSuffix="", 
    strPlotSubtitle="",
    dfImpulseResults=dfImpulseResults,
    vecRefPval=vecRefResults,
    strNameMethod2="DESeq2",
    strMode=strMode, 
    NPARAM=NPARAM)
}
########################################
# edge
vecboolImpulseFitted <- !is.na(as.numeric(matQval[,"A_edge"]))
vecboolRefBeatsQ <- as.numeric(matQval[,"A_edge"]) < Q & !is.na(as.numeric(matQval[,"A_edge"]))
vecBoolRefBeatsImpulse <- as.numeric(matQval[,"A_ImpulseDE2"]) >= Qdelta*as.numeric(matQval[,"A_edge"]) | is.na(as.numeric(matQval[,"A_DESeq2"]))
vecboolImpulseBeatsQ <- as.numeric(matQval[,"A_ImpulseDE2"]) < Q & !is.na(as.numeric(matQval[,"A_ImpulseDE2"]))
vecboolImpulseBeatsRef <- as.numeric(matQval[,"A_edge"]) >= Qdelta*as.numeric(matQval[,"A_ImpulseDE2"]) | is.na(as.numeric(matQval[,"A_ImpulseDE2"]))

vecDEgenes_Ref_only <- as.vector(matQval[vecboolRefBeatsQ & vecBoolRefBeatsImpulse & vecboolImpulseFitted,"Gene"])
vecDEgenes_Impulse_only <- as.vector(matQval[vecboolImpulseBeatsQ & vecboolImpulseBeatsRef & vecboolImpulseFitted,"Gene"])
graphics.off()
print(paste0("Number of significant DE genes at q-value of ",Q))
print(paste0("ImpulseDE only ",length(vecDEgenes_Impulse_only)))
print(paste0("Edge only ",length(vecDEgenes_Ref_only)))

vecRefResults <- as.numeric(matQval[,"A_edge"])
names(vecRefResults) <- matQval[,"Gene"]
strControlName <- NULL
strCaseName <- "case"
strMode="batch"
NPARAM <- 6

plotDEGenes(vecGeneIDs=vecDEgenes_Ref_only,
  matCountDataProc=matCountDataProc,
  lsMatTranslationFactors=lsMatTranslationFactors,
  matSizeFactors=matSizeFactors,
  dfAnnotation=dfAnnotationProc, 
  lsImpulseFits=lsImpulseFits,
  strCaseName=strCaseName, 
  strControlName=strControlName, 
  strFileNameSuffix="RNAseq_EdgeBeatsImpulseDE2", 
  strPlotTitleSuffix="", 
  strPlotSubtitle="",
  dfImpulseResults=dfImpulseResults,
  vecRefPval=vecRefResults,
  strNameMethod2="edge",
  strMode=strMode, 
  NPARAM=NPARAM)

if(length(vecDEgenes_Impulse_only)){
  plotDEGenes(vecGeneIDs=vecDEgenes_Impulse_only,
    matCountDataProc=matCountDataProc,
    lsMatTranslationFactors=lsMatTranslationFactors,
    matSizeFactors=matSizeFactors,
    dfAnnotation=dfAnnotationProc, 
    lsImpulseFits=lsImpulseFits,
    strCaseName=strCaseName, 
    strControlName=strControlName, 
    strFileNameSuffix="RNAseq_ImpulseDE2BeatsEdge", 
    strPlotTitleSuffix="", 
    strPlotSubtitle="",
    dfImpulseResults=dfImpulseResults,
    vecRefPval=vecRefResults,
    strNameMethod2="edge",
    strMode=strMode, 
    NPARAM=NPARAM)
}
