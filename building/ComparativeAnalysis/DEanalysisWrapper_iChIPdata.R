rm(list = ls())

print("Process data iChIP")
# Load Data set iChIP
# 1. Counts
dfH3K4me1 <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/H3K4me1_EryLineage.tab",
  sep="\t",header=T,colClasses=c(
    "character","character","character",
    "numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric"))
vecPeaksIDs <- apply(dfH3K4me1[,1:3],1,function(coord){paste(coord,collapse="-")})
matDataA <- dfH3K4me1[,4:dim(dfH3K4me1)[2]]
matDataA <- t(apply(matDataA,1,as.numeric))

rownames(matDataA) <- vecPeaksIDs
colnames(matDataA) <- colnames(dfH3K4me1)[4:dim(dfH3K4me1)[2]]

# 2. Annotation
dfAnnotationiChIP <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/AnnotationTable_iChIPdata_H3K4me1_EryLineage.tab",header=T)
dfAnnotationiChIP$TimeCateg <- paste0(rep("_",length(dfAnnotationiChIP$Time)),dfAnnotationiChIP$Time)
rownames(dfAnnotationiChIP) <- dfAnnotationiChIP$Sample

dfAnnotationA <- dfAnnotationiChIP

# All zero rows
vecboolNonzeroA <- apply(matDataA,1,function(gene){any(gene>0)})
matDataA <- matDataA[vecboolNonzeroA,]

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
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/building/code_files/ImpulseDE2_main.R")

# Create input data set
# Only retain non zero
matDataA_ImpulseDE2 <- matDataA

tm_ImpulseDE2A <- system.time({
  setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/ImpulseDE2/H3K4me1")
  strControlName = NULL
  strCaseName = "case"
  lsImpulseDE_resultsA <- runImpulseDE2(matCountData=matDataA_ImpulseDE2, 
    dfAnnotation=dfAnnotationA,
    strCaseName = strCaseName, 
    strControlName=strControlName, 
    strMode="batch",
    nProc=1, 
    Q_value=10^(-3),
    boolPlotting=TRUE)
  dfImpulseResultsA <- lsImpulseDE_resultsA$dfImpulseResults
  qvals_A <- dfImpulseResultsA$adj.p
  pvals_A <- dfImpulseResultsA$p
})

# Extract Results
matQval_iChIPData[rownames(dfImpulseResultsA),"A_ImpulseDE2"] <- qvals_A
matPval_iChIPData[rownames(dfImpulseResultsA),"A_ImpulseDE2"] <- pvals_A

matRunTime_iChIPData["ImpulseDE2"] <- round(tm_ImpulseDE2A["elapsed"])

lsResDEcomparison_iChIPdataH3K4me1 <- list("matQval_iChIPData"=matQval_iChIPData, 
  "matPval_iChIPData"=matPval_iChIPData, 
  "matRunTime_iChIPData"=matRunTime_iChIPData)

save(lsResDEcomparison_iChIPdataH3K4me1,file=file.path("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/lsResDEcomparison_iChIPdataH3K4me1_EryLineage.RData"))
print("Finished ImpulseDE2")

# Currently no comparative analysis
if(FALSE){
  ################################################################################
  # ImpulseDE
  print("Run ImpulseDE")
  source("/Users/davidsebastianfischer/MasterThesis/software/ImpulseDE/R/Impulse_DE_fin.R")
  
  # Create input data set
  # Only retain non zero
  matDataA_ImpulseDE <- matDataA
  
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
  })
  
  # Extract Results
  matQval_iChIPData[rownames(dds_resultsTableA),"A_DESeq2"] <- qvals_A
  matPval_iChIPData[rownames(dds_resultsTableA),"A_DESeq2"] <- pvals_A
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
  matDataA_EDGE <- matDataA
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
    #edge_odp <- odp(edge_norm, bs.its = 30, verbose=FALSE)
    # likelihood ratio test: throws error
    edge_lrt <- lrt(edge_sva)
    
    # Extract results
    qval_obj_A <- qvalueObj(edge_lrt)
    qvals_A <- qval_obj_A$qvalues
    pvals_A <- qval_obj_A$pvalues
    lfdr_A <- qval_obj_A$lfdr
    pi0_A <- qval_obj_A$pi0
  })
  
  # Extract Results
  matQval_iChIPData[rownames(edge_lrt),"A_edge"] <- qvals_A
  matPval_iChIPData[rownames(edge_lrt),"A_edge"] <- pvals_A
  matRunTime_iChIPData["edge"] <- round(tm_edgeA["elapsed"])
  
  lsResDEcomparison_iChIPdataH3K4me1 <- list("matQval_iChIPData"=matQval_iChIPData, 
    "matPval_iChIPData"=matPval_iChIPData, 
    "matRunTime_iChIPData"=matRunTime_iChIPData)
  
  save(lsResDEcomparison_iChIPdataH3K4me1,file=file.path("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/lsResDEcomparison_iChIPdataH3K4me1.RData"))
  print("Finished edge.")
  
}

# method comparison
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/ImpulseDE2/H3K4me1")
load("ImpulseDE2_matCountDataProc.RData")
load("ImpulseDE2_dfAnnotationProc.RData")
load("ImpulseDE2_dfImpulseResults.RData")
load("ImpulseDE2_vecDEGenes.RData")
load("ImpulseDE2_lsImpulseFits.RData")
load("ImpulseDE2_dfDESeq2Results.RData")
load("ImpulseDE2_lsMatTranslationFactors.RData")
load("ImpulseDE2_vecDispersions.RData")
load("ImpulseDE2_matSizeFactors.RData")
setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/building/code_files")
source("ImpulseDE2_main.R")
source("srcImpulseDE2_plotDEGenes.R")
source("srcImpulseDE2_compareDEMethods.R")
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/pdfs")

Q <- 10^(-3)
Qdelta <- 10^(2) # difference factor required to be plotted

t <- system.time({
for(i in seq(1,100)){
  a <- fitNBMean(vecCounts=matCountDataProc[i,],
      scaDispEst=vecDispersions,
      vecNormConst=matSizeFactors[1,])
}
})

# Load run data ImpulseDE2, ImpulseDE, DESeq2
matQval <- as.matrix(data.frame( Gene=rownames(dfImpulseResults), ImpulseDE2=dfImpulseResults$adj.p, DESeq2=NA))
rownames(matQval) <- rownames(dfImpulseResults)
matQval[,"DESeq2"] <- dfDESeq2Results[rownames(matQval),"padj"]

colnames(matQval) <- c("Gene", "ImpulseDE2", "DESeq2")

setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/pdfs")
compareDEMethods(matQval,
  strMethod2="DESeq2",
  Q = Q,
  Qdelta = Qdelta,
  matCountDataProc = matCountDataProc,
  lsMatTranslationFactors = lsMatTranslationFactors,
  matSizeFactors = matSizeFactors,
  dfAnnotationProc = dfAnnotationProc, 
  lsImpulseFits = lsImpulseFits,
  dfImpulseResults = dfImpulseResults,
  strCaseName="case", 
  strControlName = NULL, 
  strMode="batch",
  strDataDescriptionFilename="")

################################################################################
# DESeq2 for Gene-E
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
dfDESeq2DEgeneCounts <- matDataA[rownames(dds_resultsTableA[dds_resultsTableA$padj < 10^(-5),]),]
rownames(dfDESeq2DEgeneCounts) <- rownames(dds_resultsTableA[dds_resultsTableA$padj < 10^(-5),])
dfDESeq2DEgeneCountsNorm <- dfDESeq2DEgeneCounts/rep(sizeFactors(ddsDESeqObjectA), each = nrow(dfDESeq2DEgeneCounts))
dfDESeq2DEgeneCountsNormMeans <- do.call(cbind, lapply(unique(dfAnnotationA$Time),function(t){
  if(length(as.vector(dfAnnotationA[dfAnnotationA$Time==t,]$Sample)) > 1 ){
    apply(dfDESeq2DEgeneCountsNorm[,as.vector(dfAnnotationA[dfAnnotationA$Time==t,]$Sample)],1,mean)
  }else{
    dfDESeq2DEgeneCountsNorm[,as.vector(dfAnnotationA[dfAnnotationA$Time==t,]$Sample)]
  }
}))
colnames(dfDESeq2DEgeneCountsNormMeans) <- unique(dfAnnotationA$TimeCateg)
write.table(dfDESeq2DEgeneCountsNormMeans,"/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/matCounts_DESeq2_1e-5.tab",sep="\t",row.names=F,quote=FALSE)
# retain only high variation
dfDESeq2DEgeneCountsNormMeansFilt <- dfDESeq2DEgeneCountsNormMeans[2*apply(dfDESeq2DEgeneCountsNormMeans,1,min) <= apply(dfDESeq2DEgeneCountsNormMeans,1,max),]
write.table(dfDESeq2DEgeneCountsNormMeansFilt,"/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/2014_Weiner_ChromatinHaematopeisis/matCounts_DESeq2_1e-5_x1-5.tab",sep="\t",row.names=F,quote=FALSE)
