rm(list = ls())

print("Process data DyNB")
# Load Data set DyNB
# 1. Counts
matDataA <- as.matrix(read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE_datasets/2014_DyNB/DyNB_Th0_dataset.txt",
  sep=",",header=T,colClasses=c(
    "numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric")))
rownames(matDataA) <- paste0(rep("Gene_",dim(matDataA)[1]),seq(1,dim(matDataA)[1]),sep="")
matDataB <- as.matrix(read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE_datasets/2014_DyNB/DyNB_Th17_dataset.txt",
  sep=",",header=T,colClasses=c(
    "numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric")))
rownames(matDataB) <- paste0(rep("Gene_",dim(matDataB)[1]),seq(1,dim(matDataB)[1]),sep="")
# 2. Annotation
dfAnnotationDyNB <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE_datasets/2014_DyNB/AnnotationTable_DyNBdata.tab",header=T)
dfAnnotationDyNB$TimeCateg <- paste0(rep("_",length(dfAnnotationDyNB$Time)),dfAnnotationDyNB$Time)
rownames(dfAnnotationDyNB) <- dfAnnotationDyNB$Replicate
# 3. Rename and merge for case control comparison
dfAnnotationA <- dfAnnotationDyNB
dfAnnotationAformerge <- dfAnnotationDyNB
dfAnnotationAformerge$Sample <- paste0(rep("Case_",dim(dfAnnotationAformerge)[1]),dfAnnotationAformerge$Sample)
dfAnnotationAformerge$Replicate <- paste0(rep("Case_",dim(dfAnnotationAformerge)[1]),dfAnnotationAformerge$Replicate)
dfAnnotationAformerge$Timecourse <- paste0(rep("Case_",dim(dfAnnotationAformerge)[1]),dfAnnotationAformerge$Timecourse)

dfAnnotationB <- dfAnnotationDyNB
dfAnnotationBformerge <- dfAnnotationDyNB
dfAnnotationBformerge$Condition <- "ctrl"
dfAnnotationBformerge$Sample <- paste0(rep("Ctrl_",dim(dfAnnotationBformerge)[1]),dfAnnotationBformerge$Sample)
dfAnnotationBformerge$Replicate <- paste0(rep("Ctrl_",dim(dfAnnotationBformerge)[1]),dfAnnotationBformerge$Replicate)
dfAnnotationBformerge$Timecourse <- paste0(rep("Ctrl_",dim(dfAnnotationBformerge)[1]),dfAnnotationBformerge$Timecourse)

dfAnnotationAB <- rbind(dfAnnotationAformerge,dfAnnotationBformerge)

matDataAB <- cbind(matDataA,matDataB)
colnames(matDataAB) <- as.vector(dfAnnotationAB$Replicate)

# All zero rows
vecboolNonzeroA <- apply(matDataA,1,function(gene){any(gene>0)})
vecboolNonzeroB <- apply(matDataB,1,function(gene){any(gene>0)})
vecboolNonzeroAB <- apply(matDataAB,1,function(gene){any(gene>0)})

# Summary matrices
matQval_DyNBData <- array(NA,c(dim(matDataA)[1],4*3))
matPval_DyNBData <- array(NA,c(dim(matDataA)[1],4*3))
matRunTime_DyNBData <- array(NA,c(4,3))

rownames(matQval_DyNBData) <- paste0(rep("Gene_",dim(matQval_DyNBData)[1]),seq(1,dim(matQval_DyNBData)[1]))
colnames(matQval_DyNBData) <- c(paste0(
  c(rep("A_",4), rep("B_",4), rep("AvsB_",4)), 
  rep(c("ImpulseDE2", "ImpulseDE", "DESeq2", "edge"),3)
))

rownames(matPval_DyNBData) <- paste0(rep("Gene_",dim(matPval_DyNBData)[1]),seq(1,dim(matPval_DyNBData)[1]))
colnames(matPval_DyNBData) <- c(paste0(
  c(rep("A_",4), rep("B_",4), rep("AvsB_",4)), 
  rep(c("ImpulseDE2", "ImpulseDE", "DESeq2", "edge"),3)
))

rownames(matRunTime_DyNBData) <- c("ImpulseDE2", "ImpulseDE", "DESeq2", "edge")
colnames(matRunTime_DyNBData) <- c("A", "B", "AvsB")

################################################################################
# ImpulseDE2
print("Run ImpulseDE2")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/ImpulseDE2_main.R")

# Create input data set
# Only retain non zero
matDataA_ImpulseDE2 <- matDataA[vecboolNonzeroA,]
matDataB_ImpulseDE2 <- matDataB[vecboolNonzeroB,]
matDataAB_ImpulseDE2 <- matDataAB[vecboolNonzeroAB,]

tm_ImpulseDE2A <- system.time({
  setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE_datasets/2014_DyNB/ImpulseDE2/Th0")
  strControlName = NULL
  strCaseName = "case"
  lsImpulseDE_resultsA <- runImpulseDE2(matCountData=matDataA_ImpulseDE2, 
    dfAnnotationFull=dfAnnotationA,
    strCaseName = strCaseName, 
    strControlName=strControlName, 
    strMode="timecourses",
    nProc=3, 
    Q_value=10^(-3),
    boolPlotting=FALSE)
  dfImpulseResultsA <- lsImpulseDE_resultsA$dfImpulseResults
  qvals_A <- dfImpulseResultsA$adj.p
  pvals_A <- dfImpulseResultsA$p
})

tm_ImpulseDE2B <- system.time({
  setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE_datasets/2014_DyNB/ImpulseDE2/Th17")
  strControlName = NULL
  strCaseName = "case"
  lsImpulseDE_resultsB <- runImpulseDE2(matCountData=matDataB_ImpulseDE2, 
    dfAnnotationFull=dfAnnotationB,
    strCaseName = strCaseName, 
    strControlName=strControlName, 
    strMode="timecourses",
    nProc=3, 
    Q_value=10^(-3),
    boolPlotting=FALSE)
  dfImpulseResultsB <- lsImpulseDE_resultsB$dfImpulseResults
  qvals_B <- dfImpulseResultsB$adj.p
  pvals_B <- dfImpulseResultsB$p
})

tm_ImpulseDE2AB <- system.time({
  setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE_datasets/2014_DyNB/ImpulseDE2/Th0vs17")
  strControlName = "ctrl"
  strCaseName = "case"
  lsImpulseDE_resultsAB <- runImpulseDE2(matCountData=matDataAB_ImpulseDE2, 
    dfAnnotationFull=dfAnnotationAB,
    strCaseName = strCaseName, 
    strControlName=strControlName, 
    strMode="timecourses",
    nProc=3, 
    Q_value=10^(-3),
    boolPlotting=FALSE)
  dfImpulseResultsAB <- lsImpulseDE_resultsAB$dfImpulseResults
  qvals_AvsB <- dfImpulseResultsAB$adj.p
  pvals_AvsB <- dfImpulseResultsAB$p
})

# Extract Results
matQval_DyNBData[vecboolNonzeroA,"A_ImpulseDE2"] <- qvals_A
matQval_DyNBData[vecboolNonzeroB,"B_ImpulseDE2"] <- qvals_B
matQval_DyNBData[vecboolNonzeroAB, "AvsB_ImpulseDE2"] <- qvals_AvsB

matPval_DyNBData[vecboolNonzeroA,"A_ImpulseDE2"] <- pvals_A
matPval_DyNBData[vecboolNonzeroB,"B_ImpulseDE2"] <- pvals_B
matPval_DyNBData[vecboolNonzeroAB, "AvsB_ImpulseDE2"] <- pvals_AvsB

matRunTime_DyNBData["ImpulseDE2","A"] <- round(tm_ImpulseDE2A["elapsed"])
matRunTime_DyNBData["ImpulseDE2","B"] <- round(tm_ImpulseDE2B["elapsed"])
matRunTime_DyNBData["ImpulseDE2","AvsB"] <- round(tm_ImpulseDE2AB["elapsed"])

lsResDEcomparison_DyNBdata <- list("matQval_DyNBData"=matQval_DyNBData, 
  "matPval_DyNBData"=matPval_DyNBData, 
  "matRunTime_DyNBData"=matRunTime_DyNBData)

save(lsResDEcomparison_DyNBdata,"/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE_datasets/2014_DyNB/lsResDEcomparison_DyNBdata.RData"))
print("Finished ImpulseDE2")

################################################################################
# ImpulseDE
print("Run ImpulseDE")
source("/Users/davidsebastianfischer/MasterThesis/software/ImpulseDE/repo/Impulse_DE_fin.R")

# Create input data set
# Only retain non zero
matDataA_ImpulseDE <- matDataA[vecboolNonzeroA,]
matDataB_ImpulseDE <- matDataB[vecboolNonzeroB,]
matDataAB_ImpulseDE <- matDataAB[vecboolNonzeroAB,]

tm_ImpulseDE2AB <- system.time({
  setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE_datasets/2014_DyNB/ImpulseDE/Th0")
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

tm_ImpulseDE2AB <- system.time({
  setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE_datasets/2014_DyNB/ImpulseDE/Th17")
  strControlName = NULL
  strCaseName = "case"
  lsImpulseDE_resultsB <- impulse_DE(
    expression_table = matDataB_ImpulseDE, 
    annotation_table = dfAnnotationB,
    colname_time = "Time", colname_condition = "Condition", 
    control_timecourse = FALSE,
    control_name = strControlName, case_name = strCaseName, 
    expr_type = "Seq",
    plot_clusters = FALSE, 
    n_iter = 100, n_randoms = 50000, 
    n_process = 3,
    Q_value = 0.01)
  dfImpulseResultsB <- lsImpulseDE_resultsB$impulse_fit_results
  qvals_B <- dfImpulseResultsB$adj.p
  pvals_B <- dfImpulseResultsB$p
})

tm_ImpulseDE2AB <- system.time({
  setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE_datasets/2014_DyNB/ImpulseDE/Th0vs17")
  strControlName = "ctrl"
  strCaseName = "case"
  lsImpulseDE_resultsAB <- impulse_DE(
    expression_table = matDataAB_ImpulseDE, 
    annotation_table = dfAnnotationAB,
    colname_time = "Time", colname_condition = "Condition", 
    control_timecourse = FALSE,
    control_name = strControlName, case_name = strCaseName, 
    expr_type = "Seq",
    plot_clusters = FALSE, 
    n_iter = 100, n_randoms = 50000, 
    n_process = 3,
    Q_value = 0.01)
  dfImpulseResultsAB <- lsImpulseDE_resultsAB$impulse_fit_results
  qvals_AvsB <- dfImpulseResultsAB$adj.p
  pvals_AvsB <- dfImpulseResultsAB$p
})

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
matDataB_DESeq2 <- matDataB[vecboolNonzeroB,]
matDataAB_DESeq2 <- matDataAB[vecboolNonzeroAB,]

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

tm_DESeq2B <- system.time({
  # Create DESeq2 data object
  dds <- DESeqDataSetFromMatrix(countData = matDataB_DESeq2,
    colData = dfAnnotationB,
    design = ~TimeCateg)
  # Run DESeq2
  ddsDESeqObjectB <- DESeq(dds, test = "LRT", 
    full = ~ TimeCateg, reduced = ~ 1,
    parallel=TRUE)
  
  # DESeq results for comparison
  dds_resultsTableB <- results(ddsDESeqObjectB)
  qvals_B <- dds_resultsTableB$padj
  pvals_B <- dds_resultsTableB$pvalue
})

tm_DESeq2AB <- system.time({
  # Create DESeq2 data object
  dds <- DESeqDataSetFromMatrix(countData = matDataAB_DESeq2,
    colData = dfAnnotationAB,
    design = ~TimeCateg + Condition)
  # Run DESeq2
  ddsDESeqObjectAB <- DESeq(dds, test = "LRT", 
    full = ~ TimeCateg + Condition, reduced = ~ TimeCateg,
    parallel=TRUE)
  
  # DESeq results for comparison
  dds_resultsTableAB <- results(ddsDESeqObjectAB)
  qvals_AvsB <- dds_resultsTableAB$padj
  pvals_AvsB <- dds_resultsTableAB$pvalue
})

# Extract Results
matQval_DyNBData[vecboolNonzeroA,"A_DESeq2"] <- qvals_A
matQval_DyNBData[vecboolNonzeroB,"B_DESeq2"] <- qvals_B
matQval_DyNBData[vecboolNonzeroAB, "AvsB_DESeq2"] <- qvals_AvsB

matPval_DyNBData[vecboolNonzeroA,"A_DESeq2"] <- pvals_A
matPval_DyNBData[vecboolNonzeroB,"B_DESeq2"] <- pvals_B
matPval_DyNBData[vecboolNonzeroAB, "AvsB_DESeq2"] <- pvals_AvsB

matRunTime_DyNBData["DESeq2","A"] <- round(tm_DESeq2A["elapsed"])
matRunTime_DyNBData["DESeq2","B"] <- round(tm_DESeq2B["elapsed"])
matRunTime_DyNBData["DESeq2","AvsB"] <- round(tm_DESeq2AB["elapsed"])

lsResDEcomparison_DyNBdata <- list("matQval_DyNBData"=matQval_DyNBData, 
  "matPval_DyNBData"=matPval_DyNBData, 
  "matRunTime_DyNBData"=matRunTime_DyNBData)

save(lsResDEcomparison_DyNBdata,file=file.path("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE_datasets/2014_DyNB/lsResDEcomparison_DyNBdata.RData"))
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

matDataB_EDGE <- matDataB[vecboolNonzeroB,]
rownames(matDataB_EDGE) <- NULL
colnames(matDataB_EDGE) <- NULL

matDataAB_EDGE <- matDataAB[vecboolNonzeroAB,]
rownames(matDataAB_EDGE) <- NULL
colnames(matDataAB_EDGE) <- NULL

# (I) CASE
# A)
tm_edgeA <- system.time({
  cov <- data.frame(tc = dfAnnotationA$Timecourse, time = dfAnnotationA$Time)
  null_model <- ~tc
  full_model <- ~tc + ns(time, df=4)
  edge_obj <- build_models(data = matDataA_EDGE, cov = cov, null.model = null_model, full.model = full_model)
  
  # Normalise
  edge_norm <- apply_snm(edge_obj, int.var=1:ncol(exprs(edge_obj)), diagnose=FALSE)
  # Adjust for unmodelled variables
  edge_sva <- apply_sva(edge_norm)
  
  # optimal discovery procedure: Chose this one
  edge_odp <- odp(edge_sva, bs.its = 30, verbose=FALSE)
  # likelihood ratio test: throws error
  # edge_lrt <- lrt(edge_sva)
  
  # Extract results
  qval_obj_A <- qvalueObj(edge_odp)
  qvals_A <- qval_obj_A$qvalues
  pvals_A <- qval_obj_A$pvalues
  lfdr_A <- qval_obj_A$lfdr
  pi0_A <- qval_obj_A$pi0
})
# B)
tm_edgeB <- system.time({
  cov <- data.frame(tc = dfAnnotationB$Timecourse, time = dfAnnotationB$Time)
  null_model <- ~tc
  full_model <- ~tc + ns(time, df=4)
  edge_obj <- build_models(data = matDataB_EDGE, cov = cov, null.model = null_model, full.model = full_model)
  
  # Normalise
  edge_norm <- apply_snm(edge_obj, int.var=1:ncol(exprs(edge_obj)), diagnose=FALSE)
  # Adjust for unmodelled variables
  edge_sva <- apply_sva(edge_norm)
  
  # optimal discovery procedure: Chose this one
  edge_odp <- odp(edge_sva, bs.its = 30, verbose=FALSE)
  # likelihood ratio test: throws error
  # edge_lrt <- lrt(edge_sva)
  
  # Extract results
  qval_obj_B <- qvalueObj(edge_odp)
  qvals_B <- qval_obj_B$qvalues
  pvals_B <- qval_obj_B$pvalues
  lfdr_B <- qval_obj_B$lfdr
  pi0_B <- qval_obj_B$pi0
})
# (II) CASE CONTROL
tm_edgeAB <- system.time({
  cov <- data.frame(tc = dfAnnotationAB$Timecourse, 
    time = dfAnnotationAB$Time,
    cond = dfAnnotationAB$Condition )
  # Note: complains if tc is added - near singular model matrix
  null_model <- ~cond + ns(time, df=2, intercept = FALSE)
  full_model <- ~cond + ns(time, df=2, intercept = FALSE) + (cond):ns(time, df=2, intercept = FALSE)
  edge_obj <- build_models(data = matDataAB_EDGE, cov = cov, null.model = null_model, full.model = full_model)
  
  # Normalise
  edge_norm <- apply_snm(edge_obj, int.var=1:ncol(exprs(edge_obj)), diagnose=FALSE)
  # Adjust for unmodelled variables
  edge_sva <- apply_sva(edge_norm)
  
  # optimal discovery procedure: Chose this one
  edge_odp <- odp(edge_sva, bs.its = 30, verbose=FALSE)
  # likelihood ratio test: throws error
  # edge_lrt <- lrt(edge_sva)
  
  # Extract results
  qval_obj_AvsB <- qvalueObj(edge_odp)
  qvals_AvsB <- qval_obj_AvsB$qvalues
  pvals_AvsB <- qval_obj_AvsB$pvalues
  lfdr_AvsB <- qval_obj_AvsB$lfdr
  pi0_AvsB <- qval_obj_AvsB$pi0
})

# Extract Results
matQval_DyNBData[vecboolNonzeroA,"A_edge"] <- qvals_A
matQval_DyNBData[vecboolNonzeroB,"B_edge"] <- qvals_B
matQval_DyNBData[vecboolNonzeroAB, "AvsB_edge"] <- qvals_AvsB

matPval_DyNBData[vecboolNonzeroA,"A_edge"] <- pvals_A
matPval_DyNBData[vecboolNonzeroB,"B_edge"] <- pvals_B
matPval_DyNBData[vecboolNonzeroAB, "AvsB_edge"] <- pvals_AvsB

matRunTime_DyNBData["edge","A"] <- round(tm_edgeA["elapsed"])
matRunTime_DyNBData["edge","B"] <- round(tm_edgeB["elapsed"])
matRunTime_DyNBData["edge","AvsB"] <- round(tm_edgeAB["elapsed"])

lsResDEcomparison_DyNBdata <- list("matQval_DyNBData"=matQval_DyNBData, 
  "matPval_DyNBData"=matPval_DyNBData, 
  "matRunTime_DyNBData"=matRunTime_DyNBData)

save(lsResDEcomparison_DyNBdata,file=file.path("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE_datasets/2014_DyNB/lsResDEcomparison_DyNBdata.RData"))
print("Finished edge.")
