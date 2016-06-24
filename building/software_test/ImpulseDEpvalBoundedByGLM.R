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

matDataA <- matDataA[,as.vector(dfAnnotationA$Sample)]
colnames(matDataA) <- dfAnnotationA$Sample
#only case
matDataA <- matDataA[,dfAnnotationA$Condition=="case"]
colnames(matDataA) <- dfAnnotationA[dfAnnotationA$Condition=="case",]$Sample
dfAnnotationA <- dfAnnotationA[dfAnnotationA$Condition=="case",]
rownames(dfAnnotationA) <- dfAnnotationA$Sample


# filter
matDataA[!is.finite(matDataA)] <- NA
vecboolNonzeroA <- apply(matDataA,1,function(gene){any(gene>0 & !is.na(gene))})
matDataA <- matDataA[vecboolNonzeroA,]

# Run ImpulseDE2
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/building/code_files/ImpulseDE2_main.R")
setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/software_test_out")
strControlName = NULL
strCaseName = "case"
scaSmallRun <- 2000
lsImpulseDE_resultsA <- runImpulseDE2(matCountData=matDataA, 
  dfAnnotation=dfAnnotationA,
  strCaseName = strCaseName, 
  strControlName=strControlName, 
  strMode="batch",
  nProc=2, 
  Q_value=10^(-3),
  boolPlotting=FALSE,
  scaSmallRun=scaSmallRun)
dfImpulseResultsA <- lsImpulseDE_resultsA$dfImpulseResults
dfDESeq2ResultsA <- lsImpulseDE_resultsA$dfDESeq2Results
load("ImpulseDE2_matSizeFactors.RData")
vecSizeFactors <- matSizeFactors[1,]

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = matDataA,
  colData = dfAnnotationA,
  design = ~ TimeCateg)
# Run DESeq2
dds <- estimateSizeFactors(dds)
vecSizeFacDESeq2 <- sizeFactors(dds)
vecSizeFactors
vecSizeFacDESeq2
#ddsDESeqObject <- DESeq(dds, test = "LRT", 
#  full = ~ TimeCateg, 
#  reduced = ~ 1, 
#  parallel=TRUE)
#dds_dispersions <- 1/dispersions(ddsDESeqObject) 
#dds_resultsTable <- results(ddsDESeqObject)

vecLogLikFull <- array(NA,scaSmallRun)
vecLogLikRed <- array(NA,scaSmallRun)
vecPvalue <- array(NA,scaSmallRun)
scaDegFreedomFull <- 7
scaDegFreedomRed <- 2
scaDeltaDegFreedom <- scaDegFreedomFull - scaDegFreedomRed

matDataProc <- matDataA[,as.vector(dfAnnotationA[dfAnnotationA$Condition=="case",]$Sample)]
vecDisp <- as.numeric(as.vector(lsImpulseDE_resultsA$dfImpulseResults$size))
names(vecDisp) <- as.vector(lsImpulseDE_resultsA$dfImpulseResults$Gene)
vecDisp <- vecDisp[rownames(matDataProc[1:scaSmallRun,])]
dfAnnotationProc <- dfAnnotationA[dfAnnotationA$Condition=="case",]

vecMeans <- sapply(seq(1,scaSmallRun),function(i){
  exp(unlist(optim(
      par=log(mean(matDataProc[i,]/vecSizeFacDESeq2, na.rm=TRUE)+1),
      fn=evalLogLikNBMean_comp,
      vecCounts=matDataProc[i,],
      scaDispEst=vecDisp[i], 
      vecNormConst=vecSizeFacDESeq2,
      vecboolObserved=!is.na(matDataProc[i,]),
      method="BFGS",
      control=list(maxit=1000,fnscale=-1)
    )["par"]))
  })
matImpMeans <- lsImpulseDE_resultsA$lsImpulseFits$values_case[rownames(matDataProc)[1:scaSmallRun],]
for(i in seq(1,scaSmallRun)){
  # deterministic mean fit
  vecDetMeans <- sapply(unique(dfAnnotationProc$TimeCateg), function(tp){
    #mean((matDataProc[i,]/vecSizeFacDESeq2)[dfAnnotationProc$TimeCateg==tp], na.rm=TRUE)
    exp(unlist(optim(
      par=log(mean(matDataProc[i,dfAnnotationProc$TimeCateg==tp]/vecSizeFacDESeq2[dfAnnotationProc$TimeCateg==tp], na.rm=TRUE)+1),
      fn=evalLogLikNBMean_comp,
      vecCounts=matDataProc[i,dfAnnotationProc$TimeCateg==tp],
      scaDispEst=vecDisp[i], 
      vecNormConst=vecSizeFacDESeq2[dfAnnotationProc$TimeCateg==tp],
      vecboolObserved=!is.na(matDataProc[i,dfAnnotationProc$TimeCateg==tp]),
      method="BFGS",
      control=list(maxit=1000,fnscale=-1)
    )["par"]))
  })
  vecDetMeans <- c(vecDetMeans,vecDetMeans)
  vecLogLikFull[i] <- sum(dnbinom(matDataProc[i,],mu=vecDetMeans*vecSizeFacDESeq2,size=vecDisp[rownames(matDataProc)[i]],log=TRUE))
  vecLogLikRed[i] <- sum(dnbinom(matDataProc[i,],mu=vecMeans[i]*vecSizeFacDESeq2,size=vecDisp[rownames(matDataProc)[i]],log=TRUE))
  
  # Compute test statistic: Deviance
  scaDeviance <- 2*(vecLogLikFull[i] - vecLogLikRed[i])
  # Get p-values from Chi-square distribution (assumption about null model)
  vecPvalue[i] <- pchisq(scaDeviance,df=scaDeltaDegFreedom,lower.tail=FALSE) 
}
vecPvalueImp <- as.numeric(as.vector(dfImpulseResultsA[rownames(matDataProc)[1:scaSmallRun],"p"]))
vecPvalueDES <- dfDESeq2ResultsA[rownames(matDataProc)[1:scaSmallRun],"pvalue"]
dfLL <- data.frame(
  imppvalbiggerglm=vecPvalue<=vecPvalueImp+vecPvalueImp/10^3,
  imppvalbiggerdeseq2=vecPvalueDES<=vecPvalueImp+vecPvalueImp/10^3,
  deseq2pvalbiggerglm=vecPvalue<=vecPvalueDES+vecPvalueDES/10^3,
  ratioimppvaltoglm=round(as.numeric(as.vector(dfImpulseResultsA[rownames(matDataProc)[1:scaSmallRun],"p"]))/vecPvalue,4),
  nullLLequal=(vecLogLikRed >= 1.0000001*as.numeric(as.vector(dfImpulseResultsA[rownames(matDataProc)[1:scaSmallRun],"loglik_red"])) &
      vecLogLikRed <= 0.999999*as.numeric(as.vector(dfImpulseResultsA[rownames(matDataProc)[1:scaSmallRun],"loglik_red"]))),
  glmpval=vecPvalue,
  deseq2pval=dfDESeq2ResultsA[rownames(matDataProc)[1:scaSmallRun],"pvalue"],
  imppval=as.numeric(as.vector(dfImpulseResultsA[rownames(matDataProc)[1:scaSmallRun],"p"])),
  glmfull=vecLogLikFull,
  glmred=vecLogLikRed,
  glmmean=vecMeans,
  impfull=as.numeric(as.vector(dfImpulseResultsA[rownames(matDataProc)[1:scaSmallRun],"loglik_full"])),
  impred=as.numeric(as.vector(dfImpulseResultsA[rownames(matDataProc)[1:scaSmallRun],"loglik_red"])),
  impmean=as.numeric(as.vector(dfImpulseResultsA[rownames(matDataProc)[1:scaSmallRun],"mean"])),
  zeros=apply(matDataProc[1:scaSmallRun,],1,function(gene){sum(gene==0)}))
head(dfLL)

print("ImpulseDE2 pval bigger than optimal GLM - always - bound fulfilled")
all(dfLL$imppvalbiggerglm)
print("ImpulseDE2 pval bigger than DESeq2 - not always - DESeq2 failures?")
all(dfLL$imppvalbiggerdeseq2)
print("Failures all involve zero observations")
all(dfLL[!dfLL$imppvalbiggerdeseq2,]$zeros>0)
print("all violations have a time point completely zero - this seems to break deseq2")
min(dfLL[!dfLL$imppvalbiggerdeseq2,]$zeros)
all(apply(matDataProc[rownames(dfLL[!dfLL$imppvalbiggerdeseq2,]),],1,function(gene){
  any(sapply(unique(dfAnnotationProc$TimeCateg),function(tp){
    all(gene[dfAnnotationProc$TimeCateg==tp]==0)
  }))
}))
