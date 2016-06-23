### GENE-E clustering
# Dont do this - the data set is to confounded to look at means?
rm(list=ls())
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/clusterruns/output/ImpulseDE2")
load("ImpulseDE2_matCountDataProc.RData")
load("ImpulseDE2_dfAnnotationProc.RData")
load("ImpulseDE2_dfImpulseResults.RData")
load("ImpulseDE2_vecDEGenes.RData")
load("ImpulseDE2_lsImpulseFits.RData")
load("ImpulseDE2_dfDESeq2Results.RData")
load("ImpulseDE2_lsMatTranslationFactors.RData")
load("ImpulseDE2_matSizeFactors.RData")
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/pdfs")
# DESeq2 for Gene-E
dfDESeq2DEgeneCounts <- matCountDataProc[rownames(dfDESeq2Results[dfDESeq2Results$padj < 10^(-5) & !is.na(dfDESeq2Results$padj),]),]
dfDESeq2DEgeneCountsNorm <- dfDESeq2DEgeneCounts/rep(matSizeFactors[1,], each = nrow(dfDESeq2DEgeneCounts))
dfDESeq2DEgeneCountsNormMeans <- do.call(cbind, lapply(unique(dfAnnotationProc$Time),function(t){
  if(length(as.vector(dfAnnotationProc[dfAnnotationProc$Time==t,]$Sample)) > 1 ){
    apply(dfDESeq2DEgeneCountsNorm[,as.vector(dfAnnotationProc[dfAnnotationProc$Time==t,]$Sample)],1,mean)
  }else{
    dfDESeq2DEgeneCountsNorm[,as.vector(dfAnnotationProc[dfAnnotationProc$Time==t,]$Sample)]
  }
}))
colnames(dfDESeq2DEgeneCountsNormMeans) <- unique(dfAnnotationProc$TimeCateg)
write.table(dfDESeq2DEgeneCountsNormMeans,"/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/matCounts_DESeq2_1e-5.tab",sep="\t",row.names=F,quote=FALSE)
# retain only high variation
dfDESeq2DEgeneCountsNormMeansFilt <- dfDESeq2DEgeneCountsNormMeans[2*apply(dfDESeq2DEgeneCountsNormMeans,1,min) <= apply(dfDESeq2DEgeneCountsNormMeans,1,max),]
write.table(dfDESeq2DEgeneCountsNormMeansFilt,"/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/matCounts_DESeq2_1e-5_x1-5.tab",sep="\t",row.names=F,quote=FALSE)


########################################
# 1. Plot ImpulseDE2 DE genes
########################################
rm(list=ls())

setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/building/code_files")
source("ImpulseDE2_main.R")
source("srcImpulseDE2_plotDEGenes.R")
source("srcImpulseDE2_compareDEMethods.R")
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/clusterruns/output/ImpulseDE2")
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
# 2. Graphical comparison
########################################
## Full data set
rm(list=ls())
library(gplots)
#load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/edge/lsResDEcomparison_RNAseqdata.Rdata")
#matQval_edge <- lsResDEcomparison_RNAseqdata$matQval_RNAseqData 
#matRunTime_edge <- lsResDEcomparison_RNAseqdata$matRunTime_RNAseqData

setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE2/building/code_files")
source("ImpulseDE2_main.R")
source("srcImpulseDE2_plotDEGenes.R")
source("srcImpulseDE2_compareDEMethods.R")
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/clusterruns/output/ImpulseDE2")
load("ImpulseDE2_matCountDataProc.RData")
load("ImpulseDE2_dfAnnotationProc.RData")
load("ImpulseDE2_dfImpulseResults.RData")
load("ImpulseDE2_vecDEGenes.RData")
load("ImpulseDE2_lsImpulseFits.RData")
load("ImpulseDE2_dfDESeq2Results.RData")
load("ImpulseDE2_lsMatTranslationFactors.RData")
load("ImpulseDE2_matSizeFactors.RData")

Q <- 10^(-3)
Qdelta <- 10^(2) # difference factor required to be plotted

# Load run data ImpulseDE2, ImpulseDE, DESeq2
load("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/clusterruns/output/lsResDEcomparison_RNAseqdata.RData")
matQval <- lsResDEcomparison_RNAseqdata$matQval_RNAseqData
matRunTime <- lsResDEcomparison_RNAseqdata$matRunTime_RNAseqData

colnames(matQval) <- c("Gene", sapply(colnames(matQval),function(str){unlist(strsplit(str,"A_"))[2]})[2:5])
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/pdfs/full")

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
  strMode="longitudinal",
  strDataDescriptionFilename="")

compareDEMethods(matQval,
  strMethod2="ImpulseDE",
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
  strMode="longitudinal",
  strDataDescriptionFilename="")

strMethod1="ImpulseDE2"
strMethod2="DESeq2"
Q <- 10^(-2)
Qdelta <- 2 # difference factor required
# Note cannot plot genes which are NA in ImpulseDE2
vecboolImpulseFitted <- !is.na(as.numeric(matQval[,strMethod2]))
vecboolRefBeatsQ <- as.numeric(matQval[,strMethod2]) <= Q & !is.na(as.numeric(matQval[,strMethod2]))
vecBoolRefBeatsImpulse <- as.numeric(matQval[,strMethod1]) >= Qdelta*as.numeric(matQval[,strMethod2]) | (is.na(as.numeric(matQval[,strMethod1])) & !is.na(as.numeric(matQval[,strMethod2])))
vecboolImpulseBeatsQ <- as.numeric(matQval[,strMethod1]) <= Q & !is.na(as.numeric(matQval[,strMethod1]))
vecboolImpulseBeatsRef <- as.numeric(matQval[,strMethod2]) >= Qdelta*as.numeric(matQval[,strMethod1]) | (!is.na(as.numeric(matQval[,strMethod1])) & is.na(as.numeric(matQval[,strMethod2])))
# Look at zeros
vecTimepoints <- unique(dfAnnotationProc$TimeCateg)
scaNTP <- length(vecTimepoints)
vecTPAssign <- match(dfAnnotationProc$TimeCateg,vecTimepoints)
vecboolAllZeroSamples <- apply(matCountDataProc,1,function(gene){
  any(sapply(seq(1,scaNTP),function(tp){
    all(gene[vecTPAssign==tp & !is.na(gene)]==0)
  }))
})
vecboolNoZeros <- apply(matCountDataProc,1,function(gene){
  all(gene[!is.na(gene)]>0)
})
print(paste0("Number of genes: ",length(vecboolAllZeroSamples)))
print(paste0("Number of genes without all zero sample: ",sum(!vecboolAllZeroSamples)))

print(paste0("Number of significant DE genes at q-value of ",Q))
print(paste0(strMethod1," only ",sum(vecboolImpulseBeatsRef, na.rm=T)))
print(paste0(strMethod2," only ",sum(vecBoolRefBeatsImpulse, na.rm=T)))
print(paste0("both ",sum(vecboolImpulseBeatsQ & vecboolRefBeatsQ, na.rm=T)))
print("Only genes without all zero samples:")
print(paste0("Number of significant DE genes at q-value of ",Q))
print(paste0(strMethod1," only ",sum(vecboolImpulseBeatsRef & !vecboolAllZeroSamples, na.rm=T)))
print(paste0(strMethod2," only ",sum(vecBoolRefBeatsImpulse & !vecboolAllZeroSamples, na.rm=T)))
print(paste0("both ",sum(vecboolImpulseBeatsQ & vecboolRefBeatsQ & !vecboolAllZeroSamples, na.rm=T)))
print("Only genes without zero observations:")
print(paste0("Number of significant DE genes at q-value of ",Q))
print(paste0(strMethod1," only ",sum(vecboolImpulseBeatsRef & vecboolNoZeros, na.rm=T)))
print(paste0(strMethod2," only ",sum(vecBoolRefBeatsImpulse & vecboolNoZeros, na.rm=T)))
print(paste0("both ",sum(vecboolImpulseBeatsQ & vecboolRefBeatsQ & vecboolNoZeros, na.rm=T)))

# Plot CDF for 3
# 2. CDF p-values
setwd("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/pdfs")
strMethod3="ImpulseDE"
strMethod4="edge"
vecX <- seq(-25,-1,by=0.1)
vecCDF1 <- sapply(vecX, function(thres){sum(log(as.numeric(matQval[,strMethod1]))/log(10) <= thres, na.rm=TRUE)})
vecCDF2 <- sapply(vecX, function(thres){sum(log(as.numeric(matQval[,strMethod2]))/log(10) <= thres, na.rm=TRUE)})
vecCDF3 <- sapply(vecX, function(thres){sum(log(as.numeric(matQval[,strMethod3]))/log(10) <= thres, na.rm=TRUE)})
vecCDF4 <- sapply(vecX, function(thres){sum(log(as.numeric(matQval[,strMethod4]))/log(10) <= thres, na.rm=TRUE)})
pdf(paste0(strMethod1,"-",strMethod2,"-",strMethod3,"_ECDF-pvalues.pdf"),width=7,height=7)
plot(vecX,vecCDF1,
  col="black",pch=4,type="l",
  ylim=c(0,max(max(vecCDF1,na.rm=TRUE),max(vecCDF2,na.rm=TRUE))),
  xlab="-log_10(p-value)", 
  ylab=paste0("Cumulative p-value distribution"),
  main=paste0("Cumulative p-values distribution\n",strMethod1," versus ",strMethod2))
points(vecX,vecCDF2,
  col="red",pch=4,type="l")
points(vecX,vecCDF3,
  col="blue",pch=4,type="l")
points(vecX,vecCDF3,
  col="green",pch=4,type="l")
legend(x="topleft",
  legend=c(strMethod1,strMethod2,strMethod3,strMethod4),
  fill=c("black","red","blue","green"))
dev.off()
graphics.off()

# Find lowest p-value in DESeq2 which is significantly higher than ImpulseDE2

