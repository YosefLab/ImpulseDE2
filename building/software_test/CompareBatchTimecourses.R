# Compare batch against timecourses
print("Run after completion of DESeq")
rm(list = ls())
########################################
# Load data

# Load timecourses
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/R_data/fullrun_timecourses")
#setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
load("ImpulseDE2_dfAnnotationRed.RData")
# Load Impulse output
load("ImpulseDE2_dfImpulseResults.RData")

load("ImpulseDE2_arr3DCountData.RData")
load("ImpulseDE2_lsDEGenes.RData")
load("ImpulseDE2_lsImpulseFits.RData")
load("ImpulseDE2_dfDESeq2Results.RData")

dfImpulseResultsByTC <- dfImpulseResults

# Load batch
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/R_data/fullrun_batch")
load("ImpulseDE2_dfAnnotationRed.RData")
# Load Impulse output
load("ImpulseDE2_dfImpulseResults.RData")

dfImpulseResultsBatch <- dfImpulseResults

setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
# Compare against DESeq
# Make summary table to compare against plots
dfBatch_ByTC <- as.data.frame( cbind(
  "Gene"=as.character(as.vector( dfImpulseResultsBatch$Gene )),
  "Batch"=dfImpulseResultsBatch[as.character(as.vector( dfImpulseResultsBatch$Gene )),]$adj.p,
  "ByTC"=dfImpulseResultsByTC[as.character(as.vector( dfImpulseResultsBatch$Gene )),]$adj.p),
  stringsAsFactors=FALSE)
rownames(dfBatch_ByTC) <- as.character(as.vector( dfImpulseResultsBatch$Gene ))
dfBatch_ByTC$Batch <- as.numeric(dfBatch_ByTC$Batch)
dfBatch_ByTC$ByTC <- as.numeric(dfBatch_ByTC$ByTC)

########################################
# 1. Heatmap: Q-value thresholds
mat_overlap <- array(NA,c(11,11))
mat_intersect <- array(NA,c(11,11))
mat_union <- array(NA,c(11,11))
# DESeq on vertical
for(i in 0:10){
  # Impulse on horicontal
  for(j in 0:10){
    sig_batch <- dfBatch_ByTC$Batch <= 10^(-i)
    sig_bytc <- dfBatch_ByTC$ByTC <= 10^(-j)
    mat_overlap[i+1,j+1] <- sum(sig_batch & sig_bytc)/sum(sig_batch | sig_bytc)
    mat_intersect[i+1,j+1] <- sum(sig_batch & sig_bytc)
    mat_union[i+1,j+1] <- sum(sig_batch | sig_bytc)
  }
}
rownames(mat_overlap) <- 0:-10
colnames(mat_overlap) <- 0:-10
rownames(mat_intersect) <- 0:-10
colnames(mat_intersect) <- 0:-10
rownames(mat_union) <- 0:-10
colnames(mat_union) <- 0:-10
library(gplots)
graphics.off()
breaks <- seq(0,1,by=0.01)
hm.colors <- colorpanel( length(breaks)-1, "yellow", "red" )
graphics.off()
pdf(paste('DESeq-Impulse_JaccardCoeff-SignificantGenes_Heatmap.pdf',sep=''),width=7,height=7)
heatmap.2(mat_overlap, dendrogram="none", Rowv=FALSE,Colv=FALSE, 
  xlab =  paste0("log(p-value) ImpulseDE2 Batch"), ylab = "log(p-value) ImpulseDE2 By Timecourse",
  breaks=breaks,col=hm.colors, scale="none",
  trace="none",density.info="none",
  key.title = " ", key.xlab = paste0("Jaccard coefficient"), key.ylab = NULL,
  symkey=FALSE,
  cellnote=round(mat_overlap,digits=2),notecol="grey",
  lmat=rbind( c(3,4),c(2,1) ),lhei=c(1,4), lwid=c(1,4), margins=c(5,5))
dev.off()
print("Intersection")
print(mat_intersect)
print("Union")
print(mat_union)

########################################
# 2. Impulse fits, differentially called DE genes

setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files")
# Plot the impulse fits and input data
source("ImpulseDE2_main.R")
source("srcImpulseDE2_plotDEGenes.R")
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")

Q = 10^(-3)
NPARAM=6
strMode="timecourses" # for plotting style

DEgenes_both <- dfBatch_ByTC[(dfBatch_ByTC[,"Batch"] < Q & dfBatch_ByTC[,"ByTC"] < Q),"Gene"]
DEgenes_Batch_only <- dfBatch_ByTC[(dfBatch_ByTC[,"Batch"] < Q & dfBatch_ByTC[,"ByTC"] >= Q),"Gene"]
DEgenes_ByTC_only <- dfBatch_ByTC[(dfBatch_ByTC[,"Batch"] >= Q & dfBatch_ByTC[,"ByTC"] < Q),"Gene"]
graphics.off()
print(paste0("Number of significant DE genes at q-value of ",Q))
print(paste0("Batch only ",length(DEgenes_Batch_only)))
print(paste0("ByTC only ",length(DEgenes_ByTC_only)))
print(paste0("Both ",length(DEgenes_both)))

vecByTCResults <- dfBatch_ByTC$ByTC
names(vecByTCResults) <- dfBatch_ByTC$Gene

plotDEGenes(lsGeneIDs=DEgenes_both,
  arr3DCountData=arr3DCountData, dfAnnotationRed=dfAnnotationRed, 
  lsImpulseFits=lsImpulseFits,
  dfImpulseResults=dfImpulseResultsBatch,vecMethod2Results=vecByTCResults,strMode=strMode,
  strCaseName="case", strControlName=NULL, 
  strFileNameSuffix="DE_both", strPlotTitleSuffix="", strPlotSubtitle="",
  strNameMethod1="Batch", strNameMethod2="ByTC",
  NPARAM=NPARAM)

# sort DESeq2 only genes by padj of DESeq2
dfBatch_ByTC_ByTCOnly <- dfBatch_ByTC[DEgenes_ByTC_only,]
DEgenes_ByTC_onlySorted <- dfBatch_ByTC_ByTCOnly[order(dfBatch_ByTC_ByTCOnly$ByTC),]$Gene
plotDEGenes(lsGeneIDs=DEgenes_ByTC_onlySorted,
  arr3DCountData=arr3DCountData, dfAnnotationRed=dfAnnotationRed, 
  lsImpulseFits=lsImpulseFits,
  dfImpulseResults=dfImpulseResultsBatch,vecMethod2Results=vecByTCResults,strMode=strMode,
  strCaseName="case", strControlName=NULL, 
  strFileNameSuffix="DE_ByTC_only", strPlotTitleSuffix="", strPlotSubtitle="",
  strNameMethod1="Batch", strNameMethod2="ByTC",
  NPARAM=NPARAM)

plotDEGenes(lsGeneIDs=DEgenes_Batch_only,
  arr3DCountData=arr3DCountData, dfAnnotationRed=dfAnnotationRed, 
  lsImpulseFits=lsImpulseFits,
  strCaseName="case", strControlName=NULL,
  dfImpulseResults=dfImpulseResultsBatch,vecMethod2Results=vecByTCResults,strMode=strMode,
  strFileNameSuffix="DE_Batch_only", strPlotTitleSuffix="", strPlotSubtitle="",
  strNameMethod1="Batch", strNameMethod2="ByTC",
  NPARAM=NPARAM)