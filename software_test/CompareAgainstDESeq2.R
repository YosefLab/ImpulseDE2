print("Run after completion of DESeq")
rm(list = ls())
########################################
# Load data

# Load files from interior of ImpulseDE
#setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/R_data/fullrun_batch")
#strMode="batch"
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/R_data/fullrun_timecourses")
strMode="timecourses"
#setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
load("ImpulseDE2_arr2DCountData.RData")
load("ImpulseDE2_dfAnnotationFull.RData")
load("ImpulseDE2_vecNormConst.RData")
# Load Impulse output
load("ImpulseDE2_dfImpulseResults.RData")
load("ImpulseDE2_lsDEGenes.RData")
load("ImpulseDE2_lsImpulseFits.RData")
load("ImpulseDE2_dfDESeq2Results.RData")

setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
# Compare against DESeq
# Make summary table to compare against plots
dfDESeq_Impulse <- as.data.frame( cbind(
  "Gene"=as.character(as.vector( dfImpulseResults$Gene )),
  "DESeq"=dfDESeq2Results[as.character(as.vector( dfImpulseResults$Gene )),]$padj,
  "Impulse"=dfImpulseResults[as.character(as.vector( dfImpulseResults$Gene )),]$adj.p),
  stringsAsFactors=FALSE)
rownames(dfDESeq_Impulse) <- as.character(as.vector( dfImpulseResults$Gene ))
dfDESeq_Impulse$DESeq <- as.numeric(dfDESeq_Impulse$DESeq)
# Set NA as not detected:
dfDESeq_Impulse$DESeq[is.na(dfDESeq_Impulse$DESeq)] <- 1
dfDESeq_Impulse$Impulse <- as.numeric(dfDESeq_Impulse$Impulse)

# Compare inferred means
lsIndToCompare <- 1:10
apply(arr2DCountData[lsIndToCompare,],1,mean)
dfDESeq2Results[rownames(arr2DCountData[lsIndToCompare,]),]$baseMean
dfImpulseResults[rownames(arr2DCountData[lsIndToCompare,]),]$mean

# Find gene where ImpulseDE2 has lower p value
ImpulseLower <- dfDESeq_Impulse[dfDESeq_Impulse$Impulse < dfDESeq_Impulse$DESeq,]
pvalratio <- (ImpulseLower$Impulse - ImpulseLower$DESeq)
sortedratio <- sort(pvalratio, decreasing=TRUE)
sortedratio[1:50]
tail(sortedratio)
ImpulseLower[pvalratio %in% sortedratio[1:10],]
outliersIDs <- ImpulseLower[pvalratio %in% sortedratio[1:10],]$Gene
as.data.frame(dfDESeq2Results[outliersIDs,])
dfImpulseResults[outliersIDs,]

########################################
# 1. Heatmap: Q-value thresholds
mat_overlap <- array(NA,c(11,11))
mat_intersect <- array(NA,c(11,11))
mat_union <- array(NA,c(11,11))
# DESeq on vertical
for(i in 0:10){
  # Impulse on horicontal
  for(j in 0:10){
    sig_DESeq <- dfDESeq_Impulse$DESeq <= 10^(-i)
    sig_Impulse <- dfDESeq_Impulse$Impulse <= 10^(-j)
    mat_overlap[i+1,j+1] <- sum(sig_DESeq & sig_Impulse)/sum(sig_DESeq | sig_Impulse)
    mat_intersect[i+1,j+1] <- sum(sig_DESeq & sig_Impulse)
    mat_union[i+1,j+1] <- sum(sig_DESeq | sig_Impulse)
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
heatmap(mat_overlap, keep.dendro = FALSE,Rowv=NA,Colv= "Rowv",symm=FALSE,
  xlab =  paste0("ImpulseDE scores"), ylab = "DESeq scores")
breaks <- seq(0,1,by=0.01)
hm.colors <- colorpanel( length(breaks)-1, "yellow", "red" )
graphics.off()
pdf(paste('DESeq-Impulse_JaccardCoeff-SignificantGenes_Heatmap.pdf',sep=''),width=7,height=7)
heatmap.2(mat_overlap, dendrogram="none", Rowv=FALSE,Colv=FALSE, 
  xlab =  paste0("log(p-value) ImpulseDE"), ylab = "log(p-value) DESeq2",
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

Q = 10^(-5)
NPARAM=6
DEgenes_both <- dfDESeq_Impulse[(dfDESeq_Impulse[,"DESeq"] < Q & dfDESeq_Impulse[,"Impulse"] < Q),"Gene"]
DEgenes_DESeq_only <- dfDESeq_Impulse[(dfDESeq_Impulse[,"DESeq"] < Q & dfDESeq_Impulse[,"Impulse"] >= Q),"Gene"]
DEgenes_Impulse_only <- dfDESeq_Impulse[(dfDESeq_Impulse[,"DESeq"] >= Q & dfDESeq_Impulse[,"Impulse"] < Q),"Gene"]
graphics.off()
print(paste0("Number of significant DE genes at q-value of ",Q))
print(paste0("ImpulseDE only ",length(DEgenes_Impulse_only)))
print(paste0("DESeq2 only ",length(DEgenes_DESeq_only)))
print(paste0("Both ",length(DEgenes_both)))

vecDESeq2Results <- dfDESeq_Impulse$DESeq
names(vecDESeq2Results) <- dfDESeq_Impulse$Gene
strControlName <- NULL
strCaseName <- "case"
NPARAM <- 6

plotDEGenes(lsGeneIDs=DEgenes_both,
          arr2DCountData=arr2DCountData,
          vecNormConst=vecNormConst,
          dfAnnotationFull=dfAnnotationFull, 
          lsImpulseFits=lsImpulseFits,
          strCaseName=strCaseName, 
          strControlName=strControlName, 
          strFileNameSuffix="DE_ImpulseDE2andDESeq2", 
          strPlotTitleSuffix="", 
          strPlotSubtitle="",
          dfImpulseResults=dfImpulseResults,
          vecMethod2Results=vecDESeq2Results,
          strMode=strMode, 
          NPARAM=NPARAM)

# sort DESeq2 only genes by padj of DESeq2
dfDESeq_ImpulseDESeqOnly <- dfDESeq_Impulse[DEgenes_DESeq_only,]
DEgenes_DESeq_onlySorted <- dfDESeq_ImpulseDESeqOnly[order(dfDESeq_ImpulseDESeqOnly$DESeq),]$Gene
plotDEGenes(lsGeneIDs=DEgenes_DESeq_onlySorted,
          arr2DCountData=arr2DCountData,
          vecNormConst=vecNormConst,
          dfAnnotationFull=dfAnnotationFull, 
          lsImpulseFits=lsImpulseFits,
          strCaseName=strCaseName, 
          strControlName=strControlName, 
          strFileNameSuffix="DE_DESeq2_only", 
          strPlotTitleSuffix="", 
          strPlotSubtitle="",
          dfImpulseResults=dfImpulseResults,
          vecMethod2Results=vecDESeq2Results,
          strMode=strMode, 
          NPARAM=NPARAM)

plotDEGenes(lsGeneIDs=DEgenes_Impulse_only,
          arr2DCountData=arr2DCountData,
          vecNormConst=vecNormConst,
          dfAnnotationFull=dfAnnotationFull, 
          lsImpulseFits=lsImpulseFits,
          strCaseName=strCaseName, 
          strControlName=strControlName, 
          strFileNameSuffix="DE_ImpulseDE2_only", 
          strPlotTitleSuffix="", 
          strPlotSubtitle="",
          dfImpulseResults=dfImpulseResults,
          vecMethod2Results=vecDESeq2Results,
          strMode=strMode, 
          NPARAM=NPARAM)