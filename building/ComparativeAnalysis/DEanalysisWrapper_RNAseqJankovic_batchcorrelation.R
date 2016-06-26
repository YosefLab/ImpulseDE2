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

matDataAraw <- round(t(apply(dfRNA,1,as.numeric)))
vecGeneIDs <- as.vector(dfGeneIDs[,2])
vecboolidxDupIDs <- duplicated(vecGeneIDs)
if(sum(vecboolidxDupIDs)>0){
  vecGeneIDs[vecboolidxDupIDs] <- paste0(rep("DuplicatedGene_",sum(vecboolidxDupIDs)),seq(1,sum(vecboolidxDupIDs)))
}
rownames(matDataAraw) <- vecGeneIDs
colnames(matDataAraw) <- vecSamples

dfAnnotationRNA <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/AnnotationTable_RNAseqJankovic.tab",header=T)
dfAnnotationRNA$TimeCateg <- paste0(rep("_",length(dfAnnotationRNA$Time)),dfAnnotationRNA$Time)
dfAnnotationA <- dfAnnotationRNA

matDataA <- matDataAraw[,as.vector(dfAnnotationA$Sample)]
colnames(matDataA) <- dfAnnotationA$Sample
#only case
if(FALSE){
  condition="ctrl"
  matDataA <- matDataA[,dfAnnotationA$Condition==condition]
  colnames(matDataA) <- dfAnnotationA[dfAnnotationA$Condition==condition,]$Sample
  dfAnnotationA <- dfAnnotationA[dfAnnotationA$Condition==condition,]
  rownames(dfAnnotationA) <- dfAnnotationA$Sample
}

matDataARed <- matDataA[apply(matDataA,1,function(gene){!any(is.na(gene) | gene==0 )}),]
#colnames(matDataA) <- c(paste0("A",dfAnnotationA[1:6,]$TimeCateg,"h"),paste0("B",dfAnnotationA[7:12,]$TimeCateg,"h"))
#colnames(matDataA) <- c(paste0("A",dfAnnotationA[1:7,]$TimeCateg,"h"),paste0("B",dfAnnotationA[8:14,]$TimeCateg,"h"))
vecSRRid <- sapply(as.vector(dfAnnotationA$Sample),function(x){unlist(strsplit(x,split="SRR15255"))[2]})
colnames(matDataARed) <- c(paste0(vecSRRid[1:7],"_",dfAnnotationA[1:7,]$Condition,"_A",dfAnnotationA[1:7,]$TimeCateg,"h"),
  paste0(vecSRRid[8:13],"_",dfAnnotationA[8:13,]$Condition,"_A",dfAnnotationA[8:13,]$TimeCateg,"h"),
  paste0(vecSRRid[14:20],"_",dfAnnotationA[14:20,]$Condition,"_B",dfAnnotationA[14:20,]$TimeCateg,"h"),
  paste0(vecSRRid[21:26],"_",dfAnnotationA[21:26,]$Condition,"_B",dfAnnotationA[21:26,]$TimeCateg,"h"))
library(ggplot2)
library(reshape2)
q <- qplot(x=Var1, y=Var2, 
  data=melt(cor(matDataARed)), 
  fill=value, 
  geom="tile",
  xlab="",ylab="")
q + theme(axis.text.x = element_text(angle = 90, hjust = 1))

plot(log(matDataARed[,8])/log(2),log(matDataARed[,21])/log(2))

pdf("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/batcheffect_analysis/correlation_heatmap_batcheffect.pdf")
q <- qplot(x=Var1, y=Var2, 
  data=melt(cor(matDataARed)), 
  fill=value, 
  geom="tile",
  xlab="",ylab="")
q + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
graphics.off()

# batch correction

# look at example data set
library(sva)
library(pamr)
library(limma)
library(bladderbatch)
data(bladderdata)
pheno = pData(bladderEset)
edata = exprs(bladderEset)
pheno
head(edata)
batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

# on RNA data
# filter out all zero and NA observations - ImpulseDE cannot handle NA anyway
library(sva)
library(pamr)
library(limma)
matDataA_noNA <- matDataA[apply(matDataA,1,function(gene){!any(is.na(gene)) & any(gene>0)}),]
condition <- dfAnnotationA$Condition
time <- dfAnnotationA$TimeCateg
batch <- array(NA,dim(dfAnnotationA)[1])
batch[dfAnnotationA$LongitudinalSeries=="A"] <- 1
batch[dfAnnotationA$LongitudinalSeries=="B"] <- 2
mm <- model.matrix(~time + condition + batch)
pheno <- data.frame(
  sample=seq(1,dim(dfAnnotationA)[1]),
  batch=batch,
  time=dfAnnotationA$TimeCateg,
  cond=dfAnnotationA$Condition )
rownames(pheno) <- dfAnnotationA$Sample
modcombat <- model.matrix(~1, data=pheno)
matCountsNobatch <- ComBat(dat=matDataA_noNA, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

matCountsNobatchRed <- matCountsNobatch[apply(matCountsNobatch,1,function(gene){!any(is.na(gene) | gene==0 )}),]
vecSRRid <- sapply(as.vector(dfAnnotationA$Sample),function(x){unlist(strsplit(x,split="SRR15255"))[2]})
colnames(matCountsNobatchRed) <- c(paste0(vecSRRid[1:7],"_",dfAnnotationA[1:7,]$Condition,"_A",dfAnnotationA[1:7,]$TimeCateg,"h"),
  paste0(vecSRRid[8:13],"_",dfAnnotationA[8:13,]$Condition,"_A",dfAnnotationA[8:13,]$TimeCateg,"h"),
  paste0(vecSRRid[14:20],"_",dfAnnotationA[14:20,]$Condition,"_B",dfAnnotationA[14:20,]$TimeCateg,"h"),
  paste0(vecSRRid[21:26],"_",dfAnnotationA[21:26,]$Condition,"_B",dfAnnotationA[21:26,]$TimeCateg,"h"))
library(ggplot2)
library(reshape2)
q <- qplot(x=Var1, y=Var2, 
  data=melt(cor(matCountsNobatchRed)), 
  fill=value, 
  geom="tile",
  xlab="",ylab="")
q + theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/RNAseqJankovic/batcheffect_analysis/correlation_heatmap_CombatCorrected.pdf")
q <- qplot(x=Var1, y=Var2, 
  data=melt(cor(matDataARed)), 
  fill=value, 
  geom="tile",
  xlab="",ylab="")
q + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
graphics.off()