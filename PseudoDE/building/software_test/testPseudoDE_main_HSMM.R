rm(list = ls())
NCORES = 3

# Load data
dfCountsHSMM <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/rsem/rsem_readCountsTable.txt", header=FALSE, colClasses="numeric")
dfFpkmHSMM <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/rsem/rsem_fpkmTable.txt", header=FALSE, colClasses="numeric")
dfCellsHSMM <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/rsem/cell_list.txt", header=FALSE, colClasses="character")
dfGenesHSMM <- read.table("/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/rsem/gene_list.txt", header=FALSE, colClasses="character")

vecSampleNames <- apply(dfCellsHSMM, 1, function(sample){ unlist(strsplit(sample,split="/"))[2] })
vecSampleNames <- unlist(lapply(vecSampleNames, function(sample){ unlist(strsplit(sample,split="_1"))[1] }))
colnames(dfCountsHSMM) <- vecSampleNames
rownames(dfCountsHSMM) <- dfGenesHSMM[,1]
colnames(dfFpkmHSMM) <- vecSampleNames
rownames(dfFpkmHSMM) <- dfGenesHSMM[,1]

# get single cell samples
fileSampleAnnot <- "/Users/davidsebastianfischer/MasterThesis/data/ImpulseDE2_datasets/HSMM/SraRunTable_scRNAMonocle.txt"
dfSampleAnnot <- read.table(fileSampleAnnot, header=TRUE, sep="\t", colClasses="character")
vecSCSamples <- dfSampleAnnot[dfSampleAnnot$library_protocol_s == "Single-cell",]$Run_s
dfCountsHSMM_SC <- dfCountsHSMM[,vecSCSamples]
dfFpkmHSMM_SC <- dfFpkmHSMM[,vecSCSamples]

# Make monolcle annotation
HSMM_gene_annotationRAW <- data.frame( gene_short_name=dfGenesHSMM[,2], biotype="N/A" )
rownames(HSMM_gene_annotationRAW) <- dfGenesHSMM[,1]
HSMM_sample_sheetRAW <- data.frame( Library=dfSampleAnnot[dfSampleAnnot$Run_s%in% vecSCSamples,]$SRA_Study_s, 
  Media="GM",
  Hours=0,
  stringsAsFactors=FALSE)
rownames(HSMM_sample_sheetRAW) <- dfSampleAnnot[dfSampleAnnot$Run_s%in% vecSCSamples,]$Run_s
# Differentiation medium DM from growth medium GM after t=0h
HSMM_sample_sheetRAW[sapply(dfSampleAnnot[dfSampleAnnot$Run_s %in% rownames(HSMM_sample_sheetRAW),]$source_name_s, function(x){unlist(strsplit(x, split="_"))[2]}) == "Cell T0",]$Media <- "DM"
HSMM_sample_sheetRAW[sapply(dfSampleAnnot[dfSampleAnnot$Run_s %in% rownames(HSMM_sample_sheetRAW),]$source_name_s, function(x){unlist(strsplit(x, split="_"))[2]}) == "Cell T24",]$Hours <- 24
HSMM_sample_sheetRAW[sapply(dfSampleAnnot[dfSampleAnnot$Run_s %in% rownames(HSMM_sample_sheetRAW),]$source_name_s, function(x){unlist(strsplit(x, split="_"))[2]}) == "Cell T48",]$Hours <- 48
HSMM_sample_sheetRAW[sapply(dfSampleAnnot[dfSampleAnnot$Run_s %in% rownames(HSMM_sample_sheetRAW),]$source_name_s, function(x){unlist(strsplit(x, split="_"))[2]}) == "Cell T72",]$Hours <- 72

# Get pseudotime
#biocLite("monocle")
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("Biobase"))
#devtools::install_github("cole-trapnell-lab/monocle-release@monocle2")
#source("https://bioconductor.org/biocLite.R")

detach("package:monocle", unload=TRUE)
library(monocle)
#library(HSMMSingleCell)
#library(data.table)
library(Biobase)

#data(HSMM_expr_matrix)
#data(HSMM_gene_annotation)
#data(HSMM_sample_sheet)
pd <- new("AnnotatedDataFrame", data=HSMM_sample_sheetRAW)
fd <- new("AnnotatedDataFrame", data=HSMM_gene_annotationRAW)
HSMM <- newCellDataSet( cellData=as.matrix(dfFpkmHSMM_SC), phenoData=pd, featureData=fd )
print(head(fData(HSMM)))

# filer
HSMM <- detectGenes(HSMM, min_expr=0.1)
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 50))
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,], fullModelFormulaStr="expression~Media",cores=NCORES)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
ordering_genes <- intersect(ordering_genes, expressed_genes)
HSMM <- setOrderingFilter(HSMM, ordering_genes)
HSMM <- reduceDimension(HSMM, use_irlba=F)
# num_paths is number of leaves of tree
HSMM <- orderCells(HSMM, num_paths=2, reverse=F)
plot_spanning_tree(HSMM)
print(pData(HSMM))

# Investigate distribution of cells over pseudotime
vecPTpointsAll_HSMM <- as.vector(pData(HSMM)$Pseudotime)
names(vecPTpointsAll_HSMM) <- as.vector(rownames(pData(HSMM)))
vecPTpoints_HSMM <- vecPTpointsAll_HSMM[!is.na(vecPTpointsAll_HSMM)]
print(paste0("HSMM: Total cells: ",length(vecPTpointsAll_HSMM),", Non NA: ",length(vecPTpoints_HSMM)))

# Run PseudoDE
matCounts <- data.matrix(dfCountsHSMM_SC)
vecPT <- vecPTpoints_HSMM
matCounts <- matCounts[,colnames(matCounts) %in% names(vecPTpoints_HSMM)]
vecPT <- vecPT[colnames(matCounts)]

dfPseudotime <- data.frame( pseudotime=as.vector(vecPT) )
plotEDF <- ggplot() +
  geom_density(data=dfPseudotime, aes(x=pseudotime), colour="black", bw=1) +
  labs(title="Density estimation of cells in pseudotime") +
  xlab("pseudotime") +
  ylab("empirical probability density")
print(plotEDF)

matCountsRed <- matCounts
matCountsRed <- round(matCountsRed)

nProc=3
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/PseudoDE_main.R")
setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
lsDEresults <- runPseudoDE(matCounts=matCountsRed,
  vecPseudotime=vecPT,
  boolPlotZINBfits=FALSE,
  nProc=nProc)

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
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_lsZINBparam.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_matCountsProc.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_lsResultsClustering.RData")
  load("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out/PseudoDE_matProbNB.RData")
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
