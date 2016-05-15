rm(list = ls())

# Load Pseudotime data p63
load("/Users/davidsebastianfischer/MasterThesis/data/Pseudotime/p63_ordered.Rdata")
dfCountsE4 <- read.table("/Users/davidsebastianfischer/MasterThesis/data/Pseudotime/E4_Tophatfilt_counts.txt", header=TRUE, row.names=1)
# Make sure got same cells
sum(colnames(dfCountsE4) %in% metadata$sample)
length(colnames(dfCountsE4))
length(metadata$sample)

# Investigate distribution of cells over pseudotime
lsPTpointsAll_p63PT1 <- as.vector(metadata$Pseudotime.1)
names(lsPTpointsAll_p63PT1) <- as.vector(metadata$sample)
lsPTpoints_p63PT1 <- lsPTpointsAll_p63PT1[!is.na(lsPTpointsAll_p63PT1)]
print(paste0("p63PT1: Total cells: ",length(lsPTpointsAll_p63PT1),", Non NA: ",length(lsPTpoints_p63PT1)))

source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/PseudoDE_main.R")
matCounts <- data.matrix(dfCountsE4)
#matCounts <- matCounts[1:20,]
#matCounts <- round(counts)
nProc=2
vecPseudotime <- lsPTpoints_p63PT1
setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
lsDEresults <- runPseudoDE(matCounts=matCounts,vecPseudotime=vecPseudotime)
print(lsDEresults)

if(FALSE){
      vecboolNonzeroGenes <- apply(matCounts,1,
        function(gene){mean(gene)>10})
      matCounts <- matCounts[vecboolNonzeroGenes,]
      MAXITER <- 3
      
      # Fit zinb model
      lsZINBparam <- estimate_zinb(
        Y = matCounts, 
        maxiter = MAXITER, 
        verbose = TRUE)
}