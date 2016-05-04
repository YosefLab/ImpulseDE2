rm(list = ls())

# Load Pseudotime data p63
load("/Users/davidsebastianfischer/MasterThesis/data/Pseudotime/p63_ordered.Rdata")

# Investigate distribution of cells over pseudotime
lsPTpointsAll_p63PT1 <- as.vector(metadata$Pseudotime.1)
names(lsPTpointsAll_p63PT1) <- as.vector(metadata$sample)
lsPTpoints_p63PT1 <- lsPTpointsAll_p63PT1[!is.na(lsPTpointsAll_p63PT1)]
print(paste0("p63PT1: Total cells: ",length(lsPTpointsAll_p63PT1),", Non NA: ",length(lsPTpoints_p63PT1)))

source("/Users/davidsebastianfischer/MasterThesis/code/PseudoDE/building/code_files/PseudoDE_main.R")
lsDEresults <- runPseudoDE(matCounts=round(counts),vecPseudotime=lsPTpoints_p63PT1)
print(lsDEresults)