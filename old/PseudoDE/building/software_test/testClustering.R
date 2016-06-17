# Load Pseudotime data p63
load("/Users/davidsebastianfischer/MasterThesis/data/Pseudotime/p63_ordered.Rdata")

# Investigate distribution of cells over pseudotime
lsPTpointsAll_p63PT1 <- as.vector(metadata$Pseudotime.1)
lsPTpoints_p63PT1 <- lsPTpointsAll_p63PT1[!is.na(lsPTpointsAll_p63PT1)]
print(paste0("p63PT1: Total cells: ",length(lsPTpointsAll_p63PT1),", Non NA: ",length(lsPTpoints_p63PT1)))

lsSS_p63PT1 <- clusterCellsInPseudotime(lsPTpoints_p63PT1)
print(lsSS_p63PT1$Centroids)
hist(lsPTpoints_p63PT1,breaks=max(lsPTpoints_p63PT1),
  xlab="Pseudotime.1",ylab="Frequency",main="Olfactory epithelium cells: p63\nDistribution of cells over Pseudotime.1")
