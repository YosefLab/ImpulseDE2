# Load Pseudotime data p63
load("/Users/davidsebastianfischer/MasterThesis/data/Pseudotime/p63_ordered.Rdata")

# Investigate distribution of cells over pseudotime
lsPTpointsAll_p63PT1 <- as.vector(metadata$Pseudotime.1)
lsPTpoints_p63PT1 <- lsPTpointsAll_p63PT1[!is.na(lsPTpointsAll_p63PT1)]
print(paste0("p63PT1: Total cells: ",length(lsPTpointsAll_p63PT1),", Non NA: ",length(lsPTpoints_p63PT1)))

lsPTpointsAll_p63PT2 <- as.vector(metadata$Pseudotime.2)
lsPTpoints_p63PT2 <- lsPTpointsAll_p63PT2[!is.na(lsPTpointsAll_p63PT2)]
print(paste0("p63PT2: Total cells: ",length(lsPTpointsAll_p63PT2),", Non NA: ",length(lsPTpoints_p63PT2)))

# Load Pseudotime data EWT
load("/Users/davidsebastianfischer/MasterThesis/data/Pseudotime/EWT_ordered.Rdata")

# Investigate distribution of cells over pseudotime
lsPTpointsAll_EWTPT1 <- as.vector(metadata$Pseudotime.1)
lsPTpoints_EWTPT1 <- lsPTpointsAll_EWTPT1[!is.na(lsPTpointsAll_EWTPT1)]
print(paste0("EWTPT1: Total cells: ",length(lsPTpointsAll_EWTPT1),", Non NA: ",length(lsPTpoints_EWTPT1)))

lsPTpointsAll_EWTPT2 <- as.vector(metadata$Pseudotime.2)
lsPTpoints_EWTPT2 <- lsPTpointsAll_EWTPT2[!is.na(lsPTpointsAll_EWTPT2)]
print(paste0("EWTPT2: Total cells: ",length(lsPTpointsAll_EWTPT2),", Non NA: ",length(lsPTpoints_EWTPT2)))

lsPTpointsAll_EWTPT3 <- as.vector(metadata$Pseudotime.3)
lsPTpoints_EWTPT3 <- lsPTpointsAll_EWTPT3[!is.na(lsPTpointsAll_EWTPT3)]
print(paste0("EWTPT3: Total cells: ",length(lsPTpointsAll_EWTPT3),", Non NA: ",length(lsPTpoints_EWTPT3)))


graphics.off()
pdf("/Users/davidsebastianfischer/MasterThesis/data/Pseudotime/Histograms_CellDistributionOverPseudotime.pdf")

# Plot p63
hist(lsPTpoints_p63PT1,breaks=max(lsPTpoints_p63PT1),
  xlab="Pseudotime.1",ylab="Frequency",main="Olfactory epithelium cells: p63\nDistribution of cells over Pseudotime.1")
hist(lsPTpoints_p63PT2,breaks=max(lsPTpoints_p63PT2),
  xlab="Pseudotime.2",ylab="Frequency",main="Olfactory epithelium cells: p63\nDistribution of cells over Pseudotime.2")

# Plot EWT
hist(lsPTpoints_EWTPT1,breaks=max(lsPTpoints_EWTPT1),
  xlab="Pseudotime.1",ylab="Frequency",main="Olfactory epithelium cells: EWT\nDistribution of cells over Pseudotime.1")
hist(lsPTpoints_EWTPT2,breaks=max(lsPTpoints_EWTPT2),
  xlab="Pseudotime.2",ylab="Frequency",main="Olfactory epithelium cells: EWT\nDistribution of cells over Pseudotime.2")
hist(lsPTpoints_EWTPT3,breaks=max(lsPTpoints_EWTPT3),
  xlab="Pseudotime.3",ylab="Frequency",main="Olfactory epithelium cells: EWT\nDistribution of cells over Pseudotime.3")

dev.off()