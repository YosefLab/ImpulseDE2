setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files")
# Compute normalisation constants
source("srcImpulseDE2_computeNormConst.R")

# generate data
N=1000
arr3DToyData <- array(NA,c(N,5,3))
# Fill toy data with gaussian rand data with different means
matMeans <- t(array(c(4,5,6),c(3,5)))
matMeans <- t(sapply(seq(10,50,by=10),function(sample){matMeans[sample/10,]+sample}))
for(i in 1:N){
  for(j in 1:5){
    arr3DToyData[i,j,] <- rnorm(n=3,mean=matMeans[j,], sd=1)
  }
}

matMeans
apply(arr3DToyData,c(2,3),mean)

# generate normalisation constants
matMeans/mean(matMeans)
computeNormConst(arr3DCountData=arr3DToyData)
