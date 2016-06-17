rm(list=ls())
if(FALSE){
  # Load files from interior of ImpulseDE
  setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
  load("ImpulseDE2_arr3DCountData.RData")
  load("ImpulseDE2_dfAnnotationRed.RData")
  # Load Impulse output
  load("ImpulseDE2_dfImpulseResults.RData")
  load("ImpulseDE2_lsDEGenes.RData")
  load("ImpulseDE2_lsImpulseFits.RData")
  load("ImpulseDE2_dfDESeq2Results.RData")
  NPARAM <- 6
  setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files")
  source("srcImpulseDE2_CostFunctionsFit.R")
  setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
  geneID <- "G0S2"
  vecX<-c(0,2,4,6,12,24)
  vecTheta <- lsImpulseFits$parameters_case[geneID,1:6]
  vecTheta[2:4] <- log(vecTheta[2:4])
  matY <- arr3DCountData[geneID,,]
  scaDispEst <- as.numeric(as.vector(dfImpulseResults[geneID,]$size))
  scaMuEst <- log(lsImpulseFits$parameters_case[geneID,"mu"])
  vecindMuByTimecourse <- grep("muByTimecourse",colnames(lsImpulseFits$parameters_case))
  vecMuEstTimecourse <- log(lsImpulseFits$parameters_case[geneID,vecindMuByTimecourse])
}

setwd("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files")
source("srcImpulseDE2_CostFunctionsFit.R")
source("srcImpulseDE2_computeNormConst.R")
source("ImpulseDE2_main.R")
# artificial data
# generate size factors
N=1000
arr3DToyData <- array(NA,c(N,5,3))
# Fill toy data with gaussian rand data with different means
matMeans <- t(array(c(2,5,8),c(3,5)))
matMeans <- t(sapply(seq(10,50,by=10),function(sample){matMeans[sample/10,]+sample}))
for(i in 1:N){
  for(j in 1:5){
    arr3DToyData[i,j,] <- round(rnorm(n=3,mean=matMeans[j,], sd=0.5))
  }
}
matNormConst <- computeNormConst(arr3DCountData=arr3DToyData)
vecTheta <- c(1,10,40,50,2,4)
vecTheta[2:4] <- log(vecTheta[2:4])
vecX<-c(1:5)
matY <- arr3DToyData[1,,]
scaDispEst <- 2
scaMu <- mean(arr3DToyData[1,,]/matNormConst)
vecMuTimecourses <- apply(arr3DToyData[1,,]/matNormConst,2,mean)
matMuTimecourses <- matrix(vecMuTimecourses, nrow=dim(matY)[1], ncol=dim(matY)[2], byrow=TRUE)
# Compute translation factors: Normalisation factor
# to scale impulse model to indivindual time courses.
# Note: Translation factor is the same for all replicates
# in a time course.
matTranslationFactors <- matMuTimecourses/scaMu
matboolObserved <- !is.na(arr3DToyData[1,,])
loglikH1 <- evalLogLikImpulseByTC(
  vecTheta=vecTheta, 
  vecX=vecX,
  matY=matY, 
  scaDispEst=scaDispEst, 
  matNormConst=matNormConst,
  matTranslationFactors=matTranslationFactors,
  matboolObserved=matboolObserved)
print("Reference:")
matImpulseValue = matrix(calcImpulse_comp(vecTheta,vecX),
  nrow=dim(matY)[1], ncol=dim(matY)[2], byrow=FALSE)*matNormConst
vecTranslationFactors <- vecMuTimecourses/scaMu
sum(dnbinom(
  matY, 
  mu=matImpulseValue * matrix(vecTranslationFactors,nrow=dim(matImpulseValue)[1],ncol=dim(matImpulseValue)[2],byrow=TRUE), 
  size=scaDispEst, 
  log=TRUE))

print("By timecourse")
loglikH1

loglikH1 <- evalLogLikImpulseBatch(
  vecTheta=vecTheta, 
  vecX=vecX,
  matY=matY, 
  scaDispEst=scaDispEst, 
  matNormConst=matNormConst,
  matboolObserved=matboolObserved )

print("As batch")
loglikH1
print("Reference:")
matImpulseValue = matrix(calcImpulse_comp(vecTheta,vecX),
  nrow=dim(matY)[1], ncol=dim(matY)[2], byrow=FALSE)*matNormConst
sum(dnbinom(
  matY,
  mu=matImpulseValue, 
  size=scaDispEst, 
  log=TRUE))

## single cell : drop outs
# artificial data
# generate size factors
N=1000
arr3DToyData <- array(0,c(N,5,3))
# Fill toy data with gaussian rand data with different means
# drop out rate 0.8
matMeans <- t(array(c(2,5,8),c(3,5)))
matMeans <- t(sapply(seq(10,50,by=10),function(sample){matMeans[sample/10,]+sample}))
for(i in 1:N){
  for(j in 1:5){
    idxNotDropped <- rbinom(n=3,size=1,prob=0.8)==1
    arr3DToyData[i,j,idxNotDropped] <- round(rnorm(n=3,mean=matMeans[j,], sd=1))[idxNotDropped]
  }
}
# matnormconst from earlier
vecTheta <- c(1,10,40,50,2,4)
vecTheta[2:4] <- log(vecTheta[2:4])
vecX<-c(1:5)
matY <- arr3DToyData[1,,]
scaDispEst <- 2
scaDropoutEst <- 0.78
matboolZero <- matY==0 

loglikH1 <- evalLogLikImpulseSC(
  vecTheta=vecTheta,
  vecX=vecX,
  matY=matY,
  scaDispEst=scaDispEst,
  scaDropoutEst=scaDropoutEst,
  matNormConst=matNormConst,
  matboolObserved=matboolObserved,
  matboolZero=matboolZero)

print("single cell")
loglikH1
print("Reference:")
matboolObserved <- !is.na(matY)
matboolZero <- matY==0
matImpulseValue = matrix(calcImpulse_comp(vecTheta,vecX),
  nrow=dim(matY)[1], ncol=dim(matY)[2], byrow=FALSE)*matNormConst
scaLogLikZeros <- sum(log(
  (1-scaDropoutEst)*
    dnbinom(
      matY[matboolZero & matboolObserved], 
      mu=matImpulseValue[matboolZero & matboolObserved], 
      size=scaDispEst, 
      log=FALSE) +
    scaDropoutEst
))
# Likelihood of non-zero counts:
scaLogLikNonzeros <- sum(log(
  (1-scaDropoutEst)*
    dnbinom(
      matY[!matboolZero & matboolObserved], 
      mu=matImpulseValue[!matboolZero & matboolObserved], 
      size=scaDispEst, 
      log=FALSE)
))
# Compute likelihood of all data:
scaLogLik <- scaLogLikZeros + scaLogLikNonzeros
print(c(scaLogLikZeros,scaLogLikNonzeros))
print(scaLogLik)