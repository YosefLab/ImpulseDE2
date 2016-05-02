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

loglikH1 <- evalLogLikImpulseByTC(vecTheta=vecTheta, vecX=vecX,
  matY=matY, scaDispEst=scaDispEst,
  scaMuEst=scaMuEst, vecMuEstTimecourse=vecMuEstTimecourse)

loglikH0 <- evalLogLikMean(scaMuEst=vecMuEstTimecourse[1],
  matY=matY[,1], scaDispEst=scaDispEst)
loglikH0 <- loglikH0 + evalLogLikMean(scaMuEst=vecMuEstTimecourse[2],
  matY=matY[,2], scaDispEst=scaDispEst)
loglikH0 <- loglikH0 + evalLogLikMean(scaMuEst=vecMuEstTimecourse[3],
  matY=matY[,3], scaDispEst=scaDispEst)

print("By timecourse")
loglikH0
loglikH1

loglikH1 <- evalLogLikImpulseBatch(vecTheta=vecTheta, vecX=vecX,
  matY=matY, scaDispEst=scaDispEst)

loglikH0 <- evalLogLikMean(scaMuEst=scaMuEst,
  matY=matY, scaDispEst=scaDispEst)

print("As batch")
loglikH0
loglikH1

dfImpulseResults
lsImpulseFits$parameters_case
