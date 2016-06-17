### Run after completion of ImpulseDE
### Plot deviance distribution
# non-DE scores
graphics.off()
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
pdf("EmpiricalDevianceDistribution.pdf",height=6.0,width=9.0)
dImp<-density(as.numeric(as.vector(DE_results[(as.numeric(DE_results$adj.p) > Q_value),"deviance"])))
dDE<-density(as.numeric((dds_resultsTable[(as.numeric(dds_resultsTable$padj) > Q_value),"stat"])))

plot(dImp$x,dImp$y,type='l',col='black',xlim=c(0,20),ylim=c(0,0.2),
  xlab="Deviance", ylab="Probability density",
  main="DESeq2 vs ImpulseDEv2 test statistic\nProbability density of non-DE gene scores")
lines(dDE$x,dDE$y,type='l',col='red')
curve( dchisq(x, df=5), col='green',add=TRUE)
legend(x="topright",c("ImpulseDEv2 statistic","DESeq2 statistic","Chi-square (df=5)"),fill=c("black","red","green"), cex=1)

# DE scores
dImp<-density(as.numeric(as.vector(DE_results[(as.numeric(DE_results$adj.p) < Q_value),"deviance"])))
dDE<-density(as.numeric((dds_resultsTable[(as.numeric(dds_resultsTable$padj) < Q_value),"stat"])))

plot(dImp$x,dImp$y,type='l',col='black',xlim=c(0,100),ylim=c(0,0.15),
  xlab="Deviance", ylab="Probability density",
  main="DESeq2 vs ImpulseDEv2 test statistic\nProbability density of DE gene scores")
lines(dDE$x,dDE$y,type='l',col='red')
curve( dchisq(x, df=5), col='green',add=TRUE)
legend(x="topright",c("ImpulseDEv2 statistic","DESeq2 statistic","Chi-square (df=5)"),fill=c("black","red","green"), cex=1)

# All scores
dImp<-density(as.numeric(as.vector(DE_results[,"deviance"])))
dDE<-density(as.numeric((dds_resultsTable[,"stat"])))

plot(dImp$x,dImp$y,type='l',col='black',xlim=c(0,100),ylim=c(0,0.15),
  xlab="Deviance", ylab="Probability density",
  main="DESeq2 vs ImpulseDEv2 test statistic\nProbability density of all gene scores")
lines(dDE$x,dDE$y,type='l',col='red')
curve( dchisq(x, df=5), col='green',add=TRUE)
legend(x="topright",c("ImpulseDEv2 statistic","DESeq2 statistic","Chi-square (df=5)"),fill=c("black","red","green"), cex=1)
dev.off()

# only plot once if simulating under null: no DE genes detected
# Simulated one single model
if(FALSE){
  ### Plot deviance distribution
  # non-DE scores
  graphics.off()
  setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
  pdf("SimulatedNullDevianceDistribution.pdf",height=6.0,width=9.0)
  dImp<-density(as.numeric(as.vector(DE_results[,"deviance"])),from=0,to=200)
  dDE<-density(as.numeric((dds_resultsTable[,"stat"])),from=0,to=200)
  
  plot(dImp$x,dImp$y,type='l',col='black',xlim=c(0,20),ylim=c(0,0.2),
    xlab="Deviance", ylab="Probability density",
    main="DESeq2 vs ImpulseDEv2 test statistic\nProbability density of all gene scores")
  lines(dDE$x,dDE$y,type='l',col='red')
  curve( dchisq(x, df=5), col='green',add=TRUE)
  legend(x="topright",c("ImpulseDEv2 statistic","DESeq2 statistic","Chi-square (df=5)"),fill=c("black","red","green"), cex=1)
  dev.off()
}