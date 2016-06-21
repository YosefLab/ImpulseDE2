#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++     Compare DE methods    ++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compare ImpulseDE2 output against other differential expression method
#' 
#' The comparison is performed based on the adjusted p-values of differential
#' expression. This methods allows the user to explore visually how 
#' ImpulseDE2 output differs from similar methods.
#' 
#' @seealso Called separately by user.
#' 
#' @param matCountDataProc: (matrix genes x samples)
#'    Count data: Reduced version of \code{matCountData}. 
#'    For internal use.
#' @param matProbNB: (probability matrix genes x samples) 
#'    Probability of observations to come from negative binomial 
#'    component of mixture model.
#' @param strMode: (str) [Default "batch"] 
#'    {"batch","longitudinal","singlecell"}
#'    Mode of model fitting.
#' 
#' @return matSizeFactors: (numeric matrix genes x samples) 
#'    Model scaling factors for each observation which take
#'    sequencing depth into account (size factors). One size
#'    factor per sample - rows of this matrix are equal.
#' @export

library(gplots)

compareDEMethods <- function(matQval,
  strMethod1="ImpulseDE2",
  strMethod2="",
  Q = 10^(-3),
  Qdelta = 10^(2),
  matCountDataProc,
  lsMatTranslationFactors,
  matSizeFactors,
  dfAnnotationProc, 
  lsImpulseFits,
  dfImpulseResults,
  strCaseName="case", 
  strControlName=NULL, 
  strMode="batch",
  strDataDescriptionFilename=""){
  
  graphics.off()
  # 1. Heatmap: Q-value thresholds
  mat_overlap <- array(NA,c(11,11))
  for(i in 0:10){
    for(j in 0:10){
      sig_Impulse <- as.numeric(matQval[,strMethod1]) <= 10^(-j) & !is.na(as.numeric(matQval[,strMethod1]))
      sig_Ref <- as.numeric(matQval[,strMethod2]) <= 10^(-i) & !is.na(as.numeric(matQval[,strMethod2]))
      mat_overlap[i+1,j+1] <- sum(sig_Ref & sig_Impulse)/sum(sig_Ref | sig_Impulse)
    }
  }
  rownames(mat_overlap) <- 0:-10
  colnames(mat_overlap) <- 0:-10
  
  heatmap(mat_overlap, keep.dendro = FALSE,Rowv=NA,Colv= "Rowv",symm=FALSE,
    xlab =  paste0("ImpulseDE2 scores"), ylab = "edge scores")
  breaks <- seq(0,1,by=0.01)
  hm.colors <- colorpanel( length(breaks)-1, "yellow", "red" )
  graphics.off()
  pdf(paste0(strMethod1,"-",strMethod2,"_JaccardCoeff_Heatmap.pdf"),width=7,height=7)
  heatmap.2(mat_overlap, 
    dendrogram="none", 
    Rowv=FALSE,Colv=FALSE, 
    xlab =  paste0("log_10(p-value) ", strMethod1), 
    ylab =  paste0("log_10(p-value) ", strMethod2),
    breaks=breaks,
    col=hm.colors, 
    scale="none",
    trace="none",
    density.info="none",
    key.title = " ", 
    key.xlab = paste0("Jaccard coefficient"), 
    key.ylab = NULL,
    symkey=FALSE,
    cellnote=round(mat_overlap,digits=2),
    notecol="grey",
    lmat=rbind( c(3,4),c(2,1) ),lhei=c(1,4), lwid=c(1,4), margins=c(5,5))
  dev.off()
  graphics.off()
  
  # 2. CDF p-values
  vecX <- seq(1,25,by=0.1)
  vecCDF1 <- sapply(vecX, function(thres){sum(log(as.numeric(matQval[,strMethod1]))/log(10) >= -thres, na.rm=TRUE)})
  vecCDF2 <- sapply(vecX, function(thres){sum(log(as.numeric(matQval[,strMethod2]))/log(10) >= -thres, na.rm=TRUE)})
  pdf(paste0(strMethod1,"-",strMethod2,"_ECDF-pvalues.pdf"),width=7,height=7)
  plot(vecX,vecCDF1,
    col="black",pch=4,type="l",
    ylim=c(0,max(max(vecCDF1,na.rm=TRUE),max(vecCDF2,na.rm=TRUE))),
    xlab="-log_10(threshold p-value)", 
    ylab=paste0("Number of accepted null\nhypothesis at threshold"),
    main=paste0("Cumulative number of p-values larger than threshold\n",strMethod1," versus ",strMethod2))
  points(vecX,vecCDF2,
    col="red",pch=4,type="l")
  legend(x="bottomright",
    legend=c(strMethod1,strMethod2),
    fill=c("black","red"))
  dev.off()
  graphics.off()
  
  # 3. Impulse fits of differentially called DE genes
  # Set detection and differential calling threshold
  
  # Note cannot plot genes which are NA in ImpulseDE2
  vecboolImpulseFitted <- !is.na(as.numeric(matQval[,strMethod2]))
  vecboolRefBeatsQ <- as.numeric(matQval[,strMethod2]) < Q & !is.na(as.numeric(matQval[,strMethod2]))
  vecBoolRefBeatsImpulse <- as.numeric(matQval[,strMethod1]) >= Qdelta*as.numeric(matQval[,strMethod2]) | is.na(as.numeric(matQval[,strMethod2]))
  vecboolImpulseBeatsQ <- as.numeric(matQval[,strMethod1]) < Q & !is.na(as.numeric(matQval[,strMethod1]))
  vecboolImpulseBeatsRef <- as.numeric(matQval[,strMethod2]) >= Qdelta*as.numeric(matQval[,strMethod1]) | is.na(as.numeric(matQval[,strMethod1]))
  
  vecDEgenes_Ref_only <- as.vector(matQval[vecboolRefBeatsQ & vecBoolRefBeatsImpulse & vecboolImpulseFitted,"Gene"])
  print(head(vecDEgenes_Ref_only))
  print(dim(lsImpulseFits$parameters_case))
  print(rownames(lsImpulseFits$parameters_case))
  vecDEgenes_Impulse_only <- as.vector(matQval[vecboolImpulseBeatsQ & vecboolImpulseBeatsRef & vecboolImpulseFitted,"Gene"])
  graphics.off()
  print(paste0("Number of significant DE genes at q-value of ",Q))
  print(paste0(strMethod1," only ",length(vecDEgenes_Impulse_only)))
  print(paste0(strMethod2," only ",length(vecDEgenes_Ref_only)))
  
  vecRefResults <- as.numeric(matQval[,strMethod2])
  names(vecRefResults) <- matQval[,"Gene"]
  strControlName <- NULL
  strCaseName <- "case"
  strMode="batch"
  NPARAM <- 6
  
  plotDEGenes(vecGeneIDs=vecDEgenes_Ref_only,
    matCountDataProc=matCountDataProc,
    lsMatTranslationFactors=lsMatTranslationFactors,
    matSizeFactors=matSizeFactors,
    dfAnnotationProc=dfAnnotationProc, 
    lsImpulseFits=lsImpulseFits,
    strCaseName=strCaseName, 
    strControlName=strControlName, 
    strFileNameSuffix=paste0(strDataDescriptionFilename,"_",strMethod2,"Beats",strMethod1), 
    strPlotTitleSuffix="", 
    strPlotSubtitle="",
    dfImpulseResults=dfImpulseResults,
    vecRefPval=vecRefResults,
    strNameMethod2=strMethod2,
    strMode=strMode, 
    NPARAM=6)
  
  if(length(vecDEgenes_Impulse_only)){
    plotDEGenes(vecGeneIDs=vecDEgenes_Impulse_only,
      matCountDataProc=matCountDataProc,
      lsMatTranslationFactors=lsMatTranslationFactors,
      matSizeFactors=matSizeFactors,
      dfAnnotationProc=dfAnnotationProc, 
      lsImpulseFits=lsImpulseFits,
      strCaseName=strCaseName, 
      strControlName=strControlName, 
      strFileNameSuffix=paste0(strDataDescriptionFilename,"_",strMethod1,"Beats",strMethod2), 
      strPlotTitleSuffix="", 
      strPlotSubtitle="",
      dfImpulseResults=dfImpulseResults,
      vecRefPval=vecRefResults,
      strNameMethod2=strMethod2,
      strMode=strMode, 
      NPARAM=NPARAM)
  }
}