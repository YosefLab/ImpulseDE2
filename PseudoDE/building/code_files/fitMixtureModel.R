fitMixtureModel <- function(matCounts, lsResultsClustering, 
  dfAnnotationClusters){
  
  --todo: remove 0s
  # 1.Identify overdispersion factor of negative binomial
  # by gene:
  # From GLM based on categorial independent variable
  # (cluster index). Dispersion estimation is performed by 
  # DESeq2, which extends GLM fitting of overdispersion by
  # an Empririal Bayes correction.
  
  lsDESeq2ResultsClusters <- runDESeq2(
    dfAnnotationFull=dfAnnotationClusters, 
    arr2DCountData=matCounts, 
    strMode="singlecell")
  vecDESeq2Dispersions <- lsDESeq2ResultsClusters[[1]]
  dfDESeq2ResultsClusters <- lsDESeq2ResultsClusters[[2]]
  
  # 2. MLE of mean of negative binomials of each cluster
  # of each gene:
  # The MLE of the mean of the negative binomial has a closed
  # form solution: the sample average.
  
  matMeans <- array(NA,c(dim(matCounts)[1],lsResultsClustering$K))
  rownames(matMeans) <- rownames(matCounts)
  colnames(matMeans) <- c(1:lsResultsClustering$K)
  for(gene in rownames(matCounts)){
    matMeans[gene,] <- sapply(c(1:lsResultsClustering$K),
      function(cluster){ matMeans[gene,cluster] <- 
          mean( matCounts[gene,lsResultsClustering$Assignments==cluster], na.rm=TRUE )})
  }
  
  # 3. Identify drop-out rate: The drop-out rate is 
  # obtained by fitting the hurdle model (zero inflated
  # negative binomial) to the data using the previously
  # identified negative binomial component of the mixture.
  
  return(list(vecDESeq2DispersionsClusters,vecDropoutRates))
}