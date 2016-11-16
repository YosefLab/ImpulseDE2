#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++   Simulate a data set  +++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Simulate a data set for ImpulseDE2

#' Simulates a data set with genes with constant and impulse
#' expression traces. Expression strength and variation in impulse
#' like traces are parameterised and random. All temporary files
#' are saved into dirOutSimulation and only the objects necessary
#' for running ImpulseDE2 (the count matrix and the annotation table
#' are returned). The remaining objects representing hidden
#' parameters can be used to evaluate parameter estimates.
#' 
#' @seealso Called by separately by user.
#' 
#' @param vecTimePointsA: (numeric vector number of time points)
#'    Number of time points in batch A.
#' @param vecTimePointsB: (numeric vector number of time points)
#'    Number of time points in batch B.
#' @param boolCaseCtrl: (bool) [Defaul FALSE] Whether A and B are 
#'    are case and control data (used for annotation data frame).
#' @param scaNConst: (scalar) Number of constant genes in data set.
#' @param scaNImp: (scalar) Number of impulse distributed genes in data set.
#' @param scaNLin: (scalar) Number of linear distributed genes in data set.
#' @param scaNSig: (scalar) Number of sigmoid distributed genes in data set.
#' @param scaMumax: (scalar) [Default 1000]
#'    Maximum expression mean parameter to be used.
#' @param scaSDImpulseAmplitude: (scalar) [Default 1]
#'    Standard deviation of normal distribution form which the 
#'    amplitude change within an impulse trace is drawn.
#' @param scaMuSizeEffect: (numeric vector number of genes) [Default NULL]
#'    Mean of normal distribution of which scaNLing factor for 
#'    size effects per sample are drawn.
#' @param scaSDSizeEffect: (numeric vector number of genes) [Default NULL]
#'    Standard deviation of normal distribution of which scaling factor for 
#'    size effects per sample are drawn.
#' @param scaMuBatchEffect: (numeric vector number of genes) [Default NULL]
#'    Mean of normal distribution of which scaling factor for 
#'    batch effects per gene are drawn (reference is batch A).
#' @param scaSDBatchEffect: (numeric vector number of genes) [Default NULL]
#'    Standard deviation of normal distribution of which scaling factor for 
#'    batch effects per gene are drawn (reference is batch A).
#' @param vecNormConstExternal: (numeric vector number of cells)
#'    [Default NULL]
#'    Size factors for data set. Size factors are set to 1 if this is
#'    not specified (NULL).
#' @param dirOutSimulation: (str directory)
#'    Directory to which simulated parameter objects are 
#'    saved to.
#' 
#' @return list: (length 2)
#' \itemize{
#'    \item vecPT: (numerical vector length number of cells)
#'    Pseudotime coordinates (1D) of cells: One scalar per cell.
#'    \item matSampledCountsObserved: (matrix genes x cells)
#'    Sampled count data of all cells after drop-out.
#' }
#'    
#' @author David Sebastian Fischer
#' 
#' @export

simulateDataSetImpulseDE2 <- function(vecTimePointsA,
  vecTimePointsB,
  boolCaseCtrl=FALSE,
  scaNConst,
  scaNImp,
  scaNLin,
  scaNSig,
  scaMumax=1000,
  scaSDImpulseAmplitude=1,
  scaMuSizeEffect=1,
  scaSDSizeEffect=0.1,
  scaMuBatchEffect=NULL,
  scaSDBatchEffect=NULL,
  dirOutSimulation){
  
  ####
  # Internal functions
  # Evalute impulse model at time points
  evalImpulse <- function(t,beta,t1,t2,h0,h1,h2){
    return(1/h1* (h0+(h1-h0)*1/(1+exp(-beta*(t-t1))))*
        (h2+(h1-h2)*1/(1+exp(beta*(t-t2)))))
  }
  # Evalute sigmoid model at time points
  evalSigmoid <- function(t,beta,t1,h0,h1){
    return(h0+(h1-h0)/(1+exp(-beta*(t-t1))))
  }
  
  ####
  # Simulate data (case only)
  # Create annotation and sample naming
  vecSamplesA <- vecTimePointsA
  names(vecSamplesA) <- sapply(seq(1, length(vecSamplesA)), function(i){
    paste0("A_", vecSamplesA[i], "_Rep", 
      match(i, which(vecTimePointsA==vecTimePointsA[i]) ))
  })
  vecSamplesB <- vecTimePointsB
  if(!is.null(vecSamplesB)){
    names(vecSamplesB) <- sapply(seq(1, length(vecSamplesB)), function(i){
      paste0("B_", vecSamplesB[i], "_Rep", 
        match(i, which(vecSamplesB==vecSamplesB[i]) ))
    })
  }
  vecSamples <- c(vecSamplesA, vecSamplesB)
  scaNSamples <- length(vecSamples)
  vecTimePointsUnique <- unique(vecSamples)
  vecindTimePointAssign <- match(vecSamples, vecTimePointsUnique)
  vecTimePointsUniqueB <- unique(vecSamplesB)
  vecindTimePointAssignB <- match(vecSamplesB, vecTimePointsUniqueB)
  
  dfAnnotation <- data.frame(
    Sample=names(vecSamples),
    Condition=rep("case", length(vecSamples)),
    Time=vecSamples,
    LongitudinalSeries=c(rep("A", length(vecSamplesA)),
      rep("B", length(vecSamplesB))),
    stringsAsFactors=FALSE
  )
  rownames(dfAnnotation) <- dfAnnotation$Sample
  if(boolCaseCtrl){
    dfAnnotation[dfAnnotation$LongitudinalSeries=="B",]$Condition <- rep("ctrl", sum(dfAnnotation$LongitudinalSeries=="B"))
  }
  
  # 1. Create hidden data set
  vecConstIDs <- paste0(rep("_",scaNConst),c(1:scaNConst))
  vecImpulseIDs <- paste0(rep("_",scaNImp),c((scaNConst+1):(scaNConst+scaNImp)))
  vecLinIDs <- paste0(rep("_",scaNLin),c((scaNConst+scaNImp+1):(scaNConst+scaNImp+scaNLin)))
  vecSigIDs <- paste0(rep("_",scaNSig),c((scaNConst+scaNImp+scaNLin+1):(scaNConst+scaNImp+scaNLin+scaNSig)))
  
  scaNGenes <- scaNConst+scaNImp+scaNLin+scaNSig
  
  # a. Draw means from uniform (first half of genes): one mean per gene
  vecMuConstHidden <- runif(scaNConst)*scaMumax
  matMuConstHidden <- matrix(vecMuConstHidden,
    nrow=scaNConst,
    ncol=scaNSamples,
    byrow=FALSE )
  rownames(matMuConstHidden) <- vecConstIDs
  colnames(matMuConstHidden) <- names(vecSamples)
  
  # b. Draw means from impulse model
  beta <- runif(scaNImp)*2+0.5
  ta <- runif(scaNImp)*max(vecTimePointsUnique)
  tb <- runif(scaNImp)*max(vecTimePointsUnique)
  t1 <- ta
  t1[tb < ta] <- tb[tb < ta]
  t2 <- tb
  t2[tb < ta] <- ta[tb < ta]
  h0 <- runif(scaNImp)*scaMumax
  h1 <- h0*abs(rnorm(n=scaNImp, mean=0,sd=scaSDImpulseAmplitude))
  h2 <- h0*abs(rnorm(n=scaNImp, mean=0,sd=scaSDImpulseAmplitude))
  h0[h0<0.00001] <-0.00001
  h1[h1<0.00001] <-0.00001
  h2[h2<0.00001] <-0.00001
  lsMuImpulseHidden <- lapply(seq(1,scaNImp), function(gene){
    evalImpulse(t=vecTimePointsUnique,
      beta=beta[gene],
      t1=t1[gene],
      t2=t2[gene],
      h0=h0[gene],
      h1=h1[gene],
      h2=h2[gene])[vecindTimePointAssign]
  })
  matImpulseModelHidden <- cbind(beta, h0, h1, h2, t1, t2)
  matMuImpulseHidden <- do.call(rbind, lsMuImpulseHidden)
  rownames(matImpulseModelHidden) <- vecImpulseIDs
  rownames(matMuImpulseHidden) <- vecImpulseIDs
  colnames(matMuImpulseHidden) <- names(vecSamples)
  if(boolCaseCtrl){
    # Keep start point the same, this doesn't work well for impulse model
    scaNImpDE <- round(scaNImp/2)
    betaDE <- runif(scaNImpDE)*2+0.5
    ta <- runif(scaNImpDE)*max(vecTimePointsUniqueB)
    tb <- runif(scaNImpDE)*max(vecTimePointsUniqueB)
    t1DE <- ta
    t1DE[tb < ta] <- tb[tb < ta]
    t2DE <- tb
    t2DE[tb < ta] <- ta[tb < ta]
    h1DE <- h0*abs(rnorm(n=scaNImpDE, mean=0,sd=scaSDImpulseAmplitude))
    h2DE <- h0*abs(rnorm(n=scaNImpDE, mean=0,sd=scaSDImpulseAmplitude))
    h1DE[h1<0.00001] <-0.00001
    h2DE[h2<0.00001] <-0.00001
    lsMuImpulseHiddenDE <- lapply(seq(1,scaNImpDE), function(gene){
      evalImpulse(t=vecTimePointsUniqueB,
                  beta=betaDE[gene],
                  t1=t1DE[gene],
                  t2=t2DE[gene],
                  h0=h0[gene],
                  h1=h1DE[gene],
                  h2=h2DE[gene])[vecindTimePointAssignB]
    })
    matMuImpulseHiddenDE <- do.call(rbind, lsMuImpulseHiddenDE)
    rownames(matMuImpulseHiddenDE) <- vecImpulseIDs[1:scaNImpDE]
    colnames(matMuImpulseHiddenDE) <- names(vecSamplesB)
    matMuImpulseHidden[seq(1,scaNImpDE),names(vecSamplesB)] <- matMuImpulseHiddenDE
  }

  
  # c. Linear functions
  # Draw linear model parameters
  vecInitialLevel <- runif(scaNLin)*scaMumax
  vecFinalLevel <- runif(scaNLin)*scaMumax
  # Evaluate linear functions
  scaDeltaTtot <- max(vecTimePointsUnique)-min(vecTimePointsUnique)
  matMuLinHidden <- do.call(rbind, lapply(seq(1, scaNLin), function(i){
    (vecInitialLevel[i]+(vecFinalLevel[i]-vecInitialLevel[i])/scaDeltaTtot*
        (vecTimePointsUnique-min(vecTimePointsUnique)))[vecindTimePointAssign]
  }))
  rownames(matMuLinHidden) <- vecLinIDs
  colnames(matMuLinHidden) <- names(vecSamples)
  if(boolCaseCtrl){
    # Keep start point the same, this doesn't work well for impulse model
    scaNLinDE <- round(scaNLin/2)
    vecFinalLevelDE <- runif(scaNLinDE)*scaMumax
    # Evaluate linear functions
    scaDeltaTtot <- max(vecTimePointsUniqueB)-min(vecTimePointsUniqueB)
    matMuLinHiddenDE <- do.call(rbind, lapply(seq(1, scaNLinDE), function(i){
      (vecInitialLevel[i]+(vecFinalLevelDE[i]-vecInitialLevel[i])/scaDeltaTtot*
         (vecTimePointsUniqueB-min(vecTimePointsUniqueB)))[vecindTimePointAssignB]
    }))
    rownames(matMuLinHiddenDE) <- vecLinIDs[1:scaNLinDE]
    colnames(matMuLinHiddenDE) <- names(vecSamplesB)
    matMuLinHidden[seq(1,scaNImpDE),names(vecSamplesB)] <- matMuLinHiddenDE
  }
  
  # d. Sigmoid functions
  # Draw sigmoid model parameters
  vech0 <- runif(scaNSig)*scaMumax
  vech1 <- runif(scaNSig)*scaMumax
  vecBeta <- runif(scaNSig)*4+0.5
  vecT1 <- runif(scaNSig)*(max(vecTimePointsUnique)-min(vecTimePointsUnique))+min(vecTimePointsUnique)
  # Evaluate sigmoid functions
  matMuSigHidden <- do.call(rbind, lapply(seq(1, scaNSig), function(i){
    evalSigmoid(t=vecTimePointsUnique,
      beta=vecBeta[i],
      t1=vecT1[i],
      h0=vech0[i],
      h1=vech1[i])[vecindTimePointAssign]
  }))
  rownames(matMuSigHidden) <- vecSigIDs
  colnames(matMuSigHidden) <- names(vecSamples)
  if(boolCaseCtrl){
    # Keep start point the same, this doesn't work well for impulse model
    scaNSigDE <- round(scaNSig/2)
    vech1DE <- runif(scaNSigDE)*scaMumax
    vecBetaDE <- runif(scaNSigDE)*4+0.5
    vecT1DE <- runif(scaNSigDE)*(max(vecTimePointsUniqueB)-min(vecTimePointsUniqueB))+
      min(vecTimePointsUniqueB)
    # Evaluate sigmoid functions
    matMuSigHiddenDE <- do.call(rbind, lapply(seq(1, scaNSigDE), function(i){
      evalSigmoid(t=vecTimePointsUniqueB,
                  beta=vecBetaDE[i],
                  t1=vecT1DE[i],
                  h0=vech0[i],
                  h1=vech1DE[i])[vecindTimePointAssignB]
    }))
    rownames(matMuSigHiddenDE) <- vecSigIDs[1:scaNSigDE]
    colnames(matMuSigHiddenDE) <- names(vecSamplesB)
    matMuSigHidden[seq(1,scaNSigDE),names(vecSamplesB)] <- matMuSigHiddenDE
  }
  
  # e. Merge data
  matMuHidden <- do.call(rbind, list(matMuConstHidden, 
    matMuImpulseHidden,
    matMuLinHidden,
    matMuSigHidden))
  if(boolCaseCtrl){
    vecCaseCtrlDEIDs <- c(vecImpulseIDs[1:scaNImpDE],
                           vecLinIDs[1:scaNLinDE],
                           vecSigIDs[1:scaNSigDE])
  }
  
  # Add scaling factors
  # Sample size factors
  vecSizeFactorsHidden <- rnorm(n=scaNSamples, 
    mean=scaMuSizeEffect, sd=scaSDSizeEffect)
  vecSizeFactorsHidden[vecSizeFactorsHidden<0.1] <- 0.1
  vecSizeFactorsHidden[vecSizeFactorsHidden>10] <- 10
  names(vecSizeFactorsHidden) <- names(vecSamples)
  # Sample batch normalisation factors
  if(!is.null(vecSamplesB) & !boolCaseCtrl){
    vecLongitudinalFactorsHidden <- rnorm(n=scaNGenes,
                                          mean=scaMuBatchEffect, sd=scaSDBatchEffect)
    vecLongitudinalFactorsHidden[vecLongitudinalFactorsHidden<0.1] <- 0.1
    vecLongitudinalFactorsHidden[vecLongitudinalFactorsHidden>10] <- 10
    names(vecLongitudinalFactorsHidden) <- rownames(matMuHidden)
  }
  # Scale
  matMuHiddenScaled <- matMuHidden*
    matrix(vecSizeFactorsHidden,
      nrow=dim(matMuHidden)[1],
      ncol=dim(matMuHidden)[2], byrow=TRUE)
  if(!is.null(vecSamplesB) & !boolCaseCtrl){
    matMuHiddenScaled <- matMuHiddenScaled*
      cbind(
        matrix(1,
          nrow=dim(matMuHidden)[1],
          ncol=length(vecSamplesA), byrow=FALSE),
        matrix(vecLongitudinalFactorsHidden,
          nrow=dim(matMuHidden)[1],
          ncol=length(vecSamplesB), byrow=FALSE)
      )
  }
  rownames(matMuHiddenScaled) <- rownames(matMuHidden)
  colnames(matMuHiddenScaled) <- colnames(matMuHidden)
  
  # e. draw dispersions by gene
  vecDispHidden <- runif(dim(matMuHiddenScaled)[1])+10
  names(vecDispHidden) <- rownames(matMuHiddenScaled)
  
  # f. add noise - draw from negative binomial
  matObservedData <- do.call(rbind, lapply(seq(1,dim(matMuHiddenScaled)[1]), function(gene){
    sapply(vecSamples, function(sample){
      rnbinom(n=1, mu=matMuHiddenScaled[gene,sample], size=vecDispHidden[sample])
    })
  }))
  rownames(matObservedData) <- rownames(matMuHiddenScaled)
  colnames(matObservedData) <- colnames(matMuHiddenScaled)
    
  #  Counts
  matObservedCounts <- round(matObservedData)
  
  # Save simulation
  save(vecSamples,file=file.path(dirOutSimulation,"Simulation_vecPT.RData"))
  
  save(vecConstIDs,file=file.path(dirOutSimulation,"Simulation_vecConstIDs.RData"))
  save(vecImpulseIDs,file=file.path(dirOutSimulation,"Simulation_vecImpulseIDs.RData"))
  save(vecLinIDs,file=file.path(dirOutSimulation,"Simulation_vecLinIDs.RData"))
  save(vecSigIDs,file=file.path(dirOutSimulation,"Simulation_vecSigIDs.RData"))
  
  save(vecSizeFactorsHidden,file=file.path(dirOutSimulation,"Simulation_vecSizeFactorsHidden.RData"))
  if(!is.null(vecSamplesB) & !boolCaseCtrl) { save(vecLongitudinalFactorsHidden,file=file.path(dirOutSimulation,"Simulation_vecLongitudinalFactorsHidden.RData")) }
  if(boolCaseCtrl){ save(vecCaseCtrlDEIDs,file=file.path(dirOutSimulation,"Simulation_vecCaseCtrlDEIDs.RData")) }
  
  save(vecDispHidden,file=file.path(dirOutSimulation,"Simulation_vecDispHidden.RData"))
  save(matImpulseModelHidden,file=file.path(dirOutSimulation,"Simulation_matImpulseModelHidden.RData"))
  save(matMuHidden,file=file.path(dirOutSimulation,"Simulation_matMuHidden.RData"))
  save(matMuHiddenScaled,file=file.path(dirOutSimulation,"Simulation_matMuHiddenScaled.RData"))
  save(matObservedData,file=file.path(dirOutSimulation,"Simulation_matObservedData.RData"))
  save(matObservedCounts,file=file.path(dirOutSimulation,"Simulation_matObservedCounts.RData"))
  
  return(list( dfAnnotation=dfAnnotation,
    matObservedCounts=matObservedCounts ))
}