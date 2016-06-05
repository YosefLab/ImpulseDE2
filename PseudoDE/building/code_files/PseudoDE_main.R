library(BiocParallel)
library(ggplot2)
#library(scone)
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/srcSCONE_zinb.R")
library(MASS)
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files/ImpulseDE2_main.R")

source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/srcPseudoDE_clusterCellsInPseudotime.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/srcPseudoDE_createAnnotation.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/srcPseudoDE_plotZINBfits.R")
source("/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/PseudoDE/building/code_files/srcPseudoDE_plotPseudotimeClustering.R")


#' 
#' @param matCounts: (matrix genes x cells)
#'    Count data of all cells, unobserved entries are NA.
#' @param vecPseudotime: (numerical vector length number of cells)
#'    Pseudotime coordinates (1D) of cells: One scalar per cell.
#'    Has to be named: Names of elements are cell names.
#' @param boolDEAnalysisImpulseModel: (bool) [Default TRUE]
#'    Whether to perform differential expression analysis with ImpulseDE2.
#' @param boolDEAnalysisModelFree: (bool) [Default FALSE]
#'    Whether to perform model-free differential expression analysis.
#' @param nProc: (scalar) [Default 1] Number of processes for 
#'    parallelisation. The specified value is internally changed 
#'    to \code{min(detectCores() - 1, nProc)} using the 
#'    \code{detectCores} function from the package \code{parallel} 
#'    to avoid overload.
#' 
#' @return (list length 2)
#'    \itemize{
#'    \item lsImpulseDE2results: (list length 4)
#'    \itemize{
#'      \item vecDEGenes: (list number of genes) Genes IDs identified
#'        as differentially expressed by ImpulseDE2 at threshold \code{Q_value}.
#'      \item dfImpulseResults: (data frame) ImpulseDE2 results.
#'      \item lsImpulseFits: (list) List of matrices which
#'        contain parameter fits and model values for given time course for the
#'        case condition (and control and combined if control is present).
#'        Each parameter matrix is called parameter_'condition' and has the form
#'        (genes x [beta, h0, h1, h2, t1, t2, logL_H1, converge_H1, mu, logL_H0, 
#'        converge_H0]) where beta to t2 are parameters of the impulse
#'        model, mu is the single parameter of the mean model, logL are
#'        log likelihoods of full (H1) and reduced model (H0) respectively, converge
#'        is convergence status of numerical optimisation of model fitting by
#'        \code{optim} from \code{stats} of either model. Each value matrix is called
#'        value_'condition' and has the form (genes x time points) and contains the
#'        counts predicted by the impulse model at the observed time points.
#'      \item dfDESeq2Results: (NULL) DESeq2 results, DESeq2 is not run within
#'        ImpulseDE2 in singlecell mode.
#'    }
#'    \item dfModelFreeDEAnalysis: (data frame) 
#'        Summary of model-free differential expression analysis.
#'    }   
#' @export

runPseudoDE <- function(matCounts, 
  vecPseudotime,
  boolDEAnalysisImpulseModel = TRUE,
  boolDEAnalysisModelFree = FALSE,
  nProc=1){
  
  # 1. Data preprocessing
  print("1. Data preprocessing:")
  tm_preproc <- system.time({
    lsProcessedSCData <- processSCData(matCounts=matCounts,
      vecPseudotime=vecPseudotime,
      boolDEAnalysisImpulseModel=boolDEAnalysisImpulseModel,
      boolDEAnalysisModelFree=boolDEAnalysisModelFree )
    matCountsProc <- lsProcessedSCData$matCountsProc
    vecPseudotimeProc <- lsProcessedSCData$vecPseudotimeProc
  })
  save(matCountsProc,file=file.path(getwd(),"PseudoDE_matCountsProc.RData"))
  save(vecPseudotimeProc,file=file.path(getwd(),"PseudoDE_vecPseudotimeProc.RData"))
  
  # 2. Cluster cells in pseudo-time
  print("2. Clustering:")
  tm_clustering <- system.time({
    # Cluster in pseudotime
    lsResultsClustering <- clusterCellsInPseudotime(vecPseudotime=vecPseudotimeProc)
    # Plot clustering
    plotPseudotimeClustering(vecPseudotime=vecPseudotime, 
      lsResultsClustering=lsResultsClustering)
  })
  save(lsResultsClustering,file=file.path(getwd(),"PseudoDE_lsResultsClustering.RData"))
  print(paste("Time elapsed during clustering: ",round(tm_clustering["elapsed"]/60,2),
    " min",sep=""))
  
  # 3. Create annotation table
  print("3. Create annotation table")
  dfAnnotation <- createAnnotationByCluster(matCounts=matCountsProc,
    vecPseudotime=vecPseudotimeProc,
    lsResultsClustering=lsResultsClustering)
  # Fit by cell
  #dfAnnotation <- createAnnotationByCell(matCounts=matCounts,
  #  vecPseudotime=vecPseudotime)
  save(dfAnnotation,file=file.path(getwd(),"PseudoDE_dfAnnotation.RData"))
  
  # 4. Fit mixture model
  print("4. Fit mixture model:")
  tm_fitmm <- system.time({
    lsResZINBFits <- fitZINB(matCounts, lsResultsClustering, nProc)
    vecDispersions <- lsResZINBFits$vecDispersions
    matDropout <- lsResZINBFits$matDropout
    matProbNB  <- lsResZINBFits$matProbNB
    matCountsImputed  <- lsResZINBFits$matCountsImputed
    matClusterMeansFitted  <- lsResZINBFits$matClusterMeansFitted
    matConvergence <- lsResZINBFits$matConvergence
  })
  save(vecDispersions,file=file.path(getwd(),"PseudoDE_vecDispersions.RData"))
  save(matDropout,file=file.path(getwd(),"PseudoDE_matDropout.RData"))
  save(matProbNB,file=file.path(getwd(),"PseudoDE_matProbNB.RData"))
  save(matCountsImputed,file=file.path(getwd(),"PseudoDE_matCountsImputed.RData"))
  save(matClusterMeansFitted,file=file.path(getwd(),"PseudoDE_matClusterMeansFitted.RData"))
  
  save(matCountsClean,file=file.path(getwd(),"PseudoDE_matCountsClean.RData"))
  save(matClusterMeansFitted,file=file.path(getwd(),"PseudoDE_matClusterMeansFitted.RData"))
  print(paste("Time elapsed during ZINB fitting: ",round(tm_fitmm["elapsed"]/60,2),
    " min",sep=""))
  
  graphics.off()
  # 5. Plot ZINB fits to data.
  if(boolPlotZINBfits){
    print("5. Plot ZINB fits to data.")
    plotZINBfits(vecGeneIDs=rownames(matCountsClean)[1:10], 
      matCounts=matCountsClean,
      matClusterMeans=matClusterMeansFitted, 
      vecDispersions=vecDispersions,
      matProbNB=matProbNB,
      vecClusterAssignments=lsResultsClustering$Assignments,
      lsResultsClustering=lsResultsClustering,
      dfAnnotation=dfAnnotation, 
      strPDFname="PseudoDE_ZINBfits.pdf" )
  }
  
  # 6. Differential expression analysis:
  print("6. Differential expression analysis:")
  if(boolDEAnalysisImpulseModel){
    # a) Model-based: Impulse model
    print("Differential expression analysis: Model-based (Impulse model)")
    print("### Begin ImpulseDE2 output ################################")
    tm_deanalysis_impulse <- system.time({
      lsInputToImpulseDE2 <- list(matDropout, matProbNB, matCountsImputed, matClusterMeansFitted)
      names(lsInputToImpulseDE2) <- c("matDropout", "matProbNB", "matCountsImputed", "matClusterMeansFitted")
      
      lsImpulseDE2results <- runImpulseDE2(
        matCountData = matCountsClean, 
        dfAnnotation = dfAnnotation,
        strCaseName = "case", 
        strControlName = NULL, 
        strMode = "singlecell", 
        nProc = nProc, 
        Q_value = 0.01,
        boolPlotting = TRUE,
        lsPseudo = lsInputToImpulseDE2,
        vecDispersionsExternal = vecDispersions,
        boolRunDESeq2 = FALSE,
        boolSimplePlot = TRUE, 
        boolLogPlot = TRUE )
    })
    print("### End ImpulseDE2 output ##################################")
    save(lsInputToImpulseDE2,file=file.path(getwd(),"PseudoDE_lsInputToImpulseDE2.RData"))
    save(lsImpulseDE2results,file=file.path(getwd(),"PseudoDE_lsImpulseDE2results.RData"))
    print(paste("Time elapsed during differential expression analysis with ImpulseDE2: ",
      round(tm_deanalysis_impulse["elapsed"]/60,2)," min",sep=""))
  } else {
    lsImpulseDE2results <- NULL
  }
  if(boolDEAnalysisModelFree){
    # b) Model-based: Impulse model
    print("Differential expression analysis: Model-free")
    tm_deanalysis_mf <- system.time({
      dfModelFreeDEAnalysis <- runModelFreeDEAnalysis(
        matCountData = matCountsClean, 
        dfAnnotation = dfAnnotation,
        strMode = "singlecell", 
        nProc = nProc, 
        Q_value = 0.01,
        boolPlotting = TRUE,
        lsPseudo = lsInputToImpulseDE2,
        vecDispersionsExternal = vecDispersions,)
    })
    save(dfModelFreeDEAnalysis,file=file.path(getwd(),"PseudoDE_dfModelFreeDEAnalysis.RData"))
    print(paste("Time elapsed during model-free differential expression analysis: ",
      round(tm_deanalysis_mf["elapsed"]/60,2)," min",sep=""))
  } else {
    dfModelFreeDEAnalysis <- NULL
  }
  
  print("Completed PseudoDE.")
  return(list(lsImpulseDE2results=lsImpulseDE2results,
    dfModelFreeDEAnalysis=dfModelFreeDEAnalysis))
}