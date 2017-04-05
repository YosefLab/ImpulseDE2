### Estimate dispersions with DESeq2

#' Wrapper function for running DESeq2
#' 
#' Run DESeq2 and extract dispersion parameter estimates.
#' Catch and remove dispersion outlier exception on samples
#' with zero-count observations.
#' 
#' @seealso Called by \link{runImpulseDE2}.
#' 
#' @param dfAnnotationProc (data frame samples x covariates) 
#' {Sample, Condition, Time (numeric), TimeCateg (str)
#' (and confounding variables if given).}
#' Processed annotation table with covariates for each sample.
#' @param matCountDataProc (matrix genes x samples)
#' Read count data.
#' @param boolCaseCtrl (bool) 
#' Whether to perform case-control analysis. Does case-only
#' analysis if FALSE.
#' @param vecConfounders (vector of strings number of 
#' confounding variables)
#' Factors to correct for during batch correction. Have to 
#' supply dispersion factors if more than one is supplied.
#' Names refer to columns in dfAnnotationProc.
#' 		
#' @return (numeric vector length number of genes)
#' Dispersion parameter estimates for each gene.
#' In format of parameter size of \link{dnbinom}
#' which is 1/dispersion factor of DESeq2.
#'    
#' @import DESeq2
#' @importFrom S4Vectors mcols
#' 
#' @author David Sebastian Fischer
runDESeq2 <- function(
    dfAnnotationProc, 
    matCountDataProc,
    boolCaseCtrl,
    vecConfounders){
    
    # Get gene-wise dispersion estimates
    # var = mean + alpha * mean^2, alpha is dispersion
    # DESeq2 dispersion is 1/size used dnbinom (used in cost function
    # for evaluation of likelihood)
    
    dds <- NULL
    if(!boolCaseCtrl){
        # Case-only
        if(is.null(vecConfounders)){
            # No batch correction
            dds <- suppressWarnings( DESeqDataSetFromMatrix(
                countData = matCountDataProc,
                colData = dfAnnotationProc,
                design = ~ TimeCateg) )
            dds <- estimateSizeFactors(dds)
            dds <- estimateDispersions(dds)
            
        } else {
            # With batch correction
            # Catch non full-rank design matrix 
            # e.g. batches have of one confounder have mutually 
            # exclusive sets of time points
            tryCatch({
                dds <- suppressWarnings( DESeqDataSetFromMatrix(
                    countData = matCountDataProc,
                    colData = dfAnnotationProc,
                    design = ~TimeCateg + Batch  ) )
                dds <- estimateSizeFactors(dds)
                dds <- estimateDispersions(dds)
            }, error=function(strErrorMsg){
                print(strErrorMsg)
                print(paste0(
                    "WARNING: DESeq2 failed on full model ",
                    "- dispersions may be inaccurate.",
                    "Estimating dispersions on reduced model ",
                    "formulation [full = ~Batch",
                    " reduced = ~1]. Supply externally generated",
                    " dispersion parameters via ",
                    "vecDispersionsExternal if there is a more ",
                    "accurate model for your data set."))
                warning("Warning generated in dispersion factor",
                        " estimation, read stdout.")
            }, finally={
                if(is.null(dds)){
                    dds <- suppressWarnings( DESeqDataSetFromMatrix(
                        countData = matCountDataProc,
                        colData = dfAnnotationProc,
                        design = ~Batch) )
                    dds <- estimateSizeFactors(dds)
                    dds <- estimateDispersions(dds)
                }
            })
        }
    } else {
        # Case-control
        if(is.null(vecConfounders)){
            # No batch correction
            # Catch non full-rank design matrix 
            # e.g. conditions have mutually 
            #exclusive sets of timepoints
            tryCatch({
                dds <- suppressWarnings( DESeqDataSetFromMatrix(
                    countData = matCountDataProc,
                    colData = dfAnnotationProc,
                    design = ~Condition+Condition:TimeCateg) )
                dds <- estimateSizeFactors(dds)
                dds <- estimateDispersions(dds)
            }, error=function(strErrorMsg){
                print(strErrorMsg)
                print(paste0(
                    "WARNING: DESeq2 failed on full model ",
                    "- dispersions may be inaccurate.",
                    "Estimating dispersions on reduced model",
                    " formulation [full = ~Condition",
                    " reduced = ~1]. Supply externally ",
                    "generated dispersion parameters via ",
                    "vecDispersionsExternal if there is a more",
                    " accurate model for your data set."))
                warning("Warning generated in dispersion factor",
                        " estimation, read stdout.")
            }, finally={
                if(is.null(dds)){
                    dds <- suppressWarnings( DESeqDataSetFromMatrix(
                        countData = matCountDataProc,
                        colData = dfAnnotationProc,
                        design = ~Condition) )
                    dds <- estimateSizeFactors(dds)
                    dds <- estimateDispersions(dds)
                }
            })
            
        } else {
            # With batch correction
            # Catch non full-rank design matrix 
            # e.g.  a) batches have of one confounder have mutually 
            #exclusive sets of time points or
            # b) conditions have mutually exclusive sets of time points.
            # Note that full rank of design matrix of batches is 
            # checked in data processing.
            # Note that design matrix of batches and condition is full 
            # rank if batch
            # matrix is full rank due to the nested batches trick 
            # (using Condition:Batch below).
            # Therefore, error catching tailored to a) and b).
            tryCatch({
                dds <- suppressWarnings( DESeqDataSetFromMatrix(
                    countData = matCountDataProc,
                    colData = dfAnnotationProc,
                    design = 
                        ~Condition+Condition:TimeCateg+Condition:BatchNested) )
                dds <- estimateSizeFactors(dds)
                dds <- estimateDispersions(dds)
            }, error=function(strErrorMsg){
                print(strErrorMsg)
                print(paste0(
                    "WARNING: DESeq2 failed on full model ",
                    "- dispersions may be inaccurate.",
                    "Estimating dispersions on reduced model.",
                    " Supply externally generated dispersion parameters via ",
                    "vecDispersionsExternal if there is a more accurate ",
                    "model for your data set."))
                warning("Warning generated in dispersion factor estimation.",
                        " Eead stdout.")
            }, finally={
                if(is.null(dds)){
                    matModelMatrixBatches <- do.call(cbind, lapply(
                        vecConfounders, function(confounder){
                            match(dfAnnotationProc[,confounder], 
                                  unique(dfAnnotationProc[,confounder]))
                        }))
                    matModelMatrixBatchesTime <- cbind(
                        matModelMatrixBatches,
                        match(dfAnnotationProc$Time, 
                              unique(dfAnnotationProc$Time)))
                    boolFullRankBatchTime <- rankMatrix(
                        matModelMatrixBatchesTime)[1] == 
                        dim(matModelMatrixBatchesTime)[2]
                    if(!boolFullRankBatchTime){
                        paste0("Model matrix based on confounding variables {",
                               vecConfounders,
                               "} and time is not full rank: ",
                               "There are confounding variables with ",
                               " batch structures which are linear",
                               " combinations of the time points.")
                        paste0("Using reduced model formulation ",
                               "[full= ~Condition+Condition:TimeCateg, ",
                               "reduced= ~TimeCateg].")
                        dds <- suppressWarnings( DESeqDataSetFromMatrix(
                            countData = matCountDataProc,
                            colData = dfAnnotationProc,
                            design = ~Condition+Condition:TimeCateg) )
                        dds <- estimateSizeFactors(dds)
                        dds <- estimateDispersions(dds)
                    } else {
                        paste0("Found 1. or {1. and 2.}:")
                        paste0("1. Model matrix based on confounding ",
                               "variables {", vecConfounders,
                               "} and time is not full rank:",
                               " There are confounding variables with ",
                               " batch structures which are ",
                               "linear combinations of the time points.")
                        paste0(" 2. Model matrix based on condition",
                               " and time is not full rank: ",
                               "Conditions case and control have",
                               " mutually exclusive sets of timepoints.")
                        paste0("Using reduced model formulation ",
                               "[full= ~Condition+Condition:BatchNested, ",
                               "reduced= ~1].")
                        dds <- suppressWarnings( DESeqDataSetFromMatrix(
                            countData = matCountDataProc,
                            colData = dfAnnotationProc,
                            design = ~Condition+Condition:BatchNested) )
                        dds <- estimateSizeFactors(dds)
                        dds <- estimateDispersions(dds)
                    }
                }
            })
        }
    }
    
    vecDispersionsInv <- mcols(dds)$dispersion
    # Catch dispersion trend outliers at the upper boundary
    # (alpha = 20 ->large variance)
    # which contain zero measurements: 
    # The zeros throw off the dispersion estimation 
    # in DESeq2 which may converge to the upper bound 
    # even though the obesrved variance is small.
    # Avoid outlier handling and replace estimates by 
    # MAP which is more stable in these cases.
    # Note: the upper bound 20 is not exactly reached - 
    # there is numeric uncertainty here - use > rather than ==
    vecindDESeq2HighOutliesFailure <- !is.na(mcols(dds)$dispOutlier) & 
        mcols(dds)$dispOutlier==TRUE &
        mcols(dds)$dispersion>(20-10^(-5)) & 
        apply(matCountDataProc, 1, function(gene) any(gene==0) )
    vecDispersionsInv[vecindDESeq2HighOutliesFailure] <- 
        mcols(dds)$dispMAP[vecindDESeq2HighOutliesFailure]
    if(sum(vecindDESeq2HighOutliesFailure)>0){
        print(paste0("Corrected ", sum(vecindDESeq2HighOutliesFailure),
                     " DESEq2 dispersion estimates which ",
                     "to avoid variance overestimation and loss of ",
                     "discriminatory power for model selection."))
    }
    # DESeq2 uses alpha=1/phi as dispersion
    vecDispersions <- 1/vecDispersionsInv
    names(vecDispersions) <- rownames(dds)
    
    return(vecDispersions)
}