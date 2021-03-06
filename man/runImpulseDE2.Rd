% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ImpulseDE2_main.R
\name{runImpulseDE2}
\alias{runImpulseDE2}
\alias{ImpulseDE2}
\alias{wrapper}
\title{ImpulseDE2 wrapper}
\usage{
runImpulseDE2(matCountData = NULL, dfAnnotation = NULL,
  boolCaseCtrl = FALSE, vecConfounders = NULL, scaNProc = 1,
  scaQThres = NULL, vecDispersionsExternal = NULL,
  vecSizeFactorsExternal = NULL, boolIdentifyTransients = FALSE,
  boolVerbose = TRUE)
}
\arguments{
\item{matCountData}{(matrix genes x samples) [Default NULL] 
Read count data, unobserved entries are NA. 
Can be SummarizedExperiment object.}

\item{dfAnnotation}{(data frame samples x covariates) 
{Sample, Condition, Time (numeric), TimeCateg (str)
(and confounding variables if given).}
Annotation table with covariates for each sample.}

\item{boolCaseCtrl}{(bool) [Default FALSE]
Whether to perform case-control analysis. Does case-only
analysis if FALSE.}

\item{vecConfounders}{(vector of strings number of confounding variables)
Factors to correct for during batch correction. Have to 
supply dispersion factors if more than one is supplied.
Names refer to columns in dfAnnotation.}

\item{scaNProc}{(scalar) [Default 1] Number of processes for 
parallelisation.}

\item{scaQThres}{(scalar) [Default NULL] 
FDR-corrected p-value cutoff for significance.}

\item{vecDispersionsExternal}{(vector length number of
genes in matCountData) [Default NULL]
Externally generated list of gene-wise dispersion factors
which overides DESeq2 generated dispersion factors.}

\item{vecSizeFactorsExternal}{(vector length number of
cells in matCountData) [Default NULL]
Externally generated list of size factors which override
size factor computation in ImpulseDE2.}

\item{boolIdentifyTransients}{(bool) [Defaul FALSE]
Whether to identify transiently activated or deactivated 
genes. This involves an additional fitting of sigmoidal models
and hypothesis testing between constant, sigmoidal and impulse model.}

\item{boolVerbose}{(bool) [Default TRUE] Whether to print
progress to stdout.}
}
\value{
(object of class ImpulseDE2Object)
This object can be treated as a list with 2 elements:
(list length 2)
\itemize{
\item vecDEGenes (list number of genes) Genes IDs identified
as differentially expressed by ImpulseDE2 at threshold \code{scaQThres}.
\item dfDEAnalysis (data frame samples x reported characteristics) 
Summary of fitting procedure and 
differential expression results for each gene.
\itemize{
\item Gene: Gene ID.
\item p: P-value for differential expression.
\item padj: Benjamini-Hochberg false-discovery rate corrected p-value
for differential expression analysis.
\item loglik_full: Loglikelihood of full model.
\item loglik_red: Loglikelihood of reduced model.
\item df_full: Degrees of freedom of full model.
\item df_red: Degrees of freedom of reduced model
\item mean: Inferred mean parameter of constant model of first batch.
From combined samples in case-ctrl. 
\item allZero (bool) Whether there were no observed non-zero 
observations of this gene.
If TRUE, fitting and DE analsysis were skipped and entry is NA.
}
Entries only present in case-only DE analysis:
\itemize{
\item converge_impulse: Convergence status of optim for 
impulse model fit (full model).
\item converge_const: Convergence status of optim for 
constant model fit (reduced model).
}
Entries only present in case-control DE analysis:
\itemize{
\item converge_combined: Convergence status of optim for 
impulse model fit to case and control samples combined (reduced model).
\item converge_case: Convergence status of optim for 
impulse model fit to samples of case condition (full model 1/2).
\item converge_control: Convergence status of optim for 
impulse model fit to samples of control condition (full model 2/2).
}
Entries only present if boolIdentifyTransients is TRUE:
\itemize{
\item converge_sigmoid: Convergence status of optim for 
sigmoid model fit to samples of case condition.
\item impulseTOsigmoid_p: P-value of loglikelihood ratio test
impulse model fit versus sigmoidal model on samples of case condition.
\item impulseTOsigmoid_padj: Benjamini-Hochberg 
false-discovery rate corrected p-value of loglikelihood ratio test
impulse model fit versus sigmoid model on samples of case condition.
\item sigmoidTOconst_p: P-value of loglikelihood ratio test
sigmoidal model fit versus constant model on samples of case condition.
\item sigmoidTOconst_padj: Benjamini-Hochberg 
false-discovery rate corrected p-value of loglikelihood ratio test
sigmoidal model fit versus constant model on samples of case condition.
\item isTransient (bool) Whether gene is transiently
activated or deactivated and differentially expressed.
\item isMonotonous (bool) Whether gene is not transiently
activated or deactivated and differentially expressed. This scenario
corresponds to a montonous expression level increase or decrease.
}
}
}
\description{
Wrapper to run ImpulseDE2 on bulk omics count data.
This wrapper can perform the entire analysis pipeline of 
ImpulseDE2 on its own if the right parameters are supplied.
To run ImpulseDE2 on bulk omics count data, use the minimal
parameter set:
\itemize{
\item matCountData
\item dfAnnotation
\item boolCaseCtrl
\item vecConfounders
}
Additionally, you can provide:
\itemize{
\item scaNProc to set the number of processes for parallelisation.
\item scaQThres to set the cut off for your DE gene list. 
\item vecDispersionsExternal to supply external dispersion parameters
which may be necessary depending on your confounding factors (runImpulseDE2
will tell you if it is necessary).
\item vecSizeFactorsExternal to supply external size factors.
\item boolVerbose to control stdout output.
}
}
\details{
ImpulseDE2 is based on the impulse model proposed by
Chechik and Koller (Chechik and Koller, 2009).
The computational complexity of ImpulseDE2 is linear in the
number of genes and linear in the number of samples.
}
\examples{
lsSimulatedData <- simulateDataSetImpulseDE2(
vecTimePointsA   = rep(seq(1,8),3),
vecTimePointsB   = NULL,
vecBatchesA      = NULL,
vecBatchesB      = NULL,
scaNConst        = 30,
scaNImp          = 10,
scaNLin          = 10,
scaNSig          = 10)
objectImpulseDE2 <- runImpulseDE2(
matCountData    = lsSimulatedData$matObservedCounts, 
dfAnnotation    = lsSimulatedData$dfAnnotation,
boolCaseCtrl    = FALSE,
vecConfounders  = NULL,
scaNProc        = 1 )
head(objectImpulseDE2$dfImpulseDE2Results)

}
\seealso{
Calls the following functions:
\link{processData}, 
\link{runDESeq2},
\link{computeNormConst},
\link{fitModels},
\link{fitSigmoidModels},
\link{runDEAnalysis}. 
The following functions are additionally available to the user:
\link{fitSigmoidModels},
\link{plotGenes},  
\link{plotHeatmap},
\link{runDEAnalysis},   
\link{simulateDataSetImpulseDE2}.
}
\author{
David Sebastian Fischer
}
