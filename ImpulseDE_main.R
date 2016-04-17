
################################################################################
########################     Impulse DE package     ############################
################################################################################

### Version:  1.3
### Date:     2016
### Author v1.0:  Jil Sander
### Author v1.1:  David Fischer --Divide file into src files, annotation
### Author v1.2:  David Fischer --Support replicate samples, WLS fitting and residual scaling
### Author v1.3:  David Fischer --WLS fitting with CV weights
### Author v1.3:  David Fischer --NB fitting

################################################################################
# Variance extension: implement this inside functions, v1.2
# Longitudinal extension: this requires a simle mean substraction before the
#   entire procedure, implement this in main (as an extra function) v1.3
# Data imputation: do this after fitting based on inferred model. modify
#   function to deal with non set values, something like !is.nan evertime
#   values are selected v1.4
################################################################################

library(compiler)
library(amap)
library(longitudinal)
library(parallel)
library(graphics)
library(MASS)
library(DESeq2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                             Functions to call                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building")

#Prepares the annotation table for internal use
source("srcImpulseDE_annotation_preparation.R")

### Compute value of impulse function given parameters.
source("srcImpulseDE_calc_impulse.R")

### Cost function for parametric model fit: Ordinary least squares (two_impulses)
### or weighted least squares (two_impulses_WLS)
source("srcImpulseDE_CostFunctionFit.R")

## Compile fitting simple functions to make them quicker
calc_impulse_comp <- cmpfun(calc_impulse)
cost_fun_logl_comp <- cmpfun(cost_fun_logl)
cost_fun_logl_meanfit_comp <- cmpfun(cost_fun_logl_meanfit)

### Fits impulse model to a timecourse dataset
source("srcImpulseDE_impulse_fit.R")

# Detect differentially expressed genes over time
source("srcImpulseDE_DE_analysis.R")

### plots the impulse fits to timecourse data and also the control data if present
source("srcImpulseDE_plot_impulse.R")

################################################################################
################################################################################
################################################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                     Final function calling all others
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Differential expression analysis using impulse models
#'
#' Fits an impulse model to time course data and uses this model as a basis
#' to detect differentially expressed (DE) genes. If a single time course
#' data set is given, DE genes are detected over time, whereas if an
#' additional control time course data set is present, DE genes are
#' detected between both datasets.
#' @aliases impulse_DE impulseDE
#' @param expression_tables numeric matrix of expression values; genes should
#' be in rows, samples in columns. Data should be properly normalized and
#' log2-transformed as well as filtered for present or variable genes.
#' @param annotation_table table providing co-variables for the samples
#' including condition and time points. Time points must be numeric numbers.
#' @param colname_time character string specifying the column name of the
#' co-variable "\code{Time}" in \code{annotation_table}
#' @param colname_condition character string specifying the column name of
#' the co-variable "\code{Condition}" in \code{annotation_table}
#' @param control_timecourse logical indicating whether a control time
#' timecourse is part of the data set (\code{TRUE}) or not (\code{FALSE}).
#' Default is \code{FALSE}.
#' @param control_name character string specifying the name of the control
#' condition in \code{annotation_table}.
#' @param case_name character string specifying the name of the case
#' condition in \code{annotation_table}. Should be set if more than two
#' conditions are present in \code{annotation_table}.
#' @param expr_type character string with allowed values "\code{Array}" or
#' "\code{Seq}". Default is "\code{Array}".
#' @param plot_clusters logical indicating whether to plot the clusters
#' (\code{TRUE}) or not (\code{FALSE}). Default is \code{TRUE}.
#' @param n_iter numeric value specifying the number of iterations, which are
#' performed to fit the impulse model to the clusters. Default is \code{100}.
#' @param n_randoms numeric value specifying the number of generated randomized
#' background iterations, which are used for differential expression analysis.
#' Default is \code{50000} and this value should not be decreased.
#' @param n_process numeric value indicating the number of processes, which can
#' be used on the machine to run the background calculation in parallel. Default
#' is \code{4}. The specified value is internally changed to
#' \code{min(detectCores() - 1, n_process)} using the \code{detectCores}
#' function from the package \code{parallel} to avoid overload.
#' @param Q_value numeric value specifying the cutoff to call genes
#' significantly differentially expressed after FDR correction (adjusted
#' p-value). Default is \code{0.01}.
#' @return List containing the following elements:
#' \itemize{
#' \item \code{impulse_fit_results} List containing fitted values and model
#' parameters:
#' \itemize{
#' \item \code{impulse_parameters_combined} Matrix of fitted impulse model
#' parameters and sum of squared fitting errors for the combined dataset.
#' Not existing in the case of a single time course experiment.
#' \item \code{impulse_fits_combined} Matrix of impulse values calculated based
#' on the analyzed time points and the fitted model parameters for the combined
#' dataset. Not existing in the case of a single time course experiment.
#' \item \code{impulse_parameters_case} Matrix of fitted impulse model
#' parameters and sum of squared fitting errors for the case dataset.
#' \item \code{impulse_fits_case}  Matrix of impulse values calculated based
#' on the analyzed time points and the fitted model parameters for the case
#' dataset.
#' \item \code{impulse_parameters_control} Matrix of fitted impulse model
#' parameters and sum of squared fitting errors for the control dataset.
#' Not existing in the case of a single time course experiment.
#' \item \code{impulse_fits_control} Matrix of impulse values calculated based
#' on the analyzed time points and the fitted model parameters for the control
#' dataset. Not existing in the case of a single time course experiment.
#' }
#' \item \code{ImpulseDE_results} Matrix containing the names of genes being called
#' as differentially expressed according to the specified cutoff \code{Q_value}
#' together with the adjusted p-values.
#' }
#' Additionally, \code{ImpulseDE} saves the following objects and tables into
#' the working directory:
#' \itemize{
#' \item \code{prepared_annotation.RData} Object containing the internally used
#' modified version of \code{annotation_table}
#' \item \code{clus_out2.txt} Text-file saving std out from the multi-threading.
#' Can be ignored.
#' \item \code{kmeans_clus_final_combined_genes.txt} Text-file containing
#' each gene together with its corresponding cluster number for the combined
#' dataset. Not existing in the case of a single time course experiment.
#' \item \code{kmeans_clus_final_case_genes.txt} Text-file containing
#' each gene together with its corresponding cluster number for the case
#' dataset.
#' \item \code{kmeans_clus_final_control_genes.txt} Text-file containing
#' each gene together with its corresponding cluster number for the control
#' dataset. Not existing in the case of a single time course experiment.
#' \item \code{clusters_combined_genes.pdf} PDF-file containing
#' the clusters for the combined dataset. Not existing in the case of a single
#' time course experiment.
#' \item \code{clusters_case_genes.pdf} PDF-file containing
#' the clusters for the case dataset.
#' \item \code{clusters_control_genes.pdf} PDF-file containing
#' the clusters for the control dataset. Not existing in the case of a single
#' time course experiment.
#' \item \code{impulse_fit_clusters.RData} Object containing a list of the
#' fitted impulse model values and paramaters to the cluster means; structure
#' is the same as for the list element \code{impulse_fit_results} of the output
#' value.
#' \item \code{impulse_fit_genes.RData} Object containing a list of the
#' fitted impulse model values and paramaters to the genes; structure is the
#' same as for the list element \code{impulse_fit_results} of the output value.
#' \item \code{cluster_out_random.txt} Text-file saving std out from the
#' multi-threading related to the randomized data. Can be ignored.
#' \item \code{background_results.RData} F-score values generated based on the
#' fits to the randomized data as used for the differential expression analysis.
#' \item \code{impulse_DE_genes.RData} Object containing names of genes being
#' called as differentially expressed according to the specified cutoff
#' \code{Q_value} together with the adjusted p-values; same as the list element
#' \code{ImpulseDE_results} of the output value.
#' \item \code{pvals_and_flags.txt} Text-file containing all gene names
#' together with the adjusted p-values and flags for differential expression
#' according to additional tests.
#' }
#' @details \code{ImpulseDE} is based on the impulse model proposed by
#' Chechik and Koller, which reflects a two-step behavior of genes within a cell
#' responding to environmental changes (Chechik and Koller, 2009). To detect
#' differentially expressed genes, a five-step workflow is followed:
#' \enumerate{
#' \item The genes are clustered into a limited number of groups using k-means
#' clustering. If \code{plot_clusters} = \code{TRUE}, PDF documents are
#' generated, which contain plots of each cluster. Additionally, a text-file is
#' produced containing each gene together with its corresponding cluster number.
#' \item The impulse model is fitted to the mean expression profiles of the
#' clusters. The best parameter sets are then used for the next step.
#' \item The impulse model is fitted to each gene separately using the parameter
#' sets from step 2 as optimal start point guesses.
#' \item The impulse model is fitted to a randomized dataset (bootstrap), which
#' is essential to detect significantly differentially expressed genes
#' (Storey et al., 2005).
#' \item Detection of differentially expressed genes utilizing the fits to the
#' real and randomized data sets. FDR-correction is performed to obtain adjusted
#' p-values (Benjamini and Hochberg, 1995).
#' }
#' @examples
#' \dontrun{
#' #' Install package longitudinal and load it
#' library(longitudinal)
#' #' Attach datasets
#' data(tcell)
#' #' check dimension of data matrix of interest
#' dim(tcell.10)
#' #' generate a proper annotation table
#' annot <- as.data.frame(cbind("Time" =
#'    sort(rep(get.time.repeats(tcell.10)$time,10)),
#'    "Condition" = "activated"), stringsAsFactors = FALSE)
#' #' Time columns must be numeric
#' annot$Time <- as.numeric(annot$Time)
#' #' rownames of annotation table must appear in data table
#' rownames(annot) = rownames(tcell.10)
#' #' apply ImpulseDE in single time course mode
#' #' since genes must be in rows, transpose data matrix using t()
#' #' For the example, reduce random iterations to 100 and number of
#' #' used processors to 1
#' impulse_results <- impulse_DE(t(tcell.10), annot, "Time", "Condition",
#'    n_randoms = 50, n_process = 1)
#' }
#' @seealso \code{\link{plot_impulse}}, \code{\link{calc_impulse}}.
#' @author Jil Sander
#' @references Benjamini, Y. and Hochberg, Y. (1995) Controlling the false
#' discovery rate: a practical and powerful approach to multiple testing.
#' J. R. Stat. Soc. Series B Stat. Methodol., 57, 289-300.
#' @references Storey, J.D. et al. (2005) Significance analysis of time course
#' microarray experiments. Proc. Natl. Acad. Sci. USA, 102, 12837-12841.
#' @references Rangel, C., Angus, J., Ghahramani, Z., Lioumi, M., Sotheran, E.,
#' Gaiba, A., Wild, D.L., Falciani, F. (2004) Modeling T-cell activation using
#' gene expression profiling and state-space models. Bioinformatics, 20(9),
#' 1361-72.
#' @references Chechik, G. and Koller, D. (2009) Timing of Gene Expression
#' Responses to Envi-ronmental Changes. J. Comput. Biol., 16, 279-290.
#' @references Yosef, N. et al. (2013) Dynamic regulatory network controlling
#' TH17 cell differentiation. Nature, 496, 461-468.
#' @export

# INPUT
#   expression_tables: (list number of expression tables)
#       List of matrices (genes x samples) corresponding to replicates.
#       Matrices must have same format. This list will be transformed into
#       (Numeric 3D array genes x samples x replicates).
#       Contains expression values or similar locus-specific read-outs.
#       Note that this object contains conditions from CONTROL AND CASE.
#   annotation_table: (Table) providing co-variables for the samples including 
#       condition and time points. Time points must be numeric numbers.
#       Provide on 2D table to represent all replicates.
#   colname_time: (str) column name of the co-variable "Time" in 
#       annotation_table
#   colanme_condition: (str) column name of the co-variable "Condition" in 
#       annotation_table
#   control_timecourse: (bool) [Default FALSE]control time timecourse is part 
#       of the data set (TRUE) or not (FALSE).
#   control_name: (str) name of the control condition in annotation_table.
#   case_name: (str) name of the case condition in annotation_table.
#   expr_type: (str) ["Array" or "Seq"] In case of Sequencing data ("Seq") a 
#       DESeq2 test is added to the results. Default is "Array".
#   n_process: (scalar) [Default 4] number of processes, which can be used on the 
#       machine to run the background calculation in parallel
#   Q_value: (scalar) [default 0.01] cutoff to call genes significantly 
#       differentially expressed after FDR correction
# OUTPUT
#   impulse_fit_genes: (list) List containing fitted values and model parameters
#   impulse_DE_genes: names of genes being called as differentially expressed
#   ...more output saved to working directory, see documentation

impulse_DE <- function(expression_tables = NULL, annotation_table = NULL,
  colname_time = NULL, colname_condition = NULL, control_timecourse = FALSE,
  control_name = NULL, case_name = NULL, n_process = 4, Q_value = 0.01){
  
  print("Impulse v1.3 loglik based")
  NPARAM=6
  
  tm_tot <- system.time({
    
    print("Testing input")
    # Test annotation table
    replicates <- unique(annotation_table$Replicate)
    # Test that all replicates contain all timepoints
    for(iRep in 1:length(replicates)){
      if(!identical(unique(annotation_table$Time),unique(annotation_table[
        annotation_table$Replicate==replicates[iRep],]$Time))){
        print(paste0("ERROR: [Annotation table] Not all timepoints listed for replicate ",replicates[iRep]))
        print("Each replicate must contain a sample for each time point.")
        print("Aborting.")
        return(NULL)
      }
    }
    # Test that all Replicates contain the same number of samples
    for(iRep in 1:length(replicates)){
      # Establish number of samples in first replicate
      if(iRep==1){
        nSamples <- sum( annotation_table$Replicate==replicates[iRep] )
        # Check for equality of number of samples in remaining replicates
      } else {
        if(sum( annotation_table$Replicate==replicates[iRep] ) != nSamples){
          print(paste0("ERROR: [Annotation table] Replicate ",replicates[iRep],
            " contains different number of samples compared to replicate ",
            replicates[1]))
          stop("All replicates should contain the same number of samples.")
        }
      }
    }
    # Test that Replicate_name do not occur twice
    if(length(unique(annotation_table$Replicate_name)) != dim(annotation_table)[1]){
      print(paste0("ERROR: [Annotation table]",
        "Number of Replicate_name instances different to annotation table number of rows."))
      stop("Replicate_name cannot occur multiple times.")
    }
    
    # Test validity of expression table column naming with respect
    # to annotation table.
    # Test that all entries in annotation table occur in expression table
    if( sum(annotation_table$Replicate_name %in% colnames(expression_table)) < 
        length(annotation_table$Replicate_name) ){
      print(paste0("ERROR: [Expression table] Replicate_names ",
        as.character( annotation_table$Replicate_name[
          !annotation_table$Replicate_name %in% colnames(expression_table)] ),
        " do not occur in expression table."))
      stop(paste0("All Replicate_names mentioned in the annotation table",
        "must occur in the exprresion table."))
    }
    # Reduce expression table to columns contained in annotation table
    if( sum(!(colnames(expression_table) %in% annotation_table$Replicate_name)) > 0 ){
      print(paste0("WARNING: Omitting samples ",
        unlist(colnames(expression_table)[
          !(colnames(expression_table) %in% annotation_table$Replicate_name)]),
        " from expression table because they were not mentioned in annotation table."))
      print(paste0("This does not affect the algorithm.",
        "Adjust annotation table if you wish that these samples are included."))
      expression_table <- expression_table[,as.character(annotation_table$Replicate_name)]
    }
    # Reduce expression table to rows containing at least one non-zero count
    # with mean expression > 1
    #rowIdx_lowCounts <- apply(expression_table,1,function(x){all(x==0)})
    rowIdx_lowCounts <- apply(expression_table,1,function(x){mean(x) < 1})
    if(sum(rowIdx_lowCounts) > 0){
      #print(paste0("WARNING: ",sum(rowIdx_lowCounts), " out of ",
      #  dim(expression_table)[1]," genes had only zero counts in all considered samples."))
      print(paste0("WARNING: ",sum(rowIdx_lowCounts), " out of ",
        dim(expression_table)[1]," genes had a mean RNA count of < 1."))
      print("These genes are omitted in the analysis.")
      expression_table <- expression_table[!rowIdx_lowCounts,]
    }
    
    # DAVID to be deprecated
    # Shorten expression table
    if(TRUE){
      ind_toKeep <- 1:100
      print(paste0("Working on subset of data: ",length(ind_toKeep)," genes."))
      expression_table <- expression_table[ind_toKeep,]
    }
    
    ## 1. Prepare data
    print("1. Preparing data.")
    
    # Assign samples to replicates
    expression_tables <- list()
    replicates <- unique(annotation_table$Replicate)
    for(iRep in 1:length(replicates)){
      expression_tables[[iRep]] <- expression_table[,as.character( annotation_table[
        annotation_table$Replicate %in% replicates[iRep],]$Replicate_name) ]
      # Call columns after sample (not specific replicate)
      colnames(expression_tables[[iRep]]) <- as.character( 
        annotation_table[annotation_table$Replicate %in% replicates[iRep],]$Sample )
    }
    head(expression_tables[1])
    
    # Convert list of matrices into 3D array  
    expression_array=array(NA,c(dim(expression_tables[[1]]),length(expression_tables)))
    for (i in 1:length(expression_tables)){
      expression_array[,,i] <- as.matrix(expression_tables[[i]])    # make numeric
    }
    rownames(expression_array) <- rownames(expression_tables[[1]])
    colnames(expression_array) <- colnames(expression_tables[[1]])
    save(expression_array,file=file.path(getwd(),"ImpulseDE_expression_array.RData"))
    
    # 2. Run DESeq2
    print("2. Run DESeq2")
    dfCountData <- expression_table[,colnames(expression_table) %in% annotation_table$Replicate_name]
    colnames(dfCountData) <- annotation_table$Sample[match(colnames(dfCountData),annotation_table$Replicate_name)]
    dds <- DESeqDataSetFromMatrix(countData = dfCountData,
      colData = annotation_table,
      design = ~ Sample)
    ddsDESeqObject <- DESeq(dds, test = "LRT", full = ~ Sample, reduced = ~1)
    # Get gene-wise dispersion estimates
    # var = mean + alpha * mean^2, alpha is dispersion
    # DESeq2 dispersion is 1/dispersion used in negative binomial in R
    dds_dispersions <- 1/dispersions(ddsDESeqObject) 
    # DESeq results for comparison
    dds_resultsTable <- results(ddsDESeqObject)
    # Correct counts by size factors
    save(dds_dispersions,file=file.path(getwd(),"ImpulseDE_DEseq2_dispersions.RData"))
    save(dds_resultsTable,file=file.path(getwd(),"ImpulseDE_DEseq2_results.RData"))
    
    ## 2. Prepare annotation table for internal usage: chose samples from input
    print("2. Prepare annotation table for internal usage")
    print("-------------------------------------------------------------------")
    tm_annot <- system.time({
      prepared_annotation <- annotation_preparation(
        data_annotation=annotation_table,data_tables=expression_tables[[1]], 
        colname_time=colname_time, colname_condition=colname_condition, 
        control_timecourse=control_timecourse, control_name=control_name, 
        case_name=case_name)
    })
    prepared_annotation <- prepared_annotation[order(
      prepared_annotation$Condition),]
    prepared_annotation <- prepared_annotation[order(prepared_annotation$Time),]
    save(prepared_annotation, file=file.path(getwd(),"ImpulseDE_prepared_annotation.RData"))
    print("DONE")
    print("###################################################################")
    
    # Get expression values of target samples specifcied in prepared_annotation
    # DAVID unnecessary, done above
    #expression_array <- expression_array[,rownames(prepared_annotation),]
    
    # DAVID: Can my model handle this?
    # exclude genes with missing values(NAs)
    indx <- apply(expression_array,1,function(x){TRUE %in% is.na(x)})
    expression_array <- expression_array[!(indx),,]
    
    # if rownames are just 1,2,3 or if there are no rownames
    if(is.null(rownames(expression_array))){
      rownames(expression_array) <- paste("G", 1:nrow(expression_array),
        sep = "_")
    } else if(length(grep("[a-zA-Z]",rownames(expression_array))) == 0){
      rownames(expression_array) <- paste(rownames(expression_array),"G",
        sep = "_")
    }
    
    ###  3. Fit Impule model to each gene 
    print("3. Fitting Impulse model to the genes")
    print("-------------------------------------------------------------------")
    tm_imp_fit_gen <- system.time({
      impulse_fit_genes <- impulse_fit(data_input=expression_array, 
        data_annotation=prepared_annotation, 
        control_timecourse=control_timecourse, control_name=control_name, 
        n_proc = n_process,dispersion_vector=dds_dispersions,NPARAM=NPARAM)
    })
    save(impulse_fit_genes,file=file.path(getwd(),"ImpulseDE_model_fit.RData"))
    print("DONE")
    print(paste("Consumed time: ",round(tm_imp_fit_gen["elapsed"]/60,2),
      " min",sep=""))
    print("###################################################################")
    
    ### 5. Detect differentially expressed genes
    print("5. DE analysis")
    print("-------------------------------------------------------------------")
    tm_DE <- system.time({
      ImpulseDE_results <- DE_analysis(data_array=expression_array,
        data_annotation=prepared_annotation,impulse_fit_results=impulse_fit_genes,
        control_timecourse=control_timecourse,control_name=control_name,Q=Q_value,
        NPARAM=NPARAM,dispersion_vector=dds_dispersions)
      impulse_DE_genes <- as.character(as.vector( 
        ImpulseDE_results[as.numeric(ImpulseDE_results$adj.p) <= Q_value,"Gene"] ))
    })
    save(ImpulseDE_results,file=file.path(getwd(),"ImpulseDE_results.RData"))
    save(impulse_DE_genes,file=file.path(getwd(),"ImpulseDE_DE_genes.RData"))
    
    ### 6. Plot the top DE genes
    print("6. Plot top DE genes")
    plot_impulse(impulse_DE_genes,
      expression_array, prepared_annotation, impulse_fit_genes,
      control_timecourse, control_name, case_ind, file_name_part = "DE",
      title_line = "", sub_line = "",
      ImpulseDE_res=ImpulseDE_results,DESeq2_res=dds_resultsTable)
    print("DONE")
    print(paste("Consumed time: ",round(tm_DE["elapsed"]/60,2)," min",sep=""))
    print("##################################################################")
  })
  print(paste("TOTAL consumed time: ",round(tm_tot["elapsed"]/60,2),
    " min",sep=""))
  
  return(list(
    "ImpulseDE_DE_genes"=impulse_DE_genes,
    "ImpulseDE_results"=ImpulseDE_results,
    "ImpulseDE_impulse_fit_results"=impulse_fit_genes,
    "DESeq2_results"=dds_resultsTable))
}
