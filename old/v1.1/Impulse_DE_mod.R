
################################################################################
########################     Impulse DE package     ############################
################################################################################

### Version:  1.1
### Date:     2016
### Author v1.0:  Jil Sander
### Author v1.1:  David Fischer

################################################################################
# DAVID: General remarks not related to variance extension
################################################################################
# 1. Bit code compiled function two_impulse_comp is created but not used for squared error
#	and impulse prediction: ## compile simple functions to make them quicker
# 2. Impulse_fit_gene_wise what is the theta initialisation after
#     if(fit_to_clusters == TRUE & is.null(start_val)){
################################################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                             Functions to call                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



################################################################################

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Annotation preparation    +++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Prepares the annotation table for internal use
source("src_Impulse_annotation_preparation.R")

################################################################################

#' Impulse model value prediction
#'
#' Calculates impulse model values for given timepoints and predicted
#' impulse parameters.
#' @aliases calc_impulse
#' @param theta numerical vector of impulse parameters with the order
#' beta, h0, h1, h2, t1, t2.
#' @param timepoints numercial vector of time point(s).
#' @return The predicted impulse model values for the given time point(s).
#' @export
calc_impulse <- function(theta,timepoints){
  beta1 = theta[1]
  h0 = theta[2]
  h1 = theta[3]
  h2 = theta[4]
  t1 = theta[5]
  t2 = theta[6]
  res = NULL
   for (i in 1:length(timepoints)){
    res[i] = (1/h1) * (h0 + (h1-h0) * (1/(1+exp(-beta1*(timepoints[i]-t1))))) *
        (h2 + (h1-h2) * (1/(1+exp(beta1*(timepoints[i]-t2)))))
  }
  res = unlist(res)
#  names(res) <- timepoints
  return(res)
}

################################################################################

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++     Calculation of input for Optimization    ++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Cost function for parametric model fit: Ordinary least squares (two_impulses)
### or weighted least squares (two_impulses_WLS)
source("src_Impulse_CostFunctionFit.R")

################################################################################

## compile fitting simple functions to make them quicker
calc_impulse_comp <- cmpfun(calc_impulse)
two_impulses_comp <- cmpfun(two_impulses)


################################################################################

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++     Clustering    +++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Performs 2-step clustering on expression data or background data
source("src_Impulse_cluster_genes_for_impulse.R")

################################################################################

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Impulse model fit    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Fits impulse model to a timecourse dataset
source("src_Impulse_impulse_fit.R")

################################################################################

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Plot impulse fits    ++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### plots the impulse fits to timecourse data and also the control data if present
source("src_Impulse_plot_impulse.R")

#################################################################################

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Background generation   +++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### generates the background for the DE analysis
source("src_Impulse_generate_background.R")

################################################################################

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++     DE analysis   ++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# Detect differentially expressed genes over time
source("src_Impulse_DE_analysis.R")

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
#' @param expression_table numeric matrix of expression values; genes should
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
#' \item \code{DE_results} Matrix containing the names of genes being called
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
#' \code{DE_results} of the output value.
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
#   expression_table: (Numeric matrix genes x samples)
#   annotation_table: (Table) providing co-variables for the samples including 
#       condition and time points. Time points must be numeric numbers.
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
#   plot_clusters: (bool) [Default TRUE] plot the clusters
#   n_iter: (scalar) [Defaul 100] number of iterations, which are performed to 
#       fit the impulse model to the clusters.
#   n_randoms: (scalar) [Default 50000] number of generated randomized background 
#       iterations, which are used for differential expression analysis Don't decrease.
#   n_process: (scalar) [Default 4] number of processes, which can be used on the 
#       machine to run the background calculation in parallel
#   Q_value: (scalar) [default 0.01] cutoff to call genes significantly 
#       differentially expressed after FDR correction
# OUTPUT
#   impulse_fit_genes: (list) List containing fitted values and model parameters
#   impulse_DE_genes: names of genes being called as differentially expressed
#   ...more output saved to working directory, see documentation

impulse_DE <- function(expression_table = NULL, annotation_table = NULL,
    colname_time = NULL, colname_condition = NULL, control_timecourse = FALSE,
    control_name = NULL, case_name = NULL, expr_type = "Array",
    plot_clusters = TRUE, n_iter = 100, n_randoms = 50000, n_process = 4,
    Q_value = 0.01){

  tm_tot <- system.time({

    #' @import compiler
    #' @importFrom grDevices dev.off pdf
    #' @importFrom graphics abline axis legend par plot points
    #' @importFrom stats aov cor dist kmeans mad median nlminb optim p.adjust
    #' @importFrom stats runif sd
    #' @importFrom utils head write.table

    # 1. READ IN DATA
    ## prepare annotation table for internal usage: chose samples from input
    print("START: Prepare annotation table for internal usage")
    print("-------------------------------------------------------------------")
    tm_annot <- system.time({
      prepared_annotation <- annotation_preparation(annotation_table,
          expression_table,
          colname_time ,colname_condition, control_timecourse, control_name,
          case_name)
    })
    prepared_annotation <- prepared_annotation[order(
            prepared_annotation$Condition),]
    prepared_annotation <- prepared_annotation[order(prepared_annotation$Time),]
    save(prepared_annotation, file=file.path(getwd(),
            "/prepared_annotation.RData"))
    print("DONE")
    print("###################################################################")

    # Get expression values of target samples specifcied in prepared_annotation
    expression_table <- as.matrix(expression_table)    # to have numeric values
    expression_table <- expression_table[,rownames(prepared_annotation)]

    # DAVID: Can my model handle this?
    # exclude genes with missing values(NAs)
    indx <- apply(expression_table,1,function(x){TRUE %in% is.na(x)})
    expression_table <- expression_table[!(indx),]

    # if rownames are just 1,2,3 or if there are no rownames
    if(is.null(rownames(expression_table))){
       rownames(expression_table) <- paste("G", 1:nrow(expression_table),
            sep = "_")
    } else if(length(grep("[a-zA-Z]",rownames(expression_table))) == 0){
       rownames(expression_table) <- paste(rownames(expression_table),"G",
            sep = "_")
    }

    ### 2. cluster genes to reduce efforts for fitting the impulse model
    print("START: Clustering genes for Impulse model fit")
    print("-------------------------------------------------------------------")
    load(file.path(getwd(),"/prepared_annotation.RData"))
    tm_clust <- system.time({
      clustering_results <- cluster_genes_for_impulse(expression_table,
          prepared_annotation, control_timecourse, control_name, plot_clusters)
    })
    save(clustering_results, file=file.path(getwd(),"clustering_results.RData"))
    print("DONE")
    print(paste("Consumed time: ",round(tm_clust["elapsed"]/60,2)," min",sep=""))
    print("##################################################################")

    # 3. fit Impulse model to the clusters
    print("START: Fitting Impulse model to the clusters")
    print("-------------------------------------------------------------------")
    load(file.path(getwd(),"clustering_results.RData"))
    tm_imp_fit_clus <- system.time({
      impulse_fit_clusters <- impulse_fit(clustering_results,prepared_annotation,
          n_iter, control_timecourse, control_name, n_proc = n_process)
    })
    save(impulse_fit_clusters,file=file.path(getwd(),
        "impulse_fit_clusters.RData"))
    print("DONE")
    print(paste("Consumed time: ",round(tm_imp_fit_clus["elapsed"]/60,2)," min",
        sep=""))
    print("###################################################################")

    ###  4. fit Impule model to each gene by using the cluster fits as start values
    print("START: Fitting Impulse model to the genes")
    print("-------------------------------------------------------------------")
    load(file.path(getwd(),"impulse_fit_clusters.RData"))
    tm_imp_fit_gen <- system.time({
      impulse_fit_genes <- impulse_fit(expression_table, prepared_annotation,
        n_iter, control_timecourse, control_name, clustering_results,
        impulse_fit_clusters, n_proc = n_process)
    })
    save(impulse_fit_genes,file=file.path(getwd(),"impulse_fit_genes.RData"))
    print("DONE")
    print(paste("Consumed time: ",round(tm_imp_fit_gen["elapsed"]/60,2),
        " min",sep=""))
    print("###################################################################")

    ### 5. generate background for the DE analysis
    print("START: Generate background")
    print("-------------------------------------------------------------------")
    load(file.path(getwd(),"impulse_fit_genes.RData"))
    tm_bg <- system.time({
      background_results <-  generate_background(expression_table,
        prepared_annotation, n_iter, impulse_fit_genes, control_timecourse,
        control_name, clustering_results, n_randoms, n_process)
    })
    save(background_results,file=file.path(getwd(),"background_results.RData"))
    print("DONE")
    print(paste("Consumed time: ",round(tm_bg["elapsed"]/60,2)," min",sep=""))
    print("###################################################################")

    ## 6. detect differentially expressed genes
    print("START: DE analysis")
    print("-------------------------------------------------------------------")
    load(file.path(getwd(),"background_results.RData"))
    tm_DE <- system.time({
      impulse_DE_genes <- DE_analysis(expression_table,prepared_annotation,
            impulse_fit_genes, background_results, control_timecourse,
            control_name, expr_type, Q_value)
    })
   
    ## 7. plot the top DE genes
    if(control_timecourse == TRUE){
      case_ind = as.character(prepared_annotation$Condition[prepared_annotation$
            Condition != control_name][1])
    } else { case_ind = NULL}
    if(is.null(impulse_DE_genes) == FALSE & length(impulse_DE_genes) > 1 ){
      if(is.list(impulse_DE_genes) == TRUE &
         is.data.frame(impulse_DE_genes) == FALSE){
            plot_impulse(as.character(impulse_DE_genes[[1]][,
              "Gene"])[1:(min(nrow(impulse_DE_genes[[1]]),2*18))],
              expression_table, prepared_annotation,impulse_fit_genes,
              control_timecourse, control_name, case_ind, file_name_part = "DE",
              title_line = "", sub_line = "")
        } else {
            plot_impulse(as.character(impulse_DE_genes[,
              "Gene"])[1:(min(nrow(impulse_DE_genes),2*18))],
              expression_table, prepared_annotation, impulse_fit_genes,
              control_timecourse, control_name, case_ind, file_name_part = "DE",
              title_line = "", sub_line = "")
      }
    }
    save(impulse_DE_genes,file=file.path(getwd(),"impulse_DE_genes.RData"))
    print("DONE")
    print(paste("Consumed time: ",round(tm_DE["elapsed"]/60,2)," min",sep=""))
    print("##################################################################")
  })
  print(paste("TOTAL consumed time: ",round(tm_tot["elapsed"]/60,2),
        " min",sep=""))
  
  return(list("impulse_fit_results"=impulse_fit_genes,
        "DE_results" = impulse_DE_genes))
}
