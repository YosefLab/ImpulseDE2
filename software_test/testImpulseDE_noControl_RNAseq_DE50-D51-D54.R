################################################################################
########################     Impulse DE Test     ###############################
################################################################################

### Test for Version:  1.2
### Date:     2016
### Author:  David Fischer

################################################################################
### Syntax test of functions
################################################################################
### Small artifical data set to validate functionality of code inside 
### individual function
################################################################################
### Clear
rm(list = ls())

### Load functions

library(compiler)
library(amap)
library(longitudinal)
library(parallel)
library(graphics)
library(MASS)

setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building")
source("srcImpulseDE_annotation_preparation.R")
source("srcImpulseDE_calc_impulse.R")
source("srcImpulseDE_CostFunctionFit.R")
calc_impulse_comp <- cmpfun(calc_impulse)
cost_fun_OLS_comp <- cmpfun(cost_fun_OLS)
cost_fun_WLS_comp <- cmpfun(cost_fun_WLS)
cost_fun_logl_comp <- cmpfun(cost_fun_logl)
source("srcImpulseDE_cluster_genes_for_impulse.R")
source("srcImpulseDE_impulse_fit.R")
source("srcImpulseDE_generate_background.R")
source("srcImpulseDE_DE_analysis.R")
source("srcImpulseDE_plot_impulse.R")
source("srcImpulseDE_compute_stdvs.R")
source("srcImpulseDE_compute_weights.R")

################################################################################
### LOAD DATA
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")

# DAVID: take out timepoints
annotation_table <- read.table("annotation_table_RNAseq_D50-D51-D54.tab",header=T)
#annotation_table <- read.table("annotation_table_RNAseq_D50-D51-D54_no6h.tab",header=T)
rownames(annotation_table) <- annotation_table$Sample

colname_time <- "Time"
colname_condition <- "Condition"
control_name <- NULL
case_name <- "A"
# DE Genes (as of DESeq)
setwd("/Users/davidsebastianfischer/MasterThesis/data/DESeq_RNAseq/DESeq_output")
expression_table_raw <- read.table("D50-D51-D54_DEseq_RNAseq_filteredGenes_TPM_p1e-05.tab",sep="\t",header=T)
# All genes
#setwd("/Users/davidsebastianfischer/MasterThesis/data/DESeq_RNAseq/raw_RNAseq")
#expression_table_raw <- read.table("hDC_TPM_D50-D51-D54.tab",sep="\t",header=T)

## short version
if(TRUE){
  expression_table_raw <- expression_table_raw[1:8,]
  const_row1 <- as.data.frame(t(c("const1",rep(100,dim(expression_table_raw)[2]-1))))
  const_row1[3] <- 102
  const_row1[20] <- 101
  const_row2 <- as.data.frame(t(c("const2",rep(100,dim(expression_table_raw)[2]-1))))
  const_rows <- rbind(const_row1,const_row2)
  colnames(const_rows) <- colnames(expression_table_raw)
  rownames(const_rows) <- c("21","22")
  expression_table_raw<- rbind(expression_table_raw,const_rows)
}

expression_table <- expression_table_raw[,2:dim(expression_table_raw)[2]]
expression_table <- t(apply(expression_table,1,as.numeric))
#expression_table <- log(expression_table)/log(2)
expression_table[expression_table<0] <- 0

rownames(expression_table) <- expression_table_raw[,1]
colnames(expression_table) <- colnames(expression_table_raw)[2:dim(expression_table_raw)[2]]

# Skip 1h only present in 51 and 54 for now - missing data for later
expression_table1 <- expression_table[,c("D50_0h","D50_2h","D50_4h","D50_6h","D50_12h","D50_24h")]
colnames(expression_table1) <- c("D_0h","D_2h","D_4h","D_6h","D_12h","D_24h")
expression_table2 <- expression_table[,c("D51_0h","D51_2h","D51_4h","D51_6h","D51_12h","D51_24h")]
colnames(expression_table2) <- c("D_0h","D_2h","D_4h","D_6h","D_12h","D_24h")
expression_table3 <- expression_table[,c("D54_0h","D54_2h","D54_4h","D54_6h","D54_12h","D54_24h")]
colnames(expression_table3) <- c("D_0h","D_2h","D_4h","D_6h","D_12h","D_24h")

# DAVID take out time point
#expression_table1 <- expression_table1[,c("D_0h","D_2h","D_4h","D_12h","D_24h")]
#expression_table2 <- expression_table2[,c("D_0h","D_2h","D_4h","D_12h","D_24h")]
#expression_table3 <- expression_table3[,c("D_0h","D_2h","D_4h","D_12h","D_24h")]

expression_tables = list(expression_table1,expression_table2,expression_table3)
#expression_tables = list(expression_table1,expression_table1)

expression_array=array(NA,c(dim(expression_tables[[1]]),length(expression_tables)))
for (i in 1:length(expression_tables)){
  expression_array[,,i] <- as.numeric( as.matrix(expression_tables[[i]]) )   # make numeric
}
rownames(expression_array) <- rownames(expression_tables[[1]])
colnames(expression_array) <- colnames(expression_tables[[1]])
control_timecourse = FALSE

# Minimal expression required in a row
min_expr_value <- 5
expression_array <- expression_array[apply(expression_array,1,function(x){any(x>min_expr_value)}),,]

setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")

################################################################################
### 1. Annotation

prepared_annotation <- annotation_preparation(annotation_table,
  expression_tables[[1]],colname_time ,colname_condition, 
  control_timecourse, control_name, case_name)
prepared_annotation

prepared_annotation <- prepared_annotation[order(
  prepared_annotation$Condition),]
prepared_annotation <- prepared_annotation[order(prepared_annotation$Time),]

# Get expression values of target samples specifcied in prepared_annotation
expression_array <- expression_array[,rownames(prepared_annotation),]

# DAVID: Can my model handle this?
# exclude genes with missing values(NAs)
indx <- apply(expression_array,1,function(x){TRUE %in% is.na(x)})
expression_array <- expression_array[!(indx),,]

# Precompute standard deviations over replicates
### DAVID: smoothen out estimates, get rid of 0s
sigma_matrix <- compute_stdvs(expression_array)
weight_matrix <- compute_weights(expression_array)
mean_matrix <- apply(expression_array,c(1,2),mean)
#sigma_matrix[mean_matrix < 5] <- min(sigma_matrix[mean_matrix >= 5])
#weight_matrix[mean_matrix < 5] <- min(weight_matrix[mean_matrix >= 5])
# Dont do weighting
#weight_matrix[] <- array(1,dim(weight_matrix))

# if rownames are just 1,2,3 or if there are no rownames
if(is.null(rownames(expression_array))){
  rownames(expression_array) <- paste("G", 1:nrow(expression_array),
    sep = "_")
} else if(length(grep("[a-zA-Z]",rownames(expression_array))) == 0){
  rownames(expression_array) <- paste(rownames(expression_array),"G",
    sep = "_")
}

################################################################################
### 2. cluster genes to reduce efforts for fitting the impulse model

#plot_clusters = FALSE
#clustering_results <- cluster_genes_for_impulse(
#  apply(expression_array,c(1,2),mean),
#  prepared_annotation, control_timecourse, control_name, plot_clusters)
#save(clustering_results,file=file.path(getwd(),"clustering_results_RNAseq_D50-D51-D54.RData"))

################################################################################
### 3. fit Impulse model to the clusters

#n_process <- 3
#n_iterations <- 100
#impulse_fit_clusters <- impulse_fit(data_input=clustering_results,
#  data_annotation=prepared_annotation,weight_mat=weight_matrix,n_iter=n_iterations,
#  control_timecourse=control_timecourse, control_name=control_name, n_proc=n_process)
#save(impulse_fit_clusters,file=file.path(getwd(),"impulse_fit_clusters_RNAseq_D50-D51-D54.RData"))

################################################################################
###  4. Fit Impule model to each gene by using the cluster fits as start values
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building")
source("srcImpulseDE_annotation_preparation.R")
source("srcImpulseDE_calc_impulse.R")
source("srcImpulseDE_CostFunctionFit.R")
calc_impulse_comp <- cmpfun(calc_impulse)
cost_fun_OLS_comp <- cmpfun(cost_fun_OLS)
cost_fun_WLS_comp <- cmpfun(cost_fun_WLS)
source("srcImpulseDE_cluster_genes_for_impulse.R")
source("srcImpulseDE_impulse_fit.R")
source("srcImpulseDE_generate_background.R")
source("srcImpulseDE_DE_analysis.R")
source("srcImpulseDE_plot_impulse.R")

setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
n_process <- 3
n_iterations <- 100
#impulse_fit_genes <- impulse_fit(data_input=expression_array, 
#  data_annotation=prepared_annotation,weight_mat=weight_matrix,n_iter=n_iterations, 
#  control_timecourse=control_timecourse, control_name=control_name,
#  cluster_results=clustering_results, start_values=impulse_fit_clusters, 
#  n_proc = n_process)
impulse_fit_genes <- impulse_fit(data_input=expression_array, 
  data_annotation=prepared_annotation,weight_mat=weight_matrix,n_iter=n_iterations, 
  control_timecourse=control_timecourse, control_name=control_name, 
  n_proc = n_process)
save(impulse_fit_genes,file=file.path(getwd(),"impulse_fit_genes_RNAseq_D50-D51-D54.RData"))

################################################################################
### 5. Generate background for the DE analysis
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building")
source("srcImpulseDE_annotation_preparation.R")
source("srcImpulseDE_calc_impulse.R")
source("srcImpulseDE_CostFunctionFit.R")
calc_impulse_comp <- cmpfun(calc_impulse)
cost_fun_OLS_comp <- cmpfun(cost_fun_OLS)
cost_fun_WLS_comp <- cmpfun(cost_fun_WLS)
source("srcImpulseDE_cluster_genes_for_impulse.R")
source("srcImpulseDE_impulse_fit.R")
source("srcImpulseDE_generate_background.R")
source("srcImpulseDE_DE_analysis.R")
source("srcImpulseDE_plot_impulse.R")
source("srcImpulseDE_compute_stdvs.R")
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
n_randoms <- 60
n_process <- 3
n_iterations <- 100
#background_results <-  generate_background(expression_array,
#  prepared_annotation, weight_matrix, n_iterations, impulse_fit_genes, 
#  control_timecourse, control_name, clustering_results, n_randoms, n_process)
background_results <-  generate_background(expression_array,
  prepared_annotation, weight_matrix, n_iterations, impulse_fit_genes, 
  control_timecourse, control_name, n_randoms, n_process)
save(background_results,file=file.path(getwd(),"background_results_RNAseq_D50-D51-D54.RData"))

### 6. Detect differentially expressed genes
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building")
source("srcImpulseDE_annotation_preparation.R")
source("srcImpulseDE_calc_impulse.R")
source("srcImpulseDE_CostFunctionFit.R")
calc_impulse_comp <- cmpfun(calc_impulse)
cost_fun_OLS_comp <- cmpfun(cost_fun_OLS)
cost_fun_WLS_comp <- cmpfun(cost_fun_WLS)
source("srcImpulseDE_cluster_genes_for_impulse.R")
source("srcImpulseDE_impulse_fit.R")
source("srcImpulseDE_generate_background.R")
source("srcImpulseDE_DE_analysis.R")
source("srcImpulseDE_plot_impulse.R")
source("srcImpulseDE_compute_stdvs.R")
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
Q_value = 0.01
expr_type = "Array"
DE_results <- DE_analysis(expression_array,prepared_annotation,
  weight_matrix,impulse_fit_genes, background_results, control_timecourse,
  control_name, expr_type, Q_value)
impulse_DE_genes <- DE_results[[1]]
Flike_distribution <- DE_results[[2]]
save(impulse_DE_genes,file=file.path(getwd(),"impulse_DE_genes_RNAseq_D50-D51-D54.RData"))

### 7. Plot the top DE genes
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building")
source("srcImpulseDE_annotation_preparation.R")
source("srcImpulseDE_calc_impulse.R")
source("srcImpulseDE_CostFunctionFit.R")
calc_impulse_comp <- cmpfun(calc_impulse)
cost_fun_OLS_comp <- cmpfun(cost_fun_OLS)
cost_fun_WLS_comp <- cmpfun(cost_fun_WLS)
source("srcImpulseDE_cluster_genes_for_impulse.R")
source("srcImpulseDE_impulse_fit.R")
source("srcImpulseDE_generate_background.R")
source("srcImpulseDE_DE_analysis.R")
source("srcImpulseDE_plot_impulse.R")
source("srcImpulseDE_compute_stdvs.R")
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
#lpsGenes <- c("IRF1","IRF2","IRF3","IRF4","IRF7","IRF8","IRF9","STAT1","STAT2","RELA","RELB")
genesToPlot <- as.character(impulse_DE_genes[,
    "Gene"])[1:(min(nrow(impulse_DE_genes),36))]
genesToPlot <- rownames(expression_array[1:min(200,dim(expression_array)[1]),,])

plot_impulse(genesToPlot,
  expression_array, prepared_annotation, impulse_fit_genes,
  control_timecourse, control_name, case_ind, file_name_part = "DE",
  title_line = "", sub_line = "")

#impulse_fit_genes$impulse_parameters_case[genesToPlot,]
#impulse_fit_genes$impulse_fits_case[genesToPlot,]

