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

setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE")
source("srcImpulseDE_annotation_preparation.R")
source("srcImpulseDE_calc_impulse.R")
source("srcImpulseDE_CostFunctionFit.R")
calc_impulse_comp <- cmpfun(calc_impulse)
two_impulses_OLS_comp <- cmpfun(two_impulses_OLS)
two_impulses_WLS_comp <- cmpfun(two_impulses_WLS)
source("srcImpulseDE_cluster_genes_for_impulse.R")
source("srcImpulseDE_impulse_fit.R")
source("srcImpulseDE_generate_background.R")
source("srcImpulseDE_DE_analysis.R")
source("srcImpulseDE_plot_impulse.R")

################################################################################
### Create data
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test")
annotation_table <- read.table("annotation_table_smallArtificialData.tab",header=T)
rownames(annotation_table) <- annotation_table$Sample

colname_time <- "Time"
colname_condition <- "Condition"
control_name <- NULL
case_name <- "A"
expression_table1 <- read.table("expression_table_1_1.tab",header=T)
rownames(expression_table1) <- expression_table1[,1]
expression_table1 <- expression_table1[,2:dim(expression_table1)[2]]
expression_table2 <- read.table("expression_table_1_2.tab",header=T)
rownames(expression_table2) <- expression_table2[,1]
expression_table2 <- expression_table2[,2:dim(expression_table2)[2]]
expression_table3 <- read.table("expression_table_1_3.tab",header=T)
rownames(expression_table3) <- expression_table3[,1]
expression_table3 <- expression_table3[,2:dim(expression_table3)[2]]

expression_tables = list(expression_table1,expression_table2,expression_table3)

expression_array=array(NA,c(dim(expression_tables[[1]]),length(expression_tables)))
for (i in 1:length(expression_tables)){
  expression_array[,,i] <- as.matrix(expression_tables[[i]])    # make numeric
}
rownames(expression_array) <- rownames(expression_tables[[1]])
colnames(expression_array) <- colnames(expression_tables[[1]])
control_timecourse = FALSE

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
sigma_matrix <- apply(expression_array,c(1,2),sd)

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

plot_clusters = FALSE
clustering_results <- cluster_genes_for_impulse(
  apply(expression_array,c(1,2),mean),
  prepared_annotation, control_timecourse, control_name, plot_clusters)

################################################################################
### 3. fit Impulse model to the clusters
n_process <- 2
n_iter <- 100
impulse_fit_clusters <- impulse_fit(clustering_results,prepared_annotation,
  sigma_matrix,n_iter, control_timecourse, control_name, n_proc = n_process)

################################################################################
###  4. Fit Impule model to each gene by using the cluster fits as start values
impulse_fit_genes <- impulse_fit(expression_array, prepared_annotation,
  sigma_matrix,n_iter, control_timecourse, control_name, clustering_results,
  impulse_fit_clusters, n_proc = n_process)

################################################################################
### 5. Generate background for the DE analysis
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE")
source("srcImpulseDE_annotation_preparation.R")
source("srcImpulseDE_calc_impulse.R")
source("srcImpulseDE_CostFunctionFit.R")
calc_impulse_comp <- cmpfun(calc_impulse)
two_impulses_OLS_comp <- cmpfun(two_impulses_OLS)
two_impulses_WLS_comp <- cmpfun(two_impulses_WLS)
source("srcImpulseDE_cluster_genes_for_impulse.R")
source("srcImpulseDE_impulse_fit.R")
source("srcImpulseDE_generate_background.R")
source("srcImpulseDE_DE_analysis.R")
source("srcImpulseDE_plot_impulse.R")
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test")
n_randoms <- 20
n_process <- 2
background_results <-  generate_background(expression_array,
  prepared_annotation, sigma_matrix, n_iter, impulse_fit_genes, 
  control_timecourse, control_name, clustering_results, n_randoms, n_process)

### 6. Detect differentially expressed genes
Q_value = 0.01
expr_type = "Array"
impulse_DE_genes <- DE_analysis(expression_array,prepared_annotation,
  sigma_matrix,impulse_fit_genes, background_results, control_timecourse,
  control_name, expr_type, Q_value)
