### Clear
rm(list = ls())

### LOAD DATA
### Load annotation
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/software_test")
# DAVID: take out timepoints
annotation_table <- read.table("annotation_table_RNAseq_D50-D51-D54.tab",header=T)
rownames(annotation_table) <- annotation_table$Replicate_name

### Load count data
# DE Genes (as of DESeq)
#setwd("/Users/davidsebastianfischer/MasterThesis/data/DESeq_RNAseq/DESeq_output")
#expression_table_raw <- read.table("D50-D51-D54_DEseq_RNAseq_filteredGenes_TPM_p1e-05.tab",sep="\t",header=T)
# All genes
setwd("/Users/davidsebastianfischer/MasterThesis/data/DESeq_RNAseq/raw_RNAseq")
expression_table_raw <- read.table("hDC_TPM_D50-D51-D54.tab",sep="\t",header=T)

### Process count data
expression_table <- expression_table_raw[,2:dim(expression_table_raw)[2]]
expression_table <- t(apply(expression_table,1,as.numeric))

rownames(expression_table) <- expression_table_raw[,1]
colnames(expression_table) <- colnames(expression_table_raw)[2:dim(expression_table_raw)[2]]

# replace by simulation data
# expression_table_cut[,] <- rnbinom(dim(expression_table_cut)[1]*dim(expression_table_cut)[2],
#  size=2,mu=100)

setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building")
source("ImpulseDE_main.R")
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")

colname_time = "Time"
colname_condition = "Condition"
control_timecourse = FALSE
control_name = NULL
case_name = "case"
expr_type = "Array"
n_iter = 100
n_process = 3
Q_value = 0.01
lsImpulseDE_results <- impulse_DE(expression_tables = expression_table, annotation_table = annotation_table,
  colname_time=colname_time, colname_condition=colname_condition, control_timecourse=control_timecourse,
  control_name = control_name, case_name =case_name, expr_type =expr_type,
  n_iter =n_iter, n_process =n_process, Q_value =Q_value)

impulse_DE_genes <- lsImpulseDE_results$ImpulseDE_DE_genes
ImpulseDE_results <- lsImpulseDE_results$ImpulseDE_results
impulse_fit_genes <- lsImpulseDE_results$ImpulseDE_impulse_fit_results
dds_resultsTable <- lsImpulseDE_results$DESeq2_results
