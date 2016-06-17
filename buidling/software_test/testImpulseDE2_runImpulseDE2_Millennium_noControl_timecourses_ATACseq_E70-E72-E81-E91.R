### Clear
rm(list = ls())

### LOAD DATA
### Load annotation
setwd( "/data/yosef2/users/fischerd/code/ImpulseDE2")
dfAnnotationFull <- read.table("annotation_table_ATACseq_E70-E72-E81-E91.tab",header=T)
rownames(dfAnnotationFull) <- dfAnnotationFull$Replicate

### Load count data
# All genes
setwd("/data/yosef2/users/fischerd/data/ImpulseDE2_input")
expression_table_raw <- read.table("E70-E72-E81-E91_peak_counts_DESeq_E72_24h_LPS_downsampled.tab",
  sep="\t",header=T,colClasses=c(
    "character","character","character",
    "numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric",
    "numeric","numeric","numeric","numeric","numeric"))

### Process count data
expression_table <- expression_table_raw[,4:dim(expression_table_raw)[2]]
expression_table <- t(apply(expression_table,1,as.numeric))
lsPeaksIDs <- apply(expression_table_raw[,1:3],1,function(coord){paste(coord,collapse="-")})
# Rename elements for shaked's poster: hide coordinates
lsPeaksIDs <- paste0(rep("Region_",length(lsPeaksIDs)),c(1:length(lsPeaksIDs)))

rownames(expression_table) <- lsPeaksIDs
colnames(expression_table) <- colnames(expression_table_raw)[4:dim(expression_table_raw)[2]]

# Only chose high counts
#expression_table <- expression_table[apply(expression_table,1,function(gene){max(gene,na.rm=TRUE)}>500),]

setwd( "/data/yosef2/users/fischerd/code/ImpulseDE2")
source("ImpulseDE2_main.R")

# Run ImpulseDE2
setwd( "/data/yosef2/users/fischerd/data/ImpulseDE2_output/timecourses_ATACseq_E70-E72-E81-E91")

control_timecourse = FALSE
strControlName = NULL
strCaseName = "case"
n_process = 16
Q_value = 10^(-3)
strMode <- "timecourses"
lsImpulseDE_results <- runImpulseDE2(matCountData=expression_table, dfAnnotationFull=dfAnnotationFull,
  strCaseName = strCaseName, strControlName=strControlName, strMode=strMode,
  nProc=n_process, Q_value=Q_value)