#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Annotation preparation    +++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# INPUT:
#   annotation_table: (Table samples x 2[time and condition]) providing 
#       co-variables for the samples including condition and time points.
#       Time points must be numeric numbers.
#   data_table: (Numeric matrix genes x samples)
#   colname_time: (str) column name of the co-variable "Time" in annotation_table
#   colanme_condition: (str) column name of the co-variable "Condition" in 
#       annotation_table
#   control_timecourse: (bool) [Default FALSE]control time timecourse is part of 
#       the data set (TRUE) or not (FALSE).
#   control_name: (str) name of the control condition in annotation_table.
#   case_name: (str) name of the case condition in annotation_table.
# OUTPUT:
#   annot: (Table samples x 2[time and condition]) annotation_table reduced to 
#       target samples 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Prepares the annotation table for internal use
annotation_preparation <- function(data_annotation, data_table = NULL,
                                   colname_time = NULL, colname_condition = NULL, control_timecourse = FALSE,
                                   control_name = NULL, case_name = NULL){
  
  # 1. Check Input is valid
  if(is.null(data_table)){
    stop("ERROR: table of data values must be specified")
  }
  if((nrow(data_annotation) != ncol(data_table)) ||
       (FALSE %in% (sort(rownames(data_annotation)) ==
                      sort(colnames(data_table))))){
    stop("ERROR: column names of data table must be the same as
         row names of annotation table!")
  }
  if(is.null(colname_time)){
    stop("ERROR: name of the column containing the timepoints must be
         specified as character!")
  }
  if(is.null(colname_condition)){
    stop("ERROR: name of the column containing the condition(s) must be
         specified as character!")
  }
  if(is.numeric(data_annotation[,colname_time]) == FALSE){
    stop("ERROR: time variable must be numeric and is not allowed to
         contain characters!")
  }
  if(control_timecourse == FALSE &
       length(summary(as.factor(data_annotation[,colname_condition]))) >= 2){
    stop("ERROR: if you have only one time course do not provide more than
         one condition!")
  }
  if(control_timecourse == TRUE & (length(summary(as.factor(data_annotation[,
                                                                            colname_condition]))) > 2) & is.null(case_name)){
    stop("ERROR: please specify case and control names of interest
         since you provide more than three conditions!")
  }
  
  # 2. Format data for output
  # Note on double assignment: annot is filled with whatever column name annotation input is present
  annot <- data_annotation[,c(colname_time,colname_condition)]
  annot <- data_annotation[colnames(data_table),]
  colnames(annot) <- c("Time","Condition")
  annot$Condition <- as.character(annot$Condition)
  
  if(length(summary(as.factor(data_annotation[,colname_condition]))) > 2){
    print(paste("Case condition: ", case_name, sep  = ""))
  } else {
    print(paste("Case condition: ",  data_annotation[!(data_annotation[,
                                                                       colname_condition] %in% control_name),colname_condition][1], sep = ""))
  }
  if(control_timecourse == TRUE){
    print(paste("Control condition: ", control_name, sep  = ""))
  }
  if(control_timecourse == TRUE & (length(summary(as.factor(data_annotation[,
                                                                            colname_condition]))) > 2)){
    annot <- annot[annot$Condition %in% c(control_name, case_name),]
  }
  return(annot)
  }