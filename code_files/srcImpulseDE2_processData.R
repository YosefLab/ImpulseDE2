#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Annotation preparation    +++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Process annotation and count data
#' 
#' Check validity of input and process count data matrix and annotation
#' into data structures used later in \code{runImpulseDE2}.
#' \code{processData} is structure in the following way:
#' (I) Helper functions:
#'   checkData()
#'   nameGenes()
#'   reduceCountData()
#' (II) Script body
#' 
#' @seealso Called by \code{runImpulseDE2}.
#' 
#' @param matCountData (matrix genes x replicates) [Default NULL] Count data of all conditions, 
#'    unobserved entries are NA. Column labels are replicate names, row labels
#'    gene names.
#' @param dfAnnotationFull (Table) [Default NULL] Lists co-variables of replicates: 
#'    Replicate, Sample, Condition, Time. Time must be numeric.
#' @param strCaseName (str) [Default NULL] Name of the case condition in \code{dfAnnotationFull}.
#' @param strControlName: (str) [Default NULL] Name of the control condition in 
#'    \code{dfAnnotationFull}.
#'    
#' @return arr2DCountData: (2D array genes x replicates) Count data: Reduced 
#'    version of \code{matCountData}. For internal use.
#'
#' @export

processData <- function(dfAnnotationFull=NULL, matCountData=NULL,
  strCaseName=NULL, strControlName=NULL, strMode=NULL){
  
  ###############################################################
  # (I) Helper functions
  
  # Check format and presence of input data
  checkData <- function(dfAnnotationFull=NULL, arr2DCountData=NULL,
    strCaseName=NULL, strControlName=NULL, strMode=NULL){
    
    ### 1. Check that all necessary input was specified
    # dfAnnotationFull
    if(is.null(dfAnnotationFull)){
      stop("ERROR: Annotation table was not given as input: dfAnnotationFull")
    }
    # arr2DCountData
    if(is.null(arr2DCountData)){
      stop("ERROR: List of count tables was not given as input: arr2DCountData")
    }
    # strCaseName
    if(is.null(strCaseName)){
      stop("ERROR: No name for case condition given.")
    }
    if(is.null(strMode)){
      stop("ERROR: ImpulseDE2 mode (strMode) was not given as input.")
    }
    
    ### 2. Check annotation table content
    ### a) Check column names
    vecColNamesRequired <- c("Replicate","Sample","Condition","Time")
    if( !all(vecColNamesRequired %in% colnames(dfAnnotationFull)) ){
      stop(paste0("Could not find column ",
        vecColNamesRequired[!(vecColNamesRequired %in% colnames(dfAnnotationFull))],
        " in annotation table."))
    }
    ### b) Replicates
    # Check that replicate name do not occur twice
    if(length(unique(dfAnnotationFull$Replicate)) != dim(dfAnnotationFull)[1]){
      stop(paste0("ERROR: [Annotation table]",
        "Number of Replicate instances different to annotation table number of rows.",
        "Replicate cannot occur multiple times."))
    }
    ### c) Time points
    # Check that there is not multiple samples per time point in one condition
    vecTimepoints <- unique(as.vector( dfAnnotationFull$Time ))
    for(tp in vecTimepoints){
      if( length(unique(as.vector( dfAnnotationFull[dfAnnotationFull$Time %in% tp &
          dfAnnotationFull$Sample %in% strCaseName,]$Sample )))>1 ){
        stop(paste0("ERROR: In case condition, assigned multiple samples to time point ",tp))
      }
    }
    # Check that time points are numeric
    if(is.numeric(dfAnnotationFull$Time) == FALSE){
      stop("ERROR: Time variable in annotation table must be numeric.")
    }
    ### d) Conditions
    lsConditions <- unique(as.vector( dfAnnotationFull$Condition ))
    # Check that given conditions exisit in annotation table
    if(!(strCaseName %in% lsConditions)){
      stop("ERROR: Condition given for case does not occur in annotation table condition column.")
    }
    if(!is.null(strControlName)){
      if(!(strControlName %in% lsConditions)){
        stop("ERROR: Condition given for control does not occur in annotation table condition column.")
      }
    }
    
    ### 3. Check expression table content
    # Check that all entries in count data table occur in annotation table
    if( any(!(colnames(arr2DCountData) %in% dfAnnotationFull$Replicate)) ){
      print(paste0("WARNING: The column(s) ",
        paste0(as.character( colnames(arr2DCountData)[
          !(colnames(arr2DCountData) %in% dfAnnotationFull$Replicate)] ),collapse=","),
        " in the count data table do(es) not occur in annotation table and will be ignored."))
    }
    
    ### 4. Check mode
    lsAllowedModes <- c("batch", "timecourses", "singlecell")
    if(sum(strMode==lsAllowedModes) != 1){
      stop(paste0( "ERROR: ImpulseDE2 mode given as input, strMode=", strMode,
        ", is not recognised. Chose from {", paste0(lsAllowedModes, collapse=","), "}." ))
    }
    
    ### Summarise which conditions, samples, lsReplicates were found and mode
    print(paste0("Found conditions: ",paste0(lsConditions,collapse=",")))
    print(paste0("Case condition: ", strCaseName))
    if(!is.null(strControlName)){
      print(paste0("Control condition: ", strControlName))
    }
    print(paste0( "Found time points: ",
      paste( vecTimepoints,collapse=",") ))
    for(tp in vecTimepoints){
      print(paste0( "Case: Found the following replicates for sample ",
          unique(dfAnnotationFull[
            dfAnnotationFull$Time %in% tp &
              dfAnnotationFull$Condition %in% strCaseName &
              dfAnnotationFull$Replicate %in% colnames(arr2DCountData),
            ]$Sample),
        " at time point ", tp,
        ": ", paste0(dfAnnotationFull[
          (dfAnnotationFull$Time %in% tp) &
            (dfAnnotationFull$Condition %in% strCaseName) &
            (dfAnnotationFull$Replicate %in% colnames(arr2DCountData)),
          ]$Replicate,collapse=","),collapse="," ))
    }
    if(!is.null(strControlName)){
      for(tp in vecTimepoints){
        print(paste0( "Control: Found the following replicates for sample ", 
          unique(dfAnnotationFull[
            dfAnnotationFull$Time %in% tp &
              dfAnnotationFull$Condition %in% strControlName &
              dfAnnotationFull$Replicate %in% colnames(arr2DCountData),
            ]$Sample),
          " at time point ", tp,
          ":", paste(dfAnnotationFull[
            dfAnnotationFull$Time %in% tp &
              dfAnnotationFull$Condition %in% strControlName &
              dfAnnotationFull$Replicate %in% colnames(arr2DCountData),
            ]$Replicate,collapse=",") ))
      }
    }
    print(paste0( "ImpulseDE2 runs in mode: ", strMode ))
    
    return(NULL)
  }
  
  # Name genes if names are not given
  nameGenes <- function(arr2DCountData=NULL){
    if(is.null(rownames(arr2DCountData))){
      rownames(arr2DCountData) <- paste("G", 1:nrow(arr2DCountData),
        sep = "_")
    }
    
    return(arr2DCountData)
  }
  
  # Reduce count data to data which are utilised later
  reduceCountData <- function(dfAnnotationFull=NULL, arr2DCountData=NULL){
    
    ### 1. Columns (Conditions):
    # Reduce expression table to columns of considered conditions
    if(!is.null(strControlName)){
     vecReplicateNames_Case <- as.character(as.vector( 
        dfAnnotationFull[dfAnnotationFull$Condition %in% strCaseName,]$Replicate ))
     vecReplicateNames_Ctrl <- as.character(as.vector( 
        dfAnnotationFull[dfAnnotationFull$Condition %in% strControlName,]$Replicate ))
      arr2DCountData <- arr2DCountData[,c(vecReplicateNames_Case,vecReplicateNames_Ctrl)]
    } else {
     vecReplicateNames_Case <- as.character(as.vector( 
        dfAnnotationFull[dfAnnotationFull$Condition %in% strCaseName,]$Replicate ))
      arr2DCountData <- arr2DCountData[,vecReplicateNames_Case]
    }
    # Check that every replicate contains at least one observed value (not NA)
    lsNARep <- any(apply(arr2DCountData,2,function(rep){all(is.na(rep))}))
    if(any(lsNARep)){
      print(paste0("WARNING: Replicate(s) ",paste0(colnames(arr2DCountData[lsNARep,]),collapse=","),
        " only contain(s) NA values and will be removed from the analysis."))
      arr2DCountData <- arr2DCountData[,!lsNARep]
    }
    ### 2. Rows (Genes):
    # Reduce expression table to rows containing at least one non-zero count
    rowIdx_lowCounts <- apply(arr2DCountData,1,function(x){all(x==0)})
    if(sum(rowIdx_lowCounts) > 0){
      print(paste0("WARNING: ",sum(rowIdx_lowCounts), " out of ",
        dim(arr2DCountData)[1]," genes had only zero counts in all considered samples."))
      print("These genes are omitted in the analysis.")
      arr2DCountData <- arr2DCountData[!rowIdx_lowCounts,]
    }
    # DAVID to be deprecated
    # Shorten expression table
    if(TRUE){
      ind_toKeep <- 100
      print(paste0("Working on subset of data: ",min(ind_toKeep,dim(arr2DCountData)[1])," genes."))
      arr2DCountData <- arr2DCountData[1:min(ind_toKeep,dim(arr2DCountData)[1]),]
    }
    # Exclude genes with only missing values (NAs)
    indx_NA <- apply(arr2DCountData,1,function(x){all(is.na(x))})
    if(any(indx_NA)){
      print(paste0("WARNING: Excluded ",sum(indx_NA)," genes because no real valued samples were given."))
    }
    arr2DCountData <- arr2DCountData[!(indx_NA),]
    
    # Sort 2D array by annotation table
    vecTimepoints <- as.character(sort(unique( as.numeric(as.vector(dfAnnotationFull$Time)) )))
    arr2DCountData <- arr2DCountData[,match(colnames(arr2DCountData),as.vector(dfAnnotationFull$Replicate))]
    
    return(arr2DCountData)
  }
    
  ###############################################################
  # (II) Main body of function
  
  checkData(
    dfAnnotationFull=dfAnnotationFull,
    arr2DCountData=matCountData,
    strControlName=strControlName,
    strCaseName=strCaseName,
    strMode=strMode)
  
  arr2DCountData <- nameGenes(arr2DCountData=matCountData)
  
  arr2DCountData <- reduceCountData(
    dfAnnotationFull=dfAnnotationFull, 
    arr2DCountData=arr2DCountData)
  
  return( arr2DCountData )
}