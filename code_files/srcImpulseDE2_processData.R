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
#'   reshapeCountData()
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
#' @return (list length 3) with the following elements:
#' \itemize{
#'    \item \code{arr2DCountData} (2D array genes x replicates) Count data: Reduced 
#'        version of \code{matCountData}. For internal use.
#'    \item \code{arr3DCountData} (3D array genes x samples x replicates) Count data: 
#'        \code{arr2DCountData} reshaped into a 3D array. For internal use.
#'    \item \code{dfAnnotationRed} (data frame) Reduced version of 
#'        \code{dfAnnotationFull}. Lists co-variables of samples: 
#'        Sample, Condition, Time. Time must be numeric. For internal use.
#' }
#' @export

processData <- function(dfAnnotationFull=NULL, matCountData=NULL,
  strCaseName=NULL, strControlName=NULL){
  
  ###############################################################
  # (I) Helper functions
  
  # Check format and presence of input data
  checkData <- function(dfAnnotationFull=NULL, arr2DCountData=NULL,
    strCaseName=NULL, strControlName=NULL){
    
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
    
    ### 2. Check annotation table content
    ### a) Check column names
    lsColNamesRequired <- c("Replicate","Sample","Condition","Time")
    if( !all(lsColNamesRequired %in% colnames(dfAnnotationFull)) ){
      stop(paste0("Could not find column ",
        lsColNamesRequired[!(lsColNamesRequired %in% colnames(dfAnnotationFull))],
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
    lsTimepoints <- unique(as.vector( dfAnnotationFull$Time ))
    for(tp in lsTimepoints){
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
    
    ### Summarise which conditions, samples, lsReplicates were found
    print(paste0("Found conditions: ",lsConditions))
    print(paste0("Case condition: ", strCaseName))
    if(!is.null(strControlName)){
      print(paste0("Control condition: ", strControlName))
    }
    print(paste0( "Found time points: ",
      paste( lsTimepoints,collapse=",") ))
    for(tp in lsTimepoints){
      print(paste0( "Case: Found the following replicates for time points: ", tp,
        ": ", paste0(dfAnnotationFull[
          (dfAnnotationFull$Time %in% tp) &
            (dfAnnotationFull$Condition %in% strCaseName) &
            (dfAnnotationFull$Replicate %in% colnames(arr2DCountData)),
          ]$Replicate,collapse=","),collapse="," ))
    }
    if(!is.null(strControlName)){
      for(tp in lsTimepoints){
        print(paste0( "Control: Found the following replicates for time points: ", tp,
          ":", paste(dfAnnotationFull[
            dfAnnotationFull$Time %in% tp &
              dfAnnotationFull$Condition %in% strControlName &
              dfAnnotationFull$Replicate %in% colnames(arr2DCountData),
            ]$Replicate,collapse=",") ))
      }
    }
    
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
      lsReplicateNames_Case <- as.character(as.vector( 
        dfAnnotationFull[dfAnnotationFull$Condition %in% strCaseName,]$Replicate ))
      lsReplicateNames_Ctrl <- as.character(as.vector( 
        dfAnnotationFull[dfAnnotationFull$Condition %in% strControlName,]$Replicate ))
      arr2DCountData <- arr2DCountData[,c(lsReplicateNames_Case,lsReplicateNames_Ctrl)]
    } else {
      lsReplicateNames_Case <- as.character(as.vector( 
        dfAnnotationFull[dfAnnotationFull$Condition %in% strCaseName,]$Replicate ))
      arr2DCountData <- arr2DCountData[,lsReplicateNames_Case]
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
    # with mean expression > 1
    #rowIdx_lowCounts <- apply(arr2DCountData,1,function(x){mean(x) < 1})
    if(sum(rowIdx_lowCounts) > 0){
      print(paste0("WARNING: ",sum(rowIdx_lowCounts), " out of ",
        dim(arr2DCountData)[1]," genes had only zero counts in all considered samples."))
      print("These genes are omitted in the analysis.")
      arr2DCountData <- arr2DCountData[!rowIdx_lowCounts,]
    }
    # DAVID to be deprecated
    # Shorten expression table
    if(TRUE){
      ind_toKeep <- 1:100
      print(paste0("Working on subset of data: ",length(ind_toKeep)," genes."))
      arr2DCountData <- arr2DCountData[ind_toKeep,]
    }
    # Exclude genes with only missing values (NAs)
    indx_NA <- apply(arr2DCountData,1,function(x){all(is.na(x))})
    if(any(indx_NA)){
      print(paste0("WARNING: Excluded ",sum(indx_NA)," genes because no real valued samples were given."))
    }
    arr2DCountData <- arr2DCountData[!(indx_NA),]
    
    return(arr2DCountData)
  }
  
  # Prepare 3D count data array
  reshapeCountData <- function(dfAnnotationFull=NULL, arr2DCountData=NULL,
    strCaseName = NULL, strControlName=NULL){
    
    lsSamples <- unique(as.vector( dfAnnotationFull$Sample ))
    #lsTimepoints <- unique(as.vector( dfAnnotationFull$Time ))
    nMaxRep <- max(unlist( lapply(lsSamples,function(sample){sum(
      dfAnnotationFull$Sample %in% sample)}) ))
     
    # Convert list of matrices into 3D array
    # Initialise 3D array
    nGenes <- dim(arr2DCountData)[1]
    nSamples <- length(lsSamples)
    nRep <- nMaxRep
    arr3DCountData=array(NA,c(nGenes,nSamples,nRep))
    rownames(arr3DCountData) <- rownames(arr2DCountData)
    colnames(arr3DCountData) <- lsSamples
    
    # Write 3D array
    for (sample in lsSamples){
      lsReplicates <- dfAnnotationFull[dfAnnotationFull$Sample %in% sample,]$Replicate
      for (replicate in lsReplicates){
        arr3DCountData[,sample,match(replicate,lsReplicates)] <- arr2DCountData[ , replicate ]
      }
    }
    
    return(arr3DCountData)
  }
  
  # Create smaller annotation table with one entry per sample.
  # This is the annotation table for the 3D count data array which has
  # sample names as column names.
  reduceAnnotation <- function(dfAnnotationFull=NULL){
    # Find first occurrence of each sample in annotation table:
    # Note: Condition and time are the same for all replicates 
    # of a given sample.
    indSamples <- match( unique(as.vector(dfAnnotationFull$Sample)), 
      as.vector(dfAnnotationFull$Sample) )
    # Take one row per sample and columns [Sample,Condition,Time]
    # From now on, the replicates are anonymous in the 3D array.
    dfAnnotationRed <- dfAnnotationFull[indSamples,
      c("Sample","Condition","Time")]
  }
  
  ###############################################################
  # (II) Main body of function
  
  checkData(dfAnnotationFull=dfAnnotationFull, arr2DCountData=matCountData,
    strControlName=strControlName, strCaseName=strCaseName)
  
  arr2DCountData <- nameGenes(arr2DCountData=matCountData)
  
  arr2DCountData <- reduceCountData(dfAnnotationFull=dfAnnotationFull, 
    arr2DCountData=arr2DCountData)
  
  arr3DCountData <- reshapeCountData(
    dfAnnotationFull=dfAnnotationFull, arr2DCountData=arr2DCountData,
    strControlName=strControlName, strCaseName=strCaseName)

  dfAnnotationRed <- reduceAnnotation(dfAnnotationFull)
  
  return(list(arr2DCountData,arr3DCountData,dfAnnotationRed))
}