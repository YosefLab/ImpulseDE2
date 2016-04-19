#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Annotation preparation    +++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Check validity of input and process count data matrix into arrays used
### later in runImpulseDE.

### Structure:
### (I) Helper functions:
###   checkData()
###   nameGenes()
###   reduceCountData()
###   reshapeCountData()
### (II) Script body

# INPUT:
#   dfAnnotationFull: (Table samples x 2 [time and condition]) providing 
#       co-variables for the samples including condition and time points.
#       Time points must be numeric numbers.
#   data_tables: (Numeric 3D array genes x samples x lsReplicates).
#       Contains expression values or similar locus-specific read-outs.
#   control_timecourse: (bool) [Default FALSE]control time timecourse is part of 
#       the data set (TRUE) or not (FALSE).
#   control_name: (str) name of the control condition in dfAnnotationFull.
#   strCaseName: (str) name of the case condition in dfAnnotationFull.
# OUTPUT:
#   annot: (Table samples x 2[time and condition]) dfAnnotationFull reduced to 
#       target samples 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

processData <- function(dfAnnotationFull=NULL, matCountData=NULL,
  strCaseName=NULL, strControlName=NULL){
  
  ###############################################################
  # (I) Helper functions
  
  # Check format and presence of input data
  checkData <- function(dfAnnotationFull=NULL, matCountData=NULL,
    strCaseName=NULL, strControlName=NULL){
    
    ### 1. Check that all necessary input was specified
    # dfAnnotationFull
    if(is.null(dfAnnotationFull)){
      stop("ERROR: Annotation table was not given as input: dfAnnotationFull")
    }
    # matCountData
    if(is.null(matCountData)){
      stop("ERROR: List of count tables was not given as input: matCountData")
    }
    # strCaseName
    if(is.null(strCaseName)){
      stop("ERROR: No name for case condition given.")
    }
    
    ### 2. Check annotation table content
    ### a) Check column names
    lsColNamesRequired <- c("Replicate_name","Replicate","Sample","Condition","Time")
    if( !all(lsColNamesRequired %in% colnames(dfAnnotationFull)) ){
      stop(paste0("Could not find column ",
        lsColNamesRequired[!(lsColNamesRequired %in% colnames(dfAnnotationFull))],
        " in annotation table."))
    }
    ### b) Replicates
    # Check that all replicates contain all timepoints
    lsReplicates <- unique(as.vector( dfAnnotationFull$Replicate ))
    for(iRep in 1:length(lsReplicates)){
      if(!identical(unique(dfAnnotationFull$Time),unique(dfAnnotationFull[
        dfAnnotationFull$Replicate==lsReplicates[iRep],]$Time))){
        stop(paste0("ERROR: [Annotation table] Not all timepoints listed for each replicate ",lsReplicates[iRep]))
      }
    }
    # Check that all replicates contain the same number of samples
    for(iRep in 1:length(lsReplicates)){
      # Establish number of samples in first replicate
      if(iRep==1){
        nSamples <- sum( dfAnnotationFull$Replicate==lsReplicates[iRep] )
        # Check for equality of number of samples in remaining replicates
      } else {
        if(sum( dfAnnotationFull$Replicate==lsReplicates[iRep] ) != nSamples){
          stop(paste0("ERROR: [Annotation table] Replicate ",lsReplicates[iRep],
            " contains different number of samples compared to replicate ",
            lsReplicates[1]))
        }
      }
    }
    # Check that replicate name do not occur twice
    if(length(unique(dfAnnotationFull$Replicate_name)) != dim(dfAnnotationFull)[1]){
      stop(paste0("ERROR: [Annotation table]",
        "Number of Replicate_name instances different to annotation table number of rows.",
        "Replicate_name cannot occur multiple times."))
    }
    # Check that each replicate is associated with one condition only
    for(replicate in lsReplicates){
      if( length(unique( dfAnnotationFull[
        dfAnnotationFull$Replicate %in% replicate,"Replicate_name"] )) > 1 ){
        
      }
    }
    ### c) Time points
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
    if(!(strControlName %in% lsConditions)){
      stop("ERROR: Condition given for control does not occur in annotation table condition column.")
    }
    
    ### 3. Check expression table content
    # Test that all entries in count data table occur in annotation table
    if( any(!(colnames(arr2DCountData) %in% dfAnnotationFull$Replicate_name)) ){
      stop(paste0("ERROR: The lsReplicates ",
        as.character( colnames(arr2DCountData)[
          !(colnames(arr2DCountData) %in% dfAnnotationFull$Replicate_name)] ),
        " indicated in count data table do not occur in annotation table."))
    }
    
    ### Summarise which conditions, samples, lsReplicates were found
    print(paste0("Found conditions:",lsConditions))
    print(paste0("Case condition: ", strCaseName))
    if(!is.null(strControlName)){
      print(paste0("Control condition: ", strControlName))
    }
    print(paste0("Found time points: ",
      as.numeric(as.vector( unique(dfAnnotationFull$Time)) )) )
    for(replicate in lsReplicates){
      print(paste0("Found the following samples for replicate ",replicate,
        ":",dfAnnotationFull[dfAnnotationFull$Replicate %in% replicate] ))
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
  reduceCountData <- function(dfAnnotationFull=NULL, matCountData=NULL){
    
    ### 1. Columns (Conditions):
    # Reduce expression table to columns of considered conditions
    if(!is.null(strControlName)){
      lsReplicateNames_Case <- as.character(as.vector( 
        dfAnnotationFull[dfAnnotationFull$Condition %in% strCaseName,]$Replicate_name ))
      lsReplicateNames_Ctrl <- as.character(as.vector( 
        dfAnnotationFull[dfAnnotationFull$Condition %in% strControlName,]$Replicate_name ))
      arr2DCountData <- matCountData[,c(lsReplicateNames_Case,lsReplicateNames_Ctrl)]
    } else {
      lsReplicateNames_Case <- as.character(as.vector( 
        dfAnnotationFull[dfAnnotationFull$Condition %in% strCaseName,]$Replicate_name ))
      arr2DCountData <- matCountData[,lsReplicateNames_Case]
    }
    ### 2. Rows (Genes):
    # Reduce expression table to rows containing at least one non-zero count
    rowIdx_lowCounts <- apply(arr2DCountData,1,function(x){all(x==0)})
    # with mean expression > 1
    #rowIdx_lowCounts <- apply(arr2DCountData,1,function(x){mean(x) < 1})
    if(sum(rowIdx_lowCounts) > 0){
      print(paste0("WARNING: ",sum(rowIdx_lowCounts), " out of ",
        dim(arr2DCountData)[1]," genes had only zero counts in all considered samples."))
      #print(paste0("WARNING: ",sum(rowIdx_lowCounts), " out of ",
      #  dim(arr2DCountData)[1]," genes had a mean RNA count of < 1."))
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
    
    return(arr2DCountData)
  }
  
  # Prepare 3D count data array
  reshapeCountData <- function(dfAnnotationFull=NULL, matCountData=NULL,
    strCaseName = NULL, strControlName=NULL){
    
    # Assign samples to lsReplicates
    arr2DCountDatas <- list()
    lsReplicates <- unique(dfAnnotationFull$Replicate)
    for(iRep in 1:length(lsReplicates)){
      arr2DCountDatas[[iRep]] <- arr2DCountData[,as.character( dfAnnotationFull[
        dfAnnotationFull$Replicate %in% lsReplicates[iRep],]$Replicate_name) ]
      # Call columns after sample (not specific replicate)
      colnames(arr2DCountDatas[[iRep]]) <- as.character( 
        dfAnnotationFull[dfAnnotationFull$Replicate %in% lsReplicates[iRep],]$Sample )
    }
    
    # Convert list of matrices into 3D array  
    arr3DCountData=array(NA,c(dim(arr2DCountDatas[[1]]),length(arr2DCountDatas)))
    for (i in 1:length(arr2DCountDatas)){
      arr3DCountData[,,i] <- as.matrix(arr2DCountDatas[[i]])
    }
    rownames(arr3DCountData) <- rownames(arr2DCountDatas[[1]])
    colnames(arr3DCountData) <- colnames(arr2DCountDatas[[1]])
    
    # DAVID: Can my model handle this?
    # Exclude genes with missing values(NAs)
    indx_NA <- apply(arr3DCountData,1,function(x){TRUE %in% is.na(x)})
    arr3DCountData <- arr3DCountData[!(indx_NA),,]
    
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
  
  checkData(dfAnnotationFull=dfAnnotationFull, matCountData=matCountData,
    strControlName=strControlName, strCaseName=strCaseName)
  
  arr2DCountData <- nameGenes(matCountData=matCountData)
  
  arr2DCountData <- reduceCountData(dfAnnotationFull=dfAnnotationFull, 
    matCountData=matCountData)
  
  arr3DCountData <- reshapeCountData(
    dfAnnotationFull=dfAnnotationFull, matCountData=matCountData,
    strControlName=strControlName, strCaseName=strCaseName)
  
  dfAnnotationRed <- reduceAnnotation(dfAnnotationFull)
  
  return(list(arr2DCountData,arr3DCountData,dfAnnotationRed))
}