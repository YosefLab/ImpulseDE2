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
    lsColNamesRequired <- c("Replicate_name","Replicate","Sample","Condition","Time")
    if( !all(lsColNamesRequired %in% colnames(dfAnnotationFull)) ){
      stop(paste0("Could not find column ",
        lsColNamesRequired[!(lsColNamesRequired %in% colnames(dfAnnotationFull))],
        " in annotation table."))
    }
    ### b) Replicates
    lsReplicates <- unique(as.vector( dfAnnotationFull$Replicate ))
    # Check that all replicates contain all timepoints -DAVID -deprecated
    #for(iRep in 1:length(lsReplicates)){
    #  if(!identical(unique(dfAnnotationFull$Time),unique(dfAnnotationFull[
    #    dfAnnotationFull$Replicate==lsReplicates[iRep],]$Time))){
    #    stop(paste0("ERROR: [Annotation table] Not all timepoints listed for each replicate ",lsReplicates[iRep]))
    #  }
    #}
    # Check that all replicates contain the same number of samples -DAVID - deprecated
    #for(iRep in 1:length(lsReplicates)){
    #  # Establish number of samples in first replicate
    #  if(iRep==1){
    #    nSamples <- sum( dfAnnotationFull$Replicate==lsReplicates[iRep] )
    #    # Check for equality of number of samples in remaining replicates
    #  } else {
    #    if(sum( dfAnnotationFull$Replicate==lsReplicates[iRep] ) != nSamples){
    #      stop(paste0("ERROR: [Annotation table] Replicate ",lsReplicates[iRep],
    #        " contains different number of samples compared to replicate ",
    #        lsReplicates[1]))
    #    }
    #  }
    #}
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
    if(!is.null(strControlName)){
      if(!(strControlName %in% lsConditions)){
        stop("ERROR: Condition given for control does not occur in annotation table condition column.")
      }
    }
    
    ### 3. Check expression table content
    # Test that all entries in count data table occur in annotation table
    if( any(!(colnames(arr2DCountData) %in% dfAnnotationFull$Replicate_name)) ){
      print(paste0("WARNING: The column ",
        as.character( colnames(arr2DCountData)[
          !(colnames(arr2DCountData) %in% dfAnnotationFull$Replicate_name)] ),
        " in the count data table do not occur in annotation table and will be ignored."))
    }
    
    ### Summarise which conditions, samples, lsReplicates were found
    print(paste0("Found conditions:",lsConditions))
    print(paste0("Case condition: ", strCaseName))
    if(!is.null(strControlName)){
      print(paste0("Control condition: ", strControlName))
    }
    print(paste0( "Found time points: ",
      paste( dfAnnotationFull$Time,collapse=",") ))
    for(replicate in lsReplicates){
      print(paste0( "Found the following time points for replicate ", replicate,
        ":", paste(dfAnnotationFull[dfAnnotationFull$Replicate %in% replicate,]$Time,collapse=",") ))
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
        dfAnnotationFull[dfAnnotationFull$Condition %in% strCaseName,]$Replicate_name ))
      lsReplicateNames_Ctrl <- as.character(as.vector( 
        dfAnnotationFull[dfAnnotationFull$Condition %in% strControlName,]$Replicate_name ))
      arr2DCountData <- arr2DCountData[,c(lsReplicateNames_Case,lsReplicateNames_Ctrl)]
    } else {
      lsReplicateNames_Case <- as.character(as.vector( 
        dfAnnotationFull[dfAnnotationFull$Condition %in% strCaseName,]$Replicate_name ))
      arr2DCountData <- arr2DCountData[,lsReplicateNames_Case]
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
    # DAVID: Can my model handle this? - yes - deprecated
    # Exclude genes with missing values(NAs)
    #indx_NA <- apply(arr2DCountData,1,function(x){any(is.na(x))})
    #arr2DCountData <- arr2DCountData[!(indx_NA),]
    
    return(arr2DCountData)
  }
  
  # Prepare 3D count data array
  reshapeCountData <- function(dfAnnotationFull=NULL, arr2DCountData=NULL,
    strCaseName = NULL, strControlName=NULL){
    
    lsReplicates <- unique(as.vector( dfAnnotationFull$Replicate ))
    lsSamples <- unique(as.vector( dfAnnotationFull$Sample ))
    
    # Convert list of matrices into 3D array
    # Initialise 3D array
    nGenes <- dim(arr2DCountData)[1]
    nSamples <- length(lsSamples)
    nRep <- length(lsReplicates)
    arr3DCountData=array(NA,c(nGenes,nSamples,nRep))
    rownames(arr3DCountData) <- rownames(arr2DCountData)
    colnames(arr3DCountData) <- colnames(lsSamples)
    dimnames(arr3DCountData)[[3]] <- lsReplicates
    
    # Write 3D array
    for (sample in lsSamples){
      for (replicate in lsReplicates){
        # Substitute NA column by observations if this combination was observed
        # I.e. if this sample is included in the given replicate
        # Using more rebust boolean indexing here
        if(any((dfAnnotationFull$Sample %in% sample) & 
            (dfAnnotationFull$Replicate %in% replicate))){
          arr3DCountData[,colnames(arr3DCountData) %in% sample,
            dimnames(arr3DCountData)[[3]] %in% replicate] <- as.vector(
              arr2DCountData[ , dfAnnotationFull[(dfAnnotationFull$Sample %in% sample) &
                  (dfAnnotationFull$Replicate %in% replicate),]$Replicate_name ] )
        }
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