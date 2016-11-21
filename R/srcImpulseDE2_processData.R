#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++     Annotation preparation    +++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Process annotation and count data
#' 
#' Check validity of input and process count data matrix and annotation
#' into data structures used later in \code{runImpulseDE2}.
#' \code{processData} is structure in the following way:
#' (I) Subhelper functions:
#'    checkNull() Check whether object was supplied (is not NULL).
#'    checkDimMatch() Checks whether dimensions of matrices agree.
#'    checkElementMatch() Checks whether vectors are identical.
#'    checkNumeric() Checks whether elements are numeric.
#'    checkProbability() Checks whether elements are probabilities.
#'    checkCounts() Checks whether elements are count data.
#' (II) Helper functions:
#'    checkData() Check format and presence of input data.
#'    nameGenes() Name genes if names are not given.
#'    procAnnotation() Add categorial time variable to annotation table.
#'    reduceCountData() Reduce count data to data which are utilised later.
#' (III) Script body
#' 
#' @seealso Called by \code{runImpulseDE2}.
#' 
#' @param matCountData: (matrix genes x samples) [Default NULL] 
#'    Count data of all conditions, unobserved entries are NA. 
#' @param dfAnnotation: (Table) [Default NULL] 
#'    Annotation table. Lists co-variables of samples: 
#'    Sample, Condition, Time (and Batch). 
#'    Time must be numeric.
#' @param strCaseName (str) [Default NULL] 
#'    Name of the case condition in \code{dfAnnotation}.
#' @param strControlName: (str) [Default NULL] 
#'    Name of the control condition in \code{dfAnnotation}.
#' @param strMode: (str) [Default "singelbatch"] 
#'    {"singelbatch","batcheffects"}
#'    Batch model.
#' @param vecDispersionsExternal: (vector length number of
#'    genes in matCountData) [Default NULL]
#'    Externally generated list of gene-wise dispersion factors
#'    which overides DESeq2 generated dispersion factors.
#' @param vecSizeFactorsExternal: (vector length number of
#'    cells in matCountData) [Default NULL]
#'    Externally generated list of size factors which override
#'    size factor computation in ImpulseDE2.
#' @param matBatchFactorsExternal: (numeric matrix genes x cells) 
#'    [Default NULL] Externally generated list of translation factors 
#'    which override translation factor computation in ImpulseDE2.
#' @param boolRunDESeq2: (bool) [Default TRUE]
#'    Whether to run DESeq2.
#'    
#' @return (list length 3) with the following elements:
#' \itemize{
#'  \item \code{matCountDataProc}: (count matrix  genes x samples) 
#'      Count data: Reduced version of \code{matCountData}. 
#' }
#'
#' @export

processData <- function(dfAnnotation=NULL, 
  matCountData=NULL,
  scaSmallRun=NULL,
  strCaseName=NULL, 
  strControlName=NULL, 
  strMode=NULL,
  vecDispersionsExternal=NULL,
  vecSizeFactorsExternal=NULL,
  matBatchFactorsExternal=NULL,
  boolRunDESeq2=NULL){
  
  ###############################################################
  # (I) Helper functions
  
  # Check whether object was supplied (is not NULL).
  checkNull <- function(objectInput,strObjectInput){
    if(is.null(objectInput)){
      stop(paste0( "ERROR: ", strObjectInput," was not given as input." ))
    }
  }
  # Check whether object does not have NA elements.
  checkNA <- function(objectInput,strObjectInput){
    if(is.na(objectInput)){
      stop(paste0( "ERROR: ", strObjectInput," is NA and needs to be specifief." ))
    }
  }
  # Checks whether dimensions of matrices agree.
  checkDimMatch <-function(matInput1, matInput2, strMatInput1, strMatInput2){
    if(any(dim(matInput1)!=dim(matInput2))){
      stop(paste0( "ERROR: ", strMatInput1, " does not have the dimensions as ", strMatInput2 , "." ))
    }
  }
  # Checks whether vectors are identical.
  checkElementMatch <- function(vec1, vec2, strVec1, strVec2){
    if(!any(vec1==vec2)){
      stop(paste0( "ERROR: ",strVec1 ," do not agree with ", strVec2, "." ))
    }
  }
  # Checks whether elements are numeric
  checkNumeric <- function(matInput, strMatInput){
    if(any(!is.numeric(matInput))){
      stop(paste0( "ERROR: ", strMatInput, " contains non-numeric elements. Requires numeric data." ))
    }
  }
  # Checks whether elements are probabilities (in [0,1]).
  checkProbability <- function(matInput, strMatInput){
    checkNumeric(matInput, strMatInput)
    if(any(matInput < 0 | matInput > 1 | is.na(matInput))){
      stop(paste0( "ERROR: ", strMatInput, " contains elements outside of interval [0,1]." ))
    }
  }
  # Checks whether elements are count data: non-negative integer finite numeric elements.
  # Note that NA are allowed. Can be used to check whether element is integer if NA 
  # is checked separately.
  checkCounts <- function(matInput, strMatInput){
    checkNumeric(matInput, strMatInput)
    if(any(matInput[!is.na(matInput)] %% 1 != 0)){
      stop(paste0( "ERROR: ", strMatInput, " contains non-integer elements. Requires count data." ))
    }
    if(any(!is.finite(matInput[!is.na(matInput)]))){
      stop(paste0( "ERROR: ", strMatInput, " contains infinite elements. Requires count data." ))
    }
    if(any(matInput[!is.na(matInput)]<0)){
      stop(paste0( "ERROR: ", strMatInput, " contains negative elements. Requires count data." ))
    }
  }
  
  # Check format and presence of input data.
  checkData <- function(dfAnnotation=NULL, 
    matCountData=NULL,
    scaSmallRun=NULL,
    strCaseName=NULL, 
    strControlName=NULL, 
    strMode=NULL,
    vecDispersionsExternal=NULL,
    vecSizeFactorsExternal=NULL,
    matBatchFactorsExternal=NULL,
    boolRunDESeq2=NULL ){
    
    ### 1. Check that all necessary input was specified
    checkNull(dfAnnotation,"dfAnnotation")
    checkNull(matCountData,"matCountData")
    checkNull(strCaseName,"strCaseName")
    checkNull(strMode,"strMode")
    
    ### 2. Check mode
    lsAllowedModes <- c("singlebatch", "batcheffects")
    if(!(strMode %in% lsAllowedModes)){
      stop(paste0( "ERROR: ImpulseDE2 mode given as input, strMode=", strMode,
        ", is not recognised. Chose from {", paste0(lsAllowedModes, collapse=","), "}." ))
    }
    
    ### 3. Check annotation table content
    ### a) Check column names
    if(strMode=="singlebatch"){
      vecColNamesRequired <- c("Sample","Condition","Time")
    } else if(strMode=="batcheffects"){
      vecColNamesRequired <- c("Sample","Condition","Time","Batch")
    } else {
      stop(paste0("ERROR: Unrecognised strMode in processData::checkData(): ",strMode))
    }
    if( !all(vecColNamesRequired %in% colnames(dfAnnotation)) ){
      stop(paste0("Could not find column ",
        vecColNamesRequired[!(vecColNamesRequired %in% colnames(dfAnnotation))],
        " in annotation table."))
    }
    ### b) Samples
    # Check that sample name do not occur twice
    if(any(duplicated(dfAnnotation$Sample))){
      stop(paste0("ERROR: [Annotation table] ",
        "Sample names must be unique: Sample(s) ",
        paste0((dfAnnotation$Sample)[duplicated(dfAnnotation$Sample)],collapse=","),
        " is/are duplicated."))
    }
    ### c) Time points
    vecTimepoints <- unique(as.vector( dfAnnotation$Time ))
    # Check that time points are numeric
    checkNumeric(dfAnnotation$Time, "dfAnnotation$Time")
    ### d) Conditions
    lsConditions <- unique( dfAnnotation$Condition )
    # Check that given conditions exisit in annotation table
    if(!(strCaseName %in% lsConditions)){
      stop("ERROR: Condition given for case does not occur in annotation table condition column.")
    }
    if(!is.null(strControlName)){
      if(!(strControlName %in% lsConditions)){
        stop("ERROR: Condition given for control does not occur in annotation table condition column.")
      }
    }
    ### e) Batch
    if(strMode=="batcheffects"){
      # Check that number of time courses given is > 1
      if(length(unique( dfAnnotation$Batch ))==1){
        stop("ERROR: Only one time course given in annotation table in mode batcheffects. Use singlebatch mode.")
      }
    }
    
    ### 4. Check expression table content
    # Check that all entries in count data table occur in annotation table
    if( any(!(colnames(matCountData) %in% dfAnnotation$Sample)) ){
      print(paste0("WARNING: The column(s) ",
        paste0(as.character( colnames(matCountData)[
          !(colnames(matCountData) %in% dfAnnotation$Sample)] ),collapse=","),
        " in the count data table do(es) not occur in annotation table and will be ignored."))
    }
    checkNull(rownames(matCountData),"[Rownames of matCountData]")
    checkCounts(matCountData, "matCountData")
    # Check that gene IDs don't start with "_" -TODO
    
    ### 5. Check supplied dispersion vector
    if(!is.null(vecDispersionsExternal)){
      # Check that dispersions were named
      if(is.null(names(vecDispersionsExternal))){
        stop("ERROR: vecDispersionsExternal was not named. Name according to rownames of matCountData.")
      }
      # Check that one dispersion was supplied per gene
      if(any( !(names(vecDispersionsExternal) %in% rownames(matCountData)) ) |
          any( !(rownames(matCountData) %in% names(vecDispersionsExternal)) ) ){
        stop("ERROR: vecDispersionsExternal supplied but names do not agree with rownames of matCountData.")
      }
      # Check that dispersion vector is numeric
      checkNumeric(vecDispersionsExternal, "vecDispersionsExternal")
      # Check that dispersions are positive (should not have sub-poissonian noise in count data)
      if(any(vecDispersionsExternal < 0)){
        warning(paste0( "WARNING: vecDispersionsExternal contains negative elements which corresponds to sub-poissonian noise.",
          "These elements are kept as they are in the following." ))
      }
    }
    
    ### 6. Check supplied size facotrs
    if(!is.null(vecSizeFactorsExternal)){
      # Check that size factors were named
      if(is.null(names(vecSizeFactorsExternal))){
        stop("ERROR: vecSizeFactorsExternal was not named. Name according to colnames of matCountData.")
      }
      # Check that one size factors was supplied per cell
      if(any( !(names(vecSizeFactorsExternal) %in% colnames(matCountData)) ) |
          any( !(colnames(matCountData) %in% names(vecSizeFactorsExternal)) ) ){
        stop("ERROR: vecSizeFactorsExternal supplied but names do not agree with colnames of matCountData.")
      }
      # Check that size factors vector is numeric
      checkNumeric(vecSizeFactorsExternal, "vecSizeFactorsExternal")
      # Check that size factors are positive
      if(any(vecSizeFactorsExternal <= 0)){
        stop(paste0( "WARNING: vecSizeFactorsExternal contains negative or zero elements which leads.",
          "Size factors must be positive, remove samples if size factor is supposed to be zero." ))
      }
    }
    
    ### 7. Check DESeq2 settings
    if(is.null(vecDispersionsExternal) & !boolRunDESeq2){
      stop(paste0( "ERROR: vecDispersionsExternal not supplied and boolRunDESeq2 is FALSE.",
        "Dispersions have to be computed by DESeq2 or provided externally." ))
    }
    
    ### 8. Check scaSmallRun
    if(!is.null(scaSmallRun)){
      checkCounts(scaSmallRun, "scaSmallRun")
      if(scaSmallRun > dim(matCountData)[1]){
        stop(paste0( "ERROR: scaSmallRun (",scaSmallRun,
          ") larger then data set (",dim(matCountData)[1]," genes)."))
      }
    }
    
    ### Summarise which mode, conditions, samples and
    ### batch were found
    print(paste0("Found conditions: ",paste0(lsConditions,collapse=",")))
    print(paste0("Case condition: '", strCaseName, "'"))
    if(!is.null(strControlName)){
      print(paste0("Control condition: '", strControlName, "'"))
    }
    print(paste0( "ImpulseDE2 runs in mode: ", strMode ))
    if(strMode=="singlebatch" | strMode=="batcheffects"){
      print(paste0( "Found time points: ",
        paste( vecTimepoints, collapse=",") ))
      for(tp in vecTimepoints){
        print(paste0( "Case: Found the samples at time point ", 
          tp,": ", 
          paste0(dfAnnotation[
            (dfAnnotation$Time %in% tp) &
              (dfAnnotation$Condition %in% strCaseName) &
              (dfAnnotation$Sample %in% colnames(matCountData)),
            ]$Sample,collapse=","),collapse="," ))
      }
      if(!is.null(strControlName)){
        for(tp in vecTimepoints){
          print(paste0( "Control: Found the following samples at time point ", 
            tp, ":", 
            paste0(dfAnnotation[
              dfAnnotation$Time %in% tp &
                dfAnnotation$Condition %in% strControlName &
                dfAnnotation$Sample %in% colnames(matCountData),
              ]$Sample,collapse=","),collapse="," ))
        }
      }
    }
    if(strMode=="batcheffects"){
      for(batch in unique( dfAnnotation$Batch )){
        print(paste0( "Found the following samples for batches ",
          batch, ": ",
          paste0( dfAnnotation[
            dfAnnotation$Batch %in% batch &
              dfAnnotation$Sample %in% colnames(matCountData),
            ]$Sample,collapse=",") ))
      }
    }
    
    return(NULL)
  }
  
  # Add categorial time variable to annotation table which
  # differentiates case and control time points (given to DESeq2).
  # Add column with time scalars with underscore prefix.
  procAnnotation <- function(dfAnnotation,
    matCountData,
    strMode,
    strCaseName, 
    strControlName){
    
    if(is.null(strControlName)){
      dfAnnotationProc <- dfAnnotation[dfAnnotation$Condition==strCaseName,]
    } else {
      dfAnnotationProc <- dfAnnotation[dfAnnotation$Condition==strCaseName | 
          dfAnnotation$Condition==strControlName,]
      # The following is used to build full-rank model matrices
      # for DESeq2 even though conditions are a linear combination of
      # longitudinal time series.
      vecBatchCase <- unique(dfAnnotationProc[dfAnnotation$Condition==strCaseName,]$Batch)
      vecBatchCtrl <- unique(dfAnnotationProc[dfAnnotation$Condition==strControlName,]$Batch)
      vecindBatchNested <- array(NA, dim(dfAnnotationProc)[1])
      vecindBatchNested[dfAnnotation$Condition==strCaseName] <- 
        match(dfAnnotationProc$Batch, vecBatchCase)[dfAnnotation$Condition==strCaseName]
      vecindBatchNested[dfAnnotation$Condition==strControlName] <- 
        match(dfAnnotationProc$Batch, vecBatchCtrl)[dfAnnotation$Condition==strControlName]
      dfAnnotationProc$BatchNested <- paste0("_", vecindBatchNested)
    }
    dfAnnotationProc$TimeCateg <- paste0(
      "_", dfAnnotationProc$Time)
    
    # Reduce to samples which occur in count table
    # This allows to represent entire library of samples in annotation
    # file even if not all samples were measured yet. Not measured samples
    # which are not mentioned in count table are ignored. Thereby, the same
    # full annotation table can be used throughout sample collection on 
    # incomplete data sets.
    dfAnnotationProc <- dfAnnotationProc[dfAnnotationProc$Sample %in% colnames(matCountData),]
    
    return(dfAnnotationProc)
  }
  
  # Name genes if names are not given.
  nameGenes <- function(matCountDataProc=NULL){
    if(is.null(rownames(matCountDataProc))){
      rownames(matCountDataProc) <- paste0("Region_", seq(1,nrow(matCountDataProc)))
    }
    return(matCountDataProc)
  }
  
  # Reduce count data to data which are utilised later
  reduceCountData <- function(dfAnnotation=NULL, matCountDataProc=NULL){
    
    print(paste0("Input contained ",dim(matCountDataProc)[1]," genes/regions."))
    ### 1. Columns (Conditions):
    # Reduce expression table to columns of considered conditions
    if(!is.null(strControlName)){
      vecSampleNames_Case <- as.character(as.vector( 
        dfAnnotation[ 
          dfAnnotation$Condition %in% strCaseName &
          dfAnnotation$Sample %in% colnames(matCountDataProc),
          ]$Sample ))
      vecSampleNames_Ctrl <- as.character(as.vector( 
        dfAnnotation[ 
          dfAnnotation$Condition %in% strControlName &
          dfAnnotation$Sample %in% colnames(matCountDataProc),
          ]$Sample ))
      matCountDataProc <- matCountDataProc[,c(vecSampleNames_Case,vecSampleNames_Ctrl)]
    } else {
      vecSampleNames_Case <- as.character(as.vector( 
        dfAnnotation[ 
          dfAnnotation$Condition %in% strCaseName &
          dfAnnotation$Sample %in% colnames(matCountDataProc),
          ]$Sample ))
      matCountDataProc <- matCountDataProc[,vecSampleNames_Case]
    }
    # Check that every sample contains at least one observed value (not NA)
    vecNARep <- any(apply(matCountDataProc,2,function(rep){all(is.na(rep))}))
    if(any(vecNARep)){
      print(paste0( "WARNING: Sample(s) ",
        paste0(colnames(matCountDataProc[vecNARep,]), collapse=","),
        " only contain(s) NA values and will be removed from the analysis."))
      matCountDataProc <- matCountDataProc[,!vecNARep]
    }
    ### 2. Rows (Genes):
    # Exclude genes with only missing values (NAs)
    indx_NA <- apply(matCountDataProc,1,function(x){all(is.na(x))})
    if(any(indx_NA)){
      print(paste0("WARNING: Excluded ",sum(indx_NA)," genes because no real valued samples were given."))
    }
    matCountDataProc <- matCountDataProc[!(indx_NA),]
    # Reduce expression table to rows containing at least one non-zero count
    rowIdx_lowCounts <- apply(matCountDataProc,1,function(x){all(x[!is.na(x)]==0)})
    if(sum(rowIdx_lowCounts) > 0){
      print(paste0("WARNING: ",sum(rowIdx_lowCounts), " out of ",
        dim(matCountDataProc)[1],
        " genes had only zero counts in all considered samples and are omitted."))
      matCountDataProc <- matCountDataProc[!rowIdx_lowCounts,]
    }
    
    # Sort count matrix column by annotation table
    matCountDataProc <- matCountDataProc[,match(as.vector(dfAnnotation$Sample),colnames(matCountDataProc))]
    
    print(paste0("Selected ",dim(matCountDataProc)[1]," genes/regions for analysis."))
    return(matCountDataProc)
  }
  
  # Reduce data set to small run size if required.
  reduceDataRows <- function(matCountDataProc,
    scaSmallRun=NULL){
    if(!is.null(scaSmallRun)){
      scaNRows <- min(scaSmallRun,dim(matCountDataProc)[1])
      print(paste0("Working on subset of data: ",scaNRows," genes."))
      matCountDataProc <- matCountDataProc[1:scaNRows,]
    }
    return(matCountDataProc)
  }
  
  ###############################################################
  # (II) Main body of function
  
  # Check validity of input
  checkData(
    dfAnnotation=dfAnnotation,
    matCountData=matCountData,
    scaSmallRun=scaSmallRun,
    strControlName=strControlName,
    strCaseName=strCaseName,
    strMode=strMode,
    vecDispersionsExternal=vecDispersionsExternal,
    vecSizeFactorsExternal=vecSizeFactorsExternal,
    matBatchFactorsExternal=matBatchFactorsExternal,
    boolRunDESeq2=boolRunDESeq2 )
  
  # Process annotation table
  dfAnnotationProc <- procAnnotation(dfAnnotation=dfAnnotation,
    matCountData=matCountData,
    strMode=strMode,
    strCaseName=strCaseName, 
    strControlName=strControlName )
  
  # Process raw counts
  matCountDataProc <- nameGenes(matCountDataProc=matCountData)
  matCountDataProc <- reduceCountData(
    dfAnnotation=dfAnnotationProc, 
    matCountDataProc=matCountDataProc)
  # Keep full data set for size factor estimation.
  matCountDataProcFull <- matCountDataProc
  matCountDataProc <- reduceDataRows(matCountDataProc=matCountDataProc,
    scaSmallRun=scaSmallRun)
  
  # Reduce externally provided parameters according to reduced data set
  # and reorder according to given data set.
  if(!is.null(vecDispersionsExternal)){
    vecDispersionsExternal <- vecDispersionsExternal[rownames(matCountDataProc)]
  }
  if(!is.null(vecSizeFactorsExternal)){
    vecSizeFactorsExternal <- vecSizeFactorsExternal[colnames(matCountDataProc)]
  }

  return( list(matCountDataProc=matCountDataProc,
    matCountDataProcFull=matCountDataProcFull,
    dfAnnotationProc=dfAnnotationProc) )
}