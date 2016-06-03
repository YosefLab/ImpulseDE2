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
#' @param matCountData: (matrix genes x samples) [Default NULL] 
#'    Count data of all conditions, unobserved entries are NA. 
#' @param dfAnnotation: (Table) [Default NULL] 
#'    Annotation table. Lists co-variables of samples: 
#'    Sample, Condition, Time (and LongitudinalSeries). 
#'    Time must be numeric.
#' @param strCaseName (str) [Default NULL] 
#'    Name of the case condition in \code{dfAnnotation}.
#' @param strControlName: (str) [Default NULL] 
#'    Name of the control condition in \code{dfAnnotation}.
#' @param strMode: (str) [Default "batch"] 
#'    {"batch","longitudinal","singlecell"}
#'    Mode of model fitting.
#'    
#' @return (list length 3) with the following elements:
#' \itemize{
#'  \item \code{matCountDataProc}: (count matrix  genes x samples) 
#'      Count data: Reduced version of \code{matCountData}. 
#'  \item \code{matProbNB} (probability matrix genes x samples)
#'      Probability of each observation to originate from the negative
#'      binomial component in the zero inflated negative binomial
#'      mixture model. All entries are 1 if not operating in 
#'      strMode=="singlecell".
#' }
#'
#' @export

processData <- function(dfAnnotation=NULL, matCountData=NULL,
  strCaseName=NULL, strControlName=NULL, strMode=NULL,
  lsPseudoDE=NULL, vecDispersionsExternal=NULL, boolRunDESeq2=NULL){
  
  ###############################################################
  # (I) Helper functions
  
  # Check whether object was supplied (is not NULL).
  checkNull <- function(objectInput){
    if(is.null(objectInput)){
      stop(paste0( "ERROR: ", objectInput," was not given as input." ))
    }
  }
  # Checks whether dimensions of matrices agree.
  checkDimMatch <-function(matInput1, matInput2){
    if(any(dim(matInput1)!=dim(matInput2))){
      stop(paste0( "ERROR: ", matInput1, " does not have the dimensions as ", matInput2 , "." ))
    }
  }
  # Checks whether vectors are identical.
  checkElementMatch <- function(vec1, vec2){
    if(!any(vec1==vec2)){
      stop(paste0( "ERROR: ",vec1 ," do not agree with ", vec2, "." ))
    }
  }
  # Checks whether elements are numeric
  checkNumeric <- function(matInput){
    if(any(!is.numeric(matInput))){
      stop(paste0( "ERROR: ", matInput, " contains non-numeric elements. Requires count data." ))
    }
  }
  # Checks whether elements are probabilities (in [0,1]).
  checkProbability <- function(matInput){
    checkNumeric(matInput)
    if(any(matInput < 0 | matInput > 1 | is.na(matInput))){
      stop(paste0( "ERROR: ", matInput, " contains elements outside of interval [0,1]." ))
    }
  }
  # Checks whether elements are count data: non-negative integer finite numeric elements.
  # Note that NA are allowed.
  checkCounts <- function(matInput){
    checkNumeric(matInput)
    if(any(matInput %% 1 != 0)){
      stop(paste0( "ERROR: ", matInput, " contains non-integer elements. Requires count data." ))
    }
    if(any(!is.finite(matInput))){
      stop(paste0( "ERROR: ", matInput, " contains infinite elements. Requires count data." ))
    }
    if(any(matInput<0)){
      stop(paste0( "ERROR: ", matInput, " contains negative elements. Requires count data." ))
    }
  }
  
  # Check format and presence of input data
  checkData <- function(dfAnnotation=NULL, matCountData=NULL,
    strCaseName=NULL, strControlName=NULL, strMode=NULL,
    lsPseudoDE=NULL, vecDispersionsExternal=NULL,
    boolRunDESeq2=NULL ){
    
    ### 1. Check that all necessary input was specified
    checkNull(dfAnnotation)
    checkNull(matCountData)
    checkNull(strCaseName)
    checkNull(strMode)
    
    ### 2. Check mode
    lsAllowedModes <- c("batch", "longitudinal", "singlecell")
    if(sum(strMode==lsAllowedModes) != 1){
      stop(paste0( "ERROR: ImpulseDE2 mode given as input, strMode=", strMode,
        ", is not recognised. Chose from {", paste0(lsAllowedModes, collapse=","), "}." ))
    }
    
    ### 3. Check annotation table content
    ### a) Check column names
    if(strMode=="batch" | strMode=="singlecell"){
      vecColNamesRequired <- c("Sample","Condition","Time")
    } else if(strMode=="longitudinal"){
      vecColNamesRequired <- c("Sample","Condition","Time","LongitudinalSeries")
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
    if(length(unique(dfAnnotation$Sample)) != length(dfAnnotation$Sample)){
      stop(paste0("ERROR: [Annotation table] ",
        "Number of samples different to number of unique sample names. ",
        "Sample names must be unique."))
    }
    ### c) Time points
    vecTimepoints <- unique(as.vector( dfAnnotation$Time ))
    # Check that time points are numeric
    checkNumeric(dfAnnotation$Time)
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
    ### e) LongitudinalSeries
    if(strMode=="longitudinal"){
      # Check that number of time courses given is > 1
      if(length(unique( dfAnnotation$LongitudinalSeries ))==1){
        stop("ERROR: Only one time course given in annotation table in mode longitudinal. Use batch mode.")
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
    checkNull(rownames(matCountData))
    checkCounts(matCountData)
    
    ### 5. Check PseudoDE objects
    if(strMode=="singlecell"){
      # Check that PseudoDE object was supplied
      checkNull(lsPseudoDE)
      ### a) Dropout rates
      checkNull(lsPseudoDE$matDropout)
      checkDimMatch(lsPseudoDE$matDropout,matCountData)
      checkElementMatch(rownames(matCountData), rownames(lsPseudoDE$matDropout))
      checkElementMatch(colnames(matCountData), colnames(lsPseudoDE$matDropout))
      checkProbability(lsPseudoDE$matDropout)
      
      ### b) Posterior of observation belonging to negative binomial
      ### component in mixture model.
      checkNull(lsPseudoDE$matProbNB)
      checkDimMatch(lsPseudoDE$matProbNB,matCountData)
      checkElementMatch(rownames(matCountData), rownames(lsPseudoDE$matProbNB))
      checkElementMatch(colnames(matCountData), colnames(lsPseudoDE$matProbNB))
      checkProbability(lsPseudoDE$matProbNB)
      
      ### c) Imputed counts
      checkNull(lsPseudoDE$matCountsImputed)
      checkDimMatch(lsPseudoDE$matCountsImputed,matCountData)
      checkElementMatch(rownames(matCountData), rownames(lsPseudoDE$matCountsImputed))
      checkElementMatch(colnames(matCountData), colnames(lsPseudoDE$matCountsImputed))
      checkCounts(lsPseudoDE$matCountsImputed)
    }
    
    ### 6. Check supplied dispersion vector
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
      checkNumeric(vecDispersionsExternal)
      # Check that dispersions are positive (should not have sub-poissonian noise in count data)
      if(any(vecDispersionsExternal < 0)){
        warning(paste0( "WARNING: vecDispersionsExternal contains negative elements which corresponds to sub-poissonian noise.",
          "These elements are kept as they are in the following." ))
      }
    }
    
    ### 7. Check DESeq2 settings
    if(is.null(vecDispersionsExternal) & !boolRunDESeq2){
      stop(paste0( "ERROR: vecDispersionsExternal not supplied and boolRunDESeq2 is FALSE.",
        "Dispersions have to be computed by DESeq2 or provided externally." ))
    }
    
    ### Summarise which mode, conditions, samples and
    ### longitudinal series were found
    print(paste0("Found conditions: ",paste0(lsConditions,collapse=",")))
    print(paste0("Case condition: '", strCaseName, "'"))
    if(!is.null(strControlName)){
      print(paste0("Control condition: '", strControlName, "'"))
    }
    print(paste0( "ImpulseDE2 runs in mode: ", strMode ))
    if(strMode=="batch" | strMode=="longitudinal"){
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
    } else if(strMode=="singlecell"){
      # Shorten output as there are potentially many single cells
      print(paste0( "Case: Number of cells found is ",
        length( dfAnnotation[
          dfAnnotation$Condition %in% strCaseName &
            dfAnnotation$Sample %in% colnames(matCountData),
          ]$Sample),collapse="," ))
      if(!is.null(strControlName)){
        print(paste0( "Control: Number of cells found is ",
          length( dfAnnotation[
            dfAnnotation$Condition %in% strControlName &
              dfAnnotation$Sample %in% colnames(matCountData),
            ]$Sample),collapse="," ))
      }
      print(paste0( "Found ", length(vecTimepoints), " time points."))
    }
    if(strMode=="longitudinal"){
      for(longser in unique( dfAnnotation$LongitudinalSeries )){
        print(paste0( "Found the following samples for longitudinal series ",
          longser, ": ",
          paste0( dfAnnotation[
            dfAnnotation$LongitudinalSeries %in% longser &
              dfAnnotation$Sample %in% colnames(matCountData),
            ]$Sample,collapse=",") ))
      }
    }
    
    return(NULL)
  }
  
  # Add categorial time variable:
  # Add column with time scalars with underscore prefix
  procAnnotation <- function(dfAnnotation=NULL){
    dfAnnotationProc <- dfAnnotation
    dfAnnotationProc$TimeCateg <- paste0(
      rep("_",length(dfAnnotation$Time)), 
      dfAnnotation$Time )
    return(dfAnnotationProc)
  }
  
  # Name genes if names are not given
  nameGenes <- function(matCountDataProc=NULL){
    if(is.null(rownames(matCountDataProc))){
      rownames(matCountDataProc) <- paste("G", 1:nrow(matCountDataProc),
        sep = "_")
    }
    
    return(matCountDataProc)
  }
  
  # Reduce count data to data which are utilised later
  reduceCountData <- function(dfAnnotation=NULL, matCountDataProc=NULL){
    
    ### 1. Columns (Conditions):
    # Reduce expression table to columns of considered conditions
    if(!is.null(strControlName)){
      vecSampleNames_Case <- as.character(as.vector( 
        dfAnnotation[dfAnnotation$Condition %in% strCaseName,]$Sample ))
      vecSampleNames_Ctrl <- as.character(as.vector( 
        dfAnnotation[dfAnnotation$Condition %in% strControlName,]$Sample ))
      matCountDataProc <- matCountDataProc[,c(vecSampleNames_Case,vecSampleNames_Ctrl)]
    } else {
      vecSampleNames_Case <- as.character(as.vector( 
        dfAnnotation[dfAnnotation$Condition %in% strCaseName,]$Sample ))
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
    # Reduce expression table to rows containing at least one non-zero count
    rowIdx_lowCounts <- apply(matCountDataProc,1,function(x){all(x==0)})
    if(sum(rowIdx_lowCounts) > 0){
      print(paste0("WARNING: ",sum(rowIdx_lowCounts), " out of ",
        dim(matCountDataProc)[1]," genes had only zero counts in all considered samples."))
      print("These genes are omitted in the analysis.")
      matCountDataProc <- matCountDataProc[!rowIdx_lowCounts,]
    }
    # DAVID to be deprecated
    # Shorten expression table
    if(TRUE){
      ind_toKeep <- 100
      print(paste0("Working on subset of data: ",min(ind_toKeep,dim(matCountDataProc)[1])," genes."))
      matCountDataProc <- matCountDataProc[1:min(ind_toKeep,dim(matCountDataProc)[1]),]
    }
    # Exclude genes with only missing values (NAs)
    indx_NA <- apply(matCountDataProc,1,function(x){all(is.na(x))})
    if(any(indx_NA)){
      print(paste0("WARNING: Excluded ",sum(indx_NA)," genes because no real valued samples were given."))
    }
    matCountDataProc <- matCountDataProc[!(indx_NA),]
    
    # Sort copunt matrix column by annotation table
    matCountDataProc <- matCountDataProc[,match(as.vector(dfAnnotation$Sample),colnames(matCountDataProc))]
    
    return(matCountDataProc)
  }
  
  ###############################################################
  # (II) Main body of function
  
  # Check validity of input
  checkData(
    dfAnnotation=dfAnnotation,
    matCountData=matCountData,
    strControlName=strControlName,
    strCaseName=strCaseName,
    strMode=strMode,
    lsPseudoDE=lsPseudoDE,
    vecDispersionsExternal=vecDispersionsExternal,
    boolRunDESeq2=boolRunDESeq2 )
  
  # Process annotation table
  dfAnnotationProc <- procAnnotation(dfAnnotation=dfAnnotation)
  
  # Process raw counts
  matCountDataProc <- nameGenes(matCountDataProc=matCountData)
  matCountDataProc <- reduceCountData(
    dfAnnotation=dfAnnotation, 
    matCountDataProc=matCountDataProc)

  # Process single cell hyperparameters
  if(strMode=="singlecell"){
    matProbNB <- lsPseudoDE$matProbNB
    matDropout <- lsPseudoDE$matDropout
    matClusterMeansFitted <- lsPseudoDE$matClusterMeansFitted
  } else {
    matProbNB <- NULL
    matDropout <- NULL
    matClusterMeansFitted <- NULL
  }
  
  lsProcessedData <- list(matCountDataProc, 
    dfAnnotationProc,
    matProbNB, 
    matDropout, 
    matClusterMeansFitted)
  names(lsProcessedData) <- c("matCountDataProc",
    "dfAnnotationProc",
    "matProbNB", 
    "matDropout", 
    "matClusterMeansFitted")
  return( lsProcessedData )
}