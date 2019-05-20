### Check and process input

#' Check and process input to runImpulseDE2()
#' 
#' Check validity of input and process count data matrix and annotation
#' into data structures used later in \link{runImpulseDE2}.
#' \link{processData} is structure in the following way:
#' \itemize{
#'    \item Subhelper functions:
#'    \itemize{
#'      \item checkNull() Check whether object was supplied (is not NULL).
#'      \item checkDimMatch() Checks whether dimensions of matrices agree.
#'      \item checkElementMatch() Checks whether vectors are identical.
#'      \item checkNumeric() Checks whether elements are numeric.
#'      \item checkProbability() Checks whether elements are probabilities.
#'      \item checkCounts() Checks whether elements are count data.
#'    }
#'    \item Helper functions:
#'    \itemize{
#'      \item checkData() Check format and presence of input data.
#'      \item nameGenes() Name genes if names are not given.
#'      \item procAnnotation() Add categorial time variable to 
#'      annotation table.
#'      Add nested batch column if necessary.
#'      Reduce to samples used.
#'      \item reduceCountData() Reduce count data to data which are 
#'      utilised later.
#'    }
#'    \item Script body
#' }
#' 
#' @seealso Called by \link{runImpulseDE2}.
#' 
#' @param matCountData (matrix genes x samples) [Default NULL] 
#'    Read count data, unobserved entries are NA.
#' @param dfAnnotation (data frame samples x covariates) 
#'    {Sample, Condition, Time (numeric), TimeCateg (str)
#'    (and confounding variables if given).}
#'    Annotation table with covariates for each sample.
#' @param boolCaseCtrl (bool) 
#' 		Whether to perform case-control analysis. Does case-only
#' 		analysis if FALSE.
#' @param vecConfounders (vector of strings number of confounding variables)
#' 		Factors to correct for during batch correction. Have to 
#' 		supply dispersion factors if more than one is supplied.
#' 		Names refer to columns in dfAnnotation.
#' @param vecDispersionsExternal (vector length number of
#'    genes in matCountData) [Default NULL]
#'    Externally generated list of gene-wise dispersion factors
#'    which overides DESeq2 generated dispersion factors.
#' @param vecSizeFactorsExternal (vector length number of
#'    cells in matCountData) [Default NULL]
#'    Externally generated list of size factors which override
#'    size factor computation in ImpulseDE2.
#'    
#' @return (list length 4)
#' \itemize{
#'    \item matCountDataProc (matrix genes x samples)
#'    Read count data.
#'    \item dfAnnotationProc (data frame samples x covariates) 
#'    {Sample, Condition, Time (numeric), TimeCateg (str)
#'    (and confounding variables if given).}
#'    Processed annotation table with covariates for each sample.
#'    \item vecSizeFactorsExternalProc (numeric vector number of samples) 
#'    Model scaling factors for each sample which take
#'    sequencing depth into account (size factors).
#'    \item vecDispersionsExternalProc (vector number of genes) Gene-wise 
#'    negative binomial dispersion hyper-parameter.
#'    \item strReportProcessing (str) String of stdout of processData().
#' }
#' 
#' @author David Sebastian Fischer
processData <- function(
    dfAnnotation, 
    matCountData,
    boolCaseCtrl,
    vecConfounders,
    vecDispersionsExternal,
    vecSizeFactorsExternal){
    
    ###############################################################
    # (I) Helper functions
    
    # Check whether object was supplied (is not NULL).
    checkNull <- function(objectInput,strObjectInput){
        if(is.null(objectInput)){
            stop(paste0( "ERROR: ", strObjectInput,
                         " was not given as input." ))
        }
    }
    # Check whether object does not have NA elements.
    checkNA <- function(objectInput,strObjectInput){
        if(is.na(objectInput)){
            stop(paste0( "ERROR: ", strObjectInput,
                         " is NA and needs to be specifief." ))
        }
    }
    # Checks whether dimensions of matrices agree.
    checkDimMatch <-function(
        matInput1, matInput2, strMatInput1, strMatInput2){
        if(any(dim(matInput1)!=dim(matInput2))){
            stop(paste0( "ERROR: ", strMatInput1, 
                         " does not have the dimensions as ", 
                         strMatInput2 , "." ))
        }
    }
    # Checks whether vectors are identical.
    checkElementMatch <- function(vec1, vec2, strVec1, strVec2){
        if(!any(vec1==vec2)){
            stop(paste0( "ERROR: ",strVec1 ," do not agree with ", 
                         strVec2, "." ))
        }
    }
    # Checks whether elements are numeric
    checkNumeric <- function(matInput, strMatInput){
        if(any(!is.numeric(matInput))){
            stop(paste0( "ERROR: ", strMatInput, 
                         " contains non-numeric elements. ",
                         "Requires numeric data." ))
        }
    }
    # Checks whether elements are probabilities (in [0,1]).
    checkProbability <- function(matInput, strMatInput){
        checkNumeric(matInput, strMatInput)
        if(any(matInput < 0 | matInput > 1 | is.na(matInput))){
            stop(paste0( "ERROR: ", strMatInput, 
                         " contains elements outside",
                         " of interval [0,1]." ))
        }
    }
    # Checks whether elements are count data: 
    # non-negative integer finite numeric elements.
    # Note that NA are allowed. 
    # Can be used to check whether element is integer if NA 
    # is checked separately.
    checkCounts <- function(matInput, strMatInput){
        checkNumeric(matInput, strMatInput)
        if(any(matInput[!is.na(matInput)] %% 1 != 0)){
            stop(paste0( "ERROR: ", strMatInput, 
                         " contains non-integer elements.",
                         " Requires count data." ))
        }
        if(any(!is.finite(matInput[!is.na(matInput)]))){
            stop(paste0( "ERROR: ", strMatInput, 
                         " contains infinite elements.",
                         " Requires count data." ))
        }
        if(any(matInput[!is.na(matInput)]<0)){
            stop(paste0( "ERROR: ", strMatInput, 
                         " contains negative elements.",
                         " Requires count data." ))
        }
    }
    
    # Check format and presence of input data.
    checkData <- function(
        dfAnnotation, 
        matCountData,
        boolCaseCtrl,
        vecConfounders,
        vecDispersionsExternal,
        vecSizeFactorsExternal,
        strReportProcessing){
        
        strReportProcessing <- "Processing Details:"
        ### 1. Check that all necessary input was specified
        checkNull(dfAnnotation,"dfAnnotation")
        checkNull(matCountData,"matCountData")
        checkNull(boolCaseCtrl,"boolCaseCtrl")
        
        
        ### 2. Check annotation table content
        ### a) Check column names
        vecColNamesRequired <- c("Sample","Condition","Time",
                                 vecConfounders)
        if( !all(vecColNamesRequired %in% 
                 colnames(dfAnnotation)) ){
            stop(paste0(
                "Could not find column ",
                vecColNamesRequired[
                    !(vecColNamesRequired %in% colnames(dfAnnotation))],
                " in annotation table."))
        }
        ### b) Samples
        # Check that sample name do not occur twice
        if(any(duplicated(dfAnnotation$Sample))){
            stop(paste0(
                "ERROR: [Annotation table] ",
                "Sample names must be unique: Sample(s) ",
                paste0((dfAnnotation$Sample)[
                    duplicated(dfAnnotation$Sample)],collapse=","),
                " is/are duplicated."))
        }
        ### c) Time points
        vecTimepoints <- unique(as.vector( dfAnnotation$Time ))
        # Check that time points are numeric
        checkNumeric(dfAnnotation$Time, "dfAnnotation$Time")
        ### d) Conditions
        lsConditions <- unique( dfAnnotation$Condition )
        # Check that given conditions exisit in annotation table
        if(!("case" %in% lsConditions)){
            stop(paste0(
                "ERROR: Condition \"case\" does",
                " not occur in annotation table condition column."))
        }
        if(boolCaseCtrl){
            if(!("control" %in% lsConditions)){
                stop(paste0(
                    "ERROR: Condition \"control\" does",
                    " not occur in annotation table condition column."))
            }
        }
        ### e) Batch
        if(!is.null(vecConfounders)){
            # Check that model matrix based on confounding variables is full rank
            # 1. Dummy check: Check that number of batches 
            # for each confounding variable is > 1
            for(confounder in vecConfounders){
                if(length(unique( dfAnnotation[,confounder] ))==1){
                    stop(paste0(
                        "ERROR: Model matrix based on confounding ",
                        "variables {", vecConfounders,
                        "} is not full rank: Only one batch",
                        " given for confounder ", confounder, 
                        ". Remove from vecConfounders or correct",
                        " dfAnnotation."))
                }
            }
            # 2. More detailed full rank check
            matModelMatrix <- do.call(cbind, lapply(
                vecConfounders,
                function(confounder){
                    if(is.numeric(dfAnnotation[,confounder])){
                      dfAnnotation[,confounder]
                    } else {
                      match(dfAnnotation[,confounder],
                          unique(dfAnnotation[,confounder]))
                    }
                }))
            if(rankMatrix(matModelMatrix)[1] != dim(matModelMatrix)[2]){
                stop(paste0(
                    "Model matrix based on confounding variables {", 
                    vecConfounders,
                    "} is not full rank. ",
                    "Correct the confounding variables. ",
                    " Note that it is not possible to model",
                    " NESTED confounding variables: ",
                    " Any confounding variables cannot",
                    " be a linear combination of the other ",
                    " confounding variables."))
            }
        }
        
        ### 3. Check expression table content
        # Check that all entries in count data table occur in annotation table
        if( any(!(colnames(matCountData) %in% dfAnnotation$Sample)) ){
            strReportProcessing <- paste0(
                strReportProcessing, "\n","WARNING: The column(s) ",
                paste0(as.character( colnames(matCountData)[
                    !(colnames(matCountData) %in% dfAnnotation$Sample)] ),
                    collapse=","),
                " in the count data table do(es) not occur",
                " in annotation table and will be ignored.")
        }
        checkNull(rownames(matCountData),"[Rownames of matCountData]")
        checkCounts(matCountData, "matCountData")
        
        ### 4. Check supplied dispersion vector
        if( is.null(vecDispersionsExternal) & length(vecConfounders)>1 ){
            stop(paste0(
                "DESeq2 cannot be run in an automated fashion for multiple ",
                "confounding variables. Run DESeq separately and supply ",
                "dispersion parameters through vecDispersionsExternal."))
        }
        if(!is.null(vecDispersionsExternal)){
            # Check that dispersion parameters were named
            if(is.null(names(vecDispersionsExternal))){
                stop(paste0(
                    "ERROR: vecDispersionsExternal was not named.",
                    " Name according to rownames of matCountData."))
            }
            # Check that one dispersion parameter factors was supplied per gene
            if(any( !(names(vecDispersionsExternal) %in% 
                      rownames(matCountData)) ) |
               any( !(rownames(matCountData) %in% 
                      names(vecDispersionsExternal)) ) ) {
                stop(paste0(
                    "ERROR: vecDispersionsExternal supplied but names",
                    " do not agree with rownames of matCountData."))
            }
            # Check that dipsersion parameter vector is numeric
            checkNumeric(vecDispersionsExternal, "vecDispersionsExternal")
            # Check that dispersion parameters are positive
            if(any(vecDispersionsExternal <= 0)){
                stop(paste0(
                    "WARNING: vecDispersionsExternal contains negative",
                    " or zero elements which leads.",
                    "Dispersion parameters must be positive." ))
            }
        }
        
        ### 5. Check supplied size facotrs
        if(!is.null(vecSizeFactorsExternal)){
            # Check that size factors were named
            if(is.null(names(vecSizeFactorsExternal))){
                stop(paste0("ERROR: vecSizeFactorsExternal was not named.",
                            " Name according to colnames of matCountData."))
            }
            # Check that one size factors was supplied per cell
            if(any( !(names(vecSizeFactorsExternal) %in% 
                      colnames(matCountData)) ) |
               any( !(colnames(matCountData) %in% 
                      names(vecSizeFactorsExternal)) ) ){
                stop(paste0(
                    "ERROR: vecSizeFactorsExternal supplied but",
                    " names do not agree with colnames of matCountData."))
            }
            # Check that size factors vector is numeric
            checkNumeric(vecSizeFactorsExternal, "vecSizeFactorsExternal")
            # Check that size factors are positive
            if(any(vecSizeFactorsExternal <= 0)){
                stop(paste0(
                    "WARNING: vecSizeFactorsExternal contains negative ",
                    "or zero elements which leads.",
                    "Size factors must be positive, remove samples if size ",
                    "factor is supposed to be zero." ))
            }
        }
        
        ### Summarise which mode, conditions, samples and
        ### batch were found
        if(boolCaseCtrl) { 
            strReportProcessing <- paste0(
                strReportProcessing,"\n","ImpulseDE2 runs in case-ctrl mode.") 
        } else { 
            strReportProcessing <- paste0(
                strReportProcessing, "\n","ImpulseDE2 runs in case-only mode.") 
        }
        strReportProcessing <- paste0(
            strReportProcessing, "\n",
            paste0( "Found time points: ",
                    paste( vecTimepoints, collapse=",") ))
        for(tp in vecTimepoints){
            strReportProcessing <- 
                paste0(strReportProcessing, "\n",
                       paste0(
                           "Case: Found the samples at time point ", 
                           tp,": ", paste0(dfAnnotation[
                               (dfAnnotation$Time %in% tp) &
                                   (dfAnnotation$Condition %in% "case") &
                                   (dfAnnotation$Sample %in% 
                                        colnames(matCountData)),
                               ]$Sample,collapse=","),collapse="," ) )
        }
        if(boolCaseCtrl){
            for(tp in vecTimepoints){
                strReportProcessing <- paste0(
                    strReportProcessing, "\n", paste0(
                        "Control: Found the following samples at time point ", 
                        tp, ":", paste0(dfAnnotation[
                            dfAnnotation$Time %in% tp &
                                dfAnnotation$Condition %in% "control" &
                                dfAnnotation$Sample %in% colnames(matCountData),
                            ]$Sample,collapse=","),collapse="," ) )
            }
        }
        if(!is.null(vecConfounders)){
            for(confounder in vecConfounders){
                for(batch in unique( dfAnnotation[,confounder] )){
                    strReportProcessing <- paste0(
                        strReportProcessing, "\n", 
                        "Found the following samples for confounder ", 
                        confounder," and batch ", batch, ": ",
                        paste0( dfAnnotation[
                            dfAnnotation[,confounder] %in% batch &
                                dfAnnotation$Sample %in% colnames(matCountData),
                            ]$Sample, collapse=",") )
                }
            }
        }
        
        return(strReportProcessing)
    }
    
    # Add categorial time variable to annotation table which
    # differentiates case and control time points (given to DESeq2).
    # Add column with time scalars with underscore prefix.
    # Add nested batch column if necssary for running DESeq2.
    # Reduce to samples uesed.
    procAnnotation <- function(dfAnnotation,
                               matCountData,
                               boolCaseCtrl,
                               vecConfounders){
        
        # Make sure all columns are not factors
        for(col in seq(1,dim(dfAnnotation)[2])) dfAnnotation[,col] <- 
                as.vector(dfAnnotation[,col])
        
        if(!boolCaseCtrl){
            dfAnnotationProc <- 
                dfAnnotation[dfAnnotation$Condition=="case",]
        } else {
            dfAnnotationProc <- 
                dfAnnotation[dfAnnotation$Condition=="case" |
                                 dfAnnotation$Condition=="control",]
        }
        
        # Reduce to samples which occur in count table
        # This allows to represent entire library of samples in annotation
        # file even if not all samples were measured yet. Not measured samples
        # which are not mentioned in count table are ignored. Thereby, the same
        # full annotation table can be used throughout sample collection on 
        # incomplete data sets.
        dfAnnotationProc <- dfAnnotationProc[dfAnnotationProc$Sample %in% 
                                                 colnames(matCountData),]
        
        # Take out columns which are not used
        dfAnnotationProc <- dfAnnotationProc[,c("Sample", "Time", "Condition", 
                                                vecConfounders)]
        # Add categorial time column for DESeq2
        dfAnnotationProc$TimeCateg <- paste0("_", dfAnnotationProc$Time)
        # Add nested batches for running DESeq2
        if(boolCaseCtrl){
            # Create one new column for each nested batch
            for(confounder in vecConfounders){
                vecBatchesCase <- dfAnnotationProc[
                    dfAnnotationProc$Condition=="case",confounder]
                vecNestedBatchesCase <- paste0(
                    "NestedBatch", match(vecBatchesCase, unique(vecBatchesCase)))
                vecBatchesCtrl <- dfAnnotationProc[
                    dfAnnotationProc$Condition=="control",confounder]
                vecNestedBatchesCtrl <- paste0(
                    "NestedBatch", match(vecBatchesCtrl, unique(vecBatchesCtrl)))
                vecBatchesNested <- array(NA, dim(dfAnnotationProc)[1])
                vecBatchesNested[dfAnnotationProc$Condition=="case"] <- 
                    vecNestedBatchesCase
                vecBatchesNested[dfAnnotationProc$Condition=="control"] <- 
                    vecNestedBatchesCtrl
                strNameConfounderNested <- paste0(confounder, "Nested")
                dfAnnotationProc[strNameConfounderNested] <- 
                    vecBatchesNested
            }
        }
        
        return(dfAnnotationProc)
    }
    
    # Name genes if names are not given.
    nameGenes <- function(matCountDataProc){
        if(is.null(rownames(matCountDataProc))){
            rownames(matCountDataProc) <- 
                paste0("Gene_", seq(1,nrow(matCountDataProc)))
        }
        return(matCountDataProc)
    }
    
    # Reduce count data to data which are utilised later
    reduceCountData <- function(dfAnnotation, matCountDataProc){
        
        strReportCountRed <- paste0(
            "Input contained ",dim(matCountDataProc)[1]," genes/regions.")
        ### 1. Columns (Conditions):
        # Reduce expression table to columns of considered conditions
        if(boolCaseCtrl){
            vecSampleNames_Case <- as.character(as.vector( 
                dfAnnotation[ 
                    dfAnnotation$Condition %in% "case" &
                        dfAnnotation$Sample %in% colnames(matCountDataProc),
                    ]$Sample ))
            vecSampleNames_Ctrl <- as.character(as.vector( 
                dfAnnotation[ 
                    dfAnnotation$Condition %in% "control" &
                        dfAnnotation$Sample %in% colnames(matCountDataProc),
                    ]$Sample ))
            matCountDataProc <- 
                matCountDataProc[,c(vecSampleNames_Case,vecSampleNames_Ctrl)]
        } else {
            vecSampleNames_Case <- as.character(as.vector( 
                dfAnnotation[ 
                    dfAnnotation$Condition %in% "case" &
                        dfAnnotation$Sample %in% colnames(matCountDataProc),
                    ]$Sample ))
            matCountDataProc <- matCountDataProc[,vecSampleNames_Case]
        }
        # Check that every sample contains at least one observed value (not NA)
        vecboolNASample <- apply(matCountDataProc,2,function(sample)
            all(is.na(sample)) )
        if(any(vecboolNASample)){
            strReportCountRed <- paste0(
                strReportCountRed,  "\n","WARNING: Sample(s) ",
                paste0(colnames(matCountDataProc)[vecboolNASample], 
                       collapse=","),
                " only contain(s) NA values and will be ",
                "removed from the analysis.")
            matCountDataProc <- matCountDataProc[,!vecboolNASample]
        }
        # Trigger warning if all zero sample is encountered - these are kept!
        vecboolAllZeroSample <- apply(matCountDataProc,2,function(sample) 
            all(sample[!is.na(sample)]==0) )
        if(any(vecboolAllZeroSample)){
            strReportCountRed <- paste0(
                strReportCountRed,  "\n", "WARNING: Sample(s) ",
                paste0(colnames(matCountDataProc)[vecboolAllZeroSample], 
                       collapse=","),
                " only contain(s) zeros (and NAs).",
                " These samples are kept for analysis.")
        }
        
        ### 2. Rows (Genes):
        # Exclude genes with only missing values (NAs)
        vecboolNonZeroGene <- apply(matCountDataProc,1,function(gene)  
            any(!is.na(gene) & gene > 0) )
        matCountDataProc <- matCountDataProc[vecboolNonZeroGene,]
        if(sum(!vecboolNonZeroGene) > 0){
            strReportCountRed <- paste0(
                strReportCountRed,  "\n","WARNING: ",
                sum(!vecboolNonZeroGene), " out of ",
                length(vecboolNonZeroGene),
                " genes do not have obserserved non-zero counts",
                " and are excluded.")
        }
        
        # Sort count matrix column by annotation table
        matCountDataProc <- matCountDataProc[
            ,match(as.vector(dfAnnotation$Sample),colnames(matCountDataProc))]
        
        strReportCountRed <- paste0(strReportCountRed,  
                                    "\n","Selected ",dim(matCountDataProc)[1],
                                    " genes/regions for analysis.")
        return(list( matCountDataProc=matCountDataProc,
                     strReportCountRed=strReportCountRed) )
    }
    
    ###############################################################
    # (II) Main body of function
    
    # Check validity of input
    strReportProcessing <- checkData(
        dfAnnotation=dfAnnotation,
        matCountData=matCountData,
        boolCaseCtrl=boolCaseCtrl,
        vecConfounders=vecConfounders,
        vecDispersionsExternal=vecDispersionsExternal,
        vecSizeFactorsExternal=vecSizeFactorsExternal )
    
    # Process annotation table
    dfAnnotationProc <- procAnnotation(dfAnnotation=dfAnnotation,
                                       matCountData=matCountData,
                                       boolCaseCtrl=boolCaseCtrl,
                                       vecConfounders=vecConfounders)
    
    # Process raw counts
    matCountDataProc <- nameGenes(matCountDataProc=matCountData)
    lsReduceCounts <- reduceCountData(
        dfAnnotation=dfAnnotationProc, 
        matCountDataProc=matCountDataProc)
    matCountDataProc <- lsReduceCounts$matCountDataProc
    strReportProcessing <- paste0(strReportProcessing, "\n",
                                  lsReduceCounts$strReportCountRed)
    
    # Reduce externally provided parameters according to reduced data set
    # and reorder according to given data set.
    if(!is.null(vecSizeFactorsExternal)){
        vecSizeFactorsExternalProc <- 
            vecSizeFactorsExternal[colnames(matCountDataProc)]
    } else { vecSizeFactorsExternalProc <-NULL }
    if(!is.null(vecDispersionsExternal)){
        vecDispersionsExternalProc <- 
            vecDispersionsExternal[rownames(matCountDataProc)]
    } else { vecDispersionsExternalProc <-NULL }
    
    return( list(matCountDataProc           = matCountDataProc,
                 dfAnnotationProc           = dfAnnotationProc,
                 vecSizeFactorsExternalProc = vecSizeFactorsExternalProc,
                 vecDispersionsExternalProc = vecDispersionsExternalProc,
                 strReportProcessing        = strReportProcessing ) )
}