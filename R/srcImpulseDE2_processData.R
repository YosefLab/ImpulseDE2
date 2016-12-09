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
#' @param vecDispersionsExternal: (vector length number of
#'    genes in matCountData) [Default NULL]
#'    Externally generated list of gene-wise dispersion factors
#'    which overides DESeq2 generated dispersion factors.
#' @param vecSizeFactorsExternal: (vector length number of
#'    cells in matCountData) [Default NULL]
#'    Externally generated list of size factors which override
#'    size factor computation in ImpulseDE2.
#'    
#' @return (list length 3) with the following elements:
#' \itemize{
#'  \item \code{matCountDataProc}: (count matrix  genes x samples) 
#'      Count data: Reduced version of \code{matCountData}. 
#' }
#'
#' @export

processData <- function(dfAnnotation, 
												matCountData,
												boolCaseCtrl,
												vecConfounders,
												lsvecDispersionsExternal,
												vecSizeFactorsExternal){
	
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
	checkData <- function(dfAnnotation, 
												matCountData,
												boolCaseCtrl,
												vecConfounders,
												lsvecDispersionsExternal,
												vecSizeFactorsExternal){
		
		### 1. Check that all necessary input was specified
		checkNull(dfAnnotation,"dfAnnotation")
		checkNull(matCountData,"matCountData")
		checkNull(boolCaseCtrl,"boolCaseCtrl")
		
		
		### 2. Check annotation table content
		### a) Check column names
		vecColNamesRequired <- c("Sample","Condition","Time",vecConfounders)
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
		if(!("case" %in% lsConditions)){
			stop("ERROR: Condition \"case\" does not occur in annotation table condition column.")
		}
		if(boolCaseCtrl){
			if(!("control" %in% lsConditions)){
				stop("ERROR: Condition \"control\" does not occur in annotation table condition column.")
			}
		}
		### e) Batch
		if(!is.null(vecConfounders)){
		  # Check that model matrix based on confounding variables is full rank
		  # 1. Dummy check: Check that number of batches for each confounding variable is > 1
			for(confounder in vecConfounders){
				if(length(unique( dfAnnotation[,confounder] ))==1){
				  stop(paste0("ERROR: Model matrix based on confounding variables {", vecConfounders,
				              "} is not full rank: Only one batch given for confounder ", 
				              confounder, ". Remove from vecConfounders or correct dfAnnotation."))
				}
			}
		  # 2. More detailed full rank check
		  matModelMatrix <- do.call(cbind, lapply(vecConfounders, function(confounder){
		    match(dfAnnotation[,confounder], unique(dfAnnotation[,confounder]))
		  }))
		  if(rankMatrix(matModelMatrix)[1] != dim(matModelMatrix)[2]){
		    stop(paste0("Model matrix based on confounding variables {", vecConfounders,
		                "} is not full rank. Correct the confounding variables. ",
		                " Note that it is not possible to model NESTED confounding variables: ",
		                " Any confounding variables cannot be a linear combination of the other ",
		                " confounding variables."))
		  }
		}
		
		### 3. Check expression table content
		# Check that all entries in count data table occur in annotation table
		if( any(!(colnames(matCountData) %in% dfAnnotation$Sample)) ){
			print(paste0("WARNING: The column(s) ",
									 paste0(as.character( colnames(matCountData)[
									 	!(colnames(matCountData) %in% dfAnnotation$Sample)] ),collapse=","),
									 " in the count data table do(es) not occur in annotation table and will be ignored."))
		}
		checkNull(rownames(matCountData),"[Rownames of matCountData]")
		checkCounts(matCountData, "matCountData")
		
		### 4. Check supplied dispersion vector
		if( is.null(lsvecDispersionsExternal) & length(vecConfounders>1) ){
		  stop(paste0("DESeq2 cannot be run in an automated fashion for multiple ",
		              "confounding variables. Run DESeq separately and supply ",
		              "dispersion parameters through lsvecDispersionsExternal.")
		}
		if(!is.null(lsvecDispersionsExternal)){
			# Check that dispersions were named
			if( ( !boolCaseCtrl & any(c("case")!=names(lsvecDispersionsExternal)) ) |
			    ( boolCaseCtrl & any(c("case", "control", "combined")!=names(lsvecDispersionsExternal)) ) ){
				stop(paste0("ERROR: lsvecDispersionsExternal was not named. Supply as: ",
				            "Case-only: Supply list(case=...), ",
				            "Case-control: Supply list(case=..., control=..., combined=...).",
				            "Note that elements of vectors have to be named as the rows in matCountData."))
			}
		  # Check vectors for each condition
		  if(!boolCaseCtrl){ vecLables <- c("case")
		  } else { vecLables <- c("case", "control", "combined") }
		  for(label in vecLabels){
		    # Check that one dispersion was supplied per gene
		    if( any(names(lsvecDispersionsExternal$label) != rownames(matCountData)) ){
		      stop("ERROR: Names of lsvecDispersionsExternal$" , label ," do not agree with rownames of matCountData.")
		    }
		    # Check that dispersion vector is numeric
		    checkNumeric(lsvecDispersionsExternal$label, paste0("lsvecDispersionsExternal$", label))
		    # Check that dispersions are positive (should not have sub-poissonian noise in count data)
		    if(any(lsvecDispersionsExternal$label < 0)){
		      warning(paste0( "WARNING: lsvecDispersionsExternal$" , label ,
		                      " contains negative elements which corresponds to sub-poissonian noise.",
		                      "These elements are kept as they are in the following." ))
		    }
		  }
		}
		
		### 5. Check supplied size facotrs
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
		
		### Summarise which mode, conditions, samples and
		### batch were found
		if(boolCaseCtrl){ print(paste0( "ImpulseDE2 runs in case-ctrl mode." )) }
		else{ print(paste0( "ImpulseDE2 runs in case-only mode." )) }
		print(paste0( "Found time points: ",
									paste( vecTimepoints, collapse=",") ))
		for(tp in vecTimepoints){
			print(paste0( "Case: Found the samples at time point ", 
										tp,": ", 
										paste0(dfAnnotation[
											(dfAnnotation$Time %in% tp) &
												(dfAnnotation$Condition %in% "case") &
												(dfAnnotation$Sample %in% colnames(matCountData)),
											]$Sample,collapse=","),collapse="," ))
		}
		if(boolCaseCtrl){
			for(tp in vecTimepoints){
				print(paste0( "Control: Found the following samples at time point ", 
											tp, ":", 
											paste0(dfAnnotation[
												dfAnnotation$Time %in% tp &
													dfAnnotation$Condition %in% "control" &
													dfAnnotation$Sample %in% colnames(matCountData),
												]$Sample,collapse=","),collapse="," ))
			}
		}
		if(!is.null(vecConfounders)){
			for(confounder in vecConfounders){
				for(batch in unique( dfAnnotation[,confounder] )){
					print(paste0( "Found the following samples for confounder ", confounder,
												" and batch ", batch, ": ",
												paste0( dfAnnotation[
													dfAnnotation[,confounder] %in% batch &
														dfAnnotation$Sample %in% colnames(matCountData),
													]$Sample,collapse=",") ))
				}
			}
		}
		
		return(NULL)
	}
	
	# Add categorial time variable to annotation table which
	# differentiates case and control time points (given to DESeq2).
	# Add column with time scalars with underscore prefix.
	procAnnotation <- function(dfAnnotation,
														 matCountData,
														 boolCaseCtrl,
														 vecConfounders){
		
		if(!boolCaseCtrl){
			dfAnnotationProc <- dfAnnotation[dfAnnotation$Condition=="case",]
		} else {
			dfAnnotationProc <- dfAnnotation[dfAnnotation$Condition=="case" | 
																			 	dfAnnotation$Condition=="control",]
		}

		# Reduce to samples which occur in count table
		# This allows to represent entire library of samples in annotation
		# file even if not all samples were measured yet. Not measured samples
		# which are not mentioned in count table are ignored. Thereby, the same
		# full annotation table can be used throughout sample collection on 
		# incomplete data sets.
		dfAnnotationProc <- dfAnnotationProc[dfAnnotationProc$Sample %in% colnames(matCountData),]
		
		# Take out columns which are not used
		dfAnnotationProc <- dfAnnotationProc[,c("Sample", "Time", "Condition", vecConfounders)]
		# Add categorial time column for DESeq2
		dfAnnotationProc$TimeCateg <- paste0("_", dfAnnotationProc$Time)
		
		return(dfAnnotationProc)
	}
	
	# Name genes if names are not given.
	nameGenes <- function(matCountDataProc){
		if(is.null(rownames(matCountDataProc))){
			rownames(matCountDataProc) <- paste0("Gene_", seq(1,nrow(matCountDataProc)))
		}
		return(matCountDataProc)
	}
	
	# Reduce count data to data which are utilised later
	reduceCountData <- function(dfAnnotation, matCountDataProc){
		
		print(paste0("Input contained ",dim(matCountDataProc)[1]," genes/regions."))
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
			matCountDataProc <- matCountDataProc[,c(vecSampleNames_Case,vecSampleNames_Ctrl)]
		} else {
			vecSampleNames_Case <- as.character(as.vector( 
				dfAnnotation[ 
					dfAnnotation$Condition %in% "case" &
						dfAnnotation$Sample %in% colnames(matCountDataProc),
					]$Sample ))
			matCountDataProc <- matCountDataProc[,vecSampleNames_Case]
		}
		# Check that every sample contains at least one observed value (not NA)
		vecboolNASample <- apply(matCountDataProc,2,function(sample){all(is.na(sample))})
		if(any(vecboolNASample)){
			print(paste0( "WARNING: Sample(s) ",
										paste0(colnames(matCountDataProc)[vecboolNASample], collapse=","),
										" only contain(s) NA values and will be removed from the analysis."))
			matCountDataProc <- matCountDataProc[,!vecboolNASample]
		}
		# Trigger warning if all zero sample is encountered - these are kept!
		vecboolAllZeroSample <- apply(matCountDataProc,2,function(sample){all(sample[!is.na(sample)]==0)})
		if(any(vecboolAllZeroSample)){
			print(paste0( "WARNING: Sample(s) ",
										paste0(colnames(matCountDataProc)[vecboolAllZeroSample], collapse=","),
										" only contain(s) zeros (and NAs). These samples are kept for analysis."))
		}
		### 2. Rows (Genes):
		# Exclude genes with only missing values (NAs)
		vecboolNAGene <- apply(matCountDataProc,1,function(gene){all(is.na(gene))})
		if(any(vecboolNAGene)){
			print(paste0("WARNING: Excluded ",sum(vecboolNAGene)," genes because no real valued samples were given."))
		}
		matCountDataProc <- matCountDataProc[!vecboolNAGene,]
		# Reduce expression table to rows containing at least one non-zero count
		vecboolAllZeroGene <- apply(matCountDataProc,1,function(gene){all(gene[!is.na(gene)]==0)})
		if(sum(vecboolAllZeroGene) > 0){
			print(paste0("WARNING: ",sum(vecboolAllZeroGene), " out of ",
									 dim(matCountDataProc)[1],
									 " genes had only zero counts in all considered samples and are omitted."))
			matCountDataProc <- matCountDataProc[!vecboolAllZeroGene,]
		}
		
		# Sort count matrix column by annotation table
		matCountDataProc <- matCountDataProc[,match(as.vector(dfAnnotation$Sample),colnames(matCountDataProc))]
		
		print(paste0("Selected ",dim(matCountDataProc)[1]," genes/regions for analysis."))
		return(matCountDataProc)
	}
	
	###############################################################
	# (II) Main body of function
	
	# Check validity of input
	checkData(
		dfAnnotation=dfAnnotation,
		matCountData=matCountData,
		boolCaseCtrl=boolCaseCtrl,
		vecConfounders=vecConfounders,
		lsvecDispersionsExternal=lsvecDispersionsExternal,
		vecSizeFactorsExternal=vecSizeFactorsExternal )
	
	# Process annotation table
	dfAnnotationProc <- procAnnotation(dfAnnotation=dfAnnotation,
																		 matCountData=matCountData,
																		 boolCaseCtrl=boolCaseCtrl,
																		 vecConfounders=vecConfounders)
	
	# Process raw counts
	matCountDataProc <- nameGenes(matCountDataProc=matCountData)
	matCountDataProc <- reduceCountData(
		dfAnnotation=dfAnnotationProc, 
		matCountDataProc=matCountDataProc)
	
	# Reduce externally provided parameters according to reduced data set
	# and reorder according to given data set.
	if(!is.null(vecSizeFactorsExternal)){
		vecSizeFactorsExternalProc <- vecSizeFactorsExternal[colnames(matCountDataProc)]
	} else { vecSizeFactorsExternalProc <-NULL }
	
	return( list(matCountDataProc=matCountDataProc,
							 dfAnnotationProc=dfAnnotationProc,
							 vecSizeFactorsExternalProc=vecSizeFactorsExternalProc) )
}