#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++     Compute standard deviations matrix for genes x timepoints matrix   +++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Compute standard deviations matrix for genes x timepoints matrix.
### Currently based on MLE of pointwise stdvs.
### DAVID: look into using DESeq here

# INPUT:
#   observed_data: (Numeric 3D array genes x samples x replicates).
#       Contains expression values or similar locus-specific read-outs.
# OUTPUT:
#   stdv_mat: (Numeric 2D array genes x samples) Standard deviations over
#             replicates.

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

compute_stdvs <- function(observed_data = NULL){
  
  # Check that observed data contains replicates
  ReplicatesIdentical <- TRUE
  if(!is.na(dim(observed_data)[3])){
    if(dim(observed_data)[3]>1){
      for(i in 1:dim(observed_data)[3]-1){
        if(!identical(observed_data[,,1],observed_data[,,i+1])){
          ReplicatesIdentical <- FALSE
          break
        }
      }
    }
  }
  
  # Only evaluate sd if replicates present and replicates not identical
  if(!is.na(dim(observed_data)[3]) && dim(observed_data)[3]>1 && !ReplicatesIdentical){
    stdv_mat <- apply(observed_data,c(1,2),sd)
    
    min_sd <- min(unique(c(stdv_mat))[unique(c(stdv_mat))>0])
    stdv_mat[stdv_mat==0] <- min_sd
  # Else, set all standard deviation to 1, they won't impact fitting.
  }else{
    print("WARNING: Supplied 2D data matrix or identical replicates, standard deviations set to 1.")
    stdv_mat <- matrix(1,dim(observed_data)[1],dim(observed_data)[2])
    rownames(stdv_mat) <- rownames(observed_data)
    colnames(stdv_mat) <- colnames(observed_data)
  }
  
  return(stdv_mat)
}