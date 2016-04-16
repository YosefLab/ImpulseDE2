print("Run after completion of DESeq")

# Load files from interior of ImpulseDE
setwd( "/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/software_test_out")
load("ImpulseDE_expression_array.RData")
load("ImpulseDE_prepared_annotation.RData")

# Compare against DESeq
# Make summary table to compare against plots
dfDESeq_Impulse <- as.data.frame( cbind(
  "Gene"=as.character(as.vector( ImpulseDE_results$Gene )),
  "DESeq"=(dds_resultsTable[as.character(as.vector( ImpulseDE_results$Gene )),])$padj,
  "Impulse"=(ImpulseDE_results[as.character(as.vector( ImpulseDE_results$Gene )),])$adj.p),
  stringsAsFactors=FALSE)
rownames(dfDESeq_Impulse) <- as.character(as.vector( ImpulseDE_results$Gene ))
dfDESeq_Impulse$DESeq <- as.numeric(dfDESeq_Impulse$DESeq)
# Set NA as not detected:
dfDESeq_Impulse$DESeq[is.na(dfDESeq_Impulse$DESeq)] <- 1
dfDESeq_Impulse$Impulse <- as.numeric(dfDESeq_Impulse$Impulse)

mat_overlap <- array(NA,c(10,10))
mat_intersect <- array(NA,c(10,10))
mat_union <- array(NA,c(10,10))
# DESeq on vertical
for(i in 1:10){
  # Impulse on horicontal
  for(j in 1:10){
    sig_DESeq <- dfDESeq_Impulse$DESeq <= 10^(-i)
    sig_Impulse <- dfDESeq_Impulse$Impulse <= 10^(-j)
    mat_overlap[i,j] <- sum(sig_DESeq & sig_Impulse)/sum(sig_DESeq | sig_Impulse)
    mat_intersect[i,j] <- sum(sig_DESeq & sig_Impulse)
    mat_union[i,j] <- sum(sig_DESeq | sig_Impulse)
  }
}
rownames(mat_overlap) <- -1:-10
colnames(mat_overlap) <- -1:-10
rownames(mat_intersect) <- -1:-10
colnames(mat_intersect) <- -1:-10
rownames(mat_union) <- -1:-10
colnames(mat_union) <- -1:-10
library(gplots)
graphics.off()
heatmap(mat_overlap, keep.dendro = FALSE,Rowv=NA,Colv= "Rowv",symm=FALSE,
  xlab =  paste0("ImpulseDE scores"), ylab = "DESeq scores")
breaks <- seq(0,1,by=0.01)
hm.colors <- colorpanel( length(breaks)-1, "yellow", "red" )
graphics.off()
pdf(paste('DESeq-Impulse_OverlapSigGenes_Heatmap.pdf',sep=''),width=7,height=7)
heatmap.2(mat_overlap, dendrogram="none", Rowv=FALSE,Colv=FALSE, 
  xlab =  paste0("log(p-value) ImpulseDE"), ylab = "log(p-value) DESeq2",
  breaks=breaks,col=hm.colors, scale="none",
  trace="none",density.info="none",
  key.title = " ", key.xlab = paste0("Overlap"), key.ylab = NULL,
  symkey=FALSE,
  cellnote=round(mat_overlap,digits=2),notecol="white",
  lmat=rbind( c(3,4),c(2,1) ),lhei=c(1,4), lwid=c(1,4), margins=c(5,5))
dev.off()
print("Intersection")
print(mat_intersect)
print("Union")
print(mat_union)

Q = 10^(-3)
DEgenes_both <- dfDESeq_Impulse[(dfDESeq_Impulse[,"DESeq"] < Q & dfDESeq_Impulse[,"Impulse"] < Q),"Gene"]
DEgenes_DESeq_only <- dfDESeq_Impulse[(dfDESeq_Impulse[,"DESeq"] < Q & dfDESeq_Impulse[,"Impulse"] >= Q),"Gene"]
DEgenes_Impulse_only <- dfDESeq_Impulse[(dfDESeq_Impulse[,"DESeq"] >= Q & dfDESeq_Impulse[,"Impulse"] < Q),"Gene"]
graphics.off()
plot_impulse(DEgenes_both,
  expression_array, prepared_annotation, impulse_fit_genes,
  control_timecourse, control_name, case_ind, file_name_part = "DE_DESeqAndImpulse",
  title_line = "", sub_line = "",ImpulseDE_res=ImpulseDE_results,DESeq2_res=dds_resultsTable)
plot_impulse(DEgenes_DESeq_only,
  expression_array, prepared_annotation, impulse_fit_genes,
  control_timecourse, control_name, case_ind, file_name_part = "DE_DESeq_only",
  title_line = "", sub_line = "",ImpulseDE_res=ImpulseDE_results,DESeq2_res=dds_resultsTable)
plot_impulse(DEgenes_Impulse_only,
  expression_array, prepared_annotation, impulse_fit_genes,
  control_timecourse, control_name, case_ind, file_name_part = "DE_Impulse_only",
  title_line = "", sub_line = "",ImpulseDE_res=ImpulseDE_results,DESeq2_res=dds_resultsTable)