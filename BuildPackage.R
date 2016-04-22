# 1. Build vignette
rm(list=ls())
# Move only .R code file into directory code_files, i.e. ImpulseDE2
setwd('/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/code_files')
package.skeleton(name="ImpulseDE2", path = "./..",
  code_files=c("ImpulseDE2_main.R","srcImpulseDE2_calcImpulse.R"
    ,"srcImpulseDE2_computePval.R","srcImpulseDE2_CostFunctionsFit.R",
    "srcImpulseDE2_fitImpulse.R","srcImpulseDE2_plotDEGenes.R",
    "srcImpulseDE2_processData.R","srcImpulseDE2_runDESeq2.R"))

setwd('/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/ImpulseDE2')
roxygen2::roxygenise()

# Delete _comp files from man, dont include in vignette
rm /Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/ImpulseDE2/man/evalLogLikImpulse_comp.Rd
rm /Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/ImpulseDE2/man/evalLogLikMean_comp.Rd
rm /Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/ImpulseDE2/man/calcImpulse_comp.Rd
R CMD Rd2pdf --pdf --title='ImpulseDE2 Vignette' \
-o /Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/ImpulseDE2/ImpulseDE2_vignette.pdf \
/Users/davidsebastianfischer/MasterThesis/code/ImpulseDE/building/ImpulseDE2/man/*.Rd