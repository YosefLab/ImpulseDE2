#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=16
### redirect stdout/stderr 
#PBS -e localhost:/data/yosef2/users/fischerd/code/scriptreports/160610_RNAseqComAnalysis.err
#PBS -o localhost:/data/yosef2/users/fischerd/code/scriptreports/160610_RNAseqComAnalysis.out

#PBS -m ae
#PBS -M fischerd@berkeley.edu

### ############################################################################
Rscript /data/yosef2/users/fischerd/data/RNAseq_Jovanovic/comparative_analysis/ComparativeAnalysisCluster.R