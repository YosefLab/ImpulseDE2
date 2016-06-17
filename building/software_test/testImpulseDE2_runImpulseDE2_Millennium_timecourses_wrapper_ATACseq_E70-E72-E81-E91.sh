#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=16
### redirect stdout/stderr 
#PBS -e localhost:/data/yosef/users/fischerd/code/script_reports/160511_runImpulseDE2_E70-E71-E81-E91_timecourses.err
#PBS -o localhost:/data/yosef/users/fischerd/code/script_reports/160511_runImpulseDE2_E70-E71-E81-E91_timecourses.out

#PBS -m ae
#PBS -M fischerd@berkeley.edu

### ############################################################################

Rscript /data/yosef2/users/fischerd/code/ImpulseDE2/testImpulseDE2_noControl_timecourses_ATACseq_E70-E72-E81-E91_Millennium.R