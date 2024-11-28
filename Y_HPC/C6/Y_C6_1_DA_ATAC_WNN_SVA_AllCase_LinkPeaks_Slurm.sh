#!/bin/bash 

#SBATCH -p multiple_il
#SBATCH -N2
#SBATCH --ntasks 40
#SBATCH --time 18:30:00 
#SBATCH -J Y_C6_1_DA_SVA_AllCase_LinkPeaks

singularity exec als_brain_multiome_latest.sif Rscript Y_C6_1_DA_ATAC_WNN_SVA_AllCase_LinkPeaks_Slurm.R 

