#!/bin/bash 

#SBATCH -p multiple
#SBATCH -N2
#SBATCH --ntasks 40
#SBATCH --time 04:00:00 
#SBATCH -J Y_C7_7SV_Case

singularity exec als_brain_multiome_latest.sif Rscript Y_C7_7_DA_ATAC_WNN_SVA_AllCase_LinkPeaks_Slurm.R

