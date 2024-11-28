#!/bin/bash 

#SBATCH -p multiple_il
#SBATCH -N2
#SBATCH --ntasks 40
#SBATCH --time 04:00:00 
#SBATCH -J Y_C4_8SV_Case

singularity exec als_brain_multiome_base.sif Rscript Y_C4_8_DA_ATAC_WNN_SVA_AllCase_Slurm.R

