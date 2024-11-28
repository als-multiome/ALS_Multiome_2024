#!/bin/bash 

#SBATCH -p fat
#SBATCH -N1
#SBATCH --ntasks 40
#SBATCH --time 18:30:00 
#SBATCH -J Y_C3_1_DA_SVA_AllCase

singularity exec als_brain_multiome_base.sif Rscript Y_C3_1_DA_ATAC_WNN_SVA_AllCase_Slurm.R

