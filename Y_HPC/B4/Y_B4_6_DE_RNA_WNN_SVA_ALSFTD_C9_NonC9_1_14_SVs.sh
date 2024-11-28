#!/bin/bash 

#SBATCH -p multiple_il
#SBATCH -N2
#SBATCH --ntasks 40
#SBATCH --time 24:00:00 
#SBATCH -J Y_B4_6_14_SVs_ALSFTD_C9_NonC9 

singularity exec als_brain_multiome_latest.sif Rscript Y_B4_6_DE_RNA_WNN_SVA_ALSFTD_C9_NonC9_1_14_SVs.R 

