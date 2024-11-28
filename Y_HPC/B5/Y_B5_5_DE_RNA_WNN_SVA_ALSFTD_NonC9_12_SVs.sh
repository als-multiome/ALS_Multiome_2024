#!/bin/bash 

#SBATCH -p single
#SBATCH -N1
#SBATCH --ntasks 40
#SBATCH --mem 120000 
#SBATCH --ntasks-per-node 40
#SBATCH --time 08:00:00 
#SBATCH -J Y_B5_5_12SVs_ALSFTD_NonC9   

singularity exec als_brain_multiome_latest.sif Rscript Y_B5_5_DE_RNA_WNN_SVA_ALSFTD_NonC9_12_SVs.R 
