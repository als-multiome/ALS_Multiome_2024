#!/bin/bash 

#SBATCH -p multiple_il
#SBATCH -N2
#SBATCH --ntasks 40
#SBATCH --time 06:30:00 
#SBATCH -J Y_C2_3 

singularity exec als_brain_multiome_base.sif Rscript Y_C2_3_ATAC_DA_Slurm.R

