#!/bin/bash 

#SBATCH -p multiple_il
#SBATCH -N2
#SBATCH --time 04:30:00 
#SBATCH -J Y_C2_4

singularity exec als_brain_multiome_base.sif Rscript Y_C2_4_ATAC_DA_Slurm.R

