#!/bin/bash 

#SBATCH -p multiple_il
#SBATCH -N2
#SBATCH --time 03:30:00 
#SBATCH -J Y_B2_3

singularity exec als_brain_multiome_base.sif Rscript Y_B2_3_RNA_DE_Slurm.R
