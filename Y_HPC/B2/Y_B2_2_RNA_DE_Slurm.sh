#!/bin/bash 

#SBATCH -p multiple_il
#SBATCH -N2
#SBATCH --time 04:30:00 
#SBATCH -J Y_B2_2

singularity exec als_brain_multiome_base.sif Rscript Y_B2_2_RNA_DE_Slurm.R

