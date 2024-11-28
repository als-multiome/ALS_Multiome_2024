#!/bin/bash 

#SBATCH -p multiple_il
#SBATCH -N2
#SBATCH --time 04:30:00 
#SBATCH -J B2_6_DE_ALSFTD_C9vsNonC9
#SBATCH -o Y_B2_6_RNA_DE_Slurm.out

singularity exec als_brain_multiome_base.sif Rscript Y_B2_6_RNA_DE_Slurm.R 

