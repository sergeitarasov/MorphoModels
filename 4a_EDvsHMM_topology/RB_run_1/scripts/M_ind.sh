#!/bin/bash -l
#SBATCH --job-name=Mind
#SBATCH --account=project_2002706
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=40:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10000

module load revbayes

srun rb Morpho_TA_SMM_ind-init.Rev > log_SMM_ind.log
