#!/bin/bash
#SBATCH --output=logs/%x.%3a.log
#SBATCH --error=logs/%x.%3a.err
#SBATCH --mail-type=ALL
#SBATCH --mem=0 
#SBATCH --cpus-per-task=1 
#SBATCH --array=1-626 

Rscript src/get_patient_farlik_regional_methylation.R $SLURM_ARRAY_TASK_ID

