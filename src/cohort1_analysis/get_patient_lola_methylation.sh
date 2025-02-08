#!/bin/bash
#SBATCH --mem=8000
#SBATCH --array=1-768
#SBATCH --output=logs/%x.%3a.log
#SBATCH --error=logs/%x.%3a.err
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1 

file=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/cpgs.txt)

Rscript src/get_lola_methylation.R $file results
