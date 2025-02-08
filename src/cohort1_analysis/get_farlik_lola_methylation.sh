#!/bin/bash
#SBATCH --mem=8000
#SBATCH --array=1-122
#SBATCH --output=logs/%x.%3a.log
#SBATCH --error=logs/%x.%3a.err
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1 

# INPUT: 
# `data/farlik/GSM*.cons.Cmethyl.CpG.txt` (only HSPCs)

# OUTPUT: 
# `results/GSM*.cons.Cmethyl_tfbs.csv`

file=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/farlik.txt)

Rscript src/get_lola_methylation.R $file results
