#!/bin/sh

# INPUT:
# results/${sample}_stats.txt

# OUTPUT:
# results/${sample}_stats_get_qc.txt

# Activate conda environment. 
# Environment specs to recreate the environment are available in the resources folder as "hichip_conda_environment_spec_file.txt"

eval "$(conda shell.bash hook)" 
conda activate hichip_conda_environment

## Run Dovetail genomics QC script

get_qc.py -p $SCRATCH_DIR/$SAMPLE_NAME/${SAMPLE_NAME}_stats.txt # print to screen
get_qc.py -p $SCRATCH_DIR/$SAMPLE_NAME/${SAMPLE_NAME}_stats.txt > results/${sample}_stats_get_qc.txt # write to file