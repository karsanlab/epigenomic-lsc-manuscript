#!/bin/bash

# INPUT:
# results/${sample}_mapped.PT.bam

# OUTPUT:
# results/${sample}.macs2_peaks.narrowPeak

# Activate conda environment. 
# Environment specs to recreate the environment are available in the resources folder as "hichip_conda_environment_spec_file.txt"

eval "$(conda shell.bash hook)" 
conda activate hichip_conda_environment

# Call sample peaks using macs2

samtools view -h -F 0x900 results/${sample}_mapped.PT.bam | bedtools bamtobed -i stdin > results/${sample}.primary.aln.bed
macs2 callpeak -t results/${sample}.primary.aln.bed -n results/${sample}.macs2

