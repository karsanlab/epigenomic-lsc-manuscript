#!/bin/bash

# INPUT:
# results/${sample}_mapped.PT.bam

# OUTPUT:
# results/${sample}_mapped.PT.bam.RPGC.bw

# Source conda environment

eval "$(conda shell.bash hook)" #lets me open conda environ from within shell script
conda activate hichip_conda_environment

## Normalize to 1x coverage using RPGC instead of RPKM and rerun heatmap plots

bamCoverage --bam results/${sample}_mapped.PT.bam  -o results/${sample}_mapped.PT.bam.RPGC.bw \
    --binSize 50
    --normalizeUsing RPGC
    --effectiveGenomeSize 2913022398
    --ignoreForNormalization chrX

