#!/bin/bash

# INPUT:
# results/${sample}_shControl_mapped.PT.bam.RPGC.bw
# results/${sample}_shZNF143_mapped.PT.bam.RPGC.bw
# results/CTCF_only_binding_sites.bed
# results/ZNF143_only_binding_sites.bed
# results/CTCF_ZNF143_cobinding_sites.bed

# OUTPUT:
# results/matrix.RPGC.cobind.mat.gz

# Source conda environment
eval "$(conda shell.bash hook)" 
conda activate hichip_conda_environment

## Generate matrix and plots for scaled bigwigs with binding sites divided by overlaps

computeMatrix scale-regions \
    -S results/${sample}_shControl_mapped.PT.bam.RPGC.bw \
    results/${sample}_shZNF143_mapped.PT.bam.RPGC.bw
    -R results/CTCF_ZNF143_cobinding_sites.bed results/CTCF_only_binding_sites.bed results/ZNF143_only_binding_sites.bed \
    --beforeRegionStartLength 4000 \
    --afterRegionStartLength 4000 \
    --regionBodyLength 500 \
    --skipZeros -o results/matrix.RPGC.cobind.mat.gz

