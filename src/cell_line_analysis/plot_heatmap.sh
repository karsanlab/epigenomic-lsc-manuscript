#!/bin/bash

# INPUT:
# results/${sample}_shControl_mapped.PT.bam.RPGC.bw
# results/${sample}_shZNF143_mapped.PT.bam.RPGC.bw
# results/CTCF_only_binding_sites.bed
# results/ZNF143_only_binding_sites.bed
# results/CTCF_ZNF143_cobinding_sites.bed

# OUTPUT:
# plots/Coverage_heatmap_reference_point_RPGC_cobind.png

# Source conda environment
eval "$(conda shell.bash hook)" 
conda activate hichip_conda_environment

# Compute matrix using the centre of the TFBS as the reference point
computeMatrix reference-point \
    -S results/${sample}_shControl_mapped.PT.bam.RPGC.bw \
    results/${sample}_shZNF143_mapped.PT.bam.RPGC.bw \
    -R results/CTCF_ZNF143_cobinding_sites.bed results/CTCF_only_binding_sites.bed results/ZNF143_only_binding_sites.bed \
    -b 4000 \
    -a 4000 \
    -o results/matrix.reference.point.RPGC.cobind.mat.gz \
    --missingDataAsZero \
    --referencePoint "center"

# Plot heat map 
plotHeatmap -m results/matrix.reference.point.RPGC.cobind.mat.gz \
    -out plots/Coverage_heatmap_reference_point_RPGC_cobind.png \
    --startLabel start \
    --endLabel end \
    --samplesLabel "shControl" "shZNF143" \
    --boxAroundHeatmaps no \
    --regionsLabel "ZNF143 & CTCF" "CTCF only" "ZNF143 only" \
    --yAxisLabel "Coverage (RPGC)" \
    --colorMap "Purples" \
    --xAxisLabel "binding site (bp)" \
    --refPointLabel "TFBS"

