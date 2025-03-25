#!/bin/bash

# INPUT:
# results/${sample}_contact_map.hic

# OUTPUT:
# results/${sample}_chr${chr}_eigenvector_50kb.txt
# results/${sample}_chr${chr}_eigenvector.txt

# Activate conda environment. 
# Environment specs to recreate the environment are available in the resources folder as "hichip_conda_environment_spec_file.txt"

eval "$(conda shell.bash hook)" 
conda activate hichip_conda_environment

echo $chr

java -Xmx48000m -Djava.awt.headless=true -jar juicertools.jar eigenvector -p\
    KR results/${sample}_contact_map.hic chr${chr} BP 50000 > \
    results/${sample}_chr${chr}_eigenvector_50kb.txt

java -Xmx48000m -Djava.awt.headless=true -jar juicertools.jar eigenvector \
    KR results/${sample}_contact_map.hic chr${chr} BP 1000000 > \
    results/${sample}_chr${chr}_eigenvector.txt