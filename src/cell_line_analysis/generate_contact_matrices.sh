#!/bin/bash

# INPUT:
# results/${sample}_mapped.pairs

# OUTPUT:
# results/${sample}_contact_map.hic
# results/${sample}_contact_map_5kb.hic

# Activate conda environment. 
# Environment specs to recreate the environment are available in the resources folder as "hichip_conda_environment_spec_file.txt"

eval "$(conda shell.bash hook)" 
conda activate hichip_conda_environment

# Generate hic contact map at several resolutions

java -Xmx48000m  -Djava.awt.headless=true -jar juicertools.jar pre --threads 16 \
  results/${sample}_mapped.pairs \
  results/${sample}_contact_map.hic \
  hg38.genome

# Generate single hic contact map at 5KB resolution
java -Xmx48000m  -Djava.awt.headless=true -jar juicertools.jar pre -r 5000 --threads 16 \
  results/${sample}_mapped.pairs \
  results/${sample}_contact_map_5kb.hic \
  hg38.genome

