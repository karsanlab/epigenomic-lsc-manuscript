#!/bin/bash

# INPUT: 
# - data/E_MTAB_5456/*.fastq.gz
# - resources/Homo_sapiens.GRCh38.cdna.all.fa.gz

# OUTPUT:
# - results/E_MTAB_5456/*_quant/quant.sf

salmon index -t resources/Homo_sapiens.GRCh38.cdna.all.fa.gz -i resources/hg38_cDNA_index -k 21

for fq in data/E_MTAB_5456/*.fastq.gz; do 
	sample=$(basename $fq "*.fastq.gz") 
	salmon quant -i resources/hg38_cDNA_index -l A -r $fq --validateMappings -o results/E_MTAB_5456/${sample}_quant
done
