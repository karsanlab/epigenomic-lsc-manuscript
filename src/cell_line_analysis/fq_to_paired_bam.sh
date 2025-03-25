#!/bin/bash

# INPUT:
# data/{sample}.FQ1.fastq.gz
# data/{sample}.FQ2.fastq.gz

# OUTPUT:
# results/${sample}_stats.txt
# results/${sample}_mapped.pairs
# results/${sample}_mapped.PT.bam

# Activate conda environment. 
# Environment specs to recreate the environment are available in the resources folder as "hichip_conda_environment_spec_file.txt"

eval "$(conda shell.bash hook)" 
conda activate hichip_conda_environment

# Create output directories
mkdir -p results/${sample} ${sample}/TEMP

# Align paired end reads with `bwa mem`, parse, sort, deduplicate and split alignments using `pairtools`, then sort and index using `samtools

bwa mem -5SP -T0 -t16 /resources/hg38.fa data/${sample}.FQ1.fastq.gz data/${sample}.FQ2.fastq.gz | \
	pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 \
		--nproc-out 8 --chroms-path resources/hg38.genome | 
	pairtools sort --tmpdir=${sample}/TEMP \
		--nproc 16 | pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats results/${sample}_stats.txt| \
	pairtools split --nproc-in 8 --nproc-out 8 --output-pairs results/${sample}_mapped.pairs \
		--output-sam -| samtools view -bS -@16 | samtools sort -@16 -o results/${sample}_mapped.PT.bam

samtools index results/${sample}_mapped.PT.bam
