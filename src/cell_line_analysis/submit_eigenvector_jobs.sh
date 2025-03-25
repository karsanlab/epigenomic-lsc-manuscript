#!/bin/bash

# Submit one job per chromosome

for CHROMOSOME in `seq 1 1 22`
do
sbatch --export=chr=${CHROMOSOME} generate_ab_compartments.sh
done