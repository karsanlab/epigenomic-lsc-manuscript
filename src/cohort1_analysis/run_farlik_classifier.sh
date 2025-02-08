#!/bin/bash
# INPUT: 
# `results/combinedMethMatrix.tsv`
# `results/combinedSampleAnnot.tsv`

# OUTPUT: 
# `results/farlik/prediction/predictionResults.rds`

Rscript resources/supp_data/Farlik_cell_type_prediction/cellTypePredictor.R \
	--features results/combinedMethMatrix.tsv \
	--annot results/combinedSampleAnnot.tsv \
	--out results/farlik/prediction
