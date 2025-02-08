#!/usr/bin/env python3
# https://github.com/aertslab/arboreto/issues/42
import pandas as pd
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

if __name__ == '__main__':
    # ex_matrix is a DataFrame with gene names as column names
    atac_matrix = pd.read_csv("cache/atac_matrix.tsv", sep='\t')
    
    # tf_names is read using a utility function included in Arboreto
    tf_names = load_tf_names("cache/64tfbs.tsv")
    
    atac_network = grnboost2(expression_data=atac_matrix, tf_names=tf_names)

    atac_network.to_csv('results/atac_GRNBoost2.tsv', sep='\t', index=False, header=False)
