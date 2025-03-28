# INPUT:
# `results/*.cons.Cmethyl_tfbs.csv`
# `results/combinedSampleAnnot.tsv`

# OUTPUT:
# `results/patient_449_keq2_clustering.tsv`

library(pcaMethods)
library(umap)
library(here)
library(tidyverse)

#' calculate_normalization_stats
#'
#' Calculate summary statistics for sample-level normalization of TFBS methylation
#' @param methyl_tib tibble containing per-TFBS methylation
calculate_normalization_stats <- function(methyl_tib) {
  methyl_tib %>%
    left_join(methyl_tib %>%
                drop_na() %>%
                group_by(filename) %>%
                summarise(mean = mean(pct_meth))) %>%
    mutate(pct_of_mean = pct_meth / mean,
           diff_from_mean = pct_meth - mean,
           fold_change_over_mean = (pct_meth - mean) / mean)
}

# Import single cell metadata
aza_metadata.tib <- read_tsv(here("data", "sample_metadata.tsv"))
farlik_metadata.tib <- read_tsv(here("data", "farlik_metadata.tsv"))

# Import single cell methylation data summarised per TFBS in LOLA database
aza_and_farlik_tfbs_methylation.tib <-
    list.files("data",
               "*.cons.Cmethyl_tfbs.csv", full.names = TRUE) %>%
    set_names(str_remove(., "data") %>%
                  str_remove(".cons.*") %>%
                  str_remove("^[\\d]+_")) %>%
	map(~ read_csv(.x)) %>%
    bind_rows() %>%
    calculate_normalization_stats()

aza_and_farlik_tfbs_methylation.mat <-
    aza_and_farlik_tfbs_methylation.tib %>%
        dplyr::distinct(filename, sample, fold_change_over_mean) %>%
        pivot_wider(names_from = filename, values_from = fold_change_over_mean)

# Perform UMAP
set.seed(1234)

aza_patients.umap <-
    umap(aza_and_farlik_tfbs_methylation.mat[aza_metadata.tib$cell_barcode, ])

# Perform co-projection of normal reference cells
farlik_projected.df <-
    predict(aza_patients.umap,
            aza_and_farlik_tfbs_methylation.mat[farlik_metadata.tib$sample_id, ]) %>%
        as_tibble(rownames = "plate_index") %>%
        set_names(c("plate_index", "UMAP_1", "UMAP_2"))
