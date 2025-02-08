# INPUT:
# `results/mgatk_outputs/*_mgatk/final`

# OUTPUT:
# `results/patient_A07_449_mito_metadata.tsv`
# `results/patient_A07_449_mito_vafs.tsv`


library(Seurat)
library(Signac)
library(here)
library(tidyverse)

MIN_STRAND_COR <- 0.65
MIN_VMR <- 10^-2
BISULFITE_MUTATIONS <- c("C>T", "G>A")

scwgbs_cell_metadata.tib <-
    read_tsv(here("data", "combined_sc_methylation_metadata.tsv"))

#' import_mgatk_as_seurat_obj
#'
#' Import output data from mgatk into a Seurat Object
#' @param mgatk_data_dir path to "final" subdirectory of mgatk data
import_mgatk_as_seurat_obj <- function(mgatk_data_dir) {
    mgatk_data <- ReadMGATK(mgatk_data_dir)
    assay_obj <- CreateAssayObject(counts = mgatk_data$counts)

    metadata <- mgatk_data$depth %>%
        rownames_to_column("barcode") %>%
        left_join(scwgbs_cell_metadata.tib) %>%
        column_to_rownames("barcode")
    
    seurat_obj <- CreateSeuratObject(
        counts = assay_obj,
        assay = "mito",
        meta.data = metadata
    )
    
    seurat_obj@assays$mito@key <- "mito_"
    
    return(seurat_obj)
}

scwgbs_mgatk.signacs <-
    list.files(here("results", "mgatk_outputs"), full.names = TRUE) %>%
    paste0("/final") %>%
    set_names(., str_remove(., ".*mgatk_outputs/") %>%
                  str_remove("_mgatk.*")) %>%
    map(import_mgatk_as_seurat_obj)

# Extract chrM refallele to be used as reference
MITO_REFALLELES <-
    list.files(here("results", "mgatk_outputs"), full.names = TRUE) %>%
    .[1] %>%
    paste0("/final") %>%
    ReadMGATK() %>%
    .$refallele

# Filter for QC passsing data
scwgbs_mgatk_qcpass.signacs <-
    map(scwgbs_mgatk.signacs,
        ~ subset(., subset = qc_pass))

# Identify chrM variants across all cells
variable_sites.dfs <-
    map(scwgbs_mgatk_qcpass.signacs,
        ~ IdentifyVariants(., refallele = MITO_REFALLELES))

# Select variants passing filtration thresholds
high_confidence_variants.dfs <-
    map(variable_sites.dfs,
        ~ subset(.,
                 subset = n_cells_conf_detected >= 5 &
                     !nucleotide %in% BISULFITE_MUTATIONS &
                     !position %in% HOMOPOLYMER_BLACKLIST &
                     strand_correlation >= MIN_STRAND_COR &
                     vmr > MIN_VMR))

# Calculate allele frequency for selected variants
scwgbs_mgatk_w_clonotype.signacs <-
    map2(scwgbs_mgatk_qcpass.signacs, high_confidence_variants.dfs,
         function(signac_obj, high_conf_variant_df) {
             AlleleFreq(signac_obj,
                        variants = high_conf_variant_df$variant,
                        assay = "mito")

            DefaultAssay(signac_obj) <- "alleles"

            res <- FindClonotypes(signac_obj)
            
            return(res)
         })

# Export data into tables for downstream use and visualization
scwgbs_mgatk_w_clonotype.signac$patient_A07_449@meta.data %>%
    as_tibble(rownames = "cell_barcode") %>%
    write_tsv(here("results", "patient_A07_449_mito_metadata.tsv"))

scwgbs_mgatk_w_clonotype.signac$patient_A07_449@assays$alleles@counts %>%
    as.matrix() %>%
    as_tibble(rownames = "variant") %>%
    pivot_longer(-variant, names_to = "cell_barcode", values_to = "vaf") %>%
    write_tsv(here("results", "patient_A07_449_mito_vafs.tsv"))

