# INPUT:
# `results/EMTAB5456_dge.csv`
# `results/Cohort2_FindMarkers_DEGs.rds`

# OUTPUT:
# `results/upreg_normal_response_intersect.rds`
# `results/upreg_normal_response_intersect_tfs.rds`

library(tidyverse)
library(here)

vyas_de <- read.csv(here("results", "EMTAB5456_dge.csv"), row.names = 1 ) %>%
    drop_na()

# Downloaded from https://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt 
tf_names <- readLines(here("resources", "TF_names_v_1.01.txt"))

upreg_vyas <- vyas_de %>%
    filter(log2FoldChange >= 2, padj < 0.1) %>%
    rownames()

upreg_vyas_tfs <- vyas_de %>%
    filter(log2FoldChange >= 2, padj < 0.1) %>%
    rownames() %>%
    intersect(tf_names)

downreg_vyas <- vyas_de %>%
    filter(log2FoldChange <= -2, padj < 0.1) %>%
    rownames() 

downreg_vyas_tfs <- vyas_de %>%
    filter(log2FoldChange <= -2, padj < 0.1) %>%
    rownames() %>%
    intersect(tf_names)

results <- readRDS(here("results", "Cohort2_FindMarkers_DEGs.rds"))

upreg_response <- results %>%
    filter(p_val_adj < 0.001, avg_log2FC >= 0.05) %>%
    rownames()

upreg_response_tfs <- results %>%
    filter(p_val_adj < 0.001, avg_log2FC >= 0.05) %>%
    rownames() %>%
    intersect(tf_names)

downreg_response <- results %>%
    filter(p_val_adj < 0.001, avg_log2FC <= -0.05) %>%
    rownames()

downreg_response_tfs <- results %>%
    filter(p_val_adj < 0.001, avg_log2FC <= -0.05) %>%
    rownames() %>%
    intersect(tf_names)

intersect.upreg <- intersect(upreg_vyas, upreg_response)
intersect.downreg <- intersect(downreg_vyas, downreg_response)

intersect.upreg.tfs <- intersect(upreg_vyas_tfs, upreg_response_tfs)
intersect.downreg.tfs <- intersect(downreg_vyas_tfs, downreg_response_tfs)

join_df <- function(fn, x, y, rownames, df.out = FALSE, ...) {
    if(df.out == TRUE) {
        rownames <- "tmp"
    }
    tib <- exec(fn, x = as_tibble(x, rownames = rownames), y = as_tibble(y, rownames = rownames), ...)
    
    if (df.out == TRUE) {
        df <- column_to_rownames(tib, rownames) %>%
            as.data.frame()
        return(df)
    } else {
        return(tib)
    }
}

upreg_genes_all <- join_df(fn = full_join, x = vyas_de[intersect.upreg,], y = results[intersect.upreg,], rownames = "Gene", suffix = c(".normal", ".patient")) %>%
    dplyr::select(Gene, log2FoldChange, padj, avg_log2FC, p_val_adj) %>%
    filter(!is.na(Gene)) %>%
    column_to_rownames("Gene")

downreg_genes_all <- join_df(fn = full_join, x = vyas_de[intersect.downreg,], y = results[intersect.downreg,], rownames = "Gene", suffix = c(".normal", ".patient")) %>%
    dplyr::select(Gene, log2FoldChange, padj, avg_log2FC, p_val_adj) %>%
    filter(!is.na(Gene)) %>%
    column_to_rownames("Gene")

upreg_tfs_all <- join_df(fn = full_join, x = vyas_de[intersect.upreg.tfs,], y = results[intersect.upreg.tfs,], rownames = "Gene", suffix = c(".normal", ".patient")) %>%
    dplyr::select(Gene, log2FoldChange, padj, avg_log2FC, p_val_adj) %>%
    filter(!is.na(Gene)) %>%
    column_to_rownames("Gene")

downreg_tfs_all <- join_df(fn = full_join, x = vyas_de[intersect.downreg.tfs,], y = results[intersect.downreg.tfs,], rownames = "Gene", suffix = c(".normal", ".patient")) %>%
    dplyr::select(Gene, log2FoldChange, padj, avg_log2FC, p_val_adj) %>%
    filter(!is.na(Gene)) %>%
    column_to_rownames("Gene")

saveRDS(upreg_genes_all, here("results", "upreg_normal_response_intersect.rds"))
saveRDS(downreg_genes_all, here("results", "downreg_normal_response_intersect.rds"))

saveRDS(upreg_tfs_all, here("results", "upreg_normal_response_intersect_tfs.rds"))
saveRDS(downreg_tfs_all, here("results", "downreg_normal_response_intersect_tfs.rds"))
