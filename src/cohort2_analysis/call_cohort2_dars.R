# INPUT:
# `results/integrated_cohort2_samples.rds`

# OUTPUT: 
# `results/Cohort2_FindMarkers_DARs.rds`

library(Seurat)
library(Signac)
library(pbmcapply)
library(here)
library(tidyverse)

obj <- readRDS(here("results", "integrated_cohort2_samples.rds"))

FindMarkers.peaks <- function(feature) {
    peak_region_frags <- str_subset(colnames(obj@meta.data), "peak_region_fragments")
    
    results <- FindMarkers(obj, 
                           assay = "peaks",
                           features = feature,
                           logfc.threshold = 0,
                           test.use = "LR",
                           min.pct = 0, 
                           latent.vars = peak_region_frags, 
                           ident.1 = "responder", 
                           ident.2 = "nonresponder")
    return(results)
}

DefaultAssay(obj) <- "peaks"
Idents(obj) <- obj$response

results <- pbmclapply(rownames(obj), FindMarkers.peaks, mc.cores = 125) %>%
    list_rbind()

results <- results %>%
    mutate(dar_id = str_glue("DAR_{str_pad(row_number(), str_length(nrow(results)), 'left', pad = '0')}"))

saveRDS(results, here("results", "Cohort2_FindMarkers_DARs.rds"))
