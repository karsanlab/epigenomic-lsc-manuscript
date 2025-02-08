# INPUT:
# `results/scrambled_integrated_cohort2_samples.rds`

# OUTPUT: 
# `results/Scrambled_Cohort2_FindMarkers_DARs.rds`

library(Seurat)
library(Signac)
library(pbmcapply)
library(here)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

iter <- args[1]

column_name <- str_glue("even_patient_scrambled_response_{iter}")

obj <- readRDS(here("results", "Scrambled_Cohort2_FindMarkers_DARs.rds"))

FindMarkers.peaks <- function(feature) {
    peak_region_frags <- str_subset(colnames(obj@meta.data), "peak_region_fragments")
    
    results <- FindMarkers(obj, 
                           assay = "peaks",
                           features = feature,
                           logfc.threshold = 0,
                           test.use = "LR",
                           min.pct = 0, 
                           latent.vars = peak_region_frags, 
                           ident.1 = "group1", 
                           ident.2 = "group2")
    return(results)
}

DefaultAssay(obj) <- "peaks"
Idents(obj) <- obj[[column name]]

results <- pbmclapply(rownames(obj), FindMarkers.peaks, mc.cores = 125) %>%
    list_rbind()

results <- results %>%
    mutate(dar_id = str_glue("DAR_{str_pad(row_number(), str_length(nrow(results)), 'left', pad = '0')}"))

saveRDS(results, here("results", "Scrambled_Cohort2_FindMarkers_DARs_v{iter}.rds"))
