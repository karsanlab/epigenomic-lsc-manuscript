# INPUT: 
# `results/integrated_cohort2_samples.rds`

# OUTPUT: 
# `results/Cohort2_FindMarkers_DEGs.rds`

library(Seurat)
library(Signac)
library(pbmcapply)
library(here)
library(tidyverse)

obj <- readRDS(here("results", "integrated_cohort2_samples.rds"))

FindMarkers.RNA <- function(feature) {
    results <- FindMarkers(obj, 
                           assay = "RNA",
                           features = feature,
                           logfc.threshold = 0,
                           test.use = "LR",
                           min.pct = 0, 
                           ident.1 = "Responder", 
                           ident.2 = "Non-responder")
    return(results)
}

DefaultAssay(obj) <- "RNA"
Idents(obj) <- obj$response

results <- pbmclapply(rownames(obj), FindMarkers.RNA, mc.cores = 125) %>%
    list_rbind()

results <- results %>%
    mutate(deg_id = str_glue("DEG_{str_pad(row_number(), str_length(nrow(results)), 'left', pad = '0')}"), 
           gene = rownames(.))

saveRDS(results, here("results", "Cohort2_FindMarkers_DEGs.rds"))

