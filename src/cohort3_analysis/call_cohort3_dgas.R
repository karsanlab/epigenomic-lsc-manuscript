# INPUT: 
# `results/integrated_cohort3_samples.rds`

# OUTPUT: 
# `results/Cohort3_FindMarkers_DGA.rds`

library(Seurat)
library(Signac)
library(pbmcapply)
library(here)
library(tidyverse)

obj <- readRDS(here("results", "integrated_cohort3_samples.rds"))

FindMarkers.GeneActivity <- function(feature) {
    results <- FindMarkers(obj, 
                           assay = "GeneActivity",
                           features = feature,
                           logfc.threshold = 0,
                           test.use = "LR",
                           min.pct = 0, 
                           ident.1 = "Responder", 
                           ident.2 = "Non-responder")
    return(results)
}
DefaultAssay(obj) <- "GeneActivity"
Idents(obj) <- obj$response

results <- pbmclapply(rownames(obj), FindMarkers.GeneActivity, mc.cores = 125) %>%
    list_rbind()

results <- results %>%
    mutate(dga_id = str_glue("DGA_{str_pad(row_number(), str_length(nrow(results)), 'left', pad = '0')}"), 
           gene = rownames(.))

saveRDS(results, here("results", "Cohort3_FindMarkers_DGAs.rds"))
