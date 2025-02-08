# INPUT: 
# `results/Cohort2_FindMarkers_DAR.rds`

# OUTPUT: 
# `results/Cohort2_hyper_DAR_padj_0.05.bed`

library(here)
library(rtracklayer)
library(GenomicRanges)
library(plyranges)
library(tidyverse)

dars.gr <- readRDS(here("results", "Cohort2_FindMarkers_DARs.rds")) %>%
    rownames_to_column("region") %>%
    separate(region, into = c("seqnames", "start", "end", sep = "-")) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)

dars.gr %>%
    plyranges::filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
    rtracklayer::export.bed(here("results", "Cohort2_hyper_DAR_padj_0.05.bed"))