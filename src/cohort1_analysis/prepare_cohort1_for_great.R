# INPUT: 
# `results/cohort1_dmrs.rds`

# OUTPUT: 
# `results/Cohort1_hypo_DMR_delta-0.4.bed`
# `results/Cohort1_hyper_DMR_delta+0.4.bed`

library(here)
library(rtracklayer)
library(GenomicRanges)
library(plyranges)
library(tidyverse)

dmrs.gr <- readRDS(here("results", "cohort1_dmrs_gr.rds"))

dmrs.gr %>% 
    plyranges::filter(diff.Methy < -0.4) %>%
    rtracklayer::export.bed(here("results", "Cohort1_hypo_DMR_delta-0.4.bed"))

dmrs.gr %>% 
    plyranges::filter(diff.Methy > 0.4) %>%
    rtracklayer::export.bed(here("results", "Cohort1_hyper_DMR_delta+0.4.bed"))
