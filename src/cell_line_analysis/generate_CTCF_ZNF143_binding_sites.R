#!/usr/bin/env Rscript

# INPUT:
# results/ctcf_lola_and_aml3_peaks_combined.bed
# resources/znf143_narrowpeaks.bed # LOLA ZNF143 peakset described in resources

# OUTPUT:
# results/CTCF_only_binding_sites.bed
# results/ZNF143_only_binding_sites.bed
# results/CTCF_ZNF143_cobinding_sites.bed


# Load libraries
library(here)
library(tidyverse)
library(GenomicRanges)
library(annotatr)

# Read in CTCF and ZNF143 peak sets
CTCF_combined_peaks <- read.table(here("data", "ctcf_lola_and_aml3_peaks_combined.bed"), 
                                  col.names = c("chr","start","stop","length","star","source")) %>% 
    makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

ZNF143_peaks <- read.table(here("data", "znf143_narrowpeaks.bed"), 
                           col.names = c("chr","start","stop","length","star")) %>% 
    makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

# Add names metadata to gr peak objects
mcols(CTCF_combined_peaks)$type <- "CTCF"
mcols(ZNF143_peaks)$type <- "ZNF143"

# Add to annotatr cache
annotatr_cache$set("CTCF", CTCF_combined_peaks)
annotatr_cache$set("ZNF143", ZNF143_peaks)

# Build annotations
ZNF143_annots <- build_annotations(genome = "hg38", annotations = c("ZNF143"))
CTCF_annots <- build_annotations(genome = "hg38", annotations = c("CTCF"))

# Overlap CTCF with ZNF143 annots and vice versa so I can filter to get CTCF only and ZNF143 only loops

# CTCF_ZNF143_shared_sites
CTCF_ZNF143_overlaps <- annotate_regions(
    regions = CTCF_combined_peaks,
    annotations = ZNF143_annots,
    ignore.strand = TRUE,
    quiet = FALSE)

# CTCF only
CTCF_only <- subsetByOverlaps(CTCF_combined_peaks,CTCF_ZNF143_overlaps, invert = TRUE, type = "equal") 

ZNF143_CTCF_overlaps <- annotate_regions(
    regions = ZNF143_peaks,
    annotations = CTCF_annots,
    ignore.strand = TRUE,
    quiet = FALSE)

ZNF143_only <- subsetByOverlaps(ZNF143_peaks, ZNF143_CTCF_overlaps, invert = TRUE, type = "equal")

# Export to use with generate heatmaps
export.bed(CTCF_only, here("results","CTCF_only_binding_sites.bed"))
export.bed(ZNF143_only, here("results","ZNF143_only_binding_sites.bed"))
export.bed(CTCF_ZNF143_overlaps, here("results","CTCF_ZNF143_cobinding_sites.bed"))

