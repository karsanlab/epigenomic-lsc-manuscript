#!/usr/bin/env Rscript

# INPUT:
# results/${sample}.macs2_peaks.narrowPeak
# resources/ctcf_narrowpeaks.bed # LOLA CTCF peakset described in resources

# OUTPUT:
# results/ctcf_lola_and_aml3_peaks_combined.bed
# results/ctcf_expanded_narrowpeaks.bed

library(here)
library(tidyverse)
library(GenomicRanges)

# Expand LOLA CTCF peaks as they are very narrow (~170bp)
lola_peaks <- read.table(here("resources", "ctcf_narrowpeaks.bed"))

lola_peaks %>% 
    mutate(width = V3 - V2) %>% 
    summarize(mean(width))

lola_peaks_expanded <- lola_peaks %>% 
    mutate(V2 = V2 -150,
           V3 = V3 + 150)

lola_peaks_expanded %>% 
    mutate(width = V3 - V2) %>% 
    summarize(mean(width))

# Write expanded lola_peaks to file
write.table(lola_peaks_expanded, here("results", "ctcf_expanded_narrowpeaks.bed"), 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Read in AML3_shControl peaks to check width (~650 bp)
# Don't need to expand
AML3_peaks <- read.table(here("resources", "AML3_shControl.macs2_peaks.narrowPeak"), header = FALSE)

AML3_peaks %>% 
    mutate(width = V3 - V2) %>% 
    summarize(mean(width))

# Convert both to granges and combine
AML3_peaks.gr <- makeGRangesFromDataFrame(AML3_peaks, 
                                          keep.extra.columns = TRUE, 
                                          seqnames.field = "V1",
                                          start.field = "V2",
                                          end.field = "V3")
mcols(AML3_peaks.gr)$source <- "AML3"

lola_peaks_expanded.gr <- makeGRangesFromDataFrame(lola_peaks_expanded, 
                                                   keep.extra.columns = TRUE, 
                                                   seqnames.field = "V1",
                                                   start.field = "V2",
                                                   end.field = "V3") 
mcols(lola_peaks_expanded.gr)$source <- "lola"

# Combine ranges keeping peak source info
concat_sites.gr <- plyranges::bind_ranges(aml3_peaks = AML3_peaks.gr, 
                                          lola_peaks = lola_peaks_expanded.gr, 
                                          .id = "peak_set") %>% 
    plyranges::select(source)

# Reduce ranges and add back source annotation
concat_sites.gr <- concat_sites.gr %>% 
    reduce() %>% 
    splicejam::annotateGRfromGR(., concat_sites.gr)

# Write to file
write.table(concat_sites.gr, here("results", "ctcf_lola_and_aml3_peaks_combined.bed"), 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

