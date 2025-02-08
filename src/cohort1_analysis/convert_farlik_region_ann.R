# INPUT: 
# `resources/supp_data/liftover/hg38ToHg19.over.chain`
# `resources/supp_data/Farlik_cell_type_prediction/regionAnnot.tsv`

# OUTPUT
# `cache/regionAnnot.bed`
# `cache/regionAnnot.sorted.bed`

library(rtracklayer)
library(readr)
library(GenomicRanges)
library(here)
library(tidyverse)

hg38.19.chain <- import.chain(here("resources", "supp_data", "liftover", "hg38ToHg19.over.chain"))

farlik_regions <- read_tsv(here("resources", "supp_data", "Farlik_cell_type_prediction", "regionAnnot.tsv")) %>%
  makeGRangesFromDataFrame()

farlik_regions.hg19 <- liftOver(farlik_regions, hg38.19.chain)
farlik_regions.hg19 <- unlist(farlik_regions.hg19) %>%
    sort()

export.bed(farlik_regions.hg19, here("cache", "regionAnnot.bed"))


system("bedtools sort -i cache/regionAnnot.bed > cache/regionAnnot.sorted.bed")
