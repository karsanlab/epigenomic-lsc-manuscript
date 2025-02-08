# INPUT: 
# `results/GSM*_farlik_region_meth.tsv`
# `results/*_farlik_region_meth.tsv`
# `results/combinedSampleAnnot.tsv`

# OUTPUT: 
# `results/combinedMethMatrix.tsv`

library(rtracklayer)
library(here)
library(tidyverse)

in.dir <- "results"
out.dir <- "results"
samp_annot <- here("results", "combinedSampleAnnot.tsv")

sample_annot <- read_tsv(samp_annot)

farlik_tsvs <- file.path(in.dir, list.files(in.dir, pattern = "*_farlik_region_meth.tsv") %>% str_subset("GSM[0-9]{7}_"))
farlik_matrix <- farlik_tsvs %>%
  map(read_tsv, na = c("NA", "NaN")) %>%
  list_cbind()

# reorder farlik_matrix to be in the same order as samp_annot

## Load sample tsvs and merge
sample_tsvs <- file.path(in.dir, list.files(in.dir, pattern = "*_farlik_region_meth.tsv") %>% str_subset("[0-9]{3}_PX"))
sample_matrix <- sample_tsvs %>%
  map(read_tsv, na = c("NA", "NaN")) %>%
  list_cbind()

# reorder sample_matrix to be the same order as samp_annot
stopifnot(identical(c(colnames(farlik_matrix), colnames(sample_matrix)), sample_annot$sampleName,  attrib.as.set = FALSE))

## Merge matrices
combined_matrix <- bind_cols(farlik_matrix, sample_matrix)

combined_matrix <- combined_matrix[, sample_annot$sampleName]

stopifnot(identical(colnames(combined_matrix), sample_annot$sampleName,  attrib.as.set = FALSE))

# Check for missingness, must be <98% missing

missingness_prop <- apply(combined_matrix, 2, function(x) sum(is.na(x))/length(x))
non_missing_libs <- names(missingness_prop[missingness_prop < 0.98])
combined_matrix <- combined_matrix %>%
  select(one_of(non_missing_libs))

## Write to file
write_tsv(combined_matrix, file.path(out.dir, "combinedMethMatrix.tsv"))

# remove missing samples from sample annotation file
sample_annot %>%
  filter(sampleName %in% non_missing_libs) %>%
  write_tsv(samp_annot)
