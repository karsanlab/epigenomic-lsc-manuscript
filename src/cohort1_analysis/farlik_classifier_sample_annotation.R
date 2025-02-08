# INPUT: 
# `resources/supp_data/Farlik_cell_type_prediction/sampleAnnot.tsv` (from Farlik)
# `data/full_metadata_qc.rds`

# OUTPUT: `results/combinedSampleAnnot.tsv`

library(here)
library(tidyverse)

farlik_sample_annot <- read_tsv(here("resources", "supp_data", "Farlik_cell_type_prediction", "sampleAnnot.tsv"))

## Create sample annotation file with our single cell libraries included ##
sample_data <- readRDS(here("data", "full_metadata_qc.rds")) %>%
  mutate(cellSourceCurated = "bone_marrow",
         class = "unknown",
         nominalCellNumber = 1,
         train = FALSE,
         color = "#800000") %>%
  select(sampleName = plate_sample,
         class, nominalCellNumber, cellSourceCurated,
         donorId = patient,
         train, color)

combined_sample_annot <- bind_rows(farlik_sample_annot, sample_data)

write_tsv(combined_sample_annot, here("results", "combinedSampleAnnot.tsv"))