# INPUT:
# `data/passed_qc.rds`
# `cache/hg19_novoalign_combined_methylation_global.tibs.rds`

# OUTPUT: 
# `cache/hg19_novoalign_global_methylation_qc.rds`

library(tidyverse)
library(here)

metadata <- readRDS(here("data", "full_metadata.rds"))
passed.qc <- readRDS(here("data", "passed_qc.rds"))

combined_methylation.tibs <- readRDS(here("cache", "hg19_novoalign_combined_methylation_global.tibs.rds"))

combined_methylation.tibs <-
    map(combined_methylation.tibs,
        function(tib) {
            tib %>%
                dplyr::select("chr" = seqnames,
                              "pos" = start,
                              "N" = n_total,
                              "X" = n_methylated,
                              cell)
        })

global_methylation.tib <- combined_methylation.tibs %>%
    imap(~ {
      .x %>%
        mutate(Plate_ID = str_extract(.y, "PX[0-9]{4}")) %>%
        semi_join(passed.qc, by = c("Plate_ID", "cell")) %>%
        group_by(cell) %>%
        summarize(percent_meth = sum(X) / sum(N) * 100)
    }) %>%
    list_rbind(names_to = "id") %>%
    separate(id, into = c("patient", "Plate_ID"), sep = "_", convert = F) %>%
    left_join(metadata, by = c("Plate_ID", "patient", "cell"))

saveRDS(global_methylation.tib, here("cache", "hg19_novoalign_global_methylation_qc.rds"))