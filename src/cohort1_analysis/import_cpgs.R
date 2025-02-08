# INPUT: 
# `data/scWGBS/*.cons.Cmethyl.CpG.txt`
# `data/sample_data.csv`
# `data/qc.csv`

# OUTPUT: 
# `results/hg19_novoalign/*_cpg_pool.tib.rds`
# `cache/hg19_novoalign_combined_methylation_global.tibs.rds`
# `data/full_metdata.rds`
# `data/full_metadata_qc.rds`
# `data/passed_qc.rds`

library(here)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
DATA_PATH <- args[1]

sample.data <- read.csv(here("data", "sample_data.csv"), row.names = 1)
qc.data <- read.csv(here("data", "qc.csv"), row.names = 1)

MIN_CPGS_COVERED <- 10^5

metadata <- full_join(sample.data, qc.data, by = "plate") %>%
    mutate(passed.qc = case_when(ncpg > MIN_CPGS_COVERED & group == "single-cell" ~ TRUE, 
                                 TRUE ~ FALSE),
           Class = case_when(str_detect(progression, "Responder") ~ "Responder", 
                             str_detect(progression, "Nonresponder") ~ "Non-responder") %>%
               factor(levels = c("Responder", "Non-responder")),
           patient = as.character(patient)) %>%
    dplyr::rename(Plate_ID = plate)

saveRDS(metadata, here("data", "full_metadata.rds"))

passed.qc <- metadata %>%
    filter(passed.qc == TRUE)

saveRDS(passed.qc, here("data", "full_metadata_qc.rds"))

passed.qc <- passed.qc %>%
    select(Plate_ID, cell)

saveRDS(passed.qc, here("data", "passed_qc.rds"))

plates <- select(metadata, patient, Plate_ID) %>%
    distinct() %>% 
    unite(tmp, patient, Plate_ID, sep =  "_") %>%
    pull(tmp)
names(plates) <- plates

bed_filepaths.lists <- map(plates, ~ {
    tmp <- list.files(DATA_PATH, full.names = T, pattern = .x) %>%
                str_subset(".cons.Cmethyl.CpG.txt") %>%
        as.list() %>%
        set_names(map_chr(., ~ {basename(.x) %>%
                str_split("\\.", simplify = T) %>%
                .[,1] %>%
                str_split("_", simplify = T) %>%
                .[,3]}))
})

# https://www.novocraft.com/documentation/novoalign-2/novoalign-user-guide/bisulphite-treated-reads/novomethyl-analyzing-methylation-status/
read_cpg_bed <- function(bed_path) {
    results <- read_tsv(bed_path,
                        col_types = c("ci--ii"),
                        col_names = c("seqnames", "start", "n_methylated", "n_total"))
    
}

import_library_beds <- function(bed_filepaths.list) {
    imap(bed_filepaths.list, ~ {
        message(str_glue("Loading {.y}..."))
        read_cpg_bed(.x)
    }) %>%
        bind_rows(.id = "cell") %>%
        group_by(seqnames, start, cell) %>%
        summarise(.groups = "drop",
                  n_methylated = sum(n_methylated),
                  n_total = sum(n_total))
}

dir.create(here("cache", "hg19_novoalign"), recursive = T, showWarnings = F)

combined_methylation.tibs <-
    imap(bed_filepaths.lists,
         function(bed_filepaths.list, library_id) {
             combined_methylation.tib <- import_library_beds(bed_filepaths.list)
             
             write_rds(combined_methylation.tib,
                       here::here("cache", "hg19_novoalign", 
                                  paste0("hg19_", library_id, "_cpg_pool.tib.rds")),
                       compress = "gz")
             return(combined_methylation.tib)
         })


write_rds(combined_methylation.tibs, 
          here("cache", str_glue("hg19_novoalign_combined_methylation_global.tibs.rds")), compress = "gz")