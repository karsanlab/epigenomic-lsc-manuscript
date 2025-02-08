# INPUT: 
# `results/hg19_novoalign/*_cpg_pool.tib.rds`

# OUTPUT: 
# `cache/hg19_novoalign/scramble/*_sc.rds`
# `cache/hg19_novoalign/scramble/scrambled_library_*.rds`
# `cache/hg19_novoalign/scramble/pt*.rds`

library(GenomicRanges)
library(tidyverse)
library(here)

passed.qc <- readRDS(here("data", "passed_qc.rds"))
dir.create(here("cache", "hg19_novoalign", "scramble", recursive = T, showWarnings = F))

list.files(here("cache", "hg19_novoalign"), pattern = "*_cpg_pool.tib.rds", full.names = T) %>%
    str_subset("hg19") %>%
    set_names(list.files(here("cache", "hg19_novoalign"), pattern = "*_cpg_pool.tib.rds") %>% str_subset("hg19") %>% str_split("_", simplify = TRUE) %>% .[,3]) %>%
    iwalk(function(rds, plate) {
        df <- readRDS(rds) %>%
            dplyr::select(seqnames,
                          start,
                          n_total,
                          "n_methyl" = n_methylated,
                          cell) %>%
            dplyr::mutate(end = start + 1,
                          Plate_ID = plate,
                          plate_sample = str_glue("{plate}_{cell}")) %>%
            semi_join(passed.qc, by = c("Plate_ID", "cell")) %>%
            dplyr::select(-Plate_ID, -cell)
        
        df <- group_by(df, plate_sample)
        split <- group_split(df, .keep = F)
        names(split) <- group_keys(df) %>% pull(plate_sample)
        
        grlist <- map(split, ~ {
            makeGRangesFromDataFrame(.x, keep.extra.columns = T, ignore.strand = T)
        }) %>%
            sample(length(.), replace = F) # reshuffle so it's not always the same 75 cell barcodes

        saveRDS(grlist, here("cache", "hg19_novoalign", 
                             "scramble",
                             str_glue("{plate}_sc.rds")))
    })

# pool PX0490 and PX0484 (both from patient 371) together
# this patient had two Pre-Aza plates whereas all the other patients only had one
px0490 <- readRDS(here("cache", "hg19_novoalign", "scramble", "PX0490_sc.rds"))
px0484 <- readRDS(here("cache", "hg19_novoalign", "scramble", "PX0484_sc.rds"))

# shuffle order so when we do indices < 75 that they are from both plates
px0490_px0484 <- sample(append(px0490, px0484), length(px0490) + length(px0484), replace = F)
saveRDS(px0490_px0484, here("cache", "hg19_novoalign", "scramble", "PX0490_PX0484_sc.rds"))

invisible(file.remove(here("cache", "hg19_novoalign", "scramble", "PX0484_sc.rds")))
invisible(file.remove(here("cache", "hg19_novoalign", "scramble", "PX0490_sc.rds")))

system("Rscript src/make_library_scramble.R cache/hg19_novoalign/scramble cache/hg19_novoalign/scramble/PX*_sc.rds")

data.frame(
    paste0("cache/hg19_novoalign/scramble/scrambled_library_", 1:7, ".rds"), 
    c(paste0("pt_S0", 1:4, "_R", 1:4), paste0("pt_S0", 5:7, "_NR", 1:3))
) %>%
    write.table(file = "data/manifest.tsv", quote = F, row.names = F, col.names = F)

system('parallel --col-sep " " Rscript src/calc_methylation_by_genomic_region.R {1} cache/hg19_annotation.rds {2} :::: data/manifest.tsv')

scramble_avg_meth_by_region_w_metadata.tib <- list.files(here::here("cache", "hg19_novoalign", "scramble"), "_by_region.rds", full.names = T) %>%
    set_names(str_remove(basename(.), "pt_") %>%
                  str_remove("_methylation.*")) %>%
    map(~ readRDS(.))
