# INPUT:
# `cache/hg19_novoalign/scramble/scrambled_library_*.rds`
# `cache/hg19_annotation.rds`

# OUTPUT:
# `cache/hg19_novoalign/scramble/*_methylation_by_region.rds`

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(tidyverse)
    library(annotatr)
    library(furrr)
})

# PARSE ARGUMENTS ----

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop(paste0("Path to single-cell methylation (1) and genomic regions .rds (2), ",
                "and library ID (3) must be provided."),
         call. = FALSE)
} else {
    SC_GRS     <- args[1]
    REGION_GRS <- args[2]
    LIB_ID     <- args[3]
}

# PARAMETERS ----

plan(multisession, workers = 8)

# MAIN ----

## IMPORT RDS ----

lib_methylation.grs <- read_rds(SC_GRS)

genomic_regions.grs <- read_rds(REGION_GRS)

## CALCULATE INTERSECTIONS ----

avg_methyl_by_region.tib <-
    future_map(lib_methylation.grs,
        function(lib.gr) {
            global.tib <-
                tibble(n_methyl = sum(lib.gr$n_methyl),
                       n_total = sum(lib.gr$n_total),
                       pct_methyl = n_methyl / n_total,
                       region = "global")
            regions.tib <-
                annotate_regions(lib.gr, annotations = genomic_regions.grs) %>%
                data.frame() %>%
                filter(!annot.type %in% c("hg19_genes_intergenic", "hg19_cpg_inter")) %>%
                dplyr::mutate(
                    region = case_match(annot.type,
                                        "hg19_cpg_islands" ~ "CpG islands",
                                        "hg19_genes_3UTRs" ~ "genes",
                                        "hg19_genes_5UTRs" ~ "genes",
                                        "hg19_cpg_shores" ~ "CpG shores",
                                        "hg19_cpg_shelves" ~ "CpG shelves",
                                        "hg19_genes_1to5kb" ~ "genes",
                                        "hg19_enhancers" ~ "enhancers",
                                        "hg19_genes_exons" ~ "genes",
                                        "hg19_genes_introns" ~ "genes",
                                        "hg19_genes_promoters" ~ "promoters"
                    )
                ) %>%
                group_by(region) %>%
                summarize(n_methyl = sum(n_methyl),
                          n_total = sum(n_total),
                          pct_methyl = n_methyl / n_total)
            df <- bind_rows(global.tib, regions.tib) %>%
                mutate(region = factor(region, levels = c("global", "genes", "promoters", "enhancers", "CpG islands", "CpG shores", "CpG shelves")))

            return(df)
        }) %>%
    bind_rows(.id = "cell")

## EXPORT to RDS ----
saveRDS(avg_methyl_by_region.tib,
          paste0(dirname(SC_GRS), "/", LIB_ID, "_methylation_by_region.rds"))
