# INPUT: 
# `results/hg19_novoalign/*_cpg_pool.tib.rds`

# OUTPUT: 
# `cache/hg19_annotation.rds`
# `cache/hg19_novoalign/*_annotated_CpGs.gr.rds`
# `cache/hg19_novoalign/*_annotated_CpGs.tbl.rds`

library(here) 
library(GenomicRanges)
library(annotatr)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

sample <- args[1]
rds <- args[2]
if (length(args) == 3) {
    subdir <- args[3]
} else {
    subdir <- ""
}

message(str_glue("Working on {sample}"))
methylation.tibble <- readRDS(rds) %>%
    dplyr::select(
        "chr" = seqnames,
        "pos" = start,
        "N" = n_total,
        "X" = n_methylated,
        cell
    )

enhancers.hg19 <-
  read_tsv(here::here("data", "enhancer_atlas", "CD34+_EP.txt"),
           col_names = c("metadata", "enhancer_gene_score")) %>%
  separate("metadata", into = c("seqnames", "coords", "gene_data"),
           sep = ":|_") %>%
  separate("coords", into = c("start", "end"), sep = "-") %>%
  separate("gene_data", into = c("gene_ens_id", "gene_ens_name",
                                 "gene_seqname", "gene_tss", "gene_strand"),
           sep = "\\$") %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)
enhancer_annotation.hg19 <- enhancers.hg19 %>%
    plyranges::mutate(
        id = gene_ens_id,
        gene_id = gene_ens_name,
        symbol = "enhancers",
        type = "hg19_enhancers"
    ) %>%
    dplyr::select(
        id, gene_id, symbol, type
    )

annotatr_cache$set("hg19_enhancers", enhancer_annotation.hg19)
custom_annotation.hg19 <- build_annotations(genome = "hg19",
                    annotations = c("hg19_basicgenes", "hg19_genes_promoters",
                                    "hg19_cpg_islands", "hg19_cpg_shores",
                                    "hg19_cpg_shelves", "hg19_enhancers"))

saveRDS(custom_annotation.hg19, here("cache", "hg19_annotation.rds"))

message("Converting to Genomic Ranges")
gr <- methylation.tibble %>%
    dplyr::rename(start = pos) %>%
    dplyr::mutate(end = start + 1) %>%
    dplyr::relocate(end, .after = "start") %>%
    makeGRangesFromDataFrame(ignore.strand = T,
                             keep.extra.columns = T)

message("Annotating regions")
annotated.gr <- annotate_regions(
    regions = gr,
    annotations = custom_annotation.hg19,
    ignore.strand = T,
    quiet = F
)

message("Converting annotation to data frame")
annotated.tibble <- annotated.gr %>%
    data.frame() %>%
    dplyr::select(seqnames, cell, start, end, N, X, annot.type)


message("Caching annotated dataframe")
saveRDS(annotated.tibble, here("cache", subdir, str_glue("hg19_{sample}_annotated_CpGs.gr.rds")))


message("Calculating percent methylation")
annotated_updated.tibble <- annotated.tibble %>%
    dplyr::mutate(
        annotation = case_when(
            annot.type == "hg19_genes_promoters" ~ "promoters",
            str_starts(annot.type, "hg19_genes_") ~ "genes",
            annot.type == "hg19_enhancers" ~ "enhancers",
            str_starts(annot.type, "hg19_cpg_") ~ str_remove(annot.type, "hg19_cpg_"),
            TRUE ~ annot.type
        )
    )  %>%
    group_by(cell, annotation) %>%
    summarize(percent_meth = sum(X) / sum(N) * 100)

message("Caching percent methylation")
saveRDS(annotated_updated.tibble, here("cache", subdir, str_glue("hg19_{sample}_annotated_CpGs.tbl.rds")))