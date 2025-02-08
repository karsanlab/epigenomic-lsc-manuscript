# INPUT: 
# `results/cohort1_dmrs_gr.rds`

# OUTPUT: 
# `results/cohort1_dmrs_tad_promoter_enhancer_annotated.rds`

library(LOLA)
library(annotatr)
library(AnnotationHub)
library(GenomicRanges)
library(here)
library(tidyverse)

# Functions ---- 
classify_tf <- function(tf) {
    tf <- case_when(str_detect(tf, "(?i)znf143") ~ "ZNF143",
                    str_detect(tf, "(?i)ctcf") ~ "CTCF", 
                    TRUE ~ "Other")
    return(tf)
}

# LOLA hg19 ----
lola.db.hg19 <- loadRegionDB(here("resources", "supp_data", "LOLA", "LOLACore", "hg19"),
                             collections = c('encode_tfbs', 'codex',
                                             'cistrome_cistrome'))

lola_tfbs_metadata.tib.hg19 <- lola.db.hg19$regionAnno %>%
    as_tibble()

lola_tfbs.grlist.hg19  <- lola.db.hg19$regionGRL
names(lola_tfbs.grlist.hg19) <- lola_tfbs_metadata.tib.hg19$filename

lola_tfbs.grlist.hg19 <- mapply(function(x, y, z, a, b) {
    mcols(x) <- DataFrame(filename = y, antibody = z, celltype = a, treatment = b)
    return(x)
}, lola_tfbs.grlist.hg19, names(lola_tfbs.grlist.hg19), lola_tfbs_metadata.tib.hg19$antibody, lola_tfbs_metadata.tib.hg19$cellType, lola_tfbs_metadata.tib.hg19$treatment) %>%
    GRangesList()

lola_tfbs.grs.hg19 <- unlist(lola_tfbs.grlist.hg19)

names(lola_tfbs.grs.hg19) <- NULL
lola_annotation.hg19 <- lola_tfbs.grs.hg19 %>%
    plyranges::mutate(id = 1:length(lola_tfbs.grs.hg19),
                      gene_id = antibody,
                      symbol = classify_tf(antibody),
                      type = "hg19_custom_lola_tfbs")

annotatr_cache$set("hg19_custom_lola_tfbs", lola_annotation.hg19)

# AML3 hg38 TADs ----
aml3.tad.gr.hg38 <- rtracklayer::import(here("data", "AML3_shControl_topdom_tads.bed")) %>%
    dplyr::mutate(
        id = name,
        gene_id = name,
        symbol = str_glue("AML3 TAD {name}"),
        type = "hg38_custom_aml3_tads"
    ) %>%
    dplyr::select(
        id, gene_id, symbol, type
    )

liftover <- function(gr) {
  rtracklayer::liftOver(gr, chain) %>%
    unlist() 
}

# AML3 hg19 TADs ----
aml3.tad.gr.hg19 <- liftover(aml3.tad.gr.hg38) %>%
    plyranges::mutate(type = "hg19_custom_aml3_tads")

annotatr_cache$set("hg19_custom_aml3_tads", aml3.tad.gr.hg19)

# Liftover ----
chain <- import.chain(here("resources", "supp_data", "liftover", "hg38ToHg19.over.chain"))

# Enhancers hg19 ---- 
# http://www.enhanceratlas.org/data/AllEPs/hs/CD34+_EP.txt
enhancers.gr.hg19 <-
    read_tsv(here::here("data", "CD34+_EP.txt"),
             col_names = c("metadata", "enhancer_gene_score")) %>%
    separate("metadata", into = c("seqnames", "coords", "gene_data"),
             sep = ":|_") %>%
    separate("coords", into = c("start", "end"), sep = "-") %>%
    separate("gene_data", into = c("gene_ens_id", "gene_ens_name",
                                   "gene_seqname", "gene_tss", "gene_strand"),
             sep = "\\$") %>%
    makeGRangesFromDataFrame(keep.extra.columns = T) %>%
    plyranges::mutate(
        id = gene_ens_id,
        gene_id = gene_ens_name,
        symbol = "enhancers",
        type = "hg19_custom_enhancers"
    ) %>%
    dplyr::select(
        id, gene_id, symbol, type
    )

annotatr_cache$set("hg19_custom_enhancers", enhancers.gr.hg19)

# hg19 annotation ----
hg19.annotations <- build_annotations(genome = "hg19",
                  annotations = c("hg19_genes_promoters",
                                  "hg19_custom_enhancers",
                                  "hg19_custom_lola_tfbs",
                                  "hg19_custom_aml3_tads"))

# DMRs (hg19) ----
dmrs.gr.hg19 <- readRDS(here("data", "cohort1_dmrs_gr.rds"))

dmrs.annotated.tib.hg19 <- dmrs.gr.hg19 %>%
  annotate_regions(
    annotations = hg19.annotations
  ) %>%
  data.frame() %>%
  dplyr::select(-c(annot.seqnames, 
                   nnot.start, 
                   annot.end, 
                   annot.width, 
                   annot.strand, 
                   annot.gene_id, 
                   annot.filename, 
                   annot.antibody, 
                   annot.celltype, 
                   annot.treatment, 
                   annot.id, 
                   annot.tx_id)) %>%
  pivot_wider(names_from = "annot.type", values_from = "annot.symbol") %>%
  mutate(hg19_custom_lola_tfbs_factor = map_chr(hg19_custom_lola_tfbs, ~ {
    if("CTCF" %in% .x & "ZNF143" %in% .x) {
      return("CTCF and ZNF143")
    } else if("CTCF" %in% .x) {
      return("CTCF without ZNF143")
    } else if("ZNF143" %in% .x) {
      return("ZNF143 without CTCF")
    } else if ("Other" %in% .x) {
      return("Other binding site")
    } else if (is.null(.x)) {
      return(NA)
    }
  }),
  hg19_custom_aml3_tads_factor = map_chr(hg19_custom_aml3_tads, ~ {
    if(is.null(.x)) {
      return(NA)
    } else {
      return(
      unique(.x) %>%
      paste(collapse = ", "))
    }
  }),
  hg19_genes_promoters_factor = map_chr(hg19_genes_promoters, ~ {
    if(is.null(.x)) {
      return(NA)
    } else {
      return("promoter")
    }
  }),
  hg19_custom_enhancers_factor = map_chr(hg19_custom_enhancers, ~ {
    if(is.null(.x)) {
      return(NA)
    } else {
      return("enhancer")
    }
  })
  )

saveRDS(dmrs.annotated.tib.hg19, here("results", "cohort1_dmrs_tad_promoter_enhancer_annotated.rds"))
