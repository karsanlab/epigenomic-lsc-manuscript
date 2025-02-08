# INPUT: 
# `results/Cohort2_FindMarkers_DARs.rds`

# OUTPUT: 
# `results/cohort2_dars_tad_promoter_enhancer_annotated.rds`

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

# LOLA hg38 ----
lola.db.hg38 <- LOLA::loadRegionDB(here("resources", "supp_data", "LOLA", "LOLACore", "hg38"),
                                   collections = c('encode_tfbs', 'codex', 'cistrome_cistrome'))
lola_tfbs_metadata.tib.hg38 <- lola.db.hg38$regionAnno %>%
    as_tibble()

lola_tfbs.grlist.hg38  <- lola.db.hg38$regionGRL
names(lola_tfbs.grlist.hg38) <- lola_tfbs_metadata.tib.hg38$filename

lola_tfbs.grlist.hg38 <- mapply(function(x, y, z, a, b) {
    mcols(x) <- DataFrame(filename = y, antibody = z, celltype = a, treatment = b)
    return(x)
}, lola_tfbs.grlist.hg38, names(lola_tfbs.grlist.hg38), lola_tfbs_metadata.tib.hg38$antibody, lola_tfbs_metadata.tib.hg38$cellType, lola_tfbs_metadata.tib.hg38$treatment) %>%
    GRangesList()

lola_tfbs.grs.hg38 <- unlist(lola_tfbs.grlist.hg38)

names(lola_tfbs.grs.hg38) <- NULL
lola_annotation.hg38 <- lola_tfbs.grs.hg38 %>%
    mutate(id = 1:length(lola_tfbs.grs.hg38),
           gene_id = antibody,
           symbol = classify_tf(antibody),
           type = "hg38_custom_lola_tfbs")

annotatr_cache$set("hg38_custom_lola_tfbs", lola_annotation.hg38)

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

annotatr_cache$set("hg38_custom_aml3_tads", aml3.tad.gr.hg38)


# Liftover ----
ahub.chain <- subset(AnnotationHub(), rdataclass == "ChainFile" & species == "Homo sapiens")
hg19.to.hg38.chain <- ahub.chain[ahub.chain$title == "hg19ToHg38.over.chain.gz"]
hg19.to.hg38.chain <- hg19.to.hg38.chain[[1]]


liftover <- function(gr) {
    rtracklayer::liftOver(gr, hg19.to.hg38.chain) %>%
        unlist() 
}

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

# Enhancers hg38 ----
enhancers.gr.hg38 <- liftover(enhancers.gr.hg19) %>%
    plyranges::mutate(type = "hg38_custom_enhancers")

annotatr_cache$set("hg38_custom_enhancers", enhancers.gr.hg38)

# hg38 annotation ----
hg38.annotations <- build_annotations(genome = "hg38",
                                      annotations = c("hg38_genes_promoters",
                                                      "hg38_custom_enhancers",
                                                      "hg38_custom_lola_tfbs",
                                                      "hg38_custom_aml3_tads"))

# DARS (hg38) ----
dars.gr.hg38 <- readRDS(here("data", "FindMarkers_DARs.rds")) %>%
    rownames_to_column("region") %>%
    separate(region, into = c("seqnames", "start", "end"), sep = "-") %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)

dars.annotated.tib.hg38 <- dars.gr.hg38 %>%
    annotate_regions(
        annotations = hg38.annotations
    ) %>%
    data.frame() %>%
    dplyr::select(-c(annot.seqnames, 
                     annot.start, 
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
    mutate(hg38_custom_lola_tfbs_factor = map_chr(hg38_custom_lola_tfbs, ~ {
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
    hg38_custom_aml3_tads_factor = map_chr(hg38_custom_aml3_tads, ~ {
        if(is.null(.x)) {
            return(NA)
        } else {
            return(
                unique(.x) %>%
                    paste(collapse = ", "))
        }
    }),
    hg38_genes_promoters_factor = map_chr(hg38_genes_promoters, ~ {
        if(is.null(.x)) {
            return(NA)
        } else {
            return("promoter")
        }
    }),
    hg38_custom_enhancers_factor = map_chr(hg38_custom_enhancers, ~ {
        if(is.null(.x)) {
            return(NA)
        } else {
            return("enhancer")
        }
    })
    )

saveRDS(dars.annotated.tib.hg38, here("results", "cohort2_dars_tad_promoter_enhancer_annotated.rds"))
