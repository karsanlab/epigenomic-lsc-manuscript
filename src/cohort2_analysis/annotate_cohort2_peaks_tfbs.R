# INPUT:
# `results/integrated_cohort2_samples.rds`
# `results/top_dmtfbs.tsv`

# OUTPUT:
# `cache/cohort2_tfbs_peaks.rds`
# `cache/cohort2_umap.rds`

library(LOLA)
library(annotatr)
library(Seurat)
library(Signac)
library(here)
library(tidyverse)

# Get peaks ----
obj <- readRDS(here("results", "integrated_cohort2_samples.rds"))

DefaultAssay(obj) <- "peaks"
gr <- granges(obj)

# Load LOLA ----
LOLA_PATH = here("resources", "supp_data", "LOLA", "LOLACore", "hg38")

lola.db <- LOLA::loadRegionDB(LOLA_PATH,
                              collections = c('encode_tfbs', 'codex', 'cistrome_cistrome', 'custom'))

# load LOLA annotations
lola_tfbs_metadata.tib <- lola.db$regionAnno %>%
    as_tibble()

# load actual LOLA GRanges
lola_tfbs.grlist <- lola.db$regionGRL

names(lola_tfbs.grlist) <- lola_tfbs_metadata.tib$filename

lola_tfbs.grlist <- mapply(function(x, y, z, a, b) {
    mcols(x) <- DataFrame(filename = y, antibody = z, celltype = a, treatment = b)
    return(x)
}, lola_tfbs.grlist, names(lola_tfbs.grlist), lola_tfbs_metadata.tib$antibody, lola_tfbs_metadata.tib$cellType, lola_tfbs_metadata.tib$treatment) %>%
    GRangesList()

lola.annotation <- unlist(lola_tfbs.grlist)
names(lola.annotation) <- NULL

most.sig <- read_csv(here("results", "top_dmtfbs.tsv"))

# filter for the 64 of interest
lola.annotation <- lola.annotation %>%
    plyranges::filter(filename %in% most.sig$filename)

# annotate peaks for the 64 of interest
annotated.gr <- annotate_regions(gr, lola.annotation, quiet = F, ignore.strand = T)

annotated.df <- annotated.gr %>%
    data.frame() %>%
    select(seqnames, start, end, filename = annot.filename) %>%
    left_join(select(most.sig, filename, tf), by = "filename") %>%
    select(-filename) %>%
    unite(region, seqnames, start, end, sep = "-") %>%
    distinct()

# get raw peaks count data
peaks.raw <- GetAssayData(obj, slot = "counts", assay = "peaks") %>%
    as.data.frame() %>%
    rownames_to_column("region") %>%
    inner_join(annotated.df, by = "region") %>%
    select(-region) %>%
    distinct()

# split raw peak data by the TFs of interest
peaks.raw <- group_by(peaks.raw, tf)

split <- group_split(peaks.raw, .keep = T)

names(split) <- group_keys(peaks.raw) %>%
    pull(tf)

# summarize counts
df.list <- imap(split, ~ {
    .x %>%
        pivot_longer(-tf, names_to = "cell", values_to = "count") %>%
        group_by(tf, cell) %>%
        summarize(count = sum(count)) %>%
        pivot_wider(id_cols = tf, names_from = "cell", values_from = "count")
})

# create a matrix where the rows are transcription factors, and the columns are the cells
mtx <- df.list %>%
    list_rbind() %>%
    column_to_rownames("tf") %>%
    as.matrix() %>%
    Matrix::Matrix(sparse = TRUE)

saveRDS(mtx, here("cache", "cohort2_tfbs_peaks.rds"))

# perform UMAP dimension reduction
umap <- umap(mtx)

saveRDS(umap, here("cache", "cohort2_umap.rds"))
