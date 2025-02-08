# INPUT: 
# `resources/supp_data/Farlik_2016_meth`

# OUTPUT:
# `results/farlik_dmls.rds`
# `results/farlik__dmls_gr.rds`
# `results/farlik_dmls.bed`
# `results/farlik_dmrs.rds`
# `results/farlik_dmrs_gr.rds`
# `results/farlik_dmrs.bed`

library(GenomicRanges)
library(rtracklayer)
library(DSS)
library(here)
library(parallel)
library(pbmcapply)
library(tidyverse)

# # Read in only the relevant ones for analysis
files <- list(MPP = list.files(here("resources", "supp_data", "Farlik_2016_meth"), pattern = "MPP_", recursive = FALSE, full.names = TRUE),
              MLP0 = list.files(here("resources", "supp_data", "Farlik_2016_meth"), pattern = "MLP0_", recursive = FALSE, full.names = TRUE))

load_data <- function(f) {
    patient_grp <- str_extract(f, "_D[0-9]+_") %>%
        str_remove_all("_")
    if(str_detect(f, ".gz$")) {
        gz <- gzfile(f)
    } else {
        gz <- f
    }
    temp <- read_delim(gz,
                       col_names = c("chr", "pos", "X", "N"),
                       show_col_types = FALSE) %>%
        mutate(`Patient Group` = patient_grp) %>%
        relocate(N, .before = X)
    return(temp)
}

# Loading data
dfs <- imap(files, function(file, celltype) {
    tmp <- pbmclapply(file, load_data, mc.cores = length(file), ignore.interactive = TRUE) %>%
        purrr::reduce(bind_rows) %>%
        group_by(chr, pos, `Patient Group`) %>%
        summarize(N = sum(N),
                  X = sum(X)) %>%
        ungroup() %>%
        split(f = as.factor(.$`Patient Group`), drop = TRUE) %>%
        map(~ select(.x, -`Patient Group`))
    names(tmp) <- names(tmp) %>%
        str_replace("^", str_glue("{celltype}_"))
    return(tmp)
})

chain <- import.chain(here("resources", "supp_data", "liftover", "hg38ToHg19.over.chain"))

farlik <- append(dfs$MPP, dfs$MLP0)

# Lifting over 
farlik <- map(farlik, function(df) {
    gr <- as.data.frame(df) %>%
        mutate(end = pos + 1) %>%
        makeGRangesFromDataFrame(seqnames.field = "chr", start.field = "pos", end.field = "end", keep.extra.columns = T)

    df.hg19 <- rtracklayer::liftOver(gr, chain) %>%
        unlist() %>%
        data.frame() %>%
        dplyr::select(chr = seqnames, pos = start, N, X) %>%
        group_by(chr, pos) %>%
        summarize(X = sum(X), N = sum(N))

    return(df.hg19)
    })

farlik <- map(farlik, function(df) {
    df.hg19 <- df %>%
        dplyr::select(chr = seqnames, pos = start, N, X) %>%
        group_by(chr, pos) %>%
        summarize(X = sum(X), N = sum(N))

    return(df.hg19)
})

# Creating BSseq object
BSobj <- makeBSseqData(farlik, names(farlik))

# Starting DML testing
dmlTest <- DMLtest(BSobj, group1 = str_subset(sampleNames(BSobj), "^MPP_"), 
                   group2 = str_subset(sampleNames(BSobj), "^MLP0_"), 
                   ncores = 128, smoothing = TRUE)

# Calling DMLs
dmls <- callDML(dmlTest, p.threshold = 0.001)

# Caching DMLs
saveRDS(dmls, here("results", "farlik_dmls.rds"))

# Converting DMLs to GenomicRanges
dmls.gr <- makeGRangesFromDataFrame(df = dmls,
                                        keep.extra.columns = T, seqnames.field = "chr", start.field = "pos", end.field = "pos")
mcols(dmls.gr)$dml_id <- paste0("DML_", str_pad(seq(1, nrow(dmls)), str_length(nrow(dmls)), "left", pad = "0"))

# Caching GenomicRanges
saveRDS(dmls.gr, here("results", "farlik_dmls_gr.rds"))

# Exporting GenomicRanges as BED file
rtracklayer::export.bed(dmls.gr, here::here("results", "farlik_dmls.bed"))

# Calling DMRs
dmrs <- callDMR(dmlTest, p.threshold = 0.01)
saveRDS(dmrs, here("results", "farlik_dmrs.rds"))

# Converting DMRs to GenomicRanges
dmrs.gr <- makeGRangesFromDataFrame(df = dmrs,
                                        keep.extra.columns = T, seqnames.field = "chr")

mcols(dmrs.gr)$dmr_id <- paste0("DMR_", str_pad(seq(1, nrow(dmrs)), 5, "left", pad = "0"))

# GenomicRanges
saveRDS(dmrs.gr, here("results", "farlik_dmrs_gr.rds"))

# Calling DMRs
dmrs <- callDMR(dmlTest, p.threshold = 0.001)
saveRDS(dmrs, here("results", "farlik_dmrs.rds"))

# Converting DMRs to GenomicRanges
dmrs.gr <- makeGRangesFromDataFrame(df = dmrs,
                                    keep.extra.columns = T, seqnames.field = "chr")

mcols(dmrs.gr)$dmr_id <- paste0("DMR_", str_pad(seq(1, nrow(dmrs)), str_length(nrow(dmrs)), "left", pad = "0"))

# Caching GenomicRanges
saveRDS(dmrs.gr, here("results", "farlik_dmrs_gr.rds"))