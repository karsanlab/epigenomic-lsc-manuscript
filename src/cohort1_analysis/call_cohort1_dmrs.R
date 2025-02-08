# INPUT: 
# `cache/hg19_novoalign_combined_methylation_global.tibs.rds`
# `cache/passed_qc.rds`
# `data/metadata.csv`

# OUTPUT: 
# `results/cohort1_dmls.rds`
# `results/cohort1_dmls_gr.rds`
# `results/cohort1_dmls.bed`
# `results/cohort1_dmrs.rds`
# `results/cohort1_dmrs_gr.rds`
# `results/cohort1_dmrs.bed`

suppressPackageStartupMessages(
    {
        library(GenomicRanges)
        library(DSS)
        library(here)
        library(parallel)
        library(rtracklayer)
        library(tidyverse)
    }
)

passed_qc <- readRDS(here("resources", "passed_qc.rds")) %>%
	separate(Sample_Barcode, into = c("Plate_ID", "cell"), sep = "_")

# re summarize  without cell info since that's not needed for DMR calling but needed for the % Methylation volcano plot
 combined_methylation.tibs <- read_rds(
            here::here("cache", str_glue("hg19_novoalign_combined_methylation.tibs.rds"))) %>%
 	imap(function(tib, plate_id) {
 			tib %>%
				mutate(Plate_ID = plate_id) %>%
				semi_join(passed_qc, by = c("Plate_ID", "cell")) %>%
 				group_by(chr, pos) %>%
 				summarize(N = sum(N), X = sum(X))
    })

# Pooling PX0490 and PX0484 for patient 371
combined_methylation.tibs[["PX0490_PX0484"]] <-  bind_rows(combined_methylation.tibs[["PX0490"]], combined_methylation.tibs[["PX0484"]]) %>%
    group_by(chr, pos) %>%
    summarize(N = sum(N), X = sum(X))

combined_methylation.tibs <- discard_at(combined_methylation.tibs, c("PX0490", "PX0484"))

# Loading metadata
aza_metadata.tib <- read_csv(here("data", "metadata.csv"))

RESPONDER_NOS <- aza_metadata.tib %>%
    filter(Class == "Responder") %>%
    pull(Plate_ID) %>%
    intersect(names(combined_methylation.tibs)) %>%
    c("PX0490_PX0484") # manually add this as a responder

NONRESPONDER_NOS <- aza_metadata.tib %>%
    filter(Class == "Non-responder") %>%
    pull(Plate_ID) %>%
    intersect(names(combined_methylation.tibs))
# collapse all Responders into one tibble in the list

# Creating BSseq object
aza.bsseq <-
    makeBSseqData(combined_methylation.tibs,
                  names(combined_methylation.tibs))

# Running DMLtest
aza.dml_test <- DMLtest(aza.bsseq,
                        group1 = RESPONDER_NOS,
                        group2 = NONRESPONDER_NOS,
                        smoothing = TRUE, 
                        ncores = 48)

# Calling DMLs
aza.dmls <- callDML(aza.dml_test, p.threshold = 0.001)

# Caching DMLs
saveRDS(aza.dmls,
        here::here("results", "patient_dmls.rds"))

# Converting DMLs to GenomicRanges
aza_dmls.gr <- makeGRangesFromDataFrame(df = aza.dmls,
                                        keep.extra.columns = T, seqnames.field = "chr", start.field = "pos", end.field = "pos")

mcols(aza_dmls.gr)$dml_id <- paste0("DML_", str_pad(seq(1, nrow(aza.dmls)), str_length(nrow(aza.dmls)), "left", pad = "0"))

# Caching GenomicRanges
saveRDS(aza_dmls.gr,
        here::here("results", "patient_dmls_gr.rds"))

# Exporting GenomicRanges as BED file
rtracklayer::export.bed(aza_dmls.gr, here::here("results", "patient_dmls.bed"))

# Calling DMRs
aza.dmrs <- callDMR(aza.dml_test, p.threshold = 10e-4)

# Caching DMRs
saveRDS(aza.dmrs,
        here::here("results", "patient_dmrs.rds"))

# Converting DMRs to GenomicRanges
aza_dmrs.gr <- makeGRangesFromDataFrame(df = aza.dmrs,
                                        keep.extra.columns = T, seqnames.field = "chr")

mcols(aza_dmrs.gr)$dmr_id <- paste0("DMR_", str_pad(seq(1, nrow(aza.dmrs)), str_length(nrow(aza.dmrs)), "left", pad = "0"))

# Caching GenomicRanges
saveRDS(aza_dmrs.gr,
        here::here("results", "patient_dmrs_gr.rds"))

# Exporting GenomicRanges as BED file
rtracklayer::export.bed(aza_dmrs.gr, here::here("results", "patient_dmrs.bed"))
