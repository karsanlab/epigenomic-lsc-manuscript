# INPUT: 
# `data/scWGBS/*.cons.Cmethyl.CpGs.txt`
# `data/full_metadata_qc.rds`
# `cache/regionAnnot.sorted.bed`

# OUTPUT: 
# `results/*_farlik_region_meth.bed`
# `results/*_farlik_region_meth.tsv`

library(rtracklayer)
library(tidyverse)
library(here)

## Functions
convert_novoalign_to_bed <- function(file, dir, temp_file, region_file) {
    sample <- read_tsv(file.path(dir, file), 
                          col_names = c("chr", "start", "end", "percent_meth", "n_meth", "n_total")) %>%
        transmute(chr, start, end, score = percent_meth/100) %>%
        makeGRangesFromDataFrame(keep.extra.columns = T)
    
    subset <- subsetByOverlaps(sample, region_file) %>%
      sort()
    
    export.bed(subset, temp_file)
}

in.dir <- here("data", "scWGBS")
out.dir <- "results"
sample.file <- here("data", "full_metadata_qc.rds")
query.row <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
farlik.regions <- here("cache", "regionAnnot.sorted.bed")

## Load Farlik supp data
farlik_region_gr <- import.bed(farlik.regions)
    

# make output dir if it doesn't exist
if(!file.exists(out.dir)){system(str_glue("mkdir {out.dir}"))}

## Load sample metadata ##
sample_data <- readRDS(sample.file) %>%
  select(patient, Plate_ID, cell, plate_sample)

id <- sample_data %>%
    dplyr::slice(query.row) %>%
    dplyr::select(patient, plate_sample) %>%
    dplyr::transmute(id = paste0(patient, "_", plate_sample)) %>%
    pull(id)

in_bed_file <- paste0(id, ".cons.Cmethyl.CpG.txt")

message(paste0("Loading ", in_bed_file, " and calculating farlik regional methylation."))

# Load our the bw file data for the row
tmp <- tempfile(fileext = ".bed")
tmp.sorted <- str_replace(tmp, ".bed", ".sorted.bed")

convert_novoalign_to_bed(file = in_bed_file, dir = in.dir, temp_file = tmp, 
                  region_file = farlik_region_gr)

file.name <- paste0(id, "_farlik_region_meth.bed")

# run bedtools from within R
message("running bedtools")
system(str_glue("bedtools sort -i {tmp} > {tmp.sorted}"))
system(str_glue("bedtools map -a cache/regionAnnot.sorted.bed -b {tmp.sorted} -o mean > {file.path(out.dir, file.name)}"))

invisible(file.remove(tmp))
invisible(file.remove(tmp.sorted))

message("writing tsv")
bed_file <- read_tsv(file.path(out.dir, file.name),
                     col_names = c("chrom", "start", "end", "name", 
                                  "score", "score2", "prop_meth"),
                     na = c(".", "NA")) %>%
  select(prop_meth)

colnames(bed_file) <- stringr::str_remove(id, "[0-9]{3}_")

file.name <- paste0(id, "_farlik_region_meth.tsv")

write_tsv(bed_file, file.path(out.dir, file.name))
