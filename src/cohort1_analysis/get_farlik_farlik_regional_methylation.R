# INPUT: 
# `data/farlik/GSM*.cons.Cmethyl.CpGs.txt`
# `resources/supp_data/Farlik_cell_type_prediction/sampleAnnot.tsv`
# `cache/regionAnnot.sorted.bed`

# OUTPUT:
# `results/GSM*_farlik_region_meth.bed`
# `results/GSM*_farlik_region_meth.tsv`

library(rtracklayer)
library(here)
library(tidyverse)

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

in.dir <- here("data", "farlik")
out.dir <- "results"
sample.file <- here("resources", "supp_data", "Farlik_cell_type_prediction/sampleAnnot.tsv")
query.row <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
farlik.regions <- here("cache", "regionAnnot.sorted.bed")

## Load Farlik supp data
farlik_region_gr <- import.bed(farlik.regions)

# make output dir if it doesn't exist
if(!file.exists(out.dir)){system(str_glue("mkdir {out.dir}")}

## Load sample metadata ##
id <- read_tsv(sample.file) %>%
    dplyr::slice(query.row) %>%
    pull(sampleName)
    
# id from the sample file doesn't have the GEO ID, need to use list.files to get the GEO_id
in_bed_file <- list.files(in.dir, pattern = str_glue("_{id}\\.")) %>%
    str_subset(".cons.Cmethyl.CpG.txt")

if(length(in_bed_file) != 1) {
    stop("There should only be one file that matches.")
}

# overwrite ID with the one with GSM? 
geo.id <- str_remove(in_bed_file, ".cons.Cmethyl.CpG.txt")

message(paste0("Loading ", in_bed_file, " and calculating farlik regional methylation."))

# Load our the bw file data for the row
tmp <- tempfile(fileext = ".bed")
tmp.sorted <- str_replace(tmp, ".bed", ".sorted.bed")

convert_novoalign_to_bed(file = in_bed_file, dir = in.dir, temp_file = tmp, 
                  region_file = farlik_region_gr)

file.name <- paste0(geo.id, "_farlik_region_meth.bed")

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

colnames(bed_file) <- id

file.name <- paste0(geo.id, "_farlik_region_meth.tsv")

write_tsv(bed_file, file.path(out.dir, file.name))
