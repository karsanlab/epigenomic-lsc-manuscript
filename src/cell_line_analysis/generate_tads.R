#!/usr/bin/env Rscript

# INPUT:
# results/${sample}_contact_map.hic

# OUTPUT:
# results/${sample}_topdom_binSignal.tsv
# results/${sample}_topdom_domain.tsv
# results/${sample}_topdom_tads.bed

# Use HiCDCplus to generate TADs -------------------------------------------------------------------

### Set up args ------------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
    stop("At least one argument must be supplied (input file)", call.=FALSE)
}

### Setup R environment ----------------------------------------------------------------------------
library(plyranges)
library(HiCDCPlus)
library(tidyverse)
library(here)

### Load command line arguments --------------------------------------------------------------------
sample <- args[1]
sample 

### Set paths --------------------------------------------------------------------------------------
# hicdc+ does not like paths made with the "here" package
# doing it the old fashioned paste0 way
outdir <- "./results"
hic_path <- paste0("./data/", sample, "_contact_map.hic")

## Generate ice normalized gi_list for calling tads ------------------------------------------------
# Although HiCDCplus has a chrs = null option that should run for all chrs but it does not work. 

# Instead supplying a list of all chrs. Note it doesn't work for chrM and chrY. 
sequence_names<-c("chr1", "chr2", "chr3", "chr4", "chr5","chr6", "chr7", "chr8", "chr9", "chr10",
                  "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17","chr18","chr19",
                  "chr20","chr21","chr22","chrX")
    
gi_list_ice <- hic2icenorm_gi_list(hic_path,
                                   binsize = 50e3,
                                   chrs = sequence_names,
                                   Dthreshold = 400e3,
                                   gen_ver = "hg38", 
                                   gen = "Hsapiens")

# Write normalized counts for chr17 to file as gi_list instance. 
gi_list_write(gi_list_ice,
              fname = str_glue('./results/{sample}_ice_normalized_counts.txt'),
              chrs = "chr17",
              columns = "minimal",
              rows = "all")

# Write ice normalized counts to hic for use in plotting
hicdc2hic(gi_list_ice, 
          hicfile = str_glue('./results/{sample}_ice_normalized_counts.hic'),
          gen_ver = "hg38",memory = 64, mode = "raw")

## Call tads on the ice normalized hic data --------------------------------------------------------
tads <- gi_list_topdom(gi_list_ice,
                       window.size = 10,
                       file_out = FALSE)

## Concatenate results across chromosomes and save results -----------------------------------------

tads %>%
    map(`[[`, "binSignal") %>% 
    bind_rows() %>% 
    write.table(., file = paste0("./results/", sample, "_topdom_binSignal.tsv"),
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
tads %>%
    map(`[[`, "domain") %>% 
    bind_rows() %>% 
    write.table(., file = paste0("./results/", sample, "_topdom_domain.tsv"),
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
tads %>%
    map(`[[`, "bed") %>% 
    bind_rows() %>% 
    write.table(., file = paste0("./results/", sample, "_topdom_tads.bed"),
                sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
