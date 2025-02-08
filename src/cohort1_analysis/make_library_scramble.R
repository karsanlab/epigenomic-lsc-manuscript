# INPUT:
# `cache/hg19_novoalign/scramble/*_sc.rds`

# OUTPUT
# `cache/hg19_novoalign/scramble/scrambled_library_*.rds`

library(fs)
library(GenomicRanges)
library(tidyverse)

# PARSE ARGUMENTS ----

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop(paste0("Path to output directory (1) and ",
                "at least one single-cell methylation library .rds (2..n) must be provided."), call. = FALSE)
} else {
    OUTDIR <- args[1]
    SC_RDS <- args[2:length(args)]
}

# MAIN ----

## GENERATE RANDOM INDICES TO SEPARATE INTO SCRAMBLED POOLS ----

N_LIBRARIES = length(SC_RDS)

# helper function to generate subsetting indices for n different pools of 96 single-cells
generateScrambledIndices <- function(n) {
    random_sequence <- sample(rep(c(1:n), 75), n * 75, replace = F)
    indices_for_pool <- map(seq(1:n) %>% set_names(., .),
                            function(i) {
                                which(random_sequence %in% c(i))
                            })

    return(indices_for_pool)
}

scrambled_indices <- generateScrambledIndices(n = N_LIBRARIES)

# separate indices into units of n * 96 for easy subsetting of individual input libraries
scrambled_indices_subset_ready <-
    map(seq(1:N_LIBRARIES),
        function(n) {
            map(scrambled_indices,
                function(indices_for_pool) {
                    keep(indices_for_pool, ~ (. > (n - 1) * 75) & (. <= n * 75))
                }) %>%
                map(~ . - (n - 1) * 75)  # reset each index to be from 1-75 since the lowest library has 75 cells?
        })

## SCRAMBLE LIBRARIES ----

pwalk(list(SC_RDS, scrambled_indices_subset_ready, seq(1:N_LIBRARIES)),
      function(rds_path, scrambled_indices_nth_library, nth_library) {
          print(paste0("> Scrambling library ", path_file(rds_path)))

          # import actual library at path
          library.grs <- read_rds(rds_path)

          iwalk(scrambled_indices_nth_library,
               function(scrambled_indices, nth_scramble) {
                   print(paste0(">> into new library n = ", nth_scramble))

                   library_subsetted.grs <-
                       if (nth_library > 1) {
                           read_rds(paste0(OUTDIR, "/", "scrambled_library_", nth_scramble, ".rds"))
                       } else {
                           c()
                       }

                   library_subsetted.grs <- c(library_subsetted.grs,
                                              library.grs[scrambled_indices])

                   write_rds(library_subsetted.grs,
                             paste0(OUTDIR, "/", "scrambled_library_", nth_scramble, ".rds"))

               })

      })
