# INPUT:
# `results/{sample}/epiAneufinder_results/{sample}_Karyotype.tsv`
# `results/{sample}/epiAneufinder_results/{sample}_karyotype.rds`

# OUTPUT:
# `results/{sample}/epiAneufinder_results/{sample}_Karyotype.recurrent_removed.tsv`
# `results/{sample}/epiAneufinder_results/{sample}_karyotype.recurrent_flagged.rds`

library(tidyverse)
library(here)

# pct of samples to be considered recurrent/artifact in the cohort
pct_samples <- 50

# number of samples in the cohort
n_samples <- 12

# "string" karyotype files
karyotypes.files <- list.files(here("results"), pattern = "_Karyotype.tsv", recursive = T) %>%
    set_names(list.files(here("results"), pattern = "_Karyotype.tsv", recursive = T) %>%
                  dirname() %>%
                  dirname())

karyotypes <- karyotypes.files %>%
    map(~ read_tsv(here("results", .x), col_names = c("sample", "karyotype"), col_types = "cc")) %>%
    list_rbind() 

# karyotype tibble by arm
karyotypes.detailed.files <- list.files(here("results"), pattern = "_karyotype.rds", recursive = T) %>%
    set_names(list.files(here("results"), pattern = "_karyotype.rds", recursive = T) %>%
                  dirname() %>%
                  dirname())

karyotypes.detailed <- karyotypes.detailed.files %>% 
    map(~ readRDS(here("results", .x))) %>%
    list_rbind(names_to = "sample") 

# look only at ones that are kept, and calculate the recurrence
recurrent <- karyotypes.detailed %>%
    filter(keep == TRUE) %>%
    count(chr, call) %>% # remove calls at the chromosome and not the arm level
    filter(n > pct_samples/100*n_samples)

# filter for recurrent gains
recurrent.gains <- recurrent %>%
    filter(call == "gain") %>%
    pull(chr) %>%
    str_remove("^chr")

# filter for recurrent losses
recurrent.losses <- recurrent %>%
    filter(call == "loss") %>%
    pull(chr) %>%
    str_remove("^chr")

# remove the recurrent gains from the string
if(length(recurrent.gains > 1)) {
    walk(recurrent.gains, ~ {
        karyotypes$karyotype <<- str_remove(karyotypes$karyotype, str_glue("\\+{.x}[pq]*(, )?"))
    })
}

# remove the recurrent losses from the string
if(length(recurrent.losses > 1)) {
    walk(recurrent.losses, ~ {
        # remove dels
        karyotypes$karyotype <<- str_remove(karyotypes$karyotype, str_glue("del\\({.x}[pq]?\\)(, )?"))
        karyotypes$karyotype <<- str_remove(karyotypes$karyotype, str_glue("-{.x}(, )?"))
    })
}

karyotypes$karyotype <- str_remove(karyotypes$karyotype, ", $") 

karyotypes <- karyotypes %>%
    mutate(karyotype = case_when(karyotype == "" ~ "normal", 
                                 TRUE ~ karyotype))

# update original files with new extension
iwalk(karyotypes.files, ~ {
    karyotypes %>%
        filter(sample == .y) %>%
        write_tsv(here("results", str_replace(.x, ".tsv", ".recurrent_removed.tsv")), col_names = F)
})

# update the karyotype tibble to change recurrents from keep = TRUE to keep = FALSE
recurrent.list <- select(recurrent, -n) %>% 
    tibble::deframe() %>% 
    as.list()

iwalk(recurrent.list, ~ {
    karyotypes.detailed <<- karyotypes.detailed %>%
        mutate(keep = case_when(chr == .y & call == .x ~ FALSE, 
                                TRUE ~ keep))
})

# update original files with new extension
iwalk(karyotypes.detailed.files, ~ {
    karyotypes.detailed %>%
        filter(sample == .y) %>%
        saveRDS(here("results", str_replace(.x, ".rds", ".recurrent_flagged.rds")))
})