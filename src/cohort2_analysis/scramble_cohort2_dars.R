# INPUT:
# `results/scrambled_integrated_cohort2_samples.rds`

# OUTPUT: 
# `results/Scrambled_Cohort2_FindMarkers_DARs.rds`

library(Seurat)
library(Signac)
library(tidyverse)
library(here)

scramble_patients <- function(obj, into = 12) {
    df <- FetchData(obj, "sample") %>%
        as.data.frame() %>%
        rownames_to_column("Barcode")
    
    df$scrambled_sample <- sample(factor(rep(1:into, length.out=nrow(df)), 
                                         labels=paste0("sample", 1:into))) 
    
    df <- df %>%
        column_to_rownames("Barcode") %>%
        select(scrambled_sample)
    
    obj <- AddMetaData(obj, metadata = df, col.name = "scrambled_sample")
    
    return(obj)
}

scramble_patient_response <- function(obj, sample = 6, iter = 1:10, overwrite = F) {
    
    min <- min(iter)
    max <- max(iter)
    
    i <- min 
    while(TRUE) {
        message(str_glue("Scrambling iteration #{i}..."))
        while(TRUE) {
            # split into 2 groups representing responder vs nonresponder
            group1 <- sample(unique(obj$scrambled_sample), sample, replace = FALSE)
            group2 <- setdiff(unique(obj$scrambled_sample), group1)
            
            even.scramble.metadata <- FetchData(obj, c("scrambled_sample", "response")) %>%
                mutate(scrambled_response = case_when(scrambled_sample %in% group1 ~ "group1", 
                                                      scrambled_sample %in% group2 ~ "group2")) %>%
                select(scrambled_response)
            
            counts <- even.scramble.metadata %>% 
                count(scrambled_response) %>%
                tibble::deframe() %>%
                as.list()
            
            # do not overwrite existing scrambles
            while(overwrite == FALSE & str_glue("even_patient_scrambled_response_{i}") %in% colnames(obj@meta.data)) {
                message(str_glue("Scramble iteration #{i} already exists. Incremented to #{i+1}."))
                i <- i + 1
                max <- max + 1
                message(str_glue("Scrambling iteration #{i}..."))
                
            }
            obj <- AddMetaData(obj, metadata = even.scramble.metadata, col.name = str_glue("even_patient_scrambled_response_{i}"))
            
            # check that the same scramble does not exist already
            if(i > 1) {
                bools <- map_lgl(seq(i-1), ~ {
                    message(str_glue("Checking scramble iteration #{.x} and #{i}..."))
                    !identical(obj[[str_glue("even_patient_scrambled_response_{.x}")]], 
                               obj[[str_glue("even_patient_scrambled_response_{i}")]], 
                               attrib.as.set = F)
                })
                if (all(bools)) {
                    message(str_glue("Scramble iteration #{i} is unique!"))
                    i <- i + 1
                    break
                } else {
                    identical_iter <- which(bool) + (min-1)
                    message(str_glue("Scramble iteration #{i} is the same as scramble iteration #{identical_iter}."))
                    obj[[str_glue("even_patient_scrambled_response_{i}")]] <- NULL
                    # do not increment if identical, scramble it again.
                }
            } else {
                i <- i + 1
                break
            }
        }
        
        if(i > max) {
            break
        }
    }
    return(obj)
}

obj <- readRDS(here("results", "scrambled_integrated_cohort2_samples.rds"))

obj <- scramble_patients(obj)
obj <- scramble_patient_response(obj)

saveRDS(obj, here("results", "Scrambled_Cohort2_FindMarkers_DARs.rds"))