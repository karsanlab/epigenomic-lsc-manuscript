# INPUT: 
# `cache/hg19_novoalign_global_methylation_qc.rds`
# `cache/hg19_novoalign/*_annotated_CpGs.tbl.rds`
# `data/full_metadata.rds`
# `results/patient_449_keq2_clustering.tsv`

# OUTPUT:
# `results/combined_methylation_qc.rds`

library(nlme)
library(tidyverse)
library(here)

# takes value and # of digits; outputs formatted p-value
format_pval <- function(pval, digits) {
    c = paste0("p = ", format(round(pval, digits), scientific = F))
    return(c)
}

# calculate ANOVA from feature and grouping
calculateOneAnovaPvalV1 <- function(this.feature, grouping, data) {
    f = as.formula(paste0("percent_meth ~ ", grouping))
    
    model = lme(f, random = ~ 1 | patient,
                data = data %>%
                    filter(annotation == this.feature & !is.na(percent_meth)),
                method = "REML")
    
    anova <- anova.lme(model, type = "sequential", adjustSigma = F)
    
    return(anova$`p-value`[2])
}

calculateOneAnovaPvalV2 <- function(this.feature, grouping, data) {
    f = as.formula(paste0("methylation ~ ", grouping))
    
    model = lme(f, random = ~ 1 | patient,
                data = data %>%
                    filter(feature == this.feature & !is.na(methylation)),
                method = "REML")
    
    anova <- anova.lme(model, type = "sequential", adjustSigma = F)
    
    return(anova$`p-value`[2])
}

# read in the cached object from getGlobalMethylation.R
global_methylation.tib <- readRDS(here("cache", "hg19_novoalign_global_methylation_qc.rds"))

# read in the cached objects produced by getMethylationByGenomicRegion.R
annotated.tibs <- list.files(here("cache", "hg19_novoalign"), pattern = "*.tbl.rds") %>%
    str_subset("hg19") %>%
    str_split("_", simplify = TRUE) %>%
    .[,3] %>%
    set_names(list.files(here("cache", "hg19_novoalign"), pattern = "*.tbl.rds", full.names = T) %>%
                  str_subset("hg19"), .) %>%
    map(~ {
        readRDS(.x)
    })

metadata <- readRDS(here("data", "full_metadata.rds"))
passed.qc <- readRDS(here("cache", "passed_qc.rds"))


# combine them together
plot.data <- annotated.tibs %>%
    list_rbind(names_to = "Plate_ID") %>%
    left_join(metadata, by = c("Plate_ID", "cell")) %>%
    bind_rows(
        global_methylation.tib %>%
            mutate(annotation = "global")
    ) %>%
    # need to do this again because get_regional_methylation.R is run for all 96 cells in each sample, not only what passes qc
    semi_join(passed.qc, by = c("Plate_ID", "cell")) %>%
    mutate(annotation = case_match(annotation,
                                   "islands" ~ "CpG islands", 
                                   "shores" ~ "CpG shores",
                                   "shelves" ~ "CpG shelves",
                                   .default = annotation), 
           annotation = factor(annotation, levels = c("global", "genes", "promoters", "enhancers", "CpG islands", "CpG shores", "CpG shelves")))


cluster.449 <- read_tsv("results", "patient_449_global_methylation_keq2_clustering.tsv") %>%
    separate(uid, into = c("Plate_ID", "cell"), sep = "_") %>%
    mutate(cell_response = case_match(
        response,
        "Responder" ~ "Responder", 
        "Nonresponder" ~ "Non-responder"
    ), 
    cell_response = factor(cell_response, levels = c("Responder", "Non-responder")), 
    Aza_ID = "449") %>%
    select(-cluster, -response)

# split 449 into the two sub populations, using cell_response column 
# for those that were dropped from the global UMAP, assign based on distribution
plot.data.split <- bind_rows(
    global_methylation.tib %>%
        filter(patient == "449") %>%
        left_join(cluster.449, by = c("Plate_ID", "cell", "patient")) %>%
        mutate(cell_response = case_when(
                is.na(cell_response) & percent_meth >= 72 ~ "Non-responder", 
                is.na(cell_response) & percent_meth < 72 ~ "Responder", 
                TRUE ~ cell_response
            ),
               cell_response = factor(cell_response, levels = c("Responder", "Non-responder"))) %>%
        select(-percent_meth),
    
    # Grab the other patients
    global_methylation.tib %>%
        filter(patient != "449") %>%
        mutate(cell_response = Class) %>%
        select(-percent_meth)
) %>%
    # we now have split global methylation, re-join 
    left_join(plot.data)  %>%
    mutate(annotation = case_match(as.character(annotation),
                                   "islands" ~ "CpG islands", 
                                   "shores" ~ "CpG shores",
                                   "shelves" ~ "CpG shelves",
                                   .default = as.character(annotation)), 
           annotation = factor(annotation, levels = c("global", "genes", "promoters", "enhancers", "CpG islands", "CpG shores", "CpG shelves"))) 

anova.split <- select(plot.data.split, annotation) %>% 
    unique() %>%
    mutate(overall_anova_pval = map_dbl(annotation, ~ calculateOneAnovaPvalV1(., "cell_response", plot.data.split)),
           overall_anova_pval_adj = p.adjust(p = overall_anova_pval, method = "BH"),
           overall_label = case_when(overall_anova_pval_adj < 0.0001 ~ "* p < 0.0001",
                                     overall_anova_pval_adj < 0.05 ~ 
                                         paste0("* ",
                                                format_pval(overall_anova_pval_adj, 4)),
                                     TRUE ~ format_pval(overall_anova_pval_adj, 4)),
           Class = "TBD")

scramble_avg_meth_by_region_w_metadata.tib <- list.files(here::here("cache", "hg19_novoalign", "scramble"), "_by_region.rds", full.names = T) %>%
    set_names(str_remove(basename(.), "pt_") %>%
                  str_remove("_methylation.*")) %>%
    map(~ readRDS(.)) %>%
    list_rbind(names_to = "patient") %>%
    dplyr::rename(methylation = pct_methyl, 
                  feature = region) %>%
    dplyr::mutate(overall_response = case_when(str_detect(patient, "NR") ~ "Non-responder", 
                                               TRUE ~ "Responder"))

anova.scramble <- select(scramble_avg_meth_by_region_w_metadata.tib, feature) %>% 
    unique() %>%
    mutate(overall_anova_pval = map_dbl(feature, ~ calculateOneAnovaPvalV2(., "overall_response", scramble_avg_meth_by_region_w_metadata.tib)),
           overall_anova_pval_adj = p.adjust(p = overall_anova_pval, method = "BH"),
           overall_label = case_when(overall_anova_pval_adj < 0.0001 ~ "* p < 0.0001",
                                     overall_anova_pval_adj < 0.05 ~ 
                                         paste0("* ",
                                                format_pval(overall_anova_pval_adj, 4)),
                                     TRUE ~ format_pval(overall_anova_pval_adj, 4)),
           overall_response = "TBD")

# export the data with the original metadata
plot.data.split %>%
    pivot_wider(names_from = "annotation", values_from = "percent_meth") %>%
    dplyr::relocate(all_of(c("global", "genes", "promoters", "enhancers", "CpG islands", "CpG shores", "CpG shelves")), .after = "cell") %>%
    dplyr::relocate(Class, cell_response, .before = "cell") %>%
    saveRDS(here("results", "combined_methylation_qc_metadata.rds"))
