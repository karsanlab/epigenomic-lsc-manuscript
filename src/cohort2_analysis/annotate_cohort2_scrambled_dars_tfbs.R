# INPUT: 
# `results/Scrambled_Cohort2_FindMarkers_DARs_v{1..10}.rds`

# OUTPUT: 
# `results/Scrambled_Cohort2_FindMarkers_DARs_tfbs_annotated_v{1..10}.rds`

library(annotatr)
library(LOLA)
library(here)
library(rstatix)
library(tidyverse)

args <- commandArgs(trailingOnly = T)

# functions  ----
t_test_log2fc <- function(data, features = NULL, feature.col = NULL, p.adj.threshold = NULL, logfc.threshold = NULL) {
    set.seed(1234)
    if(is.null(features)) {
        
        if(!is.null(p.adj.threshold)) {
            data <- data %>%
                filter(p_val_adj < p.adj.threshold)
        }
        
        if(!is.null(logfc.threshold)) {
            data <- data %>%
                filter(abs(avg_log2FC) > logfc.threshold)
        }        
        
        log2fc <- data %>%
            pull(avg_log2FC)
        
        
        norm_dist <- rnorm(length(log2fc), sd = sd(log2fc))
        right_tail_test <- t.test(log2fc, norm_dist, "greater") %>%
            broom::tidy() %>%
            pluck("p.value")
        left_tail_test <- t.test(log2fc, norm_dist, "less") %>%
            broom::tidy() %>%
            pluck("p.value")
        cohens_df <- data.frame(log2fc_val = c(log2fc, norm_dist), 
                                dist_group = c(rep("log2fc", 
                                                   length(log2fc)),
                                               rep("normal", length(norm_dist))) %>%
                                    factor(levels = c("log2fc", 
                                                      "normal")))
        
        cohens_effect_size <- rstatix::cohens_d(cohens_df,
                                                log2fc_val ~ dist_group)
        
        stat_df <- tibble(right_tail_p_val = right_tail_test,
                          left_tail_p_val = left_tail_test,
                          effect_size = cohens_effect_size$effsize,
                          effect_magnitude = cohens_effect_size$magnitude) %>%
            mutate(effect_magnitude = factor(effect_magnitude, 
                                             levels = c("negligible", "small", 
                                                        "moderate", "large")))
        return(stat_df)
        
    } else {
        stat_df <- map(features, function(feat) {
            
            if(!is.null(p.adj.threshold)) {
                data <- data %>%
                    filter(p_val_adj < p.adj.threshold)
            }
            
            if(!is.null(logfc.threshold)) {
                data <- data %>%
                    filter(abs(avg_log2FC) > logfc.threshold)
            }       
            
            log2fc <- data %>%
                dplyr::filter(!!sym(feature.col) == feat) %>%
                pull("avg_log2FC")
            # make variance of both distributions equal to keep t-test assumptions
            norm_dist <- rnorm(length(log2fc), sd = sd(log2fc))
            right_tail_test <- t.test(log2fc, norm_dist, "greater") %>%
                broom::tidy() %>%
                pluck("p.value")
            left_tail_test <- t.test(log2fc, norm_dist, "less") %>%
                broom::tidy() %>%
                pluck("p.value")
            cohens_df <- data.frame(log2fc_val = c(log2fc, norm_dist), 
                                    dist_group = c(rep(str_glue("{feat} log2fc"), 
                                                       length(log2fc)),
                                                   rep("normal", length(norm_dist))) %>%
                                        factor(levels = c(str_glue("{feat} log2fc"), 
                                                          "normal")))
            cohens_effect_size <- rstatix::cohens_d(cohens_df,
                                                    log2fc_val ~ dist_group)
            return(tibble(feature = feat,
                          right_tail_p_val = right_tail_test,
                          left_tail_p_val = left_tail_test,
                          effect_size = cohens_effect_size$effsize,
                          effect_magnitude = cohens_effect_size$magnitude))
        }) %>%
            list_rbind() %>%
            mutate({{feature.col}} := factor(feature, levels = features),
                   effect_magnitude = factor(effect_magnitude, 
                                             levels = c("negligible", "small", 
                                                        "moderate", "large"))) %>%
            dplyr::select(-feature)
        return(stat_df)
    }
}

annotate_with_lola <- function(input.gr) {
    annotated.lola.gr <- annotate_regions(
        regions = input.gr,
        annotations = lola_tfbs.grs,
        ignore.strand = T,
        quiet = F
    ) 
    
    annotated.lola.tib <- annotated.lola.gr %>%
        data.frame() %>%
        unite(col = "region", seqnames, start, end, sep = "-")
}

# controls picked by:
# readRDS(here("results", "tfbs_volcano.rds")) %>%
#     rename(pval = pbkr.pval,
#            delta_meth = delta.meth) %>% arrange(abs(delta_meth), desc(pval)) %>%
#     filter(!tf %in% pull(read_csv(here("results/top_dmtfbs.tsv")), tf))

# load data ----
input <- readRDS(args[1])
iter <- readRDS(args[2])

input.gr <- input %>%
    rownames_to_column("region") %>%
    separate(region, into = c("seqnames", "start", "end"), sep = "-") %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)

# LOLA ----
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

lola_tfbs.grs <- unlist(lola_tfbs.grlist)
names(lola_tfbs.grs) <- NULL


# analysis ----
annotated.lola.tib <- annotate_with_lola(input.gr)

most.sig <- read_csv(here("results", "top_dmtfbs.tsv")) %>%
    mutate(tf = factor(tf))

controls <- readRDS(here("results", "controls.rds")) %>%
    mutate(tf = as.character(tf), 
           tf = factor(tf, levels = .$tf))

annotated.lola.tib_all <- annotated.lola.tib %>%
    filter(annot.filename %in% most.sig$filename | annot.filename %in% controls$filename) %>%
    left_join(
        bind_rows(select(most.sig, filename, tf), 
                  select(controls, filename, tf)
        ), 
        by = c("annot.filename" = "filename")) %>%
    rename(tf_class = tf) %>%
    mutate(tf_class = factor(tf_class, levels = c(levels(most.sig$tf), levels(controls$tf)))) %>%
    select(dar_id, tf_class, avg_log2FC, p_val_adj) %>%
    distinct()

saveRDS(annotated.lola.tib_all, here("results", str_glue("Scrambled_Cohort2_FindMarkers_DARs_tfbs_annotated_v{iter}.rds")))
