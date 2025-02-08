# INPUT: 
# `results/Scrambled_Cohort2_FindMarkers_DARs_v{1..10}.rds`

# OUTPUT: 
# `results/Scrambled_Cohort_FindMarkers_DARs_genomic_annotated_v{1..10}.rds`

library(annotatr)
library(here)
library(rstatix)
library(broom)
library(plyranges)
library(tidyverse)

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

annotate_with_genome <- function(input.gr) {
    
    annotated.genomic.gr <- annotate_regions(
        regions = input.gr,
        annotations = genomic_annotation,
        ignore.strand = T,
        quiet = F
    )
    
    annotated.genomic.tib <- annotated.genomic.gr %>%
        data.frame() %>%
        dplyr::mutate(
            annotation = case_match(annot.type,
                                    "hg38_cpg_islands" ~ "CpG islands",
                                    "hg38_genes_3UTRs" ~ "3' UTRs",
                                    "hg38_genes_5UTRs" ~ "5' UTRs",
                                    "hg38_cpg_shores" ~ "CpG shores",
                                    "hg38_cpg_shelves" ~ "CpG shelves",
                                    "hg38_genes_1to5kb" ~ "1 to 5 kb",
                                    "hg38_genes_intergenic" ~ "intergenic",
                                    "hg38_enhancers" ~ "enhancers",
                                    "hg38_genes_exons" ~ "exons",
                                    "hg38_genes_introns" ~ "introns",
                                    "hg38_cpg_inter" ~ "inter CpG islands",
                                    "hg38_genes_promoters" ~ "promoters"
            )
        )
    
    annotated.genomic.tib <- input.gr %>%
        data.frame() %>%
        mutate(annotation = "global") %>%
        bind_rows(annotated.genomic.tib) %>%
        mutate(annotation = factor(annotation, levels = c("global", "1 to 5 kb", "3' UTRs", "5' UTRs", 
                                                          "promoters", "enhancers", "exons", "introns", 
                                                          "intergenic", "CpG shores", "CpG shelves", 
                                                          "CpG islands", "inter CpG islands"))
        ) %>%
        unite("region", seqnames, start, end, sep = "-")
}

# load data ----
args <- commandArgs(trailingOnly = T)

input <- readRDS(args[1])
iter <- args[2]

input.gr <- input %>%
    rownames_to_column("region") %>%
    separate(region, into = c("seqnames", "start", "end"), sep = "-") %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
# enhancers ----
enhancers.hg38 <- readRDS(here("resources", "hsapiens_pb_mononuclear_enhancers.rds"))
enhancer_annotation.hg38 <- enhancers.hg38 %>%
    plyranges::mutate(id = 1:length(enhancers.hg38),
                      gene_id = "enhancers",
                      symbol = "enhancers",
                      type = "hg38_enhancers")

annotatr_cache$set("hg38_enhancers", enhancer_annotation.hg38)

genomic_annotation <- build_annotations(genome = "hg38", annotations = c("hg38_basicgenes", "hg38_genes_intergenic", "hg38_enhancers", "hg38_cpgs"))

# analysis ---- 
annotated.genomic.tib <- annotate_with_genome(input.gr)

annotated.genomic.tib <- annotated.genomic.tib %>%
    mutate(annotation = fct_recode(annotation, genes = "1 to 5 kb", genes = "3' UTRs", genes = "5' UTRs", genes = "exons", genes = "introns")) %>%
    filter(str_starts(annotation, "inter", negate = TRUE)) %>%
    mutate(annotation = fct_drop(annotation),
           annotation = fct_relevel(annotation, "enhancers", after = 3),
           annotation = fct_relevel(annotation, "CpG islands", after = 4)) %>%
    select(dar_id, avg_log2FC, p_val_adj, annotation) %>%
    distinct()

genomic.stats <- t_test_log2fc(annotated.genomic.tib,
                               features = levels(annotated.genomic.tib$annotation),
                               feature.col = "annotation")

saveRDS(annotated.genomic.tib, here("results", str_glue("Scrambled_Cohort2_FindMarkers_DARs_genomic_annotated_v{iter}.rds")))
