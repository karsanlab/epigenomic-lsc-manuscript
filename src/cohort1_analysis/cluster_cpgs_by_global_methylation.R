# INPUT: 
# `data/scWGBS/*.cons.Cmethyl.CpG.txt`
# `results/combinedSampleAnnot.tsv`

# OUTPUT:
# `results/patient_449_keq2_clustering.tsv`

# Build metadata matrix for all samples

library(pcaMethods)
library(umap)
library(here)
library(tidyverse)

aza_metadata.tib <- read_tsv(here("data", "sample_metadata.tsv"))
farlik_metadata.tib <- read_tsv(here("data", "farlik_metadata.tsv"))

combined_annotation_metadata.tib <-
  bind_rows(aza_metadata.tib,
            farlik_metadata.tib)

# List methylation files
cmethyl_files.list <-
    list.files("data", "scWGBS",
               "*.cons.Cmethyl.CpG.txt", full.names = TRUE) %>%
    set_names(str_remove(., "data") %>%
                  str_remove(".cons.*") %>%
                  str_remove("^[\\d]+_"))

# Initialize matrix
common_cpg_methylation_aza.tib <-
    cmethyl_files.list[1] %>%
    read_tsv(col_types = "cdd---",
             col_names = c("seqnames", "start", "end"))

# Import all methylation data and population matrix
iwalk(cmethyl_files.list,
      function(filepath, plate_barcode) {
        methyl_data.tib <-
            read_tsv(filepath,
                     col_types = "cddd",
                     col_names = c("seqnames", "start", "end", plate_barcode))
        
        common_cpg_methylation_aza.tib <<-
            full_join(common_cpg_methylation_aza.tib,
                      methyl_data.tib,
                      by = c("seqnames", "end"))
    })

common_cpg_methylation.mat <- common_cpg_methylation_aza.tib %>%
    data.matrix() %>%
    t()

na_counts_by_col <- colSums(is.na(common_cpg_methylation.mat))

# Filter CpGs for completness across single cells
MIN_PCT_COMPLETENESS <- 40

nonzero_aza_cpgs.tib <-
    na_counts_by_col %>%
    as_tibble(rownames = "position") %>%
    mutate(nonzero_pct = (TOTAL_CELLS - value) / TOTAL_CELLS * 100)

cpgs_gt_min_complete.list <-
    nonzero_aza_cpgs.tib %>%
    filter(nonzero_pct > MIN_PCT_COMPLETENESS) %>%
    pull(position)

cpg_methylation_gt_min_complete.mat <-
  common_cpg_methylation.mat[, cpgs_gt_min_complete.list] %>%
  .[rowSums(is.na(.)) != ncol(.), ]

# Annotate cells with predicted cell type classificaiton based on methylation
aza_predicted_celltypes.tib <-
    read_tsv(here("results", "combinedSampleAnnot.tsv"))

# Perform PCA on methylation matrix

cpg_methylation_gt_min_complete.pca <-
  cpg_methylation_gt_min_complete.mat %>%
  pca(method = "ppca", nPcs = 30)

# Perform UMAP

set.seed(1234)

cpg_methylation_gt_min_complete.umap <-
  umap(cpg_methylation_gt_min_complete.pca@scores[, 1:4])

cpg_methylation_gt_min_complete_umap_annotated.tib <-
  cpg_methylation_gt_min_complete.umap$layout %>%
  as_tibble(rownames = "uid") %>%
  dplyr::select(everything(), "UMAP_1" = V1, "UMAP_2" = V2) %>%
  left_join(combined_annotation_metadata.tib %>%
              dplyr::select(uid, response, patient, origin, celltype))

# Separate biphenotypic patient cells into clusters by response using k-means

set.seed(123)

aza_449_global_methylation.kmeans <-
  cpg_methylation_gt_min_complete_umap_annotated.tib %>%
  filter(patient == "449") %>%
  dplyr::select(uid, UMAP_1, UMAP_2) %>%
  column_to_rownames("uid") %>%
  kmeans(centers = 2, nstart = 20)

aza_449_barcodes.tib <-
  tibble(uid = names(aza_449_global_methylation.kmeans$cluster),
         cluster = as.character(aza_449_global_methylation.kmeans$cluster)) %>%
  mutate(response = case_when(cluster == 2 ~ "responder",
                              cluster == 1 ~ "nonresponder"))

write_tsv(aza_449_barcodes.tib,
          here("results",
               "patient_449_global_methylation_keq2_clustering.tsv"))

