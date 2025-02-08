# INPUT: 
# `data/cohort3/*.seurat.rds`

# OUTPUT: 
# `data/cohort3/*.seurat.processed.rds`

set.seed(1234)
library(Seurat)
library(Signac)
library(scDblFinder)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

signac <- readRDS(args[1])

# General QC based on the violin plots
upper.qc <- 0.95
lower.qc <- 0.05

signac <- subset(
    x = signac,
    subset = nFeature_ATAC < quantile(signac$nFeature_ATAC, upper.qc, na.rm = TRUE) &
        nFeature_ATAC > quantile(signac$nFeature_ATAC, lower.qc, na.rm = TRUE) &
        mtATAC_depth < quantile(signac$mtATAC_depth, upper.qc, na.rm = TRUE) &
        mtATAC_depth > quantile(signac$mtATAC_depth, lower.qc, na.rm = TRUE) &
        peak_region_fragments < quantile(signac$atac_peak_region_fragments, upper.qc, na.rm = TRUE) &
        peak_region_fragments > quantile(signac$atac_peak_region_fragments, lower.qc, na.rm = TRUE) &
        nucleosome_signal < quantile(signac$nucleosome_signal, upper.qc, na.rm = TRUE) &
        TSS.enrichment > quantile(signac$TSS.enrichment, lower.qc, na.rm = TRUE)
)

sce_atac <- scDblFinder(signac[["ATAC"]]@counts, aggregateFeatures=TRUE, BPPARAM=BiocParallel::SerialParam(RNGseed = 1234)) 
doublet_status_atac <- SingleCellExperiment::colData(sce_atac)$scDblFinder.class
names(doublet_status_atac) <- rownames(SingleCellExperiment::colData(sce_atac))
signac <- AddMetaData(signac, doublet_status_atac, "atac.scDblFinder.class")

# Keep only those found to be singlets by both assays
signac <- subset(signac, subset = atac.scDblFinder.class == "singlet")

# Process the ATAC
DefaultAssay(signac) <- "ATAC"
signac <- FindTopFeatures(signac, min.cutoff = 'q25')
signac <- RunTFIDF(signac)
signac <- RunSVD(signac, reduction.name = "lsi", reduction.key = "LSI_")

peaks <- CallPeaks(signac, 
                   macs2.path = "macs2") %>%
    keepStandardChromosomes(pruning.mode = "coarse") %>%
    subsetByOverlaps(x = ., ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts <- FeatureMatrix(
    fragments = Fragments(signac),
    features = peaks.final,
    cells = colnames(signac), 
    verbose = verbosity
)

signac[["peaks"]] <- CreateChromatinAssay(
    counts = macs2_counts,
    fragments = Fragments(signac),
    annotation = Annotation(signac),
    genome = "GRCh38.p13"
)

saveRDS(signac, str_replace(args[1], ".rds", ".processed.rds"))
