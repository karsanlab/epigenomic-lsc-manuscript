# INPUT: 
# `results/cohort2/*.seurat.rds`

# OUTPUT: 
# `results/cohort2/*.seurat.processed.rds`

set.seed(1234)
library(Seurat)
library(Signac)
library(scDblFinder)
library(SoupX)
library(Azimuth)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

DOGMA <- readRDS(args[1])

# General QC based on the violin plots
upper.qc <- 0.95
lower.qc <- 0.05

DOGMA <- subset(
    x = DOGMA,
    subset = nFeature_ATAC < quantile(DOGMA$nFeature_ATAC, upper.qc, na.rm = TRUE) &
        nFeature_ATAC > quantile(DOGMA$nFeature_ATAC, lower.qc, na.rm = TRUE) &
        nFeature_RNA < quantile(DOGMA$nFeature_RNA, upper.qc, na.rm = TRUE) &
        nFeature_RNA > quantile(DOGMA$nFeature_RNA, lower.qc, na.rm = TRUE) &
        mtATAC_depth < quantile(DOGMA$mtATAC_depth, upper.qc, na.rm = TRUE) &
        mtATAC_depth > quantile(DOGMA$mtATAC_depth, lower.qc, na.rm = TRUE) &
        mtRNA_depth < quantile(DOGMA$mtRNA_depth, upper.qc, na.rm = TRUE) &
        mtRNA_depth > quantile(DOGMA$mtRNA_depth, lower.qc, na.rm = TRUE) &
        atac_peak_region_fragments < quantile(DOGMA$atac_peak_region_fragments, upper.qc, na.rm = TRUE) &
        atac_peak_region_fragments > quantile(DOGMA$atac_peak_region_fragments, lower.qc, na.rm = TRUE) &
        nucleosome_signal < quantile(DOGMA$nucleosome_signal, upper.qc, na.rm = TRUE) &
        TSS.enrichment > quantile(DOGMA$TSS.enrichment, lower.qc, na.rm = TRUE)
)

# Remove doublets using scDblFinder
sce_rna <- scDblFinder(DOGMA[["RNA"]]@counts)
doublet_status_rna <- SingleCellExperiment::colData(sce_rna)$scDblFinder.class
names(doublet_status_rna) <- rownames(SingleCellExperiment::colData(sce_rna))
DOGMA <- AddMetaData(DOGMA, doublet_status_rna, "rna.scDblFinder.class")

sce_atac <- scDblFinder(DOGMA[["ATAC"]]@counts, aggregateFeatures=TRUE, BPPARAM=BiocParallel::SerialParam(RNGseed = 1234)) 
doublet_status_atac <- SingleCellExperiment::colData(sce_atac)$scDblFinder.class
names(doublet_status_atac) <- rownames(SingleCellExperiment::colData(sce_atac))
DOGMA <- AddMetaData(DOGMA, doublet_status_atac, "atac.scDblFinder.class")

# Keep only those found to be singlets by both assays
DOGMA <- subset(DOGMA, subset = atac.scDblFinder.class == "singlet" & rna.scDblFinder.class == "singlet")

# Remove RNA contamination using SoupX
sc <- load10X(cellranger.out)
sc1 <- try(autoEstCont(sc))

if (class(sc1)[1] == "try-error") {
    cat("Error estimating contamination fraction.", sep = "\n")
    sc1 <- setContaminationFraction(sc, 0.05)
    cat("Contamination set to 5%.", sep = "\n")
    Misc(object = DOGMA, slot = 'RNA_contam') <- 0.05
} else if (sc1$fit$rhoEst > 0.2) {
    cat(str_glue("Contamination: {sc1$fit$rhoEst*100}%"), sep = "\n")
    sc1 <- setContaminationFraction(sc, 0.05)
    cat("Contamination too high: manually set to 5%.", sep = "\n")
    Misc(object = DOGMA, slot = 'RNA_contam') <- 0.05
} else {
    cat(str_glue("Contamination: {sc1$fit$rhoEst*100}%"), sep = "\n")
    Misc(object = DOGMA, slot = 'RNA_contam') <- sc1$fit$rhoEst
}

out <- adjustCounts(sc1)

srat <-  CreateSeuratObject(out)
srat <- subset(srat, cells = Cells(DOGMA))
DOGMA[["RNA"]] <- srat[["RNA"]]
rm(srat)

# Normalize RNA 
DefaultAssay(DOGMA) <- "RNA"
DOGMA <- NormalizeData(DOGMA, normalization.method = "CLR", margin = 2) %>%
    ScaleData() %>%
    FindVariableFeatures()

DOGMA <- CellCycleScoring(DOGMA, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = FALSE)
# returns these columns in metadata S.Score, G2M.Score, and Phase
colnames(DOGMA@meta.data)[str_which(colnames(DOGMA@meta.data), '^S.Score$')] <- "S.Score.RNA"
colnames(DOGMA@meta.data)[str_which(colnames(DOGMA@meta.data), '^G2M.Score$')] <- "G2M.Score.RNA"
colnames(DOGMA@meta.data)[str_which(colnames(DOGMA@meta.data), '^Phase$')] <- "Phase.RNA"

celltype.order <- c(
    "HSC/MPP",
    "GMP",
    "MEP",
    "Prog Mk",
    "Early Eryth",
    "Late Eryth",
    "BaEoMa",
    "LMPP",
    "CLP",
    "pro B",
    "pre B",
    "transitional B",
    "Naive B",
    "Memory B",
    "pre-pDC",
    "pre-mDC",
    "pDC",
    "ASDC",
    "cDC1",
    "cDC2",
    "CD14 Mono",
    "CD16 Mono",
    "Macrophage",
    "Stromal",
    "Plasma",
    "Platelet",
    "CD4 Naive",
    "CD4 Memory",
    "CD4 Effector",
    "CD8 Naive",
    "CD8 Memory",
    "CD8 Effector_1",
    "CD8 Effector_2",
    "CD8 Effector_3",
    "MAIT",
    "T Proliferating",
    "NK Proliferating",
    "NK",
    "NK CD56+",
    "ILC"
)

DOGMA <- RunAzimuth(
    DOGMA,
    reference = "bonemarrowref",
    assay = "RNA",
    umap.name = "umap.bonemarrowref.rna"
)

# Rename the columns
colnames(DOGMA@meta.data)[str_which(colnames(DOGMA@meta.data), '^predicted.celltype.l2.score$')] <- "predicted.rna.bonemarrowref.celltype.l2.score"
colnames(DOGMA@meta.data)[str_which(colnames(DOGMA@meta.data), '^predicted.celltype.l2$')] <- "predicted.rna.bonemarrowref.celltype.l2"
colnames(DOGMA@meta.data)[str_which(colnames(DOGMA@meta.data), '^predicted.celltype.l1.score$')] <- "predicted.rna.bonemarrowref.celltype.l1.score"
colnames(DOGMA@meta.data)[str_which(colnames(DOGMA@meta.data), '^predicted.celltype.l1$')] <- "predicted.rna.bonemarrowref.celltype.l1"
Key(DOGMA[["umap.bonemarrowref.rna"]]) <- "umap.bonemarrowref.rna_"
DOGMA$predicted.rna.bonemarrowref.celltype.l2 <- str_replace(DOGMA$predicted.rna.bonemarrowref.celltype.l2, "^HSC$", "HSC/MPP") %>%
    str_replace("^EMP$", "MEP")

DOGMA$predicted.rna.bonemarrowref.celltype.l2 <- factor(DOGMA$predicted.rna.bonemarrowref.celltype.l2,
                                                        levels = celltype.order)

DOGMA$Azimuth.RNA <- DOGMA$predicted.rna.bonemarrowref.celltype.l2

# Process the ATAC
DefaultAssay(DOGMA) <- "ATAC"
DOGMA <- FindTopFeatures(DOGMA, min.cutoff = 'q25')
DOGMA <- RunTFIDF(DOGMA)
DOGMA <- RunSVD(DOGMA, reduction.name = "lsi", reduction.key = "LSI_")

peaks <- CallPeaks(DOGMA, 
                   macs2.path = "macs2") %>%
    keepStandardChromosomes(pruning.mode = "coarse") %>%
    subsetByOverlaps(x = ., ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts <- FeatureMatrix(
    fragments = Fragments(DOGMA),
    features = peaks.final,
    cells = colnames(DOGMA), 
    verbose = verbosity
)

DOGMA[["peaks"]] <- CreateChromatinAssay(
    counts = macs2_counts,
    fragments = Fragments(DOGMA),
    annotation = Annotation(DOGMA),
    genome = "GRCh38.p13"
)

saveRDS(DOGMA, str_replace(args[1], ".rds", ".processed.rds"))
