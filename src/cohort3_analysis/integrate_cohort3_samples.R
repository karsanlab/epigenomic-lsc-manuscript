# INPUT: 
# `data/cohort3/*.seurat.processed.rds`

# OUTPUT: 
# `results/integrated_cohort3_samples.rds`

library(Seurat)
library(Signac)
library(harmony)
library(pbmcapply)
library(here)
library(tidyverse)

rds <- list.files(here("results", "cohort3"), pattern = ".labelled.rds", full.names = T)
batch_effects <- c("sample")
lambda <- c(1)

# load every object into a list in parallel
tryCatch({
    message(glue("[ {now()} ] Loading data..."))
    if(length(rds) > 1) {
        DOGMA <- pbmclapply(rds, readRDS, mc.cores = length(rds), ignore.interactive = TRUE)
    } else {
        DOGMA <- map(rds, ~ readRDS(.x))
    }
}, error = function(e) {
    stop("Error loading data.")
})

# function to change the default assay of every seurat object in the list
change_defaults <- function(seurat_list, new_default) {
    seurat_new_defaults <-
        map(seurat_list,
            function(signac) {
                DefaultAssay(seurat) <- new_default
                return(seurat)
            })
    
    return(seurat_new_defaults)
}

DOGMA <- change_defaults(DOGMA, "peaks")

# ATAC analysis -- merge all objects into one by finding common peaks

# find common peaks
common_peaks <- map(DOGMA, Signac::granges) %>%
    purrr::reduce(c) %>%
    Signac::reduce()

add_common_peaks_assay <- function(signac, common_peaks) {
    DefaultAssay(signac) <- "peaks"
    counts <- FeatureMatrix(
        fragments = Fragments(signac),
        features = common_peaks,
        cells = Cells(signac)
    )
    
    assay <- CreateChromatinAssay(
        counts = counts,
        min.features = 0,
        fragments = Fragments(signac),
        annotation = Annotation(signac),
        genome = genome(signac)
    )
    signac <- CreateSeuratObject(counts = assay, assay = "peaks", project = Project(signac), meta.data = signac@meta.data)
    return(signac)
}

# find the common peaks
objects_to_merge <- map(DOGMA, ~ add_common_peaks_assay(., common_peaks))

# rename cells to include sample name
objects_to_merge <- map(objects_to_merge, function(obj) {
    obj <- RenameCells(obj, add.cell.id = Project(obj))
    return(obj)
})

# merge objects
merged_object <- merge(objects_to_merge[[1]],
                       objects_to_merge[2: length(objects_to_merge)])

# normalize ATAC
merged_object <- merged_object %>%
    RunTFIDF() %>%
    FindTopFeatures(min.cutoff = 'q25')

# calculate gene activity
gene_activity <- GeneActivity(merged_object)

merged_object[["GeneActivity"]] <- CreateAssayObject(counts = gene_activity)

# normalize gene activity
merged_object <- NormalizeData(
    object = merged_object,
    assay = 'GeneActivity',
    normalization.method = 'LogNormalize',
    scale.factor = median(merged_object$nCount_GeneActivity)
)

# Run SVD
DefaultAssay(merged_object) <- "peaks"
merged_object <- merged_object %>%
    RunSVD(reduction.name = "lsi.atac")

# Run Harmony 
integrated_object <- RunHarmony(merged_object, batch_effects, reduction.use = "lsi.atac", 
                                assay.use = "peaks", project.dim = FALSE,
                                reduction.save = "harmony_ATAC", lambda = lambda)

# Run UMAP
integrated_object <- integrated_object %>%
    RunUMAP(reduction = "harmony_ATAC", dims = 2:30, reduction.key = "atacUMAP_", reduction.name = "umap.atac", return.model = TRUE) %>%
    RunUMAP(reduction = "harmony_ATAC", dims = 2:30, reduction.key = "atac3dUMAP_", reduction.name = "umap.atac.3d", n.components = 3L, return.model = TRUE)

# Find Neighbours
integrated_object <- FindNeighbors(integrated_object, reduction = "harmony_ATAC", dims = 2:30, verbose = TRUE)

# Find Clusters
integrated_object <- FindClusters(integrated_object, resolution = 1, algorithm = 3, verbose = TRUE)

# Store cluster IDs
integrated_object$idents.harmony.atac <- Idents(integrated_object)

# save combined object
saveRDS(integrated_object, here("results", "integrated_cohort3_samples.rds"))
