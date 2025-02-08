# INPUT: 
# `results/cohort2_adt/*.seurat.processed.rds`

# OUTPUT: 
# `results/integrated_cohort2_ADT_samples.rds`

library(Seurat)
library(Signac)
library(harmony)
library(pbmcapply)
library(here)
library(tidyverse)

rds <- list.files(here("results", "cohort2_adt"), pattern = ".labelled.rds", full.names = T)
batch_effects <- c("sample", "disease_state", "source", "Phase.RNA")
lambda <- c(1,1,1,1)

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

# copy before deletion
combined <- integrated_object

rm(objects_to_merge)
rm(merged_object)
rm(integrated_object)

# RNA analysis
DOGMA <- change_defaults(DOGMA, "RNA")

# rename cell IDs to include sample name
objects_to_merge <- map(DOGMA, function(obj) {
    obj <- RenameCells(obj, add.cell.id = Project(obj))
    return(obj)
})

# Create RNA-only object for each object
objects_to_merge <- map(objects_to_merge, function(obj) {
    CreateSeuratObject(
        counts = GetAssay(obj, assay = "RNA"),
        project = Project(obj),
        meta.data = obj@meta.data
    )
}) 

# Merge RNA objects
merged_object <- merge(objects_to_merge[[1]],
                       objects_to_merge[2: length(objects_to_merge)])

# Process
merged_object <- merged_object %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures()

# Run PCA
merged_object <- merged_object %>%
    RunPCA(assay = "RNA",
           reduction.name = "pca.rna")

# Run Harmony
integrated_object <- RunHarmony(merged_object, batch_effects, reduction.use = "pca.rna", 
                                assay.use = "RNA", project.dim = FALSE, 
                                reduction.save = "harmony_RNA", lambda = lambda)

# Run UMAP
integrated_object <- integrated_object %>%
    RunUMAP(reduction = "harmony_RNA", dims = 1:30, reduction.key = "rnaUMAP_", reduction.name = "umap.rna", return.model = TRUE) %>%
    RunUMAP(reduction = "harmony_RNA", dims = 1:30, reduction.key = "rna3dUMAP_", reduction.name = "umap.rna.3d", n.components = 3L, return.model = TRUE)

# Find neighbours
integrated_object <- FindNeighbors(object = integrated_object, reduction = 'harmony_RNA', dims = 1:30)

# Find clusters
integrated_object <- FindClusters(object = integrated_object, verbose = FALSE, algorithm = 3, resolution = 1)

# Store cluster IDs
integrated_object$idents.harmony.rna <- Idents(integrated_object)

# transfer RNA assays/slots/reductions into the combined object
combined[["RNA"]] <- GetAssay(integrated_object, assay = "RNA")
combined[["pca.rna"]] <- integrated_object[["pca.rna"]]
combined[["harmony_RNA"]] <- integrated_object[["harmony_RNA"]]

combined[["umap.rna"]] <- integrated_object[["umap.rna"]]
combined[["umap.rna.3d"]] <- integrated_object[["umap.rna.3d"]]
combined$idents.harmony.rna <- integrated_object$idents.harmony.rna

rm(objects_to_merge)
rm(merged_object)
rm(integrated_object)

# ADT ----
DOGMA <- change_defaults(DOGMA, "ADT")

objects_to_merge <- map(DOGMA, function(obj) {
    obj <- RenameCells(obj, add.cell.id = Project(obj))
    return(obj)
})

objects_to_merge <- map(objects_to_merge, function(obj) {
    CreateSeuratObject(
        counts = GetAssay(obj, assay = "ADT"),
        project = Project(obj),
        meta.data = obj@meta.data, 
        assay = "ADT"
    )
}) 

merged_object <- merge(objects_to_merge[[1]],
                       objects_to_merge[2: length(objects_to_merge)])

merged_object <- merged_object %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures()

adt_features <- nrow(merged_object[["ADT"]])-1
merged_object <- merged_object %>%
    RunPCA(assay = "ADT", dims = 1:adt_features,
           reduction.name = "pca.adt")

# Integrating with Harmony
integrated_object <- RunHarmony(merged_object, batch_effects, reduction.use = "pca.adt", 
                                assay.use = "ADT", project.dim = FALSE, 
                                reduction.save = "harmony_ADT", lambda = lambda)

integrated_object <- integrated_object %>%
    RunUMAP(reduction = "harmony_ADT", dims = 1:adt_features, reduction.key = "adtUMAP_", reduction.name = "umap.adt", return.model = TRUE) %>%
    RunUMAP(reduction = "harmony_ADT", dims = 1:adt_features, reduction.key = "adt3dUMAP_", reduction.name = "umap.adt.3d", n.components = 3L, return.model = TRUE)

integrated_object <- FindNeighbors(object = integrated_object, reduction = 'harmony_ADT', dims = 1:adt_features, graph.name = c("adt.nn", "adt.snn"))
integrated_object <- FindClusters(object = integrated_object, verbose = FALSE, algorithm = 3, resolution = 1, graph.name = "adt.snn")
integrated_object$idents.harmony.adt <- Idents(integrated_object)

# transfer ADT assays/slots/reductions into the combined object
combined[["ADT"]] <- GetAssay(integrated_object, assay = "ADT")
combined[["pca.adt"]] <- integrated_object[["pca.adt"]]
combined[["harmony_ADT"]] <- integrated_object[["harmony_ADT"]]

combined[["umap.adt"]] <- integrated_object[["umap.adt"]]
combined[["umap.adt.3d"]] <- integrated_object[["umap.adt.3d"]]
combined$idents.harmony.adt <- integrated_object$idents.harmony.adt

saveRDS(combined, here("results", "integrated_cohort2_ADT_samples.rds"))
