# INPUT:
# `results/integrated_cohort2_samples.rds`
# `results/top_dmtfbs.tsv`

# OUTPUT:
# `cache/atac_matrix.tsv`
# `cache/64tfbs.tsv`
# `results/atac_GRNBoost2.tsv`
# `results/ATAC/GRNBoost2/int/1.2_corrMat.Rds`
# `results/ATAC/GRNBoost2/int/1.4_GENIE3_linkList.Rds`
# `results/ATAC/GRNBoost2/int/2.6_regulons_asIncidMat.Rds`
# `results/ATAC/GRNBoost2/output/graph.rds`
# `results/ATAC/GRNBoost2/output/layout_fr.rds`

library(data.table)
library(Seurat)
library(SCENIC)
library(SCopeLoomR)
library(igraph)
library(here)
library(tidyverse)

# export the gene accessibility scores for GRNBoost2
obj <- readRDS(here("results", "integrated_cohort2_samples.rds"))

counts <- GetAssayData(obj, assay = "GeneActivity", slot = "data")

counts %>% 
    t() %>% 
    write.table(here("cache", "atac_matrix.tsv"), quote = F, sep = "\t")

tfbs <- read_tsv(here("results", "top_dmtfbs.tsv"))

tfbs %>%
    select(tf) %>%
    write_tsv(here("cache", "64tfbs.tsv"), col_names = F)

system("python3 run_GRNboost2.py")

loom <- build_loom(here("cache", "integrated_cohort2.loom"), dgem = counts)
loom <- add_cell_annotation(loom, obj@meta.data)
close_loom(loom)

loom <- open_loom(here("cache", "integrated_cohort2.loom"))
atac.exp <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)
rm(loom)

cellInfo <- data.frame(cellInfo) %>%
    select(response)

dir.create(here("results", "ATAC", "GRNBoost2", "int"), recursive = T, showWarnings = F)
saveRDS(cellInfo, file = here("results", "ATAC", "GRNBoost2", "int", "cellInfo.Rds"))

colVars <- list(response = c("responder" = "#e69f00",
                             "nonresponder" = "#2670b0"))

colVars$response <- colVars$response[intersect(names(colVars$response), cellInfo$response)]
saveRDS(colVars, file = here("results", "ATAC", "GRNBoost2","int", "colVars.Rds"))


myDatasetTitle <- "ATAC SCENIC" 

org <- "hgnc" # or hgnc, or dmel
dbDir <- here("resources") # RcisTarget databases location


dbs <- c(
    "500bp" = "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather", 
    "10kb" = "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
)

setwd(here("results", "ATAC", "GRNBoost2"))
scenicOptions.atac.grn <- initializeScenic(org=org, dbDir=dbDir, dbs = dbs, datasetTitle=myDatasetTitle, nCores=10) 

# Modify if needed
scenicOptions.atac.grn@inputDatasetInfo$cellInfo <- here("results", "ATAC", "GRNBoost2", "int", "cellInfo.Rds")

scenicOptions.atac.grn@inputDatasetInfo$colVars <- here("results", "ATAC", "GRNBoost2","int", "colVars.Rds")

scenicOptions.atac.grn@settings$verbose <- TRUE
scenicOptions.atac.grn@settings$nCores <- 10
scenicOptions.atac.grn@settings$seed <- 123

saveRDS(scenicOptions.atac.grn, file = here("results", "ATAC", "GRNBoost2", "int", "scenicOptions.Rds"))

setwd(here("results", "ATAC", "GRNBoost2"))
runCorrelation(atac.exp, scenicOptions.atac.grn)

# Modify GRNBoost output headers and save to the expected location
# https://github.com/aertslab/SCENIC/issues/40
GRNBoost_output <- read.delim(here("results", "ATAC", "GRNBoost2", "int", "atac_GRNBoost2.tsv"), header=FALSE)
colnames(GRNBoost_output) <- c("TF","Target","weight")
saveRDS(GRNBoost_output, file= here("results", "ATAC", "GRNBoost2", "int", "1.4_GENIE3_linkList.Rds"))

setwd(here("results", "ATAC", "GRNBoost2"))
scenicOptions.atac.grn <- runSCENIC_1_coexNetwork2modules(scenicOptions.atac.grn)
scenicOptions.atac.grn <- runSCENIC_2_createRegulons(scenicOptions.atac.grn)

saveRDS(scenicOptions.atac.grn, here("results", "ATAC", "GRNBoost2", "int", "scenicOptions.completed.Rds"))

edgelist <- readRDS(here("results", "ATAC", "GRNBoost2", "int", "2.6_regulons_asIncidMat.Rds")) %>%
    as.data.frame() %>%
    filter(Freq == 1) %>%
    select(gene = Var1, tf = Var2) %>%
    mutate(tf = str_remove(tf, "_extended$"), gene = str_remove(gene, "_extended$"))

# must save because will change each time you run even if you set seed.
graph <- graph_from_data_frame(d = edgelist, directed = TRUE)
saveRDS(graph, here("results", "ATAC", "GRNBoost2", "output", "graph.rds"))

# must save because will change each time you run even if you set seed.
fr <- layout_with_fr(graph, niter = 40000, start.temp = 20)
saveRDS(fr, here("results", "ATAC", "GRNBoost2", "output", "layout_fr.rds"))

plot(graph,
     layout = fr, 
     vertex.label = ifelse(V(graph)$name %in% tfbs, V(graph)$name, NA),
     vertex.color = "#e69f00",
     vertex.size = 1,
     vertex.frame.color = "black",
     vertex.shape = "circle",
     vertex.label.dist = 0,
     vertex.label.cex = 1,
     vertex.label.family = "Helvetica",
     vertex.label.color = "#2670b0",
     vertex.label.font = 2,
     edge.color = "slategrey", 
     edge.width = 1,
     edge.arrow.size = 0,
)

plot(graph,
     layout = fr, 
     vertex.label = ifelse(V(graph)$name %in% tfbs, V(graph)$name, NA),
     vertex.color = "#e69f00",
     vertex.size = 1,
     vertex.frame.color = "black",
     vertex.shape = "circle",
     vertex.label.dist = 0,
     vertex.label.cex = 1,
     vertex.label.family = "Helvetica",
     vertex.label.color = "#2670b0",
     vertex.label.font = 2,
     edge.color = "slategrey", 
     edge.width = 1,
     edge.arrow.size = 0,
     xlim = c(-0.5, 0.3), 
     ylim = c(-0.6, -0.2)
)