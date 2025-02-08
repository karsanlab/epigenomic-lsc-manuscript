# INPUT: 
# `data/fastq_full/*_quant/quant.sf` (from Salmon)

# OUTPUT: 
# `results/EMTAB5456_dge.csv`

library(DESeq2)
library(biomaRt)
library(tximport)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(edgeR)
library(reshape2)
library(sva)
library(here)
library(tidyverse)

# write out SF files
list.files(here("data", "fastq_full"), pattern = "quant.sf", full.names = T, recursive = T) %>%
    set_names(list.files(here("data", "fastq_full"), pattern = "_quant") %>%
                  str_remove("_quant$")) %>%
    iwalk(~ {
        read_tsv(.x) %>%
            mutate(Name == gsub("\\..*", "", Name)) %>%
            write.table(file = here("data", "quants", str_glue("{.y}.sf")), row.names = F, sep = "\t")
    })

metadata <- read_tsv(here("data", "fastq", "E-MTAB-5456_metadata.txt"))

# Load transcript to gene ID reference table
edb <- EnsDb.Hsapiens.v86
tx2gene <- transcripts(edb, return.type = "DataFrame", columns = c("tx_id", "gene_id"))

# read in SF files
quants <- list.files(here("data", "quants"), pattern = ".sf", full.names = T) %>%
    set_names(list.files(here("data", "quants"), pattern = ".sf") %>%
                  str_remove(".sf$")) %>%
    map(~ {
        tximport(.x, type = "salmon", tx2gene = tx2gene)
    })

# Extract counts
EMTAB5456_counts <- imap(quants, ~ {
    counts <- as.data.frame(.x$counts)
    colnames(counts)[1] <- .y
    counts <- rownames_to_column(counts, var = "Name")
    return(counts)
}) %>%
    purrr::reduce(right_join, by = "Name")

# Remove edb and tx2gene
remove(edb)
remove(tx2gene)
rm(quants)

metadata <- as.tibble(metadata)

metadata <- metadata %>%
    dplyr::select("Source Name", "Characteristics[individual]", "Description", "Comment[ENA_RUN]", "Characteristics[phenotype]")

metadata <- metadata %>%
    dplyr::rename("OX_ID" = "Source Name") %>%
    dplyr::rename("Number" = "Characteristics[individual]") %>%
    dplyr::rename("Celltype" = "Description") %>%
    dplyr::rename("ID" = "Comment[ENA_RUN]") %>%
    dplyr::rename("Immunophenotype" = "Characteristics[phenotype]")

# Get HGNC and transcript length from biomaRt
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

hgnc <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), mart = mart, useCache = FALSE)

remove(mart)

# Merge columns and create final hgnc matrix
colnames(EMTAB5456_count)[1] <- "ensembl_gene_id"

EMTAB5456_count <- merge(hgnc, EMTAB5456_count, by="ensembl_gene_id")

EMTAB5456_count[1] <- NULL

# Clean up environment
remove(hgnc)

# Summarize to unique genes for rownames
EMTAB5456_count <- EMTAB5456_count %>%
    group_by(hgnc_symbol) %>%
    summarise_all(funs(sum))

# Remove first summary row
EMTAB5456_count <- EMTAB5456_count[-c(1),]

# Prepare metadata and matrix
EMTAB5456_count <- column_to_rownames(EMTAB5456_count, var = "hgnc_symbol")
EMTAB5456_count <- as.matrix(EMTAB5456_count)
EMTAB5456_count <- round(EMTAB5456_count)

# Rearrange metadata table
metadata <- column_to_rownames(metadata, var = "ID")
metadata$Number <- as.factor(metadata$Number)

# Merge HSC and MPP counts
EMTAB5456_merge <- as.data.frame(EMTAB5456_count)

EMTAB5456_merge$ERR1816850_sum <- EMTAB5456_merge$ERR1816834 + EMTAB5456_merge$ERR1816850
EMTAB5456_merge$ERR1816851_sum <- EMTAB5456_merge$ERR1816835 + EMTAB5456_merge$ERR1816851
EMTAB5456_merge$ERR1816852_sum <- EMTAB5456_merge$ERR1816836 + EMTAB5456_merge$ERR1816852
EMTAB5456_merge$ERR1816853_sum <- EMTAB5456_merge$ERR1816837 + EMTAB5456_merge$ERR1816853

EMTAB5456_merge$ERR1816834 <- NULL
EMTAB5456_merge$ERR1816835 <- NULL
EMTAB5456_merge$ERR1816836 <- NULL
EMTAB5456_merge$ERR1816837 <- NULL
EMTAB5456_merge$ERR1816850 <- NULL
EMTAB5456_merge$ERR1816851 <- NULL
EMTAB5456_merge$ERR1816852 <- NULL
EMTAB5456_merge$ERR1816853 <- NULL

EMTAB5456_merge <- EMTAB5456_merge %>%
    dplyr::rename("ERR1816850" = "ERR1816850_sum") %>%
    dplyr::rename("ERR1816851" = "ERR1816851_sum") %>%
    dplyr::rename("ERR1816852" = "ERR1816852_sum") %>%
    dplyr::rename("ERR1816853" = "ERR1816853_sum")

EMTAB5456_merge <- as.matrix(EMTAB5456_merge)

# Create merge metadata
metadata_merge <- subset(metadata, Celltype != "HSC")

EMTAB5456_batch_counts <- ComBat_seq(
    EMTAB5456_merge,
    batch=metadata_merge$Number,
    group=metadata_merge$Celltype,
)

remove(EMTAB5456_count)
remove(EMTAB5456_merge)
remove(metadata)

cell_of_interest <- c("MPP", "LMPP")

# Subset metadata on COI
metadata_dge <- subset(metadata_merge, Celltype %in% cell_of_interest)
metadata_dge$Celltype <- as.factor(metadata_dge$Celltype)

# Subset libraries on COI
EMTAB5456_dge <- as.data.frame(t(EMTAB5456_batch_counts))
EMTAB5456_dge <- rownames_to_column(EMTAB5456_dge, var="lib")
EMTAB5456_dge <- subset(EMTAB5456_dge, EMTAB5456_dge$lib %in% rownames(metadata_dge))

# Convert back to matrix
rownames(EMTAB5456_dge) <- NULL
EMTAB5456_dge <- column_to_rownames(EMTAB5456_dge, var="lib")
EMTAB5456_dge <- as.matrix(t(EMTAB5456_dge))

dds <- DESeqDataSetFromMatrix(countData=EMTAB5456_dge,
                              colData=metadata_dge,
                              design=~Celltype)

# no filtering of low-expressed genes since many of the genes of interest are TFs which may be lowly expressed
# in addition, knowing that some of these genes are not expressed in either MPP or LMPP is potentially of interest

# Create differential expression object
dds <- DESeq(dds)

res <- results(dds)

res <- as.data.frame(res)

write.csv(res, file=here("results", "EMTAB5456_dge.csv"), row.names=TRUE)
