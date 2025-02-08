# INPUT: 
# `results/{sample}/epiAneufinder_results/results_table.tsv`

# OUTPUT:
# `results/{sample}/epiAneufinder_results/{sample}_Karyotype.tsv`
# `results/{sample}/epiAneufinder_results/{sample}_karyotype.rds`
# `results/{sample}/epiAneufinder_results/{sample}_pct_cells_per_bin.rds`
# `results/{sample}/epiAneufinder_results/{sample}_pct_bins_per_arm.rds`
# `results/{sample}/epiAneufinder_results/{sample}_pct_cells_per_bin.png`
# `results/{sample}/epiAneufinder_results/{sample}_pct_bins_per_arm.png`

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparser))

script <- basename(sub(".*=", "", commandArgs()[4]))

# using the pipe
parser <- arg_parser(description = "Calls the karyotype depending on epiAneufinder results",
                     hide.opts = TRUE, name = script) %>%
    add_argument("results", 
                 help = "Path to results from epiAneufinder") %>%
    add_argument("--sample",
                 help = "Sample name") %>%
    add_argument("--outdir", 
                 help = "Output directory") %>%
    add_argument("--cell-pct-loss", 
                 help = "% of cells indicating loss to call a loss", 
                 default = 15) %>%
    add_argument("--cell-pct-gain", 
                 help = "% of cells indicating gain to call a gain", 
                 default = 10) %>%
    add_argument("--diff-pct", 
                 help = "if both gain and loss meet thresholds, take the larger percentage given there is x percent difference", 
                 default = 5) %>%
    add_argument("--bin-pct-loss", 
                 help = "% of bins in chromosome arm indicating loss to to call a loss", 
                 default = 40) %>%
    add_argument("--bin-pct-gain", 
                 help = "% of bins in chromosome arm indicating gain to call a gain", 
                 default = 40) %>%
    add_argument("--arm-coverage", 
                 help = "% of arm that has bins to make a call",
                 default = 50) 

argv <- parse_args(parser)

# for debugging
argv <- list(
    results =  "/projects/karsanlab/dlin_dev/mds_aza/KARSANBIO-4070_Call_copy_number_on_ATAC_data/results/H20KM3840/epiAneufinder_results/results_table.tsv",
    sample = "H20KM3840",
    outdir = NA,
    cell_pct_loss = 15,
    cell_pct_gain = 10,
    bin_pct_loss = 18,
    bin_pct_gain = 18,
    arm_coverage = 50,
    diff_pct = 5
)

sample <- argv$sample
results <- read.table(argv$results)

cell_pct_loss <- argv$cell_pct_loss
cell_pct_gain <- argv$cell_pct_gain

bin_pct_loss <- argv$bin_pct_loss
bin_pct_gain <- argv$bin_pct_gain

if(is.na(argv$outdir)) {
    outdir <- dirname(argv$results)
} else {
    outdir <- argv$outdir
}

font.family <- "Helvetica"
theme_set(theme_bw(base_family = font.family,
                   base_size = 14) +
              theme(legend.box.spacing = unit(0, "pt")))

peaks <- results[,1:3]
colnames(peaks) <- c("seqnames", "start", "end")
peaks$seqnames <- factor(peaks$seqnames, levels = paste0("chr", 1:22))

# http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
cytobands <- suppressMessages(read_tsv(here("resources", "cytoBand.txt"), col_names = F))


arm_len <- suppressMessages(
    {cytobands %>%
            mutate(arm = str_extract(X4, "^[pq]"), 
                   length = X3-X2) %>%
            group_by(X1, arm) %>%
            summarize(length = sum(length)) %>%
            drop_na(arm) %>%
            dplyr::rename(chr = X1) %>%
            mutate(chr = factor(chr, paste0("chr", 1:22)))})

centromeres <- cytobands %>%
    filter(X5 == "acen", str_detect(X4, "^q"), !X1 %in% c("chrX", "chrY")) %>%
    select(chr = X1, centromere = X2) %>%
    mutate(chr = factor(chr, levels = paste0("chr", 1:22))) %>%
    arrange(chr) 

# for plotting p and q on the dendrograms
bins_rn <- peaks %>%
    arrange(seqnames, start, end) %>%
    mutate(rn = row_number())

centromeres_rn <- centromeres %>%
    dplyr::rename(seqnames = chr, start = centromere) %>%
    mutate(end = start, type = "centromere") %>%
    bind_rows(bins_rn %>%
                  mutate(type = "bin")) %>%
    arrange(seqnames, start) 

# replace all the NAs with 0.5 of the previous one to draw on grid
centromeres_rn$rn[which(is.na(centromeres_rn$rn) == TRUE)] <- centromeres_rn$rn[which(is.na(centromeres_rn$rn) == TRUE) - 1] + 0.5

centromeres_rn <- centromeres_rn %>%
    filter(type == "centromere") %>%
    select(seqnames, rn) 

centromeres.list <- centromeres %>%
    tibble::deframe() %>%
    as.list()

arms <- imap(centromeres.list, function(coord, chr) {
    
    rownames(peaks) <- NULL
    peaks %>%
        filter(seqnames == chr) %>%
        mutate(arm = case_when(start < coord & end <= coord ~ "p", 
                               start >= coord & end > coord ~ "q",
                               .default = NA))
}) %>% 
    list_rbind() %>%
    unite(col = "bin", seqnames, start, end, sep = "-")

df <- results %>%
    pivot_longer(-c(seq, start, end), names_to = "cell", values_to = "somy") %>%
    unite(col = "bin", seq, start, end, sep = "-")

cell_counts_per_bin <- df %>%
    count(bin, somy, name = "num_cells_per_bin")

total_cells_per_bin <- df %>%
    count(bin, name = "total_cells_per_bin")

pct_cells <- full_join(cell_counts_per_bin, total_cells_per_bin, by = "bin") %>%
    mutate(pct_cells = num_cells_per_bin/total_cells_per_bin*100,
           somy = case_match(somy, 0 ~ "loss", 1 ~ "normal", 2 ~ "gain")) %>%
    select(-num_cells_per_bin, -total_cells_per_bin)


pct_calls <- pct_cells %>%
    pivot_wider(names_from = "somy", values_from = "pct_cells")

walk(c("normal", "gain", "loss"), ~ {
    if(!.x %in% colnames(pct_calls)) {
        pct_calls[[.x]] <<- 0
    }
})

pct_calls <- pct_calls %>%
    replace_na(list(normal = 0, gain = 0, loss = 0)) %>%
    mutate(chr = str_extract(bin, "chr[0-9]{1,2}"), 
           chr = factor(chr, levels = paste0("chr", 1:22))) %>%
    dplyr::relocate(chr, .after = "bin") 

pct_calls <- pct_calls %>%
    mutate(call = case_when(
        loss >= cell_pct_loss & gain >= cell_pct_gain & loss > gain & loss - gain > argv$diff_pct ~ "loss", 
        loss >= cell_pct_loss & gain >= cell_pct_gain & gain > loss & gain - loss > argv$diff_pct ~ "gain", 
        loss >= cell_pct_loss & gain >= cell_pct_gain ~ "normal", 
        loss >= cell_pct_loss ~ "loss", 
        gain >= cell_pct_gain ~ "gain", 
        .default = "normal"
    ), 
    call = factor(call, levels = c("loss", "normal", "gain"))) %>%
    left_join(arms, by = "bin") %>%
    dplyr::relocate(arm, .after = "chr")
saveRDS(pct_calls, file.path(outdir, str_glue("{sample}_pct_cells_per_bin.rds")))

cell.thresholds <- tibble(
    call = c(rep("loss", 22), rep("gain", 22)), 
    chr = rep(paste0("chr", 1:22), 2),
    pct = c(rep(cell_pct_loss, 22), rep(cell_pct_gain, 22))
) %>%
    mutate(call = factor(call, levels = c("loss", "gain", "normal")), 
           chr = factor(chr, paste0("chr", 1:22)))

cell_plot <- pct_calls %>% 
    select(-call) %>%
    pivot_longer(-c(chr, bin, arm), names_to = "call", values_to = "pct") %>%
    mutate(call = factor(call, levels = c("loss", "gain", "normal"))) %>%
    ggplot(aes(x = 1, y = pct)) +
    facet_grid(call ~ chr, scales = "free") +
    geom_boxplot(outlier.shape = NA, alpha = 0) +
    ggforce::geom_sina(size = 0.1, alpha = 0.05) +
    geom_hline(data = cell.thresholds, 
               aes(yintercept = pct), linetype = "dashed") + 
    labs(y = "% Cells", title = sample) +
    theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank())
ggsave(file.path(outdir, str_glue("{sample}_pct_cells_per_bin.png")), cell_plot, width = 16, height = 9, units = "in")


# P/Q ARM ----
bins_per_arm <- pct_calls %>%
    count(chr, arm, call, name = "bins_per_arm")

total_bins_per_arm <- pct_calls %>%
    count(chr, arm, name = "total_bins_per_arm")

pct_bins_per_arm <- full_join(bins_per_arm, total_bins_per_arm, by = c("chr", "arm")) %>%
    mutate(pct_bins_per_arm = bins_per_arm/total_bins_per_arm*100) %>%
    select(-bins_per_arm, -total_bins_per_arm) %>%
    pivot_wider(names_from = "call", values_from = "pct_bins_per_arm", values_fill = 0)

walk(c("normal", "gain", "loss"), ~ {
    if(!.x %in% colnames(pct_bins_per_arm)) {
        pct_bins_per_arm[[.x]] <<- 0
    }
})
pct_bins_per_arm <- pct_bins_per_arm %>%
    replace_na(list(normal = 0, gain = 0, loss = 0))

pct_bins_per_arm <- pct_bins_per_arm %>%
    mutate(call = case_when(
        gain >= bin_pct_gain & loss >= bin_pct_loss  & loss > gain & loss - gain > argv$diff_pct ~ "loss",
        gain >= bin_pct_gain & loss >= bin_pct_loss & gain > loss & gain - loss > argv$diff_pct ~ "gain",
        gain >= bin_pct_gain & loss >= bin_pct_loss ~ "normal", 
        loss >= bin_pct_loss ~ "loss", 
        gain >= bin_pct_gain ~ "gain", 
        .default = "normal"
    ),
    chr = factor(chr, levels = paste0("chr", 1:22)))
saveRDS(pct_bins_per_arm, file.path(outdir, str_glue("{sample}_pct_bins_per_arm.rds")))

bin.thresholds <- tibble(
    call = c(rep("loss", 22), rep("gain", 22)), 
    chr = rep(paste0("chr", 1:22), 2),
    pct = c(rep(bin_pct_loss, 22), rep(bin_pct_gain, 22))
) %>%
    mutate(call = factor(call, levels = c("loss", "gain", "normal")), 
           chr = factor(chr, paste0("chr", 1:22)))

bin_plot <- pct_bins_per_arm %>%
    select(-call) %>%
    pivot_longer(-c(chr, arm), names_to = "call", values_to = "pct") %>%
    mutate(call = factor(call, levels = c("loss", "gain", "normal")))  %>%
    ggplot(aes(x = arm, y = pct)) +
    geom_point() + 
    facet_grid(call ~ chr) +
    geom_hline(data = bin.thresholds, 
               aes(yintercept = pct), linetype = "dashed") + 
    labs(y = "% Bins", x = "Chromosome Arm", title = sample) +
    ylim(c(0,100))

ggsave(file.path(outdir, str_glue("{sample}_pct_bins_per_arm.png")), bin_plot, width = 16, height = 9, units = "in")


if(all(unique(pct_bins_per_arm$call) == "normal")) {
    karyotype <- "normal"
} else {
    karyo <- pct_bins_per_arm %>%
        group_by(chr, call) %>%
        select(chr, arm, call) %>%
        mutate(arm = paste(arm, collapse = ",")) %>%
        distinct() %>%
        arrange(chr, arm) %>%
        ungroup() %>%
        mutate(karyotype = case_when(call == "gain" & arm == "p" ~ paste0("+", str_extract(chr, '[0-9]{1,2}$'), "p"),
                                     call == "gain" & arm == "q" ~ paste0("+", str_extract(chr, '[0-9]{1,2}$'), "q"),
                                     call == "gain" & arm == "p,q" ~ paste0("+", str_extract(chr, '[0-9]{1,2}$')), 
                                     call == "loss" & arm == "p" ~ paste0("del(", str_extract(chr, '[0-9]{1,2}$'), "p)"),
                                     call == "loss" & arm == "q" ~ paste0("del(", str_extract(chr, '[0-9]{1,2}$'), "q)"),
                                     call == "loss" & arm == "p,q" ~ paste0("-", str_extract(chr, '[0-9]{1,2}$'))
        )
        ) %>%
        drop_na(karyotype) %>%
        arrange(chr)
    
    bin_size <- unique(peaks$end - peaks$start + 1)
    
    # filter out those that don't have enough bin coverage to make a call 
    cov <- suppressMessages({pct_calls %>%
            mutate(bin_size = bin_size) %>%
            group_by(chr, arm) %>%
            summarize(length = sum(bin_size)) %>%
            left_join(arm_len, by = c("chr", "arm"), suffix = c(".coverage", "")) %>%
            ungroup() %>%
            mutate(coverage = length.coverage/length*100)})
    
    cov_plot <- cov %>%
        ggplot(aes(x = arm, y = coverage)) +
        geom_point() +
        facet_wrap(~ chr, nrow = 1) +
        geom_hline(yintercept = argv$arm_coverage, linetype = "dashed") +
        labs(y = "% Arm Coverage", x = "Chromosome Arm") +
        ylim(c(0,100))
    
    ggsave(file.path(outdir, str_glue("{sample}_arm_coverage.png")), cov_plot, width = 16, height = 9, units = "in")
    drop <- cov %>%
        filter(coverage < argv$arm_coverage)
    
    if(nrow(drop) > 1) {
        dropped <- karyo %>%
            semi_join(drop, by = c("chr", "arm"))
        warning(str_glue("The following calls have been dropped due to < {argv$arm_coverage}% coverage: {paste(dropped$karyotype, collapse = ', ')}"))
    }
    
    karyo <- anti_join(karyo, drop, by = c("chr", "arm"))
    karyotype <- paste(karyo$karyotype, collapse = ", ")
    
    bind_rows(
        mutate(karyo, keep = TRUE),
        mutate(dropped, keep = FALSE),
    ) %>%
        saveRDS(str_glue("{sample}_karyotype.rds"))
    
    if(karyotype == "") {
        karyotype <- "normal"
    }
}


message(str_glue("{sample}: {karyotype}"))

tibble(sample = sample, 
       karotype = karyotype) %>%
    write_tsv(file.path(outdir, str_glue("{sample}_Karyotype.tsv")), col_names = F)
