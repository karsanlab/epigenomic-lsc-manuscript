# INPUT: 
# `results/*.cons.Cmethyl_tfbs.csv`
# `results/combined_methylation_qc_metadata.rds`

# OUTPUT: 
# `results/dm_tfbs_Responder_vs_Non-responder.tsv`
# `results/top_dmtfbs.tsv`
# `results/tfbs_heatmap.rds`
# `results/tfbs_volcano.rds`

# `results/normal_dm_tfbs_MPP_vs_LMPP_t_test.tsv`
# `results/top_normal_dmtfbs.tsv`
# `results/tfbs_volcano_normal.rds`
# `results/intersect_normal_response_tfs.rds`

library(LOLA)
require(simpleCache)
library(rtracklayer)
library(tidyverse)
library(data.table)
library(pbkrtest)
library(afex)
library(here)
library(pbapply)

# Load in LOLA data:

lola.db <- loadRegionDB(here('resources', 'supp_data', 'LOLA', 'LOLACore', 'hg19'), 
                        collections = c('encode_tfbs', 'codex', 'custom', 'cistrome_cistrome'))

lola.tfbs <- lola.db$regionAnno

lola.tfbs.gr <- lola.db$regionGRL

names(lola.tfbs.gr) <- make.unique(paste(lola.tfbs$antibody, lola.tfbs$cellType))

# Load in sample data:
sample.data.all <- readRDS(here('results', 'combined_methylation_qc_metadata.rds')) %>%
    mutate(rowname = plate_sample) %>%
    column_to_rownames("rowname")

# Load in LOLA TFBS methylation data:
tfbs.tab <- list.files("results", full.names = T, pattern = "*_tfbs.csv") %>%
    str_subset("[0-9]{3}_PX[0-9]{4}_") %>%
    map(~ read.csv(.x, row.names = 1)) %>%
    list_cbind()

colnames(tfbs.tab) <- list.files("results", pattern = "*_tfbs.csv") %>% 
  basename() %>% 
  str_subset("[0-9]{3}_PX[0-9]{4}_") %>% 
  str_remove(".cons.Cmethyl_tfbs.csv") %>% 
  str_remove("^[0-9]{3}_")

# Filter by coverage:
tfbs.filt <- tfbs.tab[,intersect(colnames(tfbs.tab), rownames(sample.data.all))]
tfbs.filt <- tfbs.filt[,which(sample.data.all[intersect(colnames(tfbs.filt), rownames(sample.data.all)),'passed.qc'])]

sample.data.tfbs <- sample.data.all[colnames(tfbs.filt),]

#' function to test differential methylation over tfbs sites in the given conditions.

calculateForOneCondition <- function(condition.column, condition.1, condition.2, sample.data.this){
  message(paste0("Performing tests for ", condition.column, ": ", condition.1, " vs ", condition.2))
  
  testOne <- function(tfbs.name, condition.column, condition.1, condition.2){
    
    this.data <- sample.data.this[intersect(rownames(sample.data.all),colnames(tfbs.tab)),]
    this.data <- this.data[which(this.data[,condition.column] %in% c(condition.1, condition.2)),]
    this.data[,condition.column] <- factor(this.data[,condition.column])
    this.data$this.meth <- unlist(tfbs.tab[tfbs.name,rownames(this.data)])
    
    pbkr.model <- suppressMessages(invisible(mixed(as.formula(paste0('this.meth ~ ', condition.column, ' + (1|patient)')), data=this.data)))

    mean.meth <- this.data %>% group_by(!!as.name(condition.column)) %>% dplyr::summarise(mean.meth=mean(this.meth))
    delta.meth <- dplyr::filter(mean.meth, get(condition.column)==condition.1)$mean.meth - 
                  dplyr::filter(mean.meth, get(condition.column)==condition.2)$mean.meth

    res <- c(pbkr.model$anova_table$`Pr(>F)`, delta.meth)
    names(res) <- c('pbkr.pval', 'delta.meth')
    return(res)
  }
  
  suppressMessages(
    tfbs.pre.tests <- pbsapply(rownames(tfbs.tab), testOne, condition.column, condition.1, condition.2 )
  )
  
  tfbs.pre.tests <- as.data.frame(t(tfbs.pre.tests))
  tfbs.pre.tests$pbkr.qval <- p.adjust(tfbs.pre.tests$pbkr.pval, method='BH')
  tfbs.pre.tests$pbkr.signif <- tfbs.pre.tests$pbkr.qval < 0.05 & abs(tfbs.pre.tests$delta.meth) > 5
  
  tfbs.pre.tests$dnabp <- rep('other', nrow(tfbs.pre.tests))
  tfbs.pre.tests$dnabp[which(grepl('ctcf', tolower(rownames(tfbs.pre.tests))))] <- 'CTCF'
  tfbs.pre.tests$dnabp[which(grepl('ezh', tolower(rownames(tfbs.pre.tests))) | grepl('suz', tolower(rownames(tfbs.pre.tests))))] <- 'PRC2'
  tfbs.pre.tests$dnabp[which(grepl('tal', tolower(rownames(tfbs.pre.tests))) |
                               grepl('gata', tolower(rownames(tfbs.pre.tests))))] <- 'TAL1/GATA1/2'
  
  tfbs.pre.tests$dnabp[which(grepl('nrsf', tolower(rownames(tfbs.pre.tests))) |
                               grepl('rest', tolower(rownames(tfbs.pre.tests))))] <- 'REST'
  
  tfbs.pre.tests$dnabp[which(grepl('SMC|RAD|STAG', toupper(rownames(tfbs.pre.tests))))] <- 'Cohesin'
  
  tfbs.pre.tests$dnabp[which(grepl('jun', tolower(rownames(tfbs.pre.tests))) |
                               grepl('fos', tolower(rownames(tfbs.pre.tests))))] <- 'Jun/Fos'
  
  tfbs.pre.tests$dnabp[which(grepl('polii', tolower(rownames(tfbs.pre.tests))) |
                               grepl('pol2', tolower(rownames(tfbs.pre.tests))))] <- 'POL2'
  
  message(paste0("Writing tsv for ", condition.column, ": ", condition.1, " vs ", condition.2))
  write.table(tfbs.pre.tests,
         file.path('results', paste0('dm_tfbs_',
                                                       condition.column,'_',
                                                       condition.1, '_vs_',
                                                       condition.2, '.tsv')))
  return(tfbs.pre.tests)
}


tfbs.pre.tests <- calculateForOneCondition('Class' , 'Responder', 'Non-responder',
                                         sample.data.all)

tfbs.pre.tests <- tfbs.pre.tests %>%
    rownames_to_column("filename")

# preserve filename as rowname too 
rownames(tfbs.pre.tests) <- tfbs.pre.tests$filename

tfbs.pre.tests$pbkr.signif <- tfbs.pre.tests$pbkr.qval < 0.1 

tfbs.pre.tests$qval.cat <- rep('Q > 0.1', nrow(tfbs.pre.tests))
tfbs.pre.tests$qval.cat[which(tfbs.pre.tests$pbkr.signif)] <- 'Q < 0.1'

tfbs.pre.tests <- tfbs.pre.tests %>%
    left_join(lola.tfbs, by = c("filename" = "filename"))

encode_antibodies <- read_csv(here('resources', 'supp_data', 'ENCODE_metadata', 'encode_antibody_metadata.csv')) %>%
    select(-Lab)
encode_cell_lines <- read_csv(here('resources', 'supp_data', 'ENCODE_metadata', 'encode_cell_line_metadata.csv'))

encode_cell_lines$Tissue[which(encode_cell_lines$Label == "Jurkat")]

# Fix up antibody & tf target information
tfbs.pre.tests <- tfbs.pre.tests %>%
    left_join(encode_antibodies, by = c("antibody" = "antibody")) %>% 
    mutate(tf = ifelse(is.na(Target),
                       toupper(antibody),
                       toupper(Target)))
tfbs.pre.tests$tf[which(grepl('POL2\\(B\\)', tfbs.pre.tests$tf))] <- 'POLR2B'
tfbs.pre.tests$tf[which(grepl('POLLL', tfbs.pre.tests$tf))] <- 'POLR2A'
tfbs.pre.tests$tf[which(grepl('POLIIS', tfbs.pre.tests$tf))] <- 'POLR2A'
tfbs.pre.tests$tf[which(grepl('POLII', tfbs.pre.tests$tf))] <- 'POLR2A'

tfbs.pre.tests$tf[which(grepl('SIN3AK', tfbs.pre.tests$tf))] <- 'SIN3A'
tfbs.pre.tests$tf[which(grepl('EGFP-GATA2', tfbs.pre.tests$tf))] <- 'GATA2'
tfbs.pre.tests$tf[which(grepl('EGFP-HDAC8', tfbs.pre.tests$tf))] <- 'HDAC8'

tfbs.pre.tests$tf[which(grepl('EGFP-JUNB', tfbs.pre.tests$tf))] <- 'JUNB'
tfbs.pre.tests$tf[which(grepl('EGFP-JUND', tfbs.pre.tests$tf))] <- 'JUND'
tfbs.pre.tests$tf[which(grepl('EGFP-FOS', tfbs.pre.tests$tf))] <- 'FOS'

tfbs.pre.tests$tf[which(grepl('C_MYC', tfbs.pre.tests$tf))] <- 'MYC'

tfbs.pre.tests$tf[which(grepl('ESR1,', tfbs.pre.tests$tf))] <- 'ESR1'

tfbs.pre.tests$tf[which(grepl('ERG\\+AR', tfbs.pre.tests$tf))] <- 'ERG'

tfbs.pre.tests$tf[which(grepl('^P300', tfbs.pre.tests$tf))] <- 'EP300'

# remove extra antibody columns
tfbs.pre.tests <- tfbs.pre.tests %>%
    select(-c(`Antibody Description`, `Target Description`, `Vendor ID`, 
              Documents, Lots, `Target Link`, Label))

# clean up cell line names
tfbs.pre.tests$cellType[which(str_detect(tfbs.pre.tests$description, "Jurkat"))] <- "Jurkat"
tfbs.pre.tests$cellType[which(str_detect(tfbs.pre.tests$description, "Kasumi-1"))] <- "Kasumi-1"

tfbs.pre.tests$cellType[which(str_detect(tfbs.pre.tests$description, "K562") &
                                  str_detect(tfbs.pre.tests$cellType, "^Erythro"))] <- "K562"
tfbs.pre.tests$cellType[which(str_detect(tfbs.pre.tests$description, "K562") &
                                  str_detect(tfbs.pre.tests$cellType, "^Myeloid"))] <- "K562"
tfbs.pre.tests$cellType[which(str_detect(tfbs.pre.tests$description, "K562") &
                                  str_detect(tfbs.pre.tests$cellType, "^Chronic"))] <- "K562"

tfbs.pre.tests$cellType[which(str_detect(tfbs.pre.tests$cellType, "A594"))] <- "A549"
tfbs.pre.tests$cellType[which(str_detect(tfbs.pre.tests$cellType, "^hESC"))] <- "H1-hESC"

# add tissue types
tfbs.pre.tests <- tfbs.pre.tests %>%
    left_join(encode_cell_lines, by = c("cellType" = "cell")) %>% 
    mutate(tissue_type = ifelse(!is.na(Tissue),
                                Tissue,
                                ifelse(!is.na(tissue),
                                       tissue,
                                       cellType)))

# remove extra cell line columns
tfbs.pre.tests <- tfbs.pre.tests %>%
    select(-c(Description, Tier, Lineage, Sex, Documents, `Vendor ID`, `Term ID`, Label))

# estrogen receptors to categories
tfbs.pre.tests <- tfbs.pre.tests %>% 
    mutate(dnabp = ifelse(str_detect(tf, "^ESR"),
                          "ESR1/2",
                          dnabp))

# tfbs.pre.tests is used for the volcano plot
tfbs.volcano <- tfbs.pre.tests %>%
    filter(str_detect(treatment, "^no|^No|None") | is.na(treatment))

saveRDS(tfbs.volcano, here("results", "tfbs_volcano.rds")) 

most.sig <- tfbs.pre.tests %>%
  filter(pbkr.signif, delta.meth < -4) %>%
  select(-tissue, -Tissue) %>% 
  arrange(pbkr.pval) %>% 
  filter(str_detect(treatment, "^no|^No|None") | is.na(treatment)) %>% 
  select(tf, cellType, pbkr.pval, delta.meth, pbkr.qval, dnabp, tissue_type, treatment, filename) %>%
  arrange(tf, pbkr.pval) %>%
  group_by(tf) %>%
  mutate(rank = row_number()) %>%
  filter(rank == 1) %>%
  select(-rank)

write_tsv(most.sig, here("results", "top_dmtfbs.tsv"))

controls <- tfbs.pre.tests %>%
  arrange(abs(delta.meth), desc(pbkr.pval)) %>%
  filter(!tf %in% pull(most.sig, tf)) %>%
  group_by(tf) %>%
  arrange(abs(delta.meth), desc(pbkr.pval)) %>%
  slice(1) %>%
  arrange(abs(delta.meth)) %>%
  filter(pbkr.signif == FALSE, !tf %in% c("NELFe", "FOS", "RDBP", "RUNX1T1")) %>% # these are not found in DARs enough for statistics
  head(8) %>% 
  select(-tf_class) %>%
  mutate(
    tf = case_match(tf,
                    "MYB" ~ "Control 1 (MYB)",
                    "NOTCH1" ~ "Control 2 (NOTCH1)",
                    "BRCA1" ~ "Control 3 (BRCA1)",
                    "RBPJ" ~ "Control 4 (RBPJ)",
                    "IRF3" ~ "Control 5 (IRF3)",
                    "TBP" ~ "Control 6 (TBP)",
                    "ETS1" ~ "Control 7 (ETS1)",
                    "CDK7" ~ "Control 8 (CDK7)"
    ),
    tf = factor(tf)) %>%
  select(all_of(colnames(most.sig)[-length(colnames(most.sig))]))

saveRDS(controls, here("results", "controls.rds"))

# Heatmap Info
saveRDS(tfbs.filt, here("results", "tfbs_heatmap.rds")) 

# Normal data
normal.data <- tfbs.tab[,str_subset(colnames(tfbs.tab), "^GSM")] %>%
  select(contains(c("MPP", "MLP0"))) %>%
  rownames_to_column("filename") %>%
  pivot_longer(cols = -filename, names_to = c("Sample", "Celltype", "x", "D", "y", "z"), values_to = "pct_meth", names_sep = "_") %>%
  select(-c(Sample, x,y,z)) 

# repeat but for the normal data
t_test_one <- function(tfbs.name, condition.column, condition.1, condition.2){
  
  this.data <- normal.data %>%
    filter(filename == tfbs.name)
  
  pbkr.model <- suppressMessages(invisible(rstatix::wilcox_test(as.formula(paste0('pct_meth ~ ', condition.column)), data=this.data)))
  
  mean.meth <- this.data %>% group_by(!!as.name(condition.column)) %>% dplyr::summarise(mean.meth=mean(pct_meth, na.rm = T))
  delta.meth <- dplyr::filter(mean.meth, get(condition.column)==condition.1)$mean.meth - 
    dplyr::filter(mean.meth, get(condition.column)==condition.2)$mean.meth
  
  res <- c(pbkr.model$p, delta.meth)
  names(res) <- c('pbkr.pval', 'delta.meth')
  return(res)
}

suppressMessages(
  tfbs.pre.tests <- pbsapply(lola.tfbs$filename, t_test_one, "Celltype", "MPP", "MLP0")
)

tfbs.pre.tests <- as.data.frame(t(tfbs.pre.tests))
tfbs.pre.tests$pbkr.qval <- p.adjust(tfbs.pre.tests$pbkr.pval, method='BH')
tfbs.pre.tests$pbkr.signif <- tfbs.pre.tests$pbkr.qval < 0.8

rownames(tfbs.pre.tests) <- lola.tfbs$filename

tfbs.pre.tests <- tfbs.pre.tests %>%
  left_join(encode_antibodies, by = c("antibody" = "antibody")) %>% 
  mutate(tf = ifelse(is.na(Target),
                     toupper(antibody),
                     toupper(Target)))

tfbs.pre.tests$tf[which(grepl('POL2\\(B\\)', tfbs.pre.tests$tf))] <- 'POLR2B'
tfbs.pre.tests$tf[which(grepl('POLLL', tfbs.pre.tests$tf))] <- 'POLR2A'
tfbs.pre.tests$tf[which(grepl('POLIIS', tfbs.pre.tests$tf))] <- 'POLR2A'
tfbs.pre.tests$tf[which(grepl('POLII', tfbs.pre.tests$tf))] <- 'POLR2A'

tfbs.pre.tests$tf[which(grepl('SIN3AK', tfbs.pre.tests$tf))] <- 'SIN3A'
tfbs.pre.tests$tf[which(grepl('EGFP-GATA2', tfbs.pre.tests$tf))] <- 'GATA2'
tfbs.pre.tests$tf[which(grepl('EGFP-HDAC8', tfbs.pre.tests$tf))] <- 'HDAC8'

tfbs.pre.tests$tf[which(grepl('EGFP-JUNB', tfbs.pre.tests$tf))] <- 'JUNB'
tfbs.pre.tests$tf[which(grepl('EGFP-JUND', tfbs.pre.tests$tf))] <- 'JUND'
tfbs.pre.tests$tf[which(grepl('EGFP-FOS', tfbs.pre.tests$tf))] <- 'FOS'

tfbs.pre.tests$tf[which(grepl('C_MYC', tfbs.pre.tests$tf))] <- 'MYC'

tfbs.pre.tests$tf[which(grepl('ESR1,', tfbs.pre.tests$tf))] <- 'ESR1'

tfbs.pre.tests$tf[which(grepl('ERG\\+AR', tfbs.pre.tests$tf))] <- 'ERG'

tfbs.pre.tests$tf[which(grepl('^P300', tfbs.pre.tests$tf))] <- 'EP300'

# remove extra antibody columns
tfbs.pre.tests <- tfbs.pre.tests %>%
  select(-c(`Antibody Description`, `Target Description`, `Vendor ID`, 
            Documents, Lots, `Target Link`, Label))

# clean up cell line names
tfbs.pre.tests$cellType[which(str_detect(tfbs.pre.tests$description, "Jurkat"))] <- "Jurkat"
tfbs.pre.tests$cellType[which(str_detect(tfbs.pre.tests$description, "Kasumi-1"))] <- "Kasumi-1"

tfbs.pre.tests$cellType[which(str_detect(tfbs.pre.tests$description, "K562") &
                                str_detect(tfbs.pre.tests$cellType, "^Erythro"))] <- "K562"
tfbs.pre.tests$cellType[which(str_detect(tfbs.pre.tests$description, "K562") &
                                str_detect(tfbs.pre.tests$cellType, "^Myeloid"))] <- "K562"
tfbs.pre.tests$cellType[which(str_detect(tfbs.pre.tests$description, "K562") &
                                str_detect(tfbs.pre.tests$cellType, "^Chronic"))] <- "K562"

tfbs.pre.tests$cellType[which(str_detect(tfbs.pre.tests$cellType, "A594"))] <- "A549"
tfbs.pre.tests$cellType[which(str_detect(tfbs.pre.tests$cellType, "^hESC"))] <- "H1-hESC"

# add tissue types
tfbs.pre.tests <- tfbs.pre.tests %>%
  left_join(encode_cell_lines, by = c("cellType" = "cell")) %>% 
  mutate(tissue_type = ifelse(!is.na(Tissue),
                              Tissue,
                              ifelse(!is.na(tissue),
                                     tissue,
                                     cellType)))

# remove extra cell line columns
tfbs.pre.tests <- tfbs.pre.tests %>%
  select(-c(Description, Tier, Lineage, Sex, Documents, `Vendor ID`, `Term ID`, Label))

tfbs.pre.tests$dnabp <- rep('other', nrow(tfbs.pre.tests))
tfbs.pre.tests$dnabp[which(grepl('ctcf', tolower(rownames(tfbs.pre.tests))))] <- 'CTCF'
tfbs.pre.tests$dnabp[which(grepl('ezh', tolower(rownames(tfbs.pre.tests))) | grepl('suz', tolower(rownames(tfbs.pre.tests))))] <- 'PRC2'
tfbs.pre.tests$dnabp[which(grepl('tal', tolower(rownames(tfbs.pre.tests))) |
                             grepl('gata', tolower(rownames(tfbs.pre.tests))))] <- 'TAL1/GATA1/2'

tfbs.pre.tests$dnabp[which(grepl('nrsf', tolower(rownames(tfbs.pre.tests))) |
                             grepl('rest', tolower(rownames(tfbs.pre.tests))))] <- 'REST'

tfbs.pre.tests$dnabp[which(grepl('SMC|RAD|STAG', toupper(rownames(tfbs.pre.tests))))] <- 'Cohesin'

tfbs.pre.tests$dnabp[which(grepl('jun', tolower(rownames(tfbs.pre.tests))) |
                             grepl('fos', tolower(rownames(tfbs.pre.tests))))] <- 'Jun/Fos'

tfbs.pre.tests$dnabp[which(grepl('polii', tolower(rownames(tfbs.pre.tests))) |
                             grepl('pol2', tolower(rownames(tfbs.pre.tests))))] <- 'POL2'

# estrogen receptors to categories
tfbs.pre.tests <- tfbs.pre.tests %>% 
  mutate(dnabp = ifelse(str_detect(tf, "^ESR"),
                        "ESR1/2",
                        dnabp))

data.table::fwrite(tfbs.pre.tests,
                   here("results", "normal_dm_tfbs_MPP_vs_LMPP_t_test.tsv"),
                   row.names = TRUE)

# tfbs.pre.tests is used for the volcano plot
tfbs.volcano <- tfbs.pre.tests %>%
  filter(str_detect(treatment, "^no|^No|None") | is.na(treatment))

saveRDS(tfbs.volcano, here("results", "tfbs_volcano_normal.rds")) 

most.sig.normal <- tfbs.pre.tests %>%
  filter(pbkr.qval < 0.05, delta.meth <= -1) %>%
  select(-tissue, -Tissue) %>% 
  arrange(pbkr.pval) %>% 
  filter(str_detect(treatment, "^no|^No|None") | is.na(treatment)) %>% 
  select(tf, cellType, pbkr.pval, delta.meth, pbkr.qval, dnabp, tissue_type, treatment, filename) %>%
  arrange(tf, pbkr.pval) %>%
  group_by(tf) %>%
  mutate(rank = row_number()) %>%
  filter(rank == 1) %>%
  select(-rank)

write_tsv(most.sig.normal, here("results", "top_normal_dmtfbs.tsv"))

saveRDS(list(patient = most.sig$tf, 
             normal = most.sig.norma$tf, 
             intersect = intersect(most.sig.normal$tf, most.sig$tf)), here("results", "intersect_normal_response_tfs.rds"))
