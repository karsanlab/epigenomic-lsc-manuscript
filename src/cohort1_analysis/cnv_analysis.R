# INPUT: 
# `results/cnv/*.deduped.bam_ratio.txt`

# OUTPUT: 
# `results/cnvs.rds`

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(here)

CUTOFF_ABNORMAL = 100

sample.data.all <- readRDS(here("data", "full_metadata_qc.rds"))

rownames(sample.data.all) <- sample.data.all$plate_sample
sample.data.all$patient <- as.factor(sample.data.all$patient)
sample.data.all$sample <- factor(sample.data.all$sample, levels = names(col.sample)[which(names(col.sample) %in% unique(sample.data.all$sample) ) ])

sample.data.all$Class <- factor(sample.data.all$Class, levels = c("Responder", "Non-responder"))

cnv.files <- list.files(here("results", "cnv"), pattern = ".deduped.bam_ratio.txt",
                        recursive = T, full.names = T)

loadOneCNV <- function(file.name)
{
    tmp.frame <- read.csv(file.name, sep='\t')
    # Fix chromosome ordering:
    tmp.chr <- as.character(tmp.frame$Chromosome)
    tmp.chr.isnum <- which(!tmp.chr %in% c('M', 'X', 'Y'))
    tmp.frame[tmp.chr.isnum,] <- tmp.frame[order(as.numeric(tmp.chr[tmp.chr.isnum])),]
    tmp.frame
}

cnv.frames <- lapply(cnv.files, loadOneCNV)
names(cnv.frames) <- basename(cnv.files) %>%
    str_remove(".deduped.bam_ratio.txt") %>%
    str_to_upper()

cnv.gr <- lapply(cnv.frames, function(x){
    x$End <- x$Start + 5*10^6; 
    x$Strand <- rep('*', nrow(x));
    x$Chromosome <- paste0('chr',x$Chromosome)
    return(x)})

cnv.gr <- lapply(cnv.gr, makeGRangesFromDataFrame, keep.extra.columns=TRUE)

cnv.data <- lapply(cnv.frames, function(x){tmp <- x$CopyNumber; tmp[which(tmp>5)] <- NA; tmp;})
cnv.frame <- do.call(cbind, cnv.data)

colnames(cnv.frame) <- basename(cnv.files) %>%
    str_remove(".deduped.bam_ratio.txt") %>%
    str_to_upper()
rownames(cnv.frame) <- paste(cnv.frames[[1]]$Chromosome, cnv.frames[[1]]$Start, sep='_')

cnv.frame <- cnv.frame[,intersect(sample.data.all$plate_sample, colnames(cnv.frame))]

cnv.abnormal <- apply(cnv.frame, 2, function(x){length(which(x!=2)) > CUTOFF_ABNORMAL})

cnv.ab.frame <- t(cnv.frame[,which(cnv.abnormal)])

template.gr <- cnv.frames[[1]]
colnames(template.gr)[1:2] <- c('chr', 'start')
template.gr$end <- template.gr$start + 5*10^6
template.gr$chr <- paste0('chr', template.gr$chr)
template.gr <- makeGRangesFromDataFrame(template.gr)

replic.mask <- data.table::fread(here('resources', 'supp_data', 'replication_timing', 'plastic_masked.txt'))
colnames(replic.mask)[1:3] <- c('chr', 'start', 'end')
replic.mask.gr <- makeGRangesFromDataFrame(replic.mask)
hg18.19.chain <- rtracklayer::import.chain(here('resources', 'supp_data', 'replication_timing', 'hg18ToHg19.over.chain'))

replic.mask.gr.hg19 <- liftOver(replic.mask.gr, hg18.19.chain)
replic.mask.gr.hg19 <- unlist(replic.mask.gr.hg19)


replic.files <- dir(here('resources', 'supp_data', 'replication_timing'), pattern='bigWig')
replic.gr <- lapply(here('resources', 'supp_data', 'replication_timing' ,replic.files), import.bw)

filterOutMasked <- function(data.gr, mask.gr)
{
    mask.overlap <- findOverlaps(mask.gr, data.gr )
    return(data.gr[-c(unique(mask.overlap@to))])
}

replic.filtered <- lapply(replic.gr, filterOutMasked, replic.mask.gr.hg19)

averageOverBins <- function(data.gr, template.gr, feature.col='score')
{
    overlaps <- findOverlaps(template.gr, data.gr)
    overlaps <- as.data.frame(overlaps)
    overlaps <- cbind(overlaps, data.gr[overlaps$subjectHits,]@elementMetadata[,feature.col] )
    colnames(overlaps)[3] <- feature.col
    by.feature.score <- group_by(overlaps, queryHits) %>%
        dplyr::summarise(mean.score=mean(score), num.cpgs=length(score))
    by.feature.score <- as.data.frame(by.feature.score)
    
    scored.features <- template.gr
    scored.features$score <- rep(0, length(scored.features))
    scored.features@elementMetadata$score[by.feature.score[,1]] <- by.feature.score$mean.score
    
    return(scored.features)
}

replic.perc <- lapply(replic.filtered, averageOverBins, template.gr)

replic.tab <- do.call(cbind, lapply(replic.perc, function(x){x$score}))
colnames(replic.tab) <- sub('PctSignalRep1.bigWig', '', sub('wgEncodeUwRepliSeq', '', replic.files))
rownames(replic.tab) <- colnames(cnv.ab.frame)


# Reorder in cell cycle order:
replic.tab <- replic.tab[,c("K562G1", "K562S1", "K562S2", "K562S3", "K562S4", "K562G2")]

replic.early <- apply(replic.tab, 1, function(x){sum(x[c('K562G1', 'K562S1')])}) > 50
replic.late <- apply(replic.tab, 1, function(x){sum(x[c('K562S4', 'K562G2')])}) > 50

cnv.low.cn <- apply(cnv.ab.frame, 1, function(x){x < 2})
cnv.high.cn <- apply(cnv.ab.frame, 1, function(x){x > 2})


cnv.early.test <- apply(cnv.low.cn, 2, function(lib){if(length(table(lib))==1){return(NA)}; fisher.test(lib, !replic.early, alternative='greater')$p.value})
cnv.early.adj <- p.adjust(cnv.early.test, method='BH')

cnv.late.test <- apply(cnv.high.cn, 2, function(lib){if(length(table(lib))==1){return(NA)} ; fisher.test(lib, !replic.late, alternative='greater')$p.value})
cnv.late.adj <- p.adjust(cnv.late.test, method='BH')

cnv.late.test.2 <- apply(cnv.high.cn, 2, function(lib){if(length(table(lib))==1){return(NA)} ; fisher.test(lib, replic.early, alternative='greater')$p.value})
cnv.late.adj.2 <- p.adjust(cnv.late.test.2, method='BH')

sample.data.all$cycle.phase <- rep('G0', nrow(sample.data.all))
sample.data.all[rownames(cnv.ab.frame), 'cycle.phase'] <- rep('?', nrow(cnv.ab.frame))
sample.data.all[rownames(cnv.ab.frame), 'cycle.phase'][which(cnv.late.adj.2 < 0.05| cnv.early.adj < 0.05 )] <- 'G1/S1'
sample.data.all[rownames(cnv.ab.frame), 'cycle.phase'][which(cnv.late.adj < 0.05)] <- 'S4/G2'

# PATIENT 437 ----
this.libs <-  dplyr::filter(sample.data.all, passed.qc, patient=='437', group!='negative')$plate_sample

this.libs <- intersect(this.libs, colnames(cnv.frame))

pt.437.chr.9 <- t(cnv.frame[which(grepl('^9_',rownames(cnv.frame))),
                            this.libs])

pt.437.chr.9 <- pt.437.chr.9[order(sample.data.all[rownames(pt.437.chr.9),'Class'],
                                   sample.data.all[rownames(pt.437.chr.9),'group']),]


pt.437.gr <- cnv.gr[intersect(names(cnv.gr), dplyr::filter(sample.data.all, passed.qc, patient=='437', cycle.phase=='G0')$plate_sample)]

pt.437.gr.9 <- lapply(pt.437.gr, function(x){as.data.frame(x[which(as.character(seqnames(x))=='chr9')])})

pt.437.frame <- bind_rows(pt.437.gr.9, .id='plate_sample')

pt.437.frame <- merge(pt.437.frame, sample.data.all)
pt.437.frame$CopyNumber <- as.factor(pt.437.frame$CopyNumber)

# PATIENT 400 ----
pt.400.chr.7 <- t(cnv.frame[which(grepl('^7_',rownames(cnv.frame))), ])

pt.400.chr.7 <- pt.400.chr.7[order(sample.data.all[rownames(pt.400.chr.7),'Class'],
                                   sample.data.all[rownames(pt.400.chr.7),'group']),]

pt.400.gr <- cnv.gr[as.character(dplyr::filter(sample.data.all, passed.qc, group!='negative', patient=='400', cycle.phase=='G0')$plate_sample)]

pt.400.gr.7 <- lapply(pt.400.gr, function(x){as.data.frame(x[which(as.character(seqnames(x))=='chr7')])})

pt.400.chr7.frame <- bind_rows(pt.400.gr.7, .id='plate_sample')

pt.400.chr7.frame <- merge(pt.400.chr7.frame, sample.data.all)
pt.400.chr7.frame$CopyNumber <- as.factor(pt.400.chr7.frame$CopyNumber)

# PATIENT 348 ----
## chr2

this.libs <-  dplyr::filter(sample.data.all, passed.qc, patient=='348')$plate_sample

this.libs <- intersect(this.libs, colnames(cnv.frame))

pt.348_20.chr.2 <- t(cnv.frame[which(grepl('^2_',rownames(cnv.frame))),
                               this.libs])

pt.348_20.chr.2 <- pt.348_20.chr.2[order(sample.data.all[rownames(pt.348_20.chr.2),'Class'],
                                         sample.data.all[rownames(pt.348_20.chr.2),'group']),]


pt.348_20.gr <- cnv.gr[intersect(names(cnv.gr), sample.data.all %>%
                                     dplyr::filter(passed.qc, 
                                                   patient=='348', 
                                                   cycle.phase=='G0') %>%
                                     dplyr::pull(plate_sample))]

pt.348_20.gr.2 <- lapply(pt.348_20.gr, function(x){as.data.frame(x[which(as.character(seqnames(x))=='chr2')])})

pt.348_20.frame <- bind_rows(pt.348_20.gr.2, .id='plate_sample')

pt.348_20.frame <- merge(pt.348_20.frame, sample.data.all)
pt.348_20.frame$CopyNumber <- as.factor(pt.348_20.frame$CopyNumber)

## CHR 17

this.libs <-  dplyr::filter(sample.data.all, passed.qc, patient=='348')$plate_sample

this.libs <- intersect(this.libs, colnames(cnv.frame))

pt.348_20.chr.17 <- t(cnv.frame[which(grepl('^17_',rownames(cnv.frame))),
                                this.libs])

pt.348_20.chr.17 <- pt.348_20.chr.17[order(sample.data.all[rownames(pt.348_20.chr.17),'Class'],
                                           sample.data.all[rownames(pt.348_20.chr.17),'group']),]


pt.348_20.gr <- cnv.gr[intersect(names(cnv.gr), sample.data.all %>%
                                     dplyr::filter(passed.qc, 
                                                   patient=='348', 
                                                   cycle.phase=='G0') %>%
                                     dplyr::pull(plate_sample))]

pt.348_20.gr.17 <- lapply(pt.348_20.gr, function(x){as.data.frame(x[which(as.character(seqnames(x))=='chr17')])})

pt.348_20.frame <- bind_rows(pt.348_20.gr.17, .id='plate_sample')

pt.348_20.frame <- merge(pt.348_20.frame, sample.data.all)
pt.348_20.frame$CopyNumber <- as.factor(pt.348_20.frame$CopyNumber)

# SUMMARIZE ALL ----
cnv.gr.all <- lapply(cnv.gr, as.data.frame)

all.cnv.frame <- bind_rows(cnv.gr.all, .id='plate_sample')
all.cnv.frame <- merge(all.cnv.frame, 
                       dplyr::filter(sample.data.all[,c('plate_sample', 'patient', 'sample', 'group', 'cycle.phase', 'passed.qc', 'Class')], group=='single-cell', cycle.phase=='G0', passed.qc, !patient %in% c('278', '446')) )
all.cnv.frame$CopyNumber[which(as.numeric(as.character(all.cnv.frame$CopyNumber)) > 5) ] <- '5'
all.cnv.frame$CopyNumber <- as.factor(all.cnv.frame$CopyNumber)

all.cnv.frame$seqnames <- sub('chr', '', all.cnv.frame$seqnames)
all.cnv.frame$seqnames <- factor(all.cnv.frame$seqnames, levels=c(as.character(1:22), 'X', 'Y'))

all.cnv.frame <- all.cnv.frame %>%
    filter(cycle.phase == "G0") %>%
    mutate(patient = factor(patient, levels = c("341", "371", "388", "400", "437", "449", "348"))) 

saveRDS(all.cnv.frame, here("results", "cnvs.rds"))
