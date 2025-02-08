# INPUT: `resources/supp_data/liftover/hg38ToHg19.over.chain`, `resources/supp_data/Farlik_2016_meth`
# OUTPUT: `results/farik/GSM*.cons.Cmethyl.CpG.txt`

library(rtracklayer)
library(data.table)
library(here)
library(tidyverse)

hg38.19.chain <- import.chain(here('resources', 'supp_data', 'liftover', 'hg38ToHg19.over.chain'))

convertOneMeth <- function(meth.file, output.dir)
{
	message('Converting ', basename(meth.file))
	meth.tab <- fread(meth.file)

	# Convert to NovoMethyl format:
	colnames(meth.tab) <- c('chr', 'start', 'meth_reads', 'all_reads')
	meth.tab$end <- meth.tab$start + 1
	meth.tab$perc_meth <- meth.tab$meth_reads / meth.tab$all_reads * 100
	meth.tab <- meth.tab[,c('chr', 'start', 'end', 'perc_meth', 'meth_reads', 'all_reads')]

	# Convert to GRange:
	meth.gr <- makeGRangesFromDataFrame(meth.tab, keep.extra.columns = TRUE)

	#Liftover:
	meth.hg19 <- liftOver(meth.gr, hg38.19.chain)
	meth.hg19 <- unlist(meth.hg19)
	meth.tab.hg19 <- as.data.frame(meth.hg19)
	meth.tab.hg19 <- meth.tab.hg19[,c('seqnames', 'start', 'end', 'perc_meth', 'meth_reads', 'all_reads')]

	#Write out:

	out.file <- sub('.txt.gz', '', basename(meth.file))

	write.table(meth.tab.hg19,
			  col.names=FALSE,
			  row.names=FALSE,
			  quote=FALSE,
			  sep='\t',
			  file=file.path(output.dir, paste0(out.file, '.cons.Cmethyl.CpG.txt')))
}

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0) {
	message('Usage: convertFarlikToNovomethyl.R path_to_meth_files output_dir')
} else {
	path.to.meth <- args[1]
	output.dir <- args[2]
  
	if(!file.exists(output.dir)){system(str_glue("mkdir {output.dir}"))}
	  
	meth.files <- dir(path.to.meth, full.names=TRUE)

	for (meth.file in meth.files)
	{
		convertOneMeth(meth.file, output.dir)
	}

}
