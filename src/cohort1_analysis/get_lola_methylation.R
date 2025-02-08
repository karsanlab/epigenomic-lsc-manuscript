# INPUT: 
# `data/scWGBS/*.cons.Cmethyl.CpG.txt`

# OUTPUT: 
# `results/*.cons.Cmethyl_tfbs.csv`

library(rtracklayer)
library(reshape)
library(GenomicRanges)
library(data.table)
library(LOLA)
require(simpleCache)
library(pbapply)
library(here)
library(tidyverse)

#' readOneLibrary
#' 
#' Function to import the bed like cpg file from a single cell and rename the
#' columns, returns a dataframe
readOneLibrary <- function(cpg.file.name){
	cpg.bed <- fread(cpg.file.name)
	colnames(cpg.bed) <- c('chr', 'start', 'end', 'perc_meth', 'meth_reads', 'all_reads')
	cpg.bed
}

#' convertToGrange
#' 
#' Filters a data frame returned from readOneLibrary to methylation percentages
#' 0 or 100 and returns the filtered table as a GRanges object
convertToGrange <- function(cpg.table)
{
	filt.table <- cpg.table[which(cpg.table$perc_meth %in% c(0, 100)),]
	filt.gr <- as(filt.table[,1:4, with=FALSE], 'GRanges')
	# filt.gr <- setNames(filt.gr, 1:length(filt.gr))
	return(filt.gr)
}

# Helper functions that allow string arguments for dplyr's data modification functions like arrange, select etc.
# Author: Sebastian Kranz
# Examples are below
#' Modified version of dplyr's filter that uses string arguments
eval.string.dplyr = function(.data, .fun.name, ...) {
	args = list(...)
	args = unlist(args)
	code = paste0(.fun.name,"(.data,", paste0(args, collapse=","), ")")
	df = eval(parse(text=code,srcfile=NULL))
	df
}

#' averageByFeature
#' 
#' Author: ?
#' Average a genome-positional value over genomic features
#' Similar to binnedAverage(), but with more control
#' @param data.grange the data to be averaged (GenomicRanges object)
#' @param feature.grange the features to be averaged over (GenomicRanges)
#' @param value.column name of the column in data.grange to use
#' @param fill.NA if true, will fill in NAs for features not represented in 
#' data.grange; if false, will omit these
averageByFeature <- function(feature.grange, data.grange,  value.column, fill.NA=TRUE) {
	# For cases where chromosomes aren't covered, ensure seqlevels are the same:
	subset.features <- subsetByOverlaps(feature.grange, data.grange, type='any', ignore.strand=TRUE)
	seqlevels(subset.features) <- seqlevels(data.grange)
  
	overlaps <- findOverlaps(subset.features, data.grange, type='any', ignore.strand=TRUE)
	overlaps <- as.data.frame(overlaps)

	# Add in the scores from data.grange:
	overlaps <- cbind(overlaps, data.grange[overlaps$subjectHits]@elementMetadata[,value.column] )

	#Test that there are actually overlapping methylation; some unplaced scaffolds may have no CpGs or something:
	tryCatch({
	  if(nrow(overlaps)==0){
	    scored.features <- feature.grange[0]
	  }	else {
	    head(overlaps)
	    colnames(overlaps)[3] <- value.column
	    #overlaps[[value.column]] <- as.numeric(as.character(overlaps[[value.column]]))  #convert back from factor to numeric as it should be
	    
	    by.feature <- group_by(overlaps, queryHits)
	    
	    #Use dplyr's summarise. Getting dplyr to work with string-named columns is ... painful:
	    by.feature.score <- eval.string.dplyr(by.feature, "summarise",
	                                          paste("mean.score=mean(",value.column,
	                                                "), num.cpgs=length(", value.column, ")", sep=''))
	    by.feature.score <- data.frame(by.feature.score)
	    
	    scored.features <- feature.grange[by.feature.score[,1]]
	    scored.features@elementMetadata[[value.column]] <- by.feature.score$mean.score
	    scored.features@elementMetadata[["num.cpgs"]] <- by.feature.score$num.cpgs
	  }
	  #Add in NAs for missing features:
	  all.features <- feature.grange
	  all.features@elementMetadata[[value.column]] <- as.numeric(rep(NA, length(all.features)))
	  all.features[which(all.features %in% subset.features)]@elementMetadata[[value.column]] <- scored.features@elementMetadata[[value.column]]
	  
	  all.features@elementMetadata[["num.cpgs"]] <- as.integer(rep(NA, length(all.features)))
	  all.features[which(all.features %in% subset.features)]@elementMetadata[["num.cpgs"]] <- scored.features@elementMetadata[["num.cpgs"]]

	})
	
	if(fill.NA){
	  return(all.features)
	}	else {
	  return(scored.features)
	}
}

#' averageAllByFeature
#' 
#' Applies the averageByFeature function across a list of feature GRanges objects.
#' 
averageAllByFeature <- function(feature.grs, cpg.gr, meth.field='perc_meth') {
	message('Averaging across features')
	averaged.meth <- pblapply(feature.grs, averageByFeature, cpg.gr, meth.field)
	getAverageMethByLib <- function(lib.feature.gr){
	  averaged <- sum(lib.feature.gr$perc_meth*lib.feature.gr$num.cpgs, na.rm=TRUE) / sum(lib.feature.gr$num.cpgs, na.rm=TRUE)
	  return(averaged)
	}
	getAverageMethBySubset <- function(subset.list){
	  pbsapply(subset.list, getAverageMethByLib)
	}
	getAverageMethBySubset(averaged.meth)
}

# Load in LOLA DB:
lola.db <- loadRegionDB(here('resources', 'supp_data', 'LOLA', 'LOLACore', 'hg19',
						collections = c('encode_tfbs', 'codex', 'custom', 'cistrome_cistrome')))
subset.grs <- lola.db$regionGRL
lola.tfbs <- lola.db$regionAnno
# this gives not unique names, want to change to using full filename
names(subset.grs) <- lola.tfbs$filename


## Actual code to get methylation by subsets:
args = commandArgs(trailingOnly=TRUE)

if(length(args) == 0){
	stop('Usage: get_lola_methylation.R path_to_bed_file output_filename')
} else {
  path.to.bed <- args[1]
  output.dir <- args[2]
  message('Loading CpG methylation calls')
  this.cpg <- readOneLibrary(path.to.bed)
  message('Converting to GenomicRanges')
  this.cpg.gr <- convertToGrange(this.cpg)
  message('Averaging within subsets')
  this.subsets.meth <- averageAllByFeature(subset.grs, this.cpg.gr)
  names(this.subsets.meth) <- names(subset.grs)
  write.csv(this.subsets.meth, file = file.path(output.dir, sub('.CpG.txt', '_tfbs.csv', basename(path.to.bed))))
}
