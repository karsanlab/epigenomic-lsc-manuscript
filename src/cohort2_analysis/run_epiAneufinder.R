# INPUT:
# results/{sample}_atac_fragments.tsv.gz
# results/{sample}_atac_fragments.tsv.gz.tbi
# resources/blacklist_hg38_signac.bed 

# OUTPUT:
# results/{sample}/epiAneufinder_results/results_table.tsv

# This script was downloaded from:
# https://github.com/colomemaria/epiAneufinder/blob/main/epiAneufinder_wrapper.R
# https://github.com/colomemaria/epiAneufinder/commit/34e5bc7e8a8ce04fae446b85fbc816e0b1b9b7fc

# Modified:
# - to take the sample name as the first argument
# - to take fragments file as second argument
# - to use a custom hg38 BSgenome (U2AF1 fix)
# - to use sample name in karyogram title
# - for use with SLURM resources

args <- commandArgs(trailingOnly=TRUE)
print("Running script for epiAneufinder")

library(epiAneufinder)

#If using a genome version different from the BSgenome.Hsapiens.UCSC.hg38 the correspnding R library should be installed and loaded
#Uncomment the following line and add the corresponding R library
#library(genome-version)
library(BSgenome.Hsapiens.NCBI.b38)

library(here)
sample <- args[1]

#Location of the input data. Can be a single fragments tsv file or a folder containing multiple bam files
input <- args[2] #path to folder containing BAM files /or/ path+file_name for a TSV fragments file. If you want to test the package, sample data can be downloaded from the sample_data folder.

# Location of the output directory. If the directory does not exist it will be created
outdir <- here("results", sample)

#BSgenome to use for the analysis. The genome should be already installed in R. In this example we use the UCSC hg38 genome from Bioconductor
genome <- "BSgenome.Hsapiens.UCSC.hg38"

#Chromosomes to be excluded from the analysis. The chromosome names should follow the naming of the genome version used. In this example we have UCSC chromosome names
#Default is NULL
exclude <- c('chrX','chrY','chrM')

#Bed file with the blacklisted regions of the genome. This file is genome-version specific and it should be downloaded by the user
# ENCODE Unified GRCh38 Exclusion List
# https://stuartlab.org/signac/reference/blacklist_hg38_unified
blacklist <- here("resources", "blacklist_hg38_signac.bed") #Path and file name of the blacklisted regions in bed format. If you use hg38 genome then the blacklisted regions can be found in the sample_data folder.

#Window size for partitioning the genome. Smaller window sizes will result in longer running times. Default is 1e5
windowSize <- 1e5

#Parameter to instruct epiAneufinder to resume from a previous run. Can be set to either True or False
#If certain parameters change, for example minsize, resuming may end in error messages. In such a case change the parameter to False 
#Default in False
reuse.existing=TRUE

#Upper quantile thrshold. Default is 0.9
uq=0.9

#Lower quantile threshold. Default is at 0.1
lq=0.1

#Title for the karyotype plot. Default is NULL
title_karyo=stringr::str_glue("{sample} Karyoplot")

#Number of cores to use for the analysis. Default is 4
if(stringr::str_starts(Sys.getenv("HOSTNAME"), "n")){
    if(Sys.getenv("SLURM_CPUS_PER_TASK") == "") {
        # use all if not set
        ncores=Sys.getenv("SLURM_CPUS_ON_NODE")
    } else {
        ncores=Sys.getenv("SLURM_CPUS_PER_TASK")
    }
} else {
    ncores=48
}

#Minimum number of fragments for a cell to be included in the analysis. This parameter is only for fragnment files. Default is 20000
minFrags = 20000

#Threshold for filtering bins if the ratio of cells with zero reads is higher than the threshold. Setting it to 0 deactivates the filter. Default is 0.85 
threshold_blacklist_bins=0.85

#Parameter on how many breakpoins to use for the CNV calculation. Default is 1, all breakpoints. If higher than one, the algorithm will calculate every n breakpoints
#Setting it to higher than 1 speeds the process with lower resolution as a result
minsize=1

#Number of segments per chromosomes (2^k). Default value is 3
k=4

#Number of consecutive bins to constitute a CNV
minsizeCNV=0

dir.create(outdir)
epiAneufinder::epiAneufinder(input=input, 
                             outdir=outdir, 
                             blacklist=blacklist, 
                             windowSize=windowSize, 
                             genome=genome, 
                             exclude=exclude, 
                             reuse.existing=reuse.existing, 
                             uq=uq, lq=lq, 
                             title_karyo=title_karyo, 
                             ncores=ncores,
                             minFrags=minFrags,
                             minsize=minsize,
                             k=k,
                             threshold_blacklist_bins=threshold_blacklist_bins,
                             minsizeCNV=minsizeCNV)
