# Downloads

## Farlik 

The Farlik methylation data is a [BLUEPRINT](https://medical-epigenomics.org/papers/BLUEPRINT_methylomes/) dataset and can be downloaded from GEO accession [GSE87196](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?token=snonoyyonvqrhkp&acc=GSE87196)

- https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE87196&format=file

The Farlik classifier can be downloaded [BLUEPRINT](https://www.medical-epigenomics.org/papers/BLUEPRINT_methylomes/#prediction):

- https://www.medical-epigenomics.org/papers/BLUEPRINT_methylomes/data/cell_type_predictor_code.zip
- https://www.medical-epigenomics.org/papers/BLUEPRINT_methylomes/data/cell_type_predictor_input.zip

## Liftover

The hg38 to hg19 liftover chain can be downloaded here: 

- https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz

The hg18 to hg19 liftover chain can be downloaded here:

- https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz

## Normal 

The Vyas dataset used for normal analysis can be downloaded from the ENA.
- fastqs: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5456/
- reference: ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

## LOLA

The LOLACore database can be downloaded here:

- http://big.databio.org/regiondb/LOLACore_170206.tgz

In addition to the database above, we also used 8 custom experiments:

```
filename    cellType    treatment   antibody    description
GSM2274676_Jurkat_RUNX1_peaks.bed   Jurkat  None    RUNX1   ChIP Jurkat RUNX1
GSM2274677_Jurkat_TET2_peaks.bed    Jurkat  None    TET2    ChIP Jurkat TET2
GSM2687534_PBMC_RUNX1_peaks.bed PBMC    None    RUNX1   ChIP PBMC RUNX1
PMID_23353889_HEK293T_OGT_peaks.bed HEK293T None    OGT ChiP HEK293T OGT
PMID_23353889_HEK293T_TET2_peaks.bed    HEK293T None    TET2    ChiP HEK293T TET2
PMID_23353889_HEK293T_TET3_peaks.bed    HEK293T None    TET3    ChiP HEK293T TET3
PMID_23727019_hESC_ELK1_peaks_hg19.bed  hESC    None    ELK1    ChIP hESC ELK1
PMID_23727019_hESC_ERK2_peaks_hg19.bed  hESC    None    ERK2    ChIP hESC ERK2
```

These can be downloaded here:

- GSM2274676: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2274676&format=file&file=GSM2274676%5FJurkat%5FRUNX1%5Fpeaks%2Ebed%2Egz
- GSM2274677: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2274677&format=file&file=GSM2274677%5FJurkat%5FTET2%5Fpeaks%2Ebed%2Egz
- GSM2687534: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2687534&format=file&file=GSM2687534%5FPBMC%5FRUNX1%5Fpeaks%2Ebed%2Egz
- PMID_23353889: https://pmc.ncbi.nlm.nih.gov/articles/instance/3590984/bin/emboj2012357s3.xls
- PMID_23727019 - ELK1: https://ars.els-cdn.com/content/image/1-s2.0-S1097276513003353-mmc3.xlsx
- PMID_23727019 - ERK2: https://ars.els-cdn.com/content/image/1-s2.0-S1097276513003353-mmc2.xlsx

## Replication Timing

Repli-seq data used in our CNV analysis was from Hansen et. al, 2010. 

- https://www.pnas.org/doi/10.1073/pnas.0912402107

The replication timing dataset can be downloaded here:

- http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562G1PctSignalRep1.bigWig
- http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562G2PctSignalRep1.bigWig
- http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562S1PctSignalRep1.bigWig
- http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562S2PctSignalRep1.bigWig
- http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562S3PctSignalRep1.bigWig
- http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562S4PctSignalRep1.bigWig

The `plastic_masked.txt` file is derived from Hansen et. al, 2010's supplementary table: 

- https://www.pnas.org/doi/suppl/10.1073/pnas.0912402107/suppl_file/st02.doc

## Human TFs

The human TF gene names were downloaded from:

- https://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt

## hg38 Blacklist

The hg38 blacklist used for `epiAneufinder` was the one packaged with the Signac R package: 
https://stuartlab.org/signac/reference/blacklist_hg38_unified

It was derived from the [ENCODE Unified GRCh38 Exclusion List](https://www.encodeproject.org/files/ENCFF356LFX/).

It was exported to a bed using `export.bed()` for use in the wrapper script. 
