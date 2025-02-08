# Cohort 2: Single-cell Multiomic Analysis

![](../../Cohort2_Workflow.png)

`cellranger-arc` is run on the ATAC and RNA libraries, and then a custom `kallisto`, `bustools`, and `kite` pipeline is run on the ADT libraries, using the same methods as Mimitou et. al, 2021. 10X data is loaded into R as described by Seurat vignettes.

## QC, Predict cell cycle, Predict cell type with _Seurat_

Figures: 
- Fig. 1F: RNA clusters by HSPC celltype
- Fig. 3M: HOXB4 expression in RNA clusters
- Fig. S6C: RNA clusters by all celltypes
- Fig. S6D: HSPC cell type counts

### `analyze_cohort2.R`

Cells are filtered for quality control based on 10X metrics and Seurat recommended metrics. Clustering and downstream analysis are conducted as recommended by Seurat. Cells are classified by RNA using Azimuth. This script is run individually for each sample in Cohort 2.

## Integrate across samples with _Harmony_

### `integrate_cohort2_samples.R`

Samples are merged and re-clustered/re-analyzed, then integrated using Harmony.

## Call differentially accessible regions (DARs) with _Seurat_

Figures: 
- Fig. 2B: DAR genomic volcano plot

### `call_cohort2_dars.R`

Differentially accessible regions (DARs) were called using Seurat's `FindMarkers` function.

## Annotate DARs with genomic regions using _annotatr_

### `annotate_cohort2_dars_genomic.R`

DARs called were annotated for genomic regions using `annotatr`. hg38 enhancers downloaded from Ensembl (<http://www.ensembl.org/biomart/martresults/55?file=martquery_1211222926_911.txt.gz>).

## Annotate DARs with TFBS using _annotatr_

Figures: 
- Fig. 2E: DAR TFBS volcano plot
- Fig. S14: DAR TFBS (full 64) volcano plot

### `annotate_cohort2_DARs_TFBS.R`

DARs called were annotated for transcription factor binding sites using LOLA.


## Annotate DARs with CTCF and ZNF143 cobinding sites and promoter, enhancer, and TAD regions

Figures: 
- Fig. 3B: TAD boundaries and promoters/enhancers within hyperaccessible DARs

### `tad_promoter_enhancer_cohort2_analysis.R`

DARs were annotated with promoter/enhancer regions and TAD boundaries. hg19 enhancers were downloaded from <http://www.enhanceratlas.org/data/AllEPs/hs/CD34+_EP.txt>.

## Call differentially expressed genes (DEGs) with _Seurat_

Figures: 
- Fig. 3K: DEG Volcano

### `call_cohort2_degs.R`

Differentially expressed genes were called using Seurat's `FindMarkers` function.

## Call differential gene accessibility (DGA) with _Seurat_

Figures: 
- Fig. 3I: DGA Volcano

### `call_cohort2_dgas.R`

Differential gene accessibility called using Seurat's `FindMarkers` function.

## _GREAT_ analysis

Figures:
- Fig. 3G: Hyperaccessible DAR GREAT scatterplot
- Fig. S21C: Hyperacessible DAR GREAT barplot
- Fig. S21D: Hypoaccessible DAR GREAT barplot
- Fig. S21F: Hypoaccessible DAR GREAT scatterplot

### `prepare_cohort2_for_great.R`

DARs were filtered for adjusted p-value \< 0.001 and `avg_log2FC` \> 0.5 (hyeraccessible) as there was no threshold set during the calling process.

## Validate accessibility with ADT data

Figures: 
- Fig. S6B: ADT Validation

### `integrate_cohort2_adt_samples.R`

Integrate samples with ADT. 

## Scramble cells across patients and response as negative control

Figures:
- Fig. S10C: DAR Scramble (genomic)
- Fig. 14: DAR Scramble (TFBS)

### `scramble_cohort2_dars.R`

Cells from all patients were randomly split into 12 equal groups to simulate the patients, and then these 12 grouops are randomly split into two groups for calling. The two-group split is iterated over 10 times.

### `annotate_cohort2_scrambled_dars_tfbs.R`

The scrambled DARs are then annotated for selected 64 transcription factor binding site experiments from LOLA.

### `annotate_cohort2_scrambled_dars_genomic.R`

The scrambled DARs are then annotated function genomic features. 

## Infer gene regulatory networks with _GRNBoost2_

Figures: 
- Fig. S15A: TF Network

### `run_SCENIC.R`

Transcription factor network analysis using SCENIC and GRNBoost2. Uses helper script `run_GRNBoost2.py`.

## Annotate peaks with TFBS using _annotatr_

Figures: 
- Fig. S15B: TF UMAP

### `annotate_cohort2_peaks_tfbs.R`

Transcription factor UMAP using transcription factor binding sites of interest within peaks. 

## Call copy number gain or loss using ATAC data

Figures:
- Fig. S2: ATAC Copy Number

### `run_epiAneufinder.R`

epiAneufinder provides a wrapper script for running on the commandline: https://github.com/colomemaria/epiAneufinder/blob/main/epiAneufinder_wrapper.R

This script was minimally modified to take additional arguments and other optimizations.

### `call_karyotype.R`

Determine a _in silico_ cytogenetics karyotype call based on the percentage of cells indicating gain or loss per bin, and the percentage of bins indicating gain or loss per chromosome arm.

### `remove_recurrent_calls.R`

Remove recurrent calls in the cohort as they could be known copy number polymorphisms or artifacts.
