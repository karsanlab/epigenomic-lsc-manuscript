# Cell Line Analysis: HiChiP

![](../../CellLine_Workflow.png)

## HiChIP conda environment

The conda environment used in the HiChiP analyses can be recreated using the specifications file provided in the resources folder: 

`conda create --name hichip_conda_environment --file resources/hichip_conda_environment_spec_file.txt`

## Fastq alignment, parsing and sorting

### `fq_to_paired_bam.sh`

Aligns, parses, sorts, deduplicates and splits fastqs, followed by sorting and indexing of the resulting bam file. 

## QC

Figures:

-   Fig. S18A: Fraction of total read pairs that are unduplicated, in cis and > 1000 bp apart.
-   Fig. S18B: Fraction of total reads within 1 kb of peak center
-   Fig. S18C: Representative screenshot of HiChIP coverage enrichment by library in a 2MB window over chr17

### `qc_analysis.sh`

Performs QC analysis following [manufacturers recommendations](https://hichip.readthedocs.io/en/latest/library_qc.html) for assessing the quality of the HiChIP data. 

Makes use of the `getqc.py` function in the HiChIP git hub repository described in the resources folder. 

## Call peaks on control sample using `MACS2`

### `get_sample_peaks.sh`

Call sample peaks using `macs2`.

## Generate `hic` contact matrices

Figures:

-   Fig. 3C: Plot of balanced contact matrices on chromosome 17 at 500KB resolution

### `generate_contact_matrices.sh`

Generate HiC contact matrices at several resolutions. Makes use of the `juicertools.jar` java archive file in the HiChIP git hub repository described in the resources folder. 

## Call A/B compartments

Figures:

-   Fig. S19A: A/B compartments for each chromosome for shControl and shZNF143 samples at 1 Mb resolution.  
-   Fig. S19B: Correlation between the eigenvector value from control and ZNF143 knockdown OCI-AML3 cell line samples at 1 Mb resolution for each chromosome.
-   Fig. S19C: Same as in (b) but showing correlation across all chromosomes.

### `submit_eigenvector_jobs.sh`

Submits jobs to call eigenvectors at 50kb and 1MB resolution.

### `generate_ab_compartments.sh`

Calls A/B compartments as eigenvectors using the `juicertools.jar` java archive file in the HiChIP git hub repository described in the resources folder. 

## Call Topologically Associated Domains (TADs) 

Figures:

-	Fig. S20A: Size distribution of TADs in control and ZNF143 knockdown OCI-AML3 cells. 
- 	Fig. S20B: Venn diagram of exact overlap between TAD domains in shControl vs shZNF143 in OCI-AML3 cells at 50 kb bin size.
-	Fig. 24A: Stacked barplot showing the distribution of changes in TADs per chromosome relative to shControl upon shZNF143 knockdown in the OCI-AML3 cell line
-   Fig. 24B: Barplot displaying the significance of differences between the distribution of overlap categories for each chromosome against the mean distribution of overlap categories across all chromosomes
-   Fig. 24C: HiC plot showing genomic contact frequency in OCI-AML3 cells with and without knockdown of ZNF143

### `generate_tads.R`

Calls Topologically Associated Domains from HiC contact matrices. 

## Merge LOLA CTCF peaks with shControl sample peaks

### `merge_lola_CTCF_and_shControl_sample_peaks.R`

Combine CTCF peaks from the LOLA databases described in resources with the shControl sample peaks.

## Normalize sequencing depth using `bamCoverage` from `deeptools`

### `normalize_sequencing_depth.sh`

Normalize sequencing depth for shControl and shZNF143 sample together as reads per genomic content.

## Calculate enrichment over CTCF and ZNF143 binding sites with `computeMatrix` from `deeptools`

### `generate_CTCF_ZNF143_binding_sites.R`

Merge and overlap CTCF and ZNF143 binding sites to generate lists of CTCF only, ZNF143 only and CTCF and ZNF143 cobinding sites. 

### `calculate_enrichment_over_binding_sites.sh`

Calculate coverage enrichment over CTCF and ZNF143 binding and co-binding sites. 

## Plot heatmap using `plotHeatmap` from `deeptools`

Figures:

-   Fig. 3D: Coverage enrichment (reads per genomic content, RPGC) over ZNF143 and CTCF binding sites in OCI-AML3 cells after knockdown of ZNF143 compared to shControl.

### `plot_heatmap.sh`

Plot RPGC normalized coverage over CTCF, ZNF143 and CTCF-ZNF143 co-binding sites.

## Identify enrichment of significant chromatin loops using `FitHiChIP`

Figures:

-   Fig. 3E: Fraction of significant loops unique to shControl or shZNF143 OCI-AML3 lines overlapping with CTCF and ZNF143 binding sites following HiCHIP.
-   Fig. 5A: Significant promoter loops (FDR < 0.01) in the anterior HOXB cluster in the OCI-AML3 cell line with and without knockdown of ZNF143.
-   Fig. S20C: Total loops detected in OCI-AML3 by HiChIP CTCF pull-down after knockdown of ZNF143 compared to control at 5 kb bin size. Venn diagram indicates degree of exact overlap between shControl and shZNF143 loops. 
-   Fig. S20D: Median loop distance (bp) at 5 KB bin size in OCI-AML3 shControl vs shZNF143 in HiChIP CTCF pull-down after knockdown of ZNF143 compared to shControl at 5 KB bin size.

### `run_fithichip.sh`

Runs the fithichip software to identify significant chromatin loops using the config file provided in the resources folder:

`fithichip_config_file`

