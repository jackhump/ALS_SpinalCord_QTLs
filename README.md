## ALS Spinal Cord differential expression and molecular QTLs

Jack Humphrey 2019-2022

![deg figure](https://github.com/jackhump/ALS_SpinalCord_QTLs/raw/master/deg_figure.png)

## Manuscript

[Link to preprint on medrxiv](https://www.medrxiv.org/content/10.1101/2021.08.31.21262682v1)

## Shiny App

[Launch Shiny app](https://jackhumphrey.shinyapps.io/als_spinal_cord_browser/)

[source code](https://github.com/jackhump/ALS_SpinalCord_QTLs/blob/master/als_spinal_cord_browser/app.R)

## Data Availability

[Raw RNA-seq data - Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137810)

[Processed gene expression count matrices with de-identified metadata - Zenodo](https://zenodo.org/record/6385747) 

[An HTML vignette on downloading the count data from Zenodo and performing differential expression](https://jackhump.github.io/ALS_SpinalCord_QTLs/html/DE_Vignette.html)

Whole genome sequencing data - dbGAP (in progress). 

[Full summary statistics for expression and splicing QTLs - Zenodo](https://zenodo.org/record/5248758).

[All TWAS weight files - Zenodo](https://zenodo.org/record/5256613)


## Rmarkdown workbooks written for the project.

Most workbooks are presented as raw .Rmd (in `scripts/` and knitted HTML files (in `html/`).

### Data QC and setup

[RNA-seq QC](https://jackhump.github.io/ALS_SpinalCord_QTLs/html/ALS_covariates_setup.html)

WGS QC pipeline hosted [here](https://github.com/jackhump/WGS-QC-Pipeline)

[Building sets of marker genes](https://github.com/jackhump/ALS_SpinalCord_QTLs/blob/master/scripts/ALS_marker_lists.Rmd)
 
### Differential expression

[Running differential expression](https://jackhump.github.io/ALS_SpinalCord_QTLs/html/ALS_differential_expression.html)

[Volcano plots, GSEA, marker gene overlaps](https://jackhump.github.io/ALS_SpinalCord_QTLs/html/ALS_deg_functions.html)

[Cell-type deconvolution](https://jackhump.github.io/ALS_SpinalCord_QTLs/html/ALS_deconvolution.html)

### Co-expression networks

[Building co-expression networks, associations with traits and GO](https://jackhump.github.io/ALS_SpinalCord_QTLs/html/ALS_WGCNA.html)

### Quantitative Trait Loci (QTLs)

[QTL results per region and sharing with GTEx](https://jackhump.github.io/ALS_SpinalCord_QTLs/html/ALS_QTL_results.html)

[Colocalisation with ALS GWAS](https://jackhump.github.io/ALS_SpinalCord_QTLs/html/ALS_COLOC.html)

[Transcriptome-wide association study - TWAS](https://jackhump.github.io/ALS_SpinalCord_QTLs/html/ALS_TWAS.html)

### Gene prioritisation

[C9orf72 and ATXN3](https://jackhump.github.io/ALS_SpinalCord_QTLs/html/ATXN3_C9orf72.html)

[Transcript plots](https://jackhump.github.io/ALS_SpinalCord_QTLs/html/ALS_locus_plots.html)

[Localising genes to cell-types](https://jackhump.github.io/ALS_SpinalCord_QTLs/html/ALS_genes_celltypes.html)
