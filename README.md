# Supplementary code accompanying manuscript: A recurrent synthetic receptor point mutation at Ser51Asn in CAR+ T cell lymphoproliferation and enterocolitis after Ciltacabtagene autoleucel

This repository contains the code used to produce the results of Mitchell Kluesner, Michael Rade et al. A recurrent synthetic receptor point mutation at Ser51Asn in CAR+ T cell lymphoproliferation and enterocolitis after Ciltacabtagene autoleucel, DOI XXXXXXXXXXXX

## Singularity

Most scripts for the single-cell analysis were developed in a Singularity image with Rstudio server. [See README](singularity/) in `./singularity/`. 

## Zenodo repository

https://doi.org/10.5281/zenodo.21457937

Includes:
- Singularity images for R/RStudio server
- All R packages used for single-cell analysis
- Seurat object used for publication
- Outout of the Variant detection workflow

## Reproduction

``` sh
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Custom reference for cellranger (hg38 + cilta)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bash code/01_references/build_gex_ref_cilta_mined.sh

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# single-cell data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Rscript code/02_cohorts/hosoya_et_al/01-prefetch-SRR.sh
Rscript code/02_cohorts/hosoya_et_al/02-fastq-dump.sh
Rscript code/02_cohorts/hosoya_et_al/03_fastq_symlinks.R
Rscript code/02_cohorts/hosoya_et_al/04_cellranger.R
Rscript code/02_cohorts/hosoya_et_al/05_cellranger_to_seurat.R
Rscript code/02_cohorts/hosoya_et_al/06_qc_anno.R

Rscript code/02_cohorts/ho_et_al/01-prefetch-SRR.sh
Rscript code/02_cohorts/ho_et_al/02-fastq-dump.sh
Rscript code/02_cohorts/ho_et_al/03_aws_download.sh
Rscript code/02_cohorts/ho_et_al/04_cellranger.R
Rscript code/02_cohorts/ho_et_al/05_cellranger_to_seurat.R
Rscript code/02_cohorts/ho_et_al/06_qc_anno.R


Rscript code/02_cohorts/aleman_et_al/01_cellranger.R.R
Rscript code/02_cohorts/aleman_et_al/02_cellranger_to_seurat.R.R
Rscript code/02_cohorts/aleman_et_al/03_qc_anno.R

Rscript code/02_cohorts/braun_et_al/03_qc_anno.R

Rscript code/02_cohorts/merz_et_al/02_qc_anno.R

Rscript code/02_cohorts/oekelen_et_al/01-prefetch-SRR.sh
Rscript code/02_cohorts/oekelen_et_al/02-fastq-dump.sh
Rscript code/02_cohorts/oekelen_et_al/03_fastq_symlinks.R
Rscript code/02_cohorts/oekelen_et_al/04_cellranger.R
Rscript code/02_cohorts/oekelen_et_al/05_cellranger_to_seurat.R
Rscript code/02_cohorts/oekelen_et_al/06_qc_anno.R

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Variant detection
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Rscript code/03_variant_workflow/01_sample_filter_merz.R
Rscript code/03_variant_workflow/02_demux_bams_merz.R
Rscript code/03_variant_workflow/02_demux_bams_oekelen.R
Rscript code/03_variant_workflow/03_align_aleman.R
Rscript code/03_variant_workflow/03_align_braun.R
Rscript code/03_variant_workflow/03_align_ho.R
Rscript code/03_variant_workflow/03_align_hosoya.R
Rscript code/03_variant_workflow/03_align_merz.R
Rscript code/03_variant_workflow/03_align_oekelen.R
Rscript code/03_variant_workflow/04_detectvariants_all.R

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main Figures
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Rscript code/04_publication/figure_2_main.R
Rscript code/04_publication/figure_2_supp_1.R

```
