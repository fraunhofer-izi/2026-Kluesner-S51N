#!/bin/bash

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310992

module load SRA-Toolkit

parallel -j 30 fastq-dump --gzip --split-files -O /mnt/ribolution/user_worktmp/michael.rade/work/2026-Kluesner-FredHutch/project_s51n/data/single_cell/hosoya_et_al/fastq {} ::: $(cat SRA-manifest)
## A nice explanation of other fastq-dump options are provided by Rob Edward's group: https://edwards.sdsu.edu/research/fastq-dump/
