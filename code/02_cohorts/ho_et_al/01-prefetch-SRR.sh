#!/bin/bash

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309403

module load SRA-Toolkit

## "custom.kfg" is a copy from "/mnt/fhgfs_ribdata/tools/sratoolkit/2.5.4-1/bin/ncbi/default.kfg"
## change in line 38 the path where files should be downloaded
export VDB_CONFIG=custom.kfg

## test: vdb-config /repository/user/main

parallel -j 10 prefetch --max-size 100G {} ::: $(cat SRR_Acc_List.txt)

## The -j 1 specifies the number of threads to use. Using 1 limits to downloading one file at a time
## (simultaneous downloads may be faster, depending on your computer and network).

## The ::: $(cat SRR_Acc_List.txt) passes the contents of SRR_Acc_List.txt as arguments to the parallel
