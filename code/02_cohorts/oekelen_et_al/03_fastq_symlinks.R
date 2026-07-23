suppressMessages({
  library(dplyr)
  library(yaml)
})

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182527
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# scRNA-Seq: CellRanger count
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

fastqs = paste0(manifest$oekelen$work, "fastq/")
sratable = read.csv("code/02_cohorts/oekelen_et_al/SraRunTable.csv")
out.dir = paste0(manifest$oekelen$work, "fastq_symlinks/")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# sample paths
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fastq.files = list.files(path = fastqs, full.names = T, recursive = T)
df = data.frame(
  SAMPLE = basename(fastq.files),
  PATH = fastq.files
)
df = df[!grepl("md5$", df$SAMPLE), ]
df = df[!grepl("_3.fastq.gz", df$SAMPLE), ]

new.lbls = c(
  "SRR15541180_1.fastq.gz" = "SAPA15_CARTBT_S3_L001_R1_001.fastq.gz",
  "SRR15541180_2.fastq.gz" = "SAPA15_CARTBT_S3_L001_R2_001.fastq.gz",
  "SRR15541181_1.fastq.gz" = "SAPA15_CARTBT_S3_L002_R1_001.fastq.gz",
  "SRR15541181_2.fastq.gz" = "SAPA15_CARTBT_S3_L002_R2_001.fastq.gz"
)

df$SAMPLE_NAMES = new.lbls[match(df$SAMPLE, names(new.lbls))]

# df$SRR = gsub("_.+", "", df$SAMPLE)
# df$GSM = sratable$Sample.Name[match(df$SRR, sratable$Run)]
# df$GSM = paste0(df$GSM, gsub(".+_", "_", df$SAMPLE))
# df$GSM = gsub("_1", "_S1_R1_001", df$GSM)
# df$GSM = gsub("_2", "_S1_R2_001", df$GSM)

df$OUT = paste0(out.dir, df$SAMPLE_NAMES)

# [Sample Name]_S1_L00[Lane Number] _[Read Type]_001.fastq.gz
# T1_Lane_1_R_S1_L001_R1_001.fastq.gz

path = paste0(out.dir, "make_symlink.sh")
file.create(path)
write("#!/bin/bash", file = path, append = TRUE)
write("", file = path, append = TRUE)
write(paste0("ln -s ", df$PATH, " ", df$OUT), file = path, append = TRUE)

cmd = paste0("bash ", path)
system(cmd)
