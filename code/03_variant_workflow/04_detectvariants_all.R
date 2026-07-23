suppressMessages({
  library(dplyr)
  library(yaml)
})

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")
work.path = paste0(manifest$work, "project_s51n/DWMCT/Analysis/braun_et_al/")

seurat.pd = openxlsx::read.xlsx("data/pub_S51A/single_cell_seurat_metadata.xlsx")
seurat.pd$SAMPLE_NAME[seurat.pd$SAMPLE_NAME == "P065-PB"] = "P065-PB-D365"


bam_files = function(p = NULL, name = NULL){

  bam.files = list.files(
    path = p, full.names = T, recursive = T, pattern = "_hg38",
  )
  bam.files = bam.files[!grepl("bck", bam.files)]
  bam.files = bam.files[grepl("\\.sorted.bam$", bam.files)]
  bam.files
}

files = bam_files(c(
    paste0(manifest$work, "project_s51n/DWMCT/Analysis/braun_et_al/"),
    paste0(manifest$work, "project_s51n/DWMCT/Analysis/aleman_et_al/"),
    paste0(manifest$work, "project_s51n/DWMCT/Analysis/hosoya_et_al/"),
    paste0(manifest$work, "project_s51n/DWMCT/Analysis/ho_et_al/"),
    paste0(manifest$work, "project_s51n/DWMCT/Analysis/demultiplexOnSeurat_oekelen/"),
    paste0(manifest$work, "project_s51n/DWMCT/Analysis/demultiplexOnSeurat_merz/")

)) %>% data.frame()
colnames(files) = "BAM_FILE"
files$SAMPLE = gsub("\\.vs.+", "", basename(files$BAM_FILE))
files$SAMPLE = gsub("orig.ident__", "", files$SAMPLE)
files$SAMPLE = gsub("PATIENT_ID_TIMEPOINT__", "", files$SAMPLE)
files$STUDY = seurat.pd$STUDY[match(files$SAMPLE, seurat.pd$SAMPLE_NAME)]
files$SOURCE = seurat.pd$SOURCE[match(files$SAMPLE, seurat.pd$SAMPLE_NAME)]
files$GROUP = seurat.pd$GROUP[match(files$SAMPLE, seurat.pd$SAMPLE_NAME)]
files$SAMPLE = gsub("P065-PB-D365", "Patient065_Very_Very_Late", files$SAMPLE)
files$STUDY[is.na(files$STUDY)] = "Merz"
files$SOURCE[is.na(files$SOURCE)] = "PB"
files$GROUP[is.na(files$GROUP)] = "Control"
files$GROUP[files$SAMPLE == "Patient065_Late"] = "L-IEC"
files$GROUP[files$SAMPLE == "Patient065_Very_Late"] = "L-IEC"
files$GROUP[files$SAMPLE == "Patient065_Very_Very_Late"] = "L-IEC"

files = files[!files$SAMPLE %in% c("BLD_HD_UnStim", "BLD_HD_Stim", "AphNB"), ]

script = paste0(manifest$homes, "/code/03_variant_workflow/detectVariants.py")

out.dir = paste0(
  manifest$work, "project_s51n/DWMCT/Analysis/detectVariants_out/",
  "detectVariants_out_2026_07_21/"
)
out.dir = paste0(
  manifest$work, "project_s51n/DWMCT/Analysis/detectVariants_out/",
  "detectVariants_out_noPhred_noMAPK_2026_07_21/"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# module load Python/3.13.5-GCCcore-14.3.0
# module load R/4.5.2-gfbf-2025b-bare
# module load SAMtools/1.22.1-GCC-14.3.0

parallel::mclapply(1:nrow(files),function(x){

  bam = files[x, ]
  ref = paste0(dirname(bam$BAM_FILE), "/ciltacel_mined__hg38.fasta")
  o = paste0(
    bam$STUDY, "_xx_", bam$SAMPLE, "_xx_", bam$SOURCE, "_xx_", bam$GROUP
  )
  print(o)
  cmd = paste0(
    "python3 ", script,
    " --bam ", bam,
    " --reference ", ref,
    " --out ", paste0(out.dir, o),
    " --write-region-bams",
    " --mut G222A,G314A,G476T",
    " --min-mapq 0",
    " --min-base-quality 0",
    " --base-alpha 0.01"
  )
  system(cmd)
  ""
}, mc.cores = 40)
