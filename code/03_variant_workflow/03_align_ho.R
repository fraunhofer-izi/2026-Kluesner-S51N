suppressMessages({
  library(dplyr)
  library(yaml)
})

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")
work.path = paste0(manifest$work, "project_s51n/DWMCT/Analysis/ho_et_al/")

# cr = paste0(manifest$work, "project_s51n/data/single_cell/ho_et_al/cellranger/")
# cr.files =  list.files(path = cr, full.names = T, recursive = T)
# cr.bam = cr.files[grepl("sample_alignments.bam", cr.files)]
# for (i in cr.bam) {
#   print(basename(i))
#   t = paste0(
#     work.path,
#     gsub("multi_", "", unique(basename(dirname(i)))),
#     gsub("sample_alignments", "", basename(i))
#   )
#   cmd = paste0("cp ", i, " ", t)
#   print(cmd)
#   system(cmd)
#
# }

bam.files = list.files(path = work.path, full.names = T, recursive = F)
bam.files = bam.files[grepl("bam$", bam.files)]
names(bam.files) = gsub(".bam", "", gsub(".+__", "", basename(bam.files)))

script = paste0(manifest$homes, "/code/03_variant_workflow/alignReference.sh")
ref.fasta = paste0(
  manifest$homes , "DWMCT/collaboration/reference_genomes/ciltacel_mined.fasta "
)
ref.hg38 = paste0(manifest$work , "references/hg38.fa")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# submit to Slurm
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for (sample in names(bam.files)) {
  job_id = sample
  print(job_id)
  bam = bam.files[[job_id]]

  cmd = paste0("bash ", script, " ", bam, " -cat ", ref.fasta, ref.hg38)
  system(cmd)
}
