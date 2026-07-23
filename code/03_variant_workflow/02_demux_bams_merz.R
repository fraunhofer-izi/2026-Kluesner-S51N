suppressMessages({
  library(dplyr)
  library(yaml)
})

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# scRNA-Seq: CellRanger count
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
nodes = 1
ntasks = 1
ttime = "66:00:00"
mail = "FAIL"
mem = 110000
cpu = 12

manifest = yaml.load_file("manifest.yaml")
work.path = "/mnt/ribolution/user_worktmp/michael.rade/work/2026-BCMA-CAR-NEW-RUN/cohorts/"
run1 = paste0(work.path, "ukl_batch_1/cellranger")
run2 = paste0(work.path, "ukl_batch_2/cellranger")
run3 = paste0(work.path, "ukl_batch_3/cellranger/cilta/")
run4 = paste0(work.path, "ukl_batch_4/cellranger")
run5 = paste0(work.path, "ukl_batch_5/cellranger")
run.ltp = paste0(
  manifest$rade_fandrei$fandrei, "multi_omics/cellranger_batch1/"
)

cellranger.dirs = list.dirs(
  path = c(run1, run2, run3, run4, run5, run.ltp), full.names = T, recursive = F
)
cellranger.samples = basename(cellranger.dirs)
fltrd.dirs = paste0(
  cellranger.dirs, "/outs/per_sample_outs/", cellranger.samples,
  "/sample_alignments.bam"
)
names(fltrd.dirs) = gsub("multi_", "", cellranger.samples)

se.obj = paste0(manifest$homes, "data/pub_S51A/seurat_t_cilta_merz.Rds")
script = paste0(manifest$homes, "/code/03_variant_workflow/demultiplexOnSeurat.py")

out.dir = paste0(
  manifest$work, "/project_s51n/DWMCT/Analysis/demultiplexOnSeurat_merz/"
)


se = readRDS(se.obj)
fltrd.dirs = fltrd.dirs[
  names(fltrd.dirs) %in% unique(as.character(se$FASTQ_FILE_NAME))
]
print(length(fltrd.dirs))
print(length(unique(se$FASTQ_FILE_NAME)))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# sbatch
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
write_subscript = function(path, bam, job_id){
  file.create(path)

  write("#!/bin/bash", file = path, append = TRUE)
  write("", file = path, append = TRUE)

  write(paste("#SBATCH -J", job_id), file = path, append = TRUE)
  write(paste("#SBATCH --nodes", nodes), file = path, append = TRUE)
  write(paste("#SBATCH --ntasks", ntasks), file = path, append = TRUE)
  write(paste("#SBATCH --time", ttime), file = path, append = TRUE)
  write(paste("#SBATCH --cpus-per-task", cpu), file = path, append = TRUE)
  write(paste("#SBATCH --mem", mem), file = path, append = TRUE)
  write(paste("#SBATCH -e", paste0("demux_", job_id, ".e")), file = path, append = TRUE)
  write(paste("#SBATCH -o", paste0("demux_", job_id, ".o")), file = path, append = TRUE)
  write("#SBATCH --mail-type=END,FAIL", file = path, append = TRUE)

  write("", file = path, append = TRUE)
  write("module load Python/3.13.5-GCCcore-14.3.0", file = path, append = TRUE)
  write("module load R/4.5.2-gfbf-2025b-bare", file = path, append = TRUE)
  write("module load SAMtools/1.22.1-GCC-14.3.0", file = path, append = TRUE)
  write("", file = path, append = TRUE)

  write(
    paste0(
      "python3 -u ", script,
      " --bam ", bam,
      " --object ", se.obj,
      " --filter-on ", paste0("FASTQ_FILE_NAME==\'", job_id, "\' "),
      " --barcode-match-strategy regex ",
      " --barcode-regex '([ACGT]{16}-[0-9]+)$' ",
      " --duplicate-normalized-barcode-policy unassigned-ambiguous ",
      " --demux-on PATIENT_ID,TIMEPOINT ",
      " --out demux ",
      " --threads ", cpu
    ), file = path, append = TRUE)

  write("", file = path, append = TRUE)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# submit to Slurm
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# sample = unique(df$SAMPLE_SHORT)[[1]]
for (sample in names(fltrd.dirs)) {
  job_id = sample
  print(job_id)
  bam = fltrd.dirs[[job_id]]
  write_subscript(paste0(out.dir, "demux_", job_id, ".slurm"), bam, job_id)

  cmd = paste0("sbatch ", paste0(out.dir, "demux_", job_id, ".slurm"))
  system(cmd)
}


