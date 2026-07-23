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
mem = 250000
cpu = 30

manifest = yaml.load_file("manifest.yaml")
work.path = paste0(manifest$oekelen$work, "cellranger/")
seurat.path = paste0(manifest$oekelen$work, "seurat/")

cellranger.dirs = list.dirs(path = work.path, full.names = T, recursive = F)
cellranger.samples = basename(cellranger.dirs)
fltrd.dirs = paste0(
  cellranger.dirs, "/outs/per_sample_outs/", cellranger.samples,
  "/sample_alignments.bam"
)
names(fltrd.dirs) = gsub("multi_", "", cellranger.samples)

# se.ori = readRDS(paste0(paste0(seurat.path, "01_seurat_ori.Rds")))
# se.ori = se.ori[c("CD4", "CD8A", "CD8B", "ciltacel-mined"), ]
# se.ori@meta.data$barcodes = gsub(".+_", "", rownames(se.ori@meta.data))
# se.ori@meta.data$FASTQ_FILE_NAME = "SAPA15_CARTBT"
# se.ori@meta.data = se.ori@meta.data %>%
#   dplyr::select(orig.ident, barcodes, FASTQ_FILE_NAME)
# saveRDS(se.ori, paste0(manifest$homes, "data/pub_S51A/seurat_t_oekelen.Rds"))

se.obj = paste0(manifest$homes, "data/pub_S51A/seurat_t_oekelen.Rds")
script = paste0(manifest$homes, "/code/03_variant_workflow/demultiplexOnSeurat.py")

out.dir = paste0(
  manifest$work, "project_s51n/DWMCT/Analysis/demultiplexOnSeurat_oekelen/"
)

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
  write(paste("#SBATCH --nodelist=ribnode[003-012,016]"), file = path, append = TRUE)
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
      " --demux-on orig.ident ",
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


