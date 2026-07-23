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
mem = 190000
cpu = 24

manifest = yaml.load_file("manifest.yaml")

fastqs = manifest$braun_rerun$fastq

ref.base = manifest$refs
ref.gex.cilta = paste0(ref.base, "GRCh38-ciltacel-mined/")

# VDJ reference: https://support.10xgenomics.com/single-cell-vdj/software/downloads/latest
ref.vdj = paste0(ref.base, "refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0/")
cd
ref.features = "/homes/olymp/michael.rade/projects/2024-BCMA-CAR-Koeln/data/feature_reference.csv"
out.dir = paste0(manifest$braun_rerun$work, "cellranger/")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# sbatch
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
write_subscript = function(path, job_id, csv.out){
  file.create(path)

  write("#!/bin/bash", file = path, append = TRUE)
  write("", file = path, append = TRUE)

  write(paste("#SBATCH -J", job_id), file = path, append = TRUE)
  write(paste("#SBATCH --nodes", nodes), file = path, append = TRUE)
  write(paste("#SBATCH --ntasks", ntasks), file = path, append = TRUE)
  write(paste("#SBATCH --time", ttime), file = path, append = TRUE)
  write(paste("#SBATCH --cpus-per-task", cpu), file = path, append = TRUE)
  write(paste("#SBATCH --mem", mem), file = path, append = TRUE)
  write(paste("#SBATCH --nodelist=ribnode[003-010]"), file = path, append = TRUE)
  write(paste("#SBATCH -e", paste0(job_id, ".e")), file = path, append = TRUE)
  write(paste("#SBATCH -o", paste0(job_id, ".o")), file = path, append = TRUE)
  write("#SBATCH --mail-type=END,FAIL", file = path, append = TRUE)

  write("", file = path, append = TRUE)
  write("module load CellRanger/10.0.0", file = path, append = TRUE)
  write("", file = path, append = TRUE)

  write(
    paste0(
      "cellranger multi",
      " --id ", job_id,
      " --csv ", csv.out,
      " --localcores=", cpu,
      " --localmem=", mem/1000
  ), file = path, append = TRUE)

  write("", file = path, append = TRUE)
}

write_multi_csv = function(sample.paths, csv.out){

  rna = sample.paths[sample.paths$SOURCE == "Gene Expression", , drop = F]
  adt = sample.paths[sample.paths$SOURCE == "Antibody Capture", , drop = F]
  tcr = sample.paths[sample.paths$SOURCE == "VDJ-T", , drop = F]
  bcr = sample.paths[sample.paths$SOURCE == "VDJ-B", , drop = F]

  file.create(csv.out)

  write(paste0("[gene-expression]"), file = csv.out, append = TRUE)
  write(paste0("reference,", rna$GEX_REF), file = csv.out, append = TRUE)
  write(paste0("create-bam,true"), file = csv.out, append = TRUE)
  write(paste0("[vdj]"), file = csv.out, append = TRUE)
  write(paste0("reference,", ref.vdj), file = csv.out, append = TRUE)
  write(paste0("[feature]"), file = csv.out, append = TRUE)
  write(paste0("reference,", ref.features), file = csv.out, append = TRUE)
  write(paste0("[libraries]"), file = csv.out, append = TRUE)
  write(paste0("fastq_id,fastqs,feature_types"), file = csv.out, append = TRUE)
  write(paste0(rna$SAMPLE, ",", rna$PATH, ",", rna$SOURCE), file = csv.out, append = TRUE)
  write(paste0(adt$SAMPLE, ",", adt$PATH, ",", adt$SOURCE), file = csv.out, append = TRUE)
  write(paste0(tcr$SAMPLE, ",", tcr$PATH, ",", tcr$SOURCE), file = csv.out, append = TRUE)
  write(paste0(bcr$SAMPLE, ",", bcr$PATH, ",", bcr$SOURCE), file = csv.out, append = TRUE)

}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# sample paths
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fastq.files = list.files(path = fastqs, full.names = T, recursive = T)
df = data.frame(
  SAMPLE = basename(fastq.files),
  PATH = dirname(fastq.files)
)
df = df[!grepl("sh$", df$SAMPLE), ]
df = df[!grepl("md5$", df$SAMPLE), ]
df$SAMPLE = gsub("_S.+", "", df$SAMPLE)

df = df %>% dplyr::mutate(
  SOURCE = dplyr::case_when(
    grepl("-R$", SAMPLE) ~ "Gene Expression",
    grepl("-A$", SAMPLE) ~ "Antibody Capture",
    grepl("-T$", SAMPLE) ~ "VDJ-T",
    grepl("-B$", SAMPLE) ~ "VDJ-B"
  )
)

df$SAMPLE_SHORT = gsub("-R|-A|-T|-B", "", df$SAMPLE)
df$GEX_REF = ref.gex.cilta

df = df[!duplicated(df$SAMPLE), ]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# submit to Slurm
# There is a link with this script under this path: "out.dir".
# The script was executed there.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for (sample in unique(df$SAMPLE_SHORT)) {
  job_id = paste0("multi_", sample)
  print(job_id)
  csv.out = paste0(out.dir, job_id, ".csv")
  sample.paths = subset(df, SAMPLE_SHORT == sample)
  write_multi_csv(sample.paths, csv.out)
  write_subscript(paste0(out.dir, job_id, ".slurm"), job_id, csv.out)

  cmd = paste0("sbatch ", paste0(out.dir, job_id, ".slurm"))
  system(cmd)
}
