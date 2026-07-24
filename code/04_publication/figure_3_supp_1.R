# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Library and functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "HGNChelper", "forcats", "cowplot",
  "remotes", "patchwork", "openxlsx", "scCustomize", "ggpubr", "tidyverse",
  "ggh4x", "ggrepel", "anndata", "scico", "matrixStats", "ggalluvial",
  "ggnewscale", "khroma"
)
.bioc_packages = c("dittoSeq", "GenomicAlignments", "Rsamtools", "GenomicRanges")

# Install CRAN packages (if not already installed)
.inst = .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}

# Install bioconductor packages (if not already installed)
.inst = .bioc_packages %in% installed.packages()
if (any(!.inst)) {
  library(BiocManager)
  BiocManager::install(.bioc_packages[!.inst], ask = T)
}

list.of.packages = c(.cran_packages, .bioc_packages)

## Loading library
for (pack in list.of.packages) {
  suppressMessages(library(
    pack,
    quietly = TRUE,
    verbose = FALSE,
    character.only = TRUE
  ))
}

if (any(!"scRepertoire" %in% installed.packages())) {
  # Sys.unsetenv("GITHUB_PAT")
  remotes::install_github(c("BorchLab/immApex", "BorchLab/scRepertoire@devel"))
}
library(scRepertoire)

if (any(!"Trex" %in% installed.packages())) {
  # Sys.unsetenv("GITHUB_PAT")
  devtools::install_github("BorchLab/Trex")
}
library(Trex)

source("code//helper/styles.R")
source("code//helper/functions.R")
source("code//helper/ora.R")
source("code//helper/functions_plots.R")
base.size = 8
theme_set(mytheme(base_size = base.size))

cust.theme = theme(
  legend.position = "right",
  legend.margin = margin(l=-3),
  legend.key.spacing.y= unit(-6, "pt"),
  legend.text = element_text(margin = margin(l = 0, unit = "pt"), size = 7),
  legend.justification = c(0,.5),
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  plot.title = element_blank()
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

t.mrkr = openxlsx::read.xlsx("data/signatures/mitch_t_marker.xlsx")
seurat.pd = openxlsx::read.xlsx("data/pub_S51A/single_cell_seurat_metadata.xlsx")

se.merz = readRDS(paste0(manifest$merz$work_p65, "seurat/03_seurat_anno_t_p65.Rds"))
se.braun = readRDS(paste0(manifest$braun$work, "seurat/03_seurat_anno_t.Rds"))
se.alem = readRDS(paste0(manifest$aleman$work, "seurat/03_seurat_anno_t.Rds"))
se.ho = readRDS(paste0(manifest$ho$work, "seurat/03_seurat_anno_t.Rds"))
se.oeke = readRDS(paste0(manifest$oekelen$work, "seurat/03_seurat_anno_t.Rds"))
se.hoso = readRDS(paste0(manifest$hosoya$work, "seurat/03_seurat_anno_t.Rds"))

se.hoso$orig.ident[se.hoso$orig.ident == "Duodenum"] = "GSM9314767"
se.hoso$orig.ident[se.hoso$orig.ident == "Stomach"] = "GSM9314768"
se.hoso$orig.ident[se.hoso$orig.ident == "Ileum"] = "GSM9314769"

se.oeke = RenameCells(
  se.oeke,
  new.names = paste0(
    se.oeke$orig.ident, "_", gsub(".+_", "", rownames(se.oeke@meta.data))
  )
)

se.ho.liec = se.ho[, se.ho$orig.ident %in% c("HTS_SO188_01_251_CSF", "HTS_SO188_05_251_PBMC")]
se.ho.liec@meta.data = droplevels(se.ho.liec@meta.data)

celltype_coarse = function(obj = NULL){
  obj@meta.data = obj@meta.data %>% dplyr::mutate(
    celltype_coarse = dplyr::case_when(
      grepl("CD4 CTL", celltype) ~ "CD4 CTL",
      grepl("CD4 Memory|CD4 NaiveLike", celltype) ~ "CD4 NaiveLike",
      grepl("CD4 Tfh|CD4 Th17", celltype) ~ "CD4 Helper",
      grepl("CD4 Treg", celltype) ~ "CD4 Treg",
      grepl("CD8 CM|CD8 NaiveLike", celltype) ~ "CD8 NaiveLike/CM",
      grepl("CD8 EM$", celltype) ~ "CD8 EM",
      grepl("CD8 EMRA", celltype) ~ "CD8 TEMRA",
      grepl("CD8 TEMRA", celltype) ~ "CD8 TEMRA",
      grepl("CD8.Cycling", celltype) ~ "CD8 Cycling",
      grepl("CD8 TPEX|CD8 TEX", celltype) ~ "CD8 Exhausted",
      TRUE ~ celltype
    )
  )
  obj = obj[, obj$celltype_coarse != "dnT"]
  k = names(table(obj$celltype_coarse)[table(obj$celltype_coarse) >= 10])
  obj = obj[, obj$celltype_coarse %in% k]
  obj@meta.data = droplevels(obj@meta.data)
  obj
}
se.merz = celltype_coarse(se.merz) # Rade
se.braun = celltype_coarse(se.braun) # Braun
se.alem = celltype_coarse(se.alem) # Aleman
se.ho = celltype_coarse(se.ho) # Ho
se.ho.liec = celltype_coarse(se.ho.liec) # Ho
se.oeke = celltype_coarse(se.oeke) # Oekelen
se.hoso = celltype_coarse(se.hoso) # Hosoya

add_clin_metadata = function(obj = NULL, pd = NULL){

  obj$SAMPLE_CLUSTER = NULL
  obj$GROUP = NULL
  obj$SOURCE = NULL
  obj$GSM = NULL

  obj@meta.data$PATIENT_ID = pd$PATIENT_ID[
    match(obj@meta.data$orig.ident, pd$SAMPLE_NAME)
  ]
  obj@meta.data$PATIENT_ID_UMAP = pd$PATIENT_ID_UMAP[
    match(obj@meta.data$orig.ident, pd$SAMPLE_NAME)
  ]
  obj@meta.data$STUDY = pd$STUDY[match(obj@meta.data$orig.ident, pd$SAMPLE_NAME)]
  obj@meta.data$GROUP = pd$GROUP[match(obj@meta.data$orig.ident, pd$SAMPLE_NAME)]
  obj@meta.data$SOURCE = pd$SOURCE[match(obj@meta.data$orig.ident, pd$SAMPLE_NAME)]
  obj@meta.data$TCR_Seq = pd$TCR_Seq[match(obj@meta.data$orig.ident, pd$SAMPLE_NAME)]
  obj@meta.data$TIME_AFTER_CAR = pd$TIME_AFTER_CAR[
    match(obj@meta.data$orig.ident, pd$SAMPLE_NAME)
  ]
  tp.lvls = c(
    "Aph", "LP", "2w", "14-128d", "1mo", "2mo", "4mo", "7mo", "8mo",
    "12mo", "13mo"
  )
  obj$TIME_AFTER_CAR = factor(obj$TIME_AFTER_CAR, levels = tp.lvls)
  obj@meta.data$CONDITION = pd$CONDITION[
    match(obj@meta.data$orig.ident, pd$SAMPLE_NAME)
  ]

  obj@meta.data <- obj@meta.data %>% dplyr::mutate_if(is.character,as.factor)
  obj@meta.data = droplevels(obj@meta.data)
  obj
}

se.merz = add_clin_metadata(obj = se.merz, pd = seurat.pd)
se.braun = add_clin_metadata(obj = se.braun, pd = seurat.pd)
se.alem = add_clin_metadata(obj = se.alem, pd = seurat.pd)
se.ho = add_clin_metadata(obj = se.ho, pd = seurat.pd)
se.ho.liec = add_clin_metadata(obj = se.ho.liec, pd = seurat.pd)
se.oeke = add_clin_metadata(obj = se.oeke, pd = seurat.pd)
se.hoso = add_clin_metadata(obj = se.hoso, pd = seurat.pd)

smpls.lbls = setNames(
  paste0(
    seurat.pd$PATIENT_ID, "|", seurat.pd$SOURCE, "\n",
    seurat.pd$TIME_AFTER_CAR
  ),
  seurat.pd$SAMPLE_NAME
)
smpls.lbls = gsub(" NA", "", smpls.lbls)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# pre-process seurat objects
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pre_process = function(se.obj = NULL, parse_vdj = T){

  if(parse_vdj == T){
    se.obj = se.obj[, !is.na(se.obj$CTnt)]
    se.obj = se.obj[, !grepl("NA", se.obj$CTnt), ]
    se.obj@meta.data = se.obj@meta.data %>%
      tibble::rownames_to_column(var="row.name") %>%
      dplyr::group_by(orig.ident) %>%
      dplyr::mutate(orig.ident.nbr = sum(!is.na(CTnt))) %>%
      dplyr::group_by(orig.ident, CTnt) %>%
      dplyr::mutate(clonalFrequency = n()) %>%
      dplyr::mutate(clonalProportion = clonalFrequency  / orig.ident.nbr) %>%
      tibble::column_to_rownames(var = "row.name") %>%
      data.frame()

    df = se.obj@meta.data %>%
      dplyr::group_by(CTnt) %>%
      dplyr::count(name = "cloneFreqAll") %>%
      dplyr::arrange(desc(cloneFreqAll))
    df$pseudo_id = paste0("Cl ", seq(1:nrow(df)))
    se.obj$CLONE_PSEUDO_ID = df$pseudo_id[match(
      se.obj$CTnt, df$CTnt
    )]
  }

  se.obj@meta.data = droplevels(se.obj@meta.data)
  se.obj@meta.data$barcode = gsub(".+_", "", rownames(se.obj@meta.data))

  set.seed(1234)
  suppressWarnings({
    suppressMessages({
      se.obj = integration(
        obj = DietSeurat(se.obj, layers = "counts"),
        no.ftrs = 1000,
        threads = 30,
        .nbr.dims = 15,
        run.integration = T,
        harmony.group.vars = c("orig.ident"),
        do.cluster = F,
        min.dist = .5
      )
    })
  })
  se.obj
}

se.merz = pre_process(se.merz)
se.braun = pre_process(se.braun)
# se.alem = pre_process(se.alem)
se.ho.liec = pre_process(se.ho.liec, parse_vdj = F)
se.hoso = pre_process(se.hoso, parse_vdj = F)

s.l = list(a = se.merz, b = se.braun, d = se.ho.liec, e = se.hoso, se.alem)
se.meta = merge(s.l[[1]], y = s.l[2:length(s.l)])
se.meta@meta.data$orig.ident = factor(se.meta@meta.data$orig.ident)
se.meta[["RNA"]] <- JoinLayers(se.meta[["RNA"]])

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Get Barcodes with muatation
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
variants.path = paste0(
  manifest$work, "project_s51n/DWMCT/Analysis/detectVariants_out/",
  "detectVariants_out_2026_07_21/"
)
variants.path.2 = paste0(
  manifest$work, "project_s51n/DWMCT/Analysis/detectVariants_out/",
  "detectVariants_out_noPhred_noMAPK_2026_07_21/"
)
read_variants = function(path = NULL){
  variants_files = list.files(
    path, full.names = T, recursive = T, pattern = "cb_mutation.tsv"
  )
  res = parallel::mclapply(variants_files,function(x){
    print(x)
    df = readr::read_tsv(x, id = "path") %>% as.data.frame()
    if(nrow(df) == 0){
      v = strsplit(basename(dirname(x)), "_xx_")[[1]]
      df[1, ] = c(NA, NA, NA, NA, NA, NA)
      df$path = x
      df$Nbr_reads = 0
      df$Nbr_reads_mut = 0
    }
    v = strsplit(basename(dirname(df$path[1])), "_xx_")[[1]]
    df$STUDY = v[1]
    df$SAMPLE = v[2]
    df$SOURCE = v[3]
    df$GROUP = v[4]
    df$path = NULL
    df$SAMPLE[df$SAMPLE == "Patient065_Late"] = "UC5A58"
    df$SAMPLE[df$SAMPLE == "Patient065_Very_Late"] = "UC5SJH"
    df$SAMPLE[df$SAMPLE == "Patient065_Very_Very_Late"] = "P065-PB"
    df$CB_ID = paste0(df$SAMPLE, "_", df$CB)
    df
  }, mc.cores = 40) %>% do.call("rbind", .)
}

variants = read_variants(path = variants.path)
variants = variants[variants$mut == "G222A" | is.na(variants$mut), ]

variants.2 = read_variants(path = variants.path.2)
variants.2 = variants.2[variants.2$mut == "G222A" | is.na(variants.2$mut), ]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Add nutation info to seurat object
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
add_mut = function(obj = NULL, bm = variants){
  obj@meta.data = obj@meta.data %>%
    tibble::rownames_to_column(var = "rowname") %>%
    dplyr::mutate(
      MUTATION = dplyr::case_when(
        rowname %in% bm[bm$Mutation == "yes", ]$CB_ID ~ "CAR+ S51N+",
        rowname %in% bm[bm$Mutation == "no", ]$CB_ID ~ "CAR+ S51N-",
        CAR_BY_EXPRS == "FALSE" ~ "CAR-",
        CAR_BY_EXPRS == "TRUE" ~ "CAR+",
        TRUE ~ "bra"
      )
    ) %>%
    tibble::column_to_rownames(var = "rowname")
  print(table(obj$MUTATION))
  obj
}

se.merz = add_mut(obj = se.merz)
se.braun = add_mut(obj = se.braun)
se.alem = add_mut(obj = se.alem)
se.ho = add_mut(obj = se.ho)
se.ho.liec = add_mut(obj = se.ho.liec)
se.hoso = add_mut(obj = se.hoso)
se.oeke = add_mut(obj = se.oeke)

se.meta = add_mut(obj = se.meta)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Supp | QC | After vs. Pre-Filtering
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
df = variants.2
df$Nbr_reads2 = variants$Nbr_reads[match(df$CB_ID, variants$CB_ID)]
df$Nbr_reads2[is.na(df$Nbr_reads2)] = 0

max.v = max(c(df$Nbr_reads2, df$Nbr_reads), na.rm = T)
read.fltr.pl =
  ggplot(df, aes(Nbr_reads, Nbr_reads2)) +
  geom_point(size = .01) +
  # geom_hex(bins = 55) +
  scico::scale_fill_scico(palette = "glasgow", direction = -1, end = .9) +
  theme(
    axis.text = element_text(size = 7),
  ) +
  ggtitle("Cell level | G222A | Reads per cell at the mutation site") +
  ylab("# reads covering the mutation site\nWith Phred and MAPK filter") +
  xlab("# reads covering the mutation site\nNo Phred and MAPK filter") +
  ylim(0, max.v) + xlim(0, max.v) +
  geom_abline(slope = 1, linetype = "dashed", linewidth = .3) +
  ggtitle("Cell level")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Supp | QC | Reads covering mutation site vs. reads with mutation
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
df = variants %>%
  dplyr::group_by(SAMPLE) %>%
  dplyr::mutate(Nbr_reads = sum(Nbr_reads)) %>%
  dplyr::mutate(Nbr_reads_mut = sum(Nbr_reads_mut)) %>%
  dplyr::distinct(SAMPLE, .keep_all = T) %>% data.frame()

smpls.lbls.2 = setNames(
  paste0(
    seurat.pd$STUDY, " | ", seurat.pd$PATIENT_ID, " | ", seurat.pd$SOURCE, " | ",
    seurat.pd$TIME_AFTER_CAR
  ),
  seurat.pd$SAMPLE_NAME
)
names(smpls.lbls.2)[names(smpls.lbls.2) == "P65_1_PB_L"] = "UC5A58"
names(smpls.lbls.2)[names(smpls.lbls.2) == "P65_1_PB_VL"] = "UC5SJH"

df$LABELS = smpls.lbls.2[match(df$SAMPLE, names(smpls.lbls.2))]
df$LABELS = ifelse(df$GROUP == "L-IEC", df$LABELS, df$GROUP)
df$LABELS = forcats::fct_relevel(df$LABELS, "Control", after = 0)
df$LABELS = forcats::fct_relevel(df$LABELS, "Neurotoxicity", after = 1)

tmp = seurat.pd
tmp$SAMPLE_NAME[tmp$SAMPLE_NAME == "P65_1_PB_L"] = "UC5A58"
tmp$SAMPLE_NAME[tmp$SAMPLE_NAME == "P65_1_PB_VL"] = "UC5SJH"

df$TIME = tmp$TIME_AFTER_CAR[match(df$SAMPLE, tmp$SAMPLE_NAME)]
df$TIME = ifelse(
  is.na(df$TIME),
  gsub(".+_", "", gsub("Late", "1mo", gsub("Very_Late", "3mo", df$SAMPLE))),
  df$TIME
)
df$TIME = gsub("^1mo$", "1-2mo", df$TIME)
df$TIME = gsub("^2mo$", "1-2mo", df$TIME)
df$TIME = gsub("^3mo$", "3-4mo", df$TIME)
df$TIME = gsub("^4mo$", "3-4mo", df$TIME)
df$TIME = gsub("^7mo$", "7-8mo", df$TIME)
df$TIME = gsub("^8mo$", "7-8mo", df$TIME)
df$TIME = gsub("^12mo$", "12-13mo", df$TIME)
df$TIME = gsub("^13mo$", "12-13mo", df$TIME)
tp.lvls = c("2w", "14-128d", "1-2mo", "3-4mo", "7-8mo", "12-13mo")
df$TIME = factor(df$TIME, levels = tp.lvls)

mut.vs.mutsite.pl1 =
ggplot() +
  theme(
    # aspect.ratio = 1,
    legend.margin = margin(l=-3),
    legend.key.spacing.y= unit(-6, "pt"),
    legend.justification = c(0,.5),
    legend.text = element_text(margin = margin(l = 0, unit = "pt"), size = 7),
    legend.spacing.y = unit(3, "pt"),
    legend.title = element_text(margin = margin(b = 0)),
    axis.text = element_text(size = 7),
    legend.box.margin = margin(20, 0, 0, 0)
  ) +
  geom_point(
    data = df[df$GROUP != "L-IEC", ],
    aes(Nbr_reads, Nbr_reads_mut, color = LABELS, shape = GROUP), size = 2
  ) +
  geom_point(
    data = df[df$GROUP == "L-IEC", ],
    aes(Nbr_reads, Nbr_reads_mut, color = LABELS, shape = GROUP), size = 2
  ) +
  ggpubr::stat_cor(
    data = df, aes(Nbr_reads, Nbr_reads_mut),
    method = "spearman",  size = 3, color = "black"
  ) +
  scale_y_continuous(
    trans='sqrt',
    breaks = c(10, 100, 1000, 5000, 10000, 20000),
    labels = c(10, 100, 1000, 5000, 10000,20000)
  ) +
  scale_x_continuous(
    trans='sqrt',
    breaks = c(10, 1000, 5000, 10000, 20000, 30000),
    labels = c(10, 1000, 5000, 10000, 20000, 30000)
  ) +
  geom_function(
    fun = function(x) x * 0.01 , colour = "black", linetype = "dashed"
  ) +
  scale_color_manual(
    values = c("black", "#BBBBBB", colors.pal.20), labels = smpls.lbls.2
  ) +
  ylab("# reads with mutation (G -> A)\nsqrt-scale") +
  xlab("# reads covering the mutation site\nsqrt-scale") +
  labs(color = NULL, shape = NULL) +
  ggtitle("Sample level")


discrete_rainbow = khroma::color("discrete rainbow")
mut.vs.mutsite.pl2 =
  ggplot() +
  theme(
    # aspect.ratio = 1,
    legend.margin = margin(l=-3),
    legend.key.spacing.y= unit(-6, "pt"),
    legend.justification = c(0,.5),
    legend.text = element_text(margin = margin(l = 0, unit = "pt"), size = 7),
    legend.spacing.y = unit(3, "pt"),
    legend.title = element_text(margin = margin(b = 0)),
    axis.text = element_text(size = 7)
  ) +
  geom_point(
    data = df,
    aes(Nbr_reads, Nbr_reads_mut, color = TIME, shape = GROUP), size = 2
  ) +
  scale_y_continuous(
    trans='sqrt',
    breaks = c(10, 100, 1000, 5000, 10000, 20000),
    labels = c(10, 100, 1000, 5000, 10000,20000)
  ) +
  scale_x_continuous(
    trans='sqrt',
    breaks = c(10, 1000, 5000, 10000, 20000, 30000),
    labels = c(10, 1000, 5000, 10000, 20000, 30000)
  ) +
  geom_function(
    fun = function(x) x * 0.01 , colour = "black", linetype = "dashed"
  ) +
  scale_color_manual(
    values = c(discrete_rainbow(6)), labels = smpls.lbls.2
  ) +
  ylab("# reads with mutation (G -> A)\nsqrt-scale") +
  xlab("# reads covering the mutation site\nsqrt-scale") +
  labs(color = "Time point after\nCAR treatment", shape = NULL) +
  ggtitle("Sample level")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Supp | QC | Reads covering mutation site vs. Sequencing depth
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# paths.1 = list.files(
#   paste0(manifest$work, "project_s51n/DWMCT/Analysis/demultiplexOnSeurat_only_demux_all/demux/"),
#   "\\.bam$", full = TRUE, recursive = F
# )
# paths.1 = paths.1[!grepl("unassigned", paths.1)]
# names(paths.1) = gsub(
#   "PATIENT_ID_TIMEPOINT__", "", gsub(".bam", "", basename(paths.1))
# )
#
# paths.2 = list.files(
#   paste0(manifest$work, "project_s51n/DWMCT/Analysis/braun_et_al/"),
#   "\\.bam$", full = TRUE, recursive = F
# )
# paths.2 = paths.2[!grepl("AphNB", paths.2)]
# names(paths.2) = gsub(".bam", "", basename(paths.2))
#
# paths.3 = dir(
#   paste0(manifest$work, "project_s51n/DWMCT/Analysis/hosoya_et_al/"),
#   "\\.bam$", full = TRUE, recursive = F
# )
# names(paths.3) = gsub(".bam", "", basename(paths.3))
#
# paths.4 = dir(
#   paste0(manifest$work, "project_s51n/DWMCT/Analysis/aleman_et_al/"),
#   "\\.bam$", full = TRUE, recursive = F
# )
# names(paths.4) = gsub(".bam", "", basename(paths.4))
#
# paths.5 = dir(
#   paste0(manifest$work, "project_s51n/DWMCT/Analysis/ho_et_al/"),
#   "\\.bam$", full = TRUE, recursive = F
# )
# names(paths.5) = gsub(".bam", "", basename(paths.5))
#
# paths.6 = dir(
#   paste0(manifest$work, "project_s51n/DWMCT/Analysis/demultiplexOnSeurat_oekelen/demux/"),
#   "\\.bam$", full = TRUE, recursive = F
# )
# paths.6 = paths.6[!grepl("unassigned", paths.6)]
# names(paths.6) = gsub(
#   "orig.ident__", "", gsub(".bam", "", basename(paths.6))
# )
#
# paths.7 = dir(
#   paste0(manifest$work, "project_s51n/DWMCT/Analysis/fandrei_p65_d365/"),
#   "\\.bam$", full = TRUE, recursive = F
# )
# names(paths.7) = gsub(".bam", "", basename(paths.7))
#
#
# paths = c(paths.1, paths.2, paths.3, paths.4, paths.5, paths.6, paths.7)
#
# # x = names(paths)[[1]]
# total.reads.l = parallel::mclapply(names(paths), function(x) {
#   id = x
#   # print(id)
#   param <- ScanBamParam(what = "qname")
#   bam <- Rsamtools::countBam(file = paths[[x]])
#   data.frame(ID = id, TOTAL_READS = bam$records)
# }, mc.cores = 50)
# total.reads = do.call("rbind", total.reads.l)
# saveRDS(object = total.reads, file = "data/pub_S51A/sequencing_depth.Rds")
total.reads = readRDS("data/pub_S51A/sequencing_depth.Rds")
total.reads$ID[total.reads$ID == "Patient065_Late"] = "UC5A58"
total.reads$ID[total.reads$ID == "Patient065_Very_Late"] = "UC5SJH"
total.reads$ID[total.reads$ID == "P065-PB-D365"] = "P065-PB"

df$SEQ_DEPTH = total.reads$TOTAL_READS[match(
  df$SAMPLE, total.reads$ID
)]
df$SEQ_DEPTH = round(df$SEQ_DEPTH / 1000000, 0)

mut.vs.mutsite.pl3 =
  ggplot() +
  theme(
    # aspect.ratio = 1,
    legend.margin = margin(l=-3),
    legend.key.spacing.y= unit(-6, "pt"),
    legend.justification = c(0,.5),
    legend.text = element_text(margin = margin(l = 0, unit = "pt"), size = 7),
    legend.spacing.y = unit(3, "pt"),
    legend.title = element_text(margin = margin(b = 0)),
    axis.text = element_text(size = 7),
    legend.box.margin = margin(20, 0, 0, 0)
  ) +
  geom_point(
    data = df[df$GROUP != "L-IEC", ],
    aes(Nbr_reads, SEQ_DEPTH, color = LABELS, shape = GROUP), size = 2
  ) +
  geom_point(
    data = df[df$GROUP == "L-IEC", ],
    aes(Nbr_reads, SEQ_DEPTH, color = LABELS, shape = GROUP), size = 2
  ) +
  ggpubr::stat_cor(
    data = df, aes(Nbr_reads, SEQ_DEPTH),
    method = "spearman",  size = 3, color = "black"
  ) +
  scale_x_continuous(
    trans='sqrt',
    breaks = c(10, 1000, 5000, 10000, 20000, 30000),
    labels = c(10, 1000, 5000, 10000, 20000, 30000)
  ) +
  scale_color_manual(
    values = c("black", "#BBBBBB", colors.pal.20), labels = smpls.lbls.2
  ) +
  ylab("Sequencing depth (million)") +
  xlab("# reads covering the mutation site\nsqrt-scale") +
  labs(color = NULL, shape = NULL) +
  ggtitle("Sample level")


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Supp | Overview | Absolute UMIs |  CAR-; CAR+; CAR+ S51N+; CAR+ S51N-
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
df = rbind(
  table(se.merz$orig.ident, se.merz$MUTATION) %>% data.frame(),
  table(se.braun$orig.ident, se.braun$MUTATION) %>% data.frame(),
  table(se.alem$orig.ident, se.alem$MUTATION) %>% data.frame(),
  table(se.ho$orig.ident, se.ho$MUTATION) %>% data.frame(),
  table(se.oeke$orig.ident, se.oeke$MUTATION) %>% data.frame(),
  table(se.hoso$orig.ident, se.hoso$MUTATION) %>% data.frame()
)
colnames(df) = c("orig.ident", "MUTATION", "Freq")

df$STUDY = seurat.pd$STUDY[match(df$orig.ident, seurat.pd$SAMPLE_NAME)]
df$GROUP = seurat.pd$GROUP[match(df$orig.ident, seurat.pd$SAMPLE_NAME)]
df$PATIENT_ID = seurat.pd$PATIENT_ID[match(df$orig.ident, seurat.pd$SAMPLE_NAME)]
df$TIME_AFTER_CAR = seurat.pd$TIME_AFTER_CAR[
  match(df$orig.ident, seurat.pd$SAMPLE_NAME)
]
df$CONDITION = seurat.pd$CONDITION[match(df$orig.ident, seurat.pd$SAMPLE_NAME)]
df = df[naturalorder(df$TIME_AFTER_CAR), ]
df$orig.ident = factor(df$orig.ident, levels = unique(df$orig.ident))

df = df[df$TIME_AFTER_CAR != "Aph", ]
df = df[df$GROUP == "L-IEC", ]

smpls.lbls.3 = setNames(
  paste0(
    seurat.pd$PATIENT_ID, " | ", seurat.pd$SOURCE, " | ",
    seurat.pd$TIME_AFTER_CAR
  ),
  seurat.pd$SAMPLE_NAME
)

df$STUDY2 = paste0(df$PATIENT_ID, "\n", df$STUDY)
overview.ss.pl =
ggplot(df, aes(x = orig.ident, y = Freq, label = Freq)) +
  geom_bar(stat="identity") +
  theme(
    plot.title = element_text(
      hjust = 0.5, face = "bold", size = rel(1), margin=margin(0,0,4,0)
    ),
    legend.key.size = unit(3, "mm"),
    legend.position = "bottom",
    legend.box.spacing = unit(0, "pt"),
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1, size = 6),
    ggh4x.facet.nestline = element_line(colour = "black", linewidth = .2),
    panel.spacing = unit(.5, "lines")
  ) +
  facet_nested(MUTATION ~  STUDY2, scales = "free_x", space = "free") +
  ylab("Nbr. of cells") + xlab(NULL) + labs(fill = NULL) +
  geom_text(
    nudge_y = 1170, size = 1.8
  ) +
  guides(fill = guide_legend(title = NULL, nrow = 1)) +
  scale_x_discrete(labels = smpls.lbls.3)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Supp | Barplot | Percentage |  CAR-; CAR+; CAR+ S51N+; CAR+ S51N-
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.merz$CLONE_PSEUDO_ID_TOPS = ifelse(
  se.merz$CLONE_PSEUDO_ID %in% c("Cl 1", "Cl 2", "Cl 3"),
  se.merz$CLONE_PSEUDO_ID,
  "Other clones"
)
se.braun$CLONE_PSEUDO_ID_TOPS = ifelse(
  se.braun$CLONE_PSEUDO_ID %in% c("Cl 1", "Cl 2", "Cl 3"),
  se.braun$CLONE_PSEUDO_ID,
  "Other clones"
)

t1 = prop.table(
  table(
    se.merz@meta.data[se.merz$orig.ident != "P65_1_PB_LP", ]$CLONE_PSEUDO_ID_TOPS,
    se.merz@meta.data[se.merz$orig.ident != "P65_1_PB_LP", ]$MUTATION
  ), margin = 1
) * 100
t2 = prop.table(
  table(
    se.braun@meta.data[se.braun$orig.ident != "AphNB", ]$CLONE_PSEUDO_ID_TOPS,
    se.braun@meta.data[se.braun$orig.ident != "AphNB", ]$MUTATION
  ), margin = 1
) * 100

df = rbind(
  data.frame(t1) %>% mutate(PATIENT_ID = "P1"),
  data.frame(t2) %>% mutate(PATIENT_ID = "P2")
)
df$Var1 = gsub("Other clones", "Other\nclones", df$Var1)

brpl.1.mut = ggplot(df, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat="identity") +
  geom_text(
    aes(label = ifelse(Freq > 10, paste0(round(Freq, 1), "%"), '')),
    position = position_stack(vjust = 0.5), size = 2.1
  ) +

  ylab("% of cells") + xlab(NULL) + labs(fill = NULL) +
  scale_fill_manual(
    values = c(
      "CAR+ S51N-" = "#BBBBBB", "CAR+ S51N+" = "#4C739F",
      "CAR+" = "#555555", "CAR-" = "#44AA99"
    )
  ) +
  facet_nested(. ~ PATIENT_ID, scales = "free", space = "free") +
  scale_x_discrete(labels = smpls.lbls.2) +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 7),
    legend.key.size = unit(3,"mm"),
    legend.key.spacing.y = unit(3, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(margin = margin(l = 2)),
    legend.position = "bottom",
    panel.spacing = unit(1, "lines")
  )

pd = se.meta@meta.data
t = prop.table(table(pd$orig.ident, pd$MUTATION), margin = 1) * 100
df = data.frame(t)
df$STUDY = seurat.pd$STUDY[match(df$Var1, seurat.pd$SAMPLE_NAME)]
df$PATIENT_ID = seurat.pd$PATIENT_ID[match(df$Var1, seurat.pd$SAMPLE_NAME)]
df = df[!df$STUDY %in% c("Braun", "Merz"), ]

smpls.lbls.tmp = setNames(
  paste0(seurat.pd$SOURCE, "\n", seurat.pd$TIME_AFTER_CAR),
  seurat.pd$SAMPLE_NAME
)
smpls.lbls.tmp = gsub("Duodenum", "Duo.", smpls.lbls.tmp)


brpl.2.mut = ggplot(df, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat="identity") +
  geom_text(
    aes(label = ifelse(Freq > 10, paste0(round(Freq, 1), "%"), '')),
    position = position_stack(vjust = 0.5), size = 2.1
  ) +

  ylab("% of cells") + xlab(NULL) + labs(fill = NULL) +
  scale_fill_manual(
    values = c(
      "CAR+ S51N-" = "#BBBBBB", "CAR+ S51N+" = "#4C739F",
      "CAR+" = "#555555", "CAR-" = "#44AA99"
    )
  ) +
  facet_nested(. ~ PATIENT_ID, scales = "free", space = "free") +
  scale_x_discrete(labels = smpls.lbls.tmp) +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 7),
    legend.key.size = unit(3,"mm"),
    legend.key.spacing.y = unit(3, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(margin = margin(l = 2)),
    legend.position = "bottom",
    panel.spacing = unit(1, "lines")
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cell type compositoin
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
smpls.lbls.3 = setNames(
  paste0(seurat.pd$TIME_AFTER_CAR, seurat.pd$SAMPLE_NAME),
  seurat.pd$SAMPLE_NAME
)
smpls.lbls.4 = setNames(
  paste0(seurat.pd$SOURCE, " | ", seurat.pd$TIME_AFTER_CAR),
  smpls.lbls.3
)
smpls.lbls.4 = gsub("Duodenum", "Duo.", smpls.lbls.4)

df.comp = dittoBarPlot(
  se.meta[, se.meta$TIME_AFTER_CAR != "Aph"], "celltype_coarse",
  group.by = "orig.ident", split.by = "MUTATION", data.out = T
)[[2]]
df.comp = df.comp %>% dplyr::mutate(
  LIN = dplyr::case_when(
    grepl("CD8", label) ~ "CD8",
    grepl("CD4", label) ~ "CD4",
    TRUE ~ label
  )
)
df.comp$percent = df.comp$percent*100
df.comp$label = factor(
  df.comp$label, levels = naturalsort(unique(df.comp$label))
)
df.comp$grouping2 = smpls.lbls.3[match(df.comp$grouping, names(smpls.lbls.3))]
df.comp$grouping2 = factor(
  df.comp$grouping2, levels = naturalsort(unique(df.comp$grouping2))
)

df.comp$PATIENT = seurat.pd$PATIENT_ID[
  match(df.comp$grouping, seurat.pd$SAMPLE_NAME)
]

bpl.comp =
  ggplot(df.comp, aes(x = grouping2, y = count, fill = label)) +
  geom_bar(stat="identity", position="fill") +
  theme(
    plot.title = element_text(
      hjust = 0.5, face = "bold", size = rel(1), margin=margin(0,0,4,0)
    ),
    legend.key.size = unit(3, "mm"),
    # legend.text = element_text(size = 6),
    legend.position = "bottom",
    legend.box.spacing = unit(0, "pt"),
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1, size = 7),
    ggh4x.facet.nestline = element_line(colour = "black", linewidth = .2),
    panel.spacing = unit(.45, "lines")
  ) +
  facet_nested(. ~ MUTATION + PATIENT, scales = "free", space = "free") +
  scale_fill_manual(
    values = t.coarse.col, # labels = setNames(df.comp$label_short, df.comp$label)
  ) +
  scale_x_discrete(labels = smpls.lbls.4) +
  scale_y_continuous(breaks = c(0, .5, 1), labels = c(0, 50, 100)) +
  ylab("% of cells") + xlab(NULL) + labs(fill = NULL) +
  ggtitle("Subtype composition") +
  geom_text(
    aes(label = ifelse(percent > 10, paste0(round(percent,1)), '')),
    position = position_fill(vjust = 0.5), size = 1.7, color = "white"
  ) +
  guides(fill = guide_legend(title = NULL, nrow = 2))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

ggsave2(
  filename="figures/paper/figure_3_supps_1.png",
  plot_grid(
    plot_grid(
      plot_grid(
        plot_grid(
          read.fltr.pl, NULL, mut.vs.mutsite.pl2, nrow = 3, align = "v",
          rel_heights = c(1, .1, 1), labels = c("A", "", "C"),
          label_fontface = "bold", label_size = 11, vjust = 1.1
        ),
        NULL,
        plot_grid(
          mut.vs.mutsite.pl1, NULL, mut.vs.mutsite.pl3, nrow = 3, align = "v",
          rel_heights = c(1, .1, 1), labels = c("B", "", "D"),
          label_fontface = "bold", label_size = 11, vjust = 1.1
        ),
        ncol = 3, rel_widths = c(1, .1, 1.25)
      ),
      NULL,
      overview.ss.pl,
      ncol = 3, rel_widths = c(1, .05, .45),
      rel_heights = c(1, .1, 1), labels = c("", "", "E"),
      label_fontface = "bold", label_size = 11, vjust = 1.1
    ),
    NULL,
    plot_grid(
      brpl.1.mut, NULL, brpl.2.mut, ncol = 3, rel_widths = c(1, .1, 1.5),
      rel_heights = c(1, .1, 1), labels = c("F", "", "G"),
      label_fontface = "bold", label_size = 11, vjust = 1.1
    ),
    NULL,
    plot_grid(
      bpl.comp, rel_heights = c(1, .1, 1), labels = c("H"),
      label_fontface = "bold", label_size = 11, vjust = 1.1
    ),
    nrow = 5, rel_heights = c(2.8, .15, 1.1, .1, 1.5)
  )
  ,
  width = 180, height = 150, dpi = 300, bg = "white", units = "mm", scale = 1.6,
  device = png, type = "cairo"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Supp | DimReduc | Merz | Patient 1
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dim = "umap"

# DimReduc ---------------------------------------------------------------------
se.merz$CLONE_PSEUDO_ID_TOPS = ifelse(
  as.character(se.merz$CLONE_PSEUDO_ID) %in% c("Cl 1", "Cl 2", "Cl 3"),
  as.character(se.merz$CLONE_PSEUDO_ID),
  "Other clones"
)

pt.scale = .75

pd = get_metadata(se.merz)
pd$DIM1 = pd[[paste0(dim, "_1")]]
pd$DIM2 = pd[[paste0(dim, "_2")]]
ct.to.pl = "celltype_coarse"

celltypes.reduc.pl  =
  ggplot() +
  scattermore::geom_scattermore(
    data = pd, aes(x = DIM1, y = DIM2), pointsize = 17, color="black"
  )+
  scattermore::geom_scattermore(
    data = pd, aes(x = DIM1, y = DIM2),  pointsize = 15, color="white"
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[grepl("CD4", pd$celltype_coarse), ],
    aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]),
    shape = ".", raster.dpi = 300, scale = .5
  ) +
  guides(colour = guide_legend(
    title = NULL, ncol = 1, override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(values = t.coarse.col) +
  ggnewscale::new_scale_color() +
  ggrastr::geom_jitter_rast(
    data = pd[grepl("CD8", pd$celltype_coarse), ],
    aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]),
    shape = ".", raster.dpi = 300, scale = .5
  ) +
  guides(colour = guide_legend(
    title = NULL, ncol = 1,  override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(values = t.coarse.col) +
  cust.theme +
  theme(
    legend.spacing.y = unit(2, 'pt'),
    legend.title = element_text(margin = margin(b = -2))
  )

time.reduc.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = TIME_AFTER_CAR)) +
  scattermore::geom_scattermore(pointsize = 17, color="black")+
  scattermore::geom_scattermore(pointsize = 15, color="white") +
  ggrastr::geom_jitter_rast(shape = ".", raster.dpi = 300, scale = .5) +
  cust.theme +
  guides(colour = guide_legend(
    title = NULL, ncol = 1, override.aes = list(shape = 16, size = 3)
  )) + scale_color_manual(values = c("#6699CC", "#EECC66", "#997700", "#994455"))

cl.c = c(
  "Cl 1" = "#004488", "Cl 2" = "#0077BB", "Cl 3" = "#88CCEE", "Other clones" = "#BBBBBB"
)
top.clones.reduc.pl  =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = CLONE_PSEUDO_ID_TOPS)) +
  scattermore::geom_scattermore(pointsize = 17, color="black")+
  scattermore::geom_scattermore(pointsize = 15, color="white") +
  ggrastr::geom_jitter_rast(
    data = pd[pd$CLONE_PSEUDO_ID_TOPS == "other", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  ggrastr::geom_jitter_rast(
    data =  pd[pd$CLONE_PSEUDO_ID_TOPS != "other", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  cust.theme +
  guides(colour = guide_legend(
    title = NULL, ncol = 1, override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(values = cl.c, na.value = "#BBBBBB")

car.reduc.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = MUTATION)) +
  scattermore::geom_scattermore(pointsize = 17, color="black")+
  scattermore::geom_scattermore(pointsize = 15, color="white") +
  ggrastr::geom_jitter_rast(
    data = pd[pd$MUTATION == "CAR-", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[pd$MUTATION == "CAR+", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[pd$MUTATION == "CAR+ S51N-", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[pd$MUTATION == "CAR+ S51N+", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  cust.theme +
  guides(colour = guide_legend(
    title = NULL,  ncol = 1, override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(values = c(
    "CAR+ S51N-" = "#6699CC", "CAR+ S51N+" = "#CC6677", "CAR-" = "#BBBBBB",
    "CAR+" = "#997700"
  )) # #CC3311


#  Annotating possible epitopes
se.merz$TRA_Epitope.species = NULL
se.merz = Trex::annotateDB(se.merz, chains = "TRA", edit.distance = 2)
f = table(se.merz$TRA_Epitope.species)
se.merz$TRA_Epitope.species = ifelse(
  se.merz$TRA_Epitope.species %in% names(f[f < 50]),
  "other", se.merz$TRA_Epitope.species
)
pd = get_metadata(se.merz)
epi.reduc.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = TRA_Epitope.species)) +
  scattermore::geom_scattermore(pointsize = 17, color="black")+
  scattermore::geom_scattermore(pointsize = 15, color="white") +
  ggrastr::geom_jitter_rast(
    data = pd[is.na(pd$TRA_Epitope.species), ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[!is.na(pd$TRA_Epitope.species), ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  cust.theme +
  guides(colour = guide_legend(
    title = NULL, ncol = 1, override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(values = colors.pal.10, na.value = "#DDDDDD")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Supp | DimReduc | Braun
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# DimReduc ---------------------------------------------------------------------
se.braun$CLONE_PSEUDO_ID_TOPS = ifelse(
  as.character(se.braun$CLONE_PSEUDO_ID) %in% c("Cl 1", "Cl 2", "Cl 3"),
  as.character(se.braun$CLONE_PSEUDO_ID),
  "Other clones"
)

pt.scale = .75

pd = get_metadata(se.braun)
pd$DIM1 = pd[[paste0(dim, "_1")]]
pd$DIM2 = pd[[paste0(dim, "_2")]]
ct.to.pl = "celltype_coarse"

celltypes.reduc.braun.pl  =
  ggplot() +
  scattermore::geom_scattermore(
    data = pd, aes(x = DIM1, y = DIM2), pointsize = 17, color="black"
  )+
  scattermore::geom_scattermore(
    data = pd, aes(x = DIM1, y = DIM2),  pointsize = 15, color="white"
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[grepl("CD4", pd$celltype_coarse), ],
    aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]),
    shape = ".", raster.dpi = 300, scale = .5
  ) +
  guides(colour = guide_legend(
    title = NULL, ncol = 1, override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(values = t.coarse.col) +
  ggnewscale::new_scale_color() +
  ggrastr::geom_jitter_rast(
    data = pd[grepl("CD8", pd$celltype_coarse), ],
    aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]),
    shape = ".", raster.dpi = 300, scale = .5
  ) +
  guides(colour = guide_legend(
    title = NULL, ncol = 1,  override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(values = t.coarse.col) +
  cust.theme +
  theme(
    legend.spacing.y = unit(2, 'pt'),
    legend.title = element_text(margin = margin(b = -2))
  )

time.reduc.braun.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = TIME_AFTER_CAR)) +
  scattermore::geom_scattermore(pointsize = 17, color="black")+
  scattermore::geom_scattermore(pointsize = 15, color="white") +
  ggrastr::geom_jitter_rast(shape = ".", raster.dpi = 300, scale = .5) +
  cust.theme +
  guides(colour = guide_legend(
    title = NULL, ncol = 1, override.aes = list(shape = 16, size = 3)
  )) + scale_color_manual(values = c("#6699CC", "#EECC66", "#994455"))

cl.c = c(
  "Cl 1" = "#999933", "Cl 2" = "#004488",
  "Cl 3" = "#DDCC77", "Other clones" = "#BBBBBB"
)

top.clones.reduc.braun.pl  =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = CLONE_PSEUDO_ID_TOPS)) +
  scattermore::geom_scattermore(pointsize = 17, color="black")+
  scattermore::geom_scattermore(pointsize = 15, color="white") +
  ggrastr::geom_jitter_rast(
    data = pd[pd$CLONE_PSEUDO_ID_TOPS == "other", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  ggrastr::geom_jitter_rast(
    data =  pd[pd$CLONE_PSEUDO_ID_TOPS != "other", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  cust.theme +
  guides(colour = guide_legend(
    title = NULL, ncol = 1, override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(values = cl.c, na.value = "#BBBBBB")

car.reduc.braun.pl  =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = MUTATION)) +
  scattermore::geom_scattermore(pointsize = 17, color="black")+
  scattermore::geom_scattermore(pointsize = 15, color="white") +
  ggrastr::geom_jitter_rast(
    data = pd[pd$MUTATION == "CAR-", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[pd$MUTATION == "CAR+", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[pd$MUTATION == "CAR+ S51N-", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[pd$MUTATION == "CAR+ S51N+", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  cust.theme +
  guides(colour = guide_legend(
    title = NULL,  ncol = 1, override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(values = c(
    "CAR+ S51N-" = "#6699CC", "CAR+ S51N+" = "#CC6677", "CAR-" = "#BBBBBB",
    "CAR+" = "#997700"
  )) # #CC3311

#  Annotating possible epitopes
se.braun$TRA_Epitope.species = NULL
se.braun = Trex::annotateDB(se.braun, chains = "TRA", edit.distance = 2)
f = table(se.braun$TRA_Epitope.species)

se.braun$TRA_Epitope.species = ifelse(
  se.braun$TRA_Epitope.species %in% names(f[f < 50]),
  "other", se.braun$TRA_Epitope.species
)
pd = get_metadata(se.braun)

pd$TRA_Epitope.species = gsub(";CMV", ";\nCMV", pd$TRA_Epitope.species)

epi.reduc.braun.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = TRA_Epitope.species)) +
  scattermore::geom_scattermore(pointsize = 17, color="black")+
  scattermore::geom_scattermore(pointsize = 15, color="white") +
  ggrastr::geom_jitter_rast(
    data = pd[is.na(pd$TRA_Epitope.species), ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[!is.na(pd$TRA_Epitope.species), ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  cust.theme +
  theme(
    legend.text = element_text(lineheight = .75),
    legend.key.spacing.y= unit(-3, "pt"),
  ) +
  guides(colour = guide_legend(
    title = NULL, ncol = 1, override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(values = colors.pal.10, na.value = "#DDDDDD")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Supp | DimReduc | Hosoya
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# DimReduc ---------------------------------------------------------------------

pt.scale = .75

pd = get_metadata(se.hoso)

pd$DIM1 = pd[[paste0(dim, "_1")]]
pd$DIM2 = pd[[paste0(dim, "_2")]]
ct.to.pl = "celltype_coarse"

celltypes.reduc.hoso.pl  =
  ggplot() +
  scattermore::geom_scattermore(
    data = pd, aes(x = DIM1, y = DIM2), pointsize = 17, color="black"
  )+
  scattermore::geom_scattermore(
    data = pd, aes(x = DIM1, y = DIM2),  pointsize = 15, color="white"
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[grepl("CD4", pd$celltype_coarse), ],
    aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]),
    shape = ".", raster.dpi = 300, scale = .5
  ) +
  guides(colour = guide_legend(
    title = NULL, ncol = 1, override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(values = t.coarse.col) +
  ggnewscale::new_scale_color() +
  ggrastr::geom_jitter_rast(
    data = pd[grepl("CD8", pd$celltype_coarse), ],
    aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]),
    shape = ".", raster.dpi = 300, scale = .5
  ) +
  guides(colour = guide_legend(
    title = NULL, ncol = 1,  override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(values = t.coarse.col) +
  cust.theme +
  theme(
    legend.spacing.y = unit(2, 'pt'),
    legend.title = element_text(margin = margin(b = -2))
  )

celltypes.source.hoso.pl  =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = SOURCE)) +
  scattermore::geom_scattermore(pointsize = 17, color="black")+
  scattermore::geom_scattermore(pointsize = 15, color="white") +
  ggrastr::geom_jitter_rast(shape = ".", raster.dpi = 300, scale = .5) +
  cust.theme +
  guides(colour = guide_legend(
    title = NULL, ncol = 1, override.aes = list(shape = 16, size = 3)
  )) + scale_color_manual(values = c("#DDAA33", "#BB5566", "#004488"))

car.reduc.hoso.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = MUTATION)) +
  scattermore::geom_scattermore(pointsize = 17, color="black")+
  scattermore::geom_scattermore(pointsize = 15, color="white") +
  ggrastr::geom_jitter_rast(
    data = pd[pd$MUTATION == "CAR-", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[pd$MUTATION == "CAR+", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[pd$MUTATION == "CAR+ S51N-", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[pd$MUTATION == "CAR+ S51N+", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  cust.theme +
  guides(colour = guide_legend(
    title = NULL,  ncol = 1, override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(values = c(
    "CAR+ S51N-" = "#6699CC", "CAR+ S51N+" = "#CC6677", "CAR-" = "#BBBBBB",
    "CAR+" = "#997700"
  )) # #CC3311

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Supp | DimReduc | Ho
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# DimReduc ---------------------------------------------------------------------

pt.scale = .75

pd = get_metadata(se.ho.liec)

pd$DIM1 = pd[[paste0(dim, "_1")]]
pd$DIM2 = pd[[paste0(dim, "_2")]]
ct.to.pl = "celltype_coarse"

celltypes.reduc.ho.pl  =
  ggplot() +
  scattermore::geom_scattermore(
    data = pd, aes(x = DIM1, y = DIM2), pointsize = 17, color="black"
  )+
  scattermore::geom_scattermore(
    data = pd, aes(x = DIM1, y = DIM2),  pointsize = 15, color="white"
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[grepl("CD4", pd$celltype_coarse), ],
    aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]),
    shape = ".", raster.dpi = 300, scale = .5
  ) +
  guides(colour = guide_legend(
    title = NULL, ncol = 1, override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(values = t.coarse.col) +
  ggnewscale::new_scale_color() +
  ggrastr::geom_jitter_rast(
    data = pd[grepl("CD8", pd$celltype_coarse), ],
    aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]),
    shape = ".", raster.dpi = 300, scale = .5
  ) +
  guides(colour = guide_legend(
    title = NULL, ncol = 1,  override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(values = t.coarse.col) +
  cust.theme +
  theme(
    legend.spacing.y = unit(2, 'pt'),
    legend.title = element_text(margin = margin(b = -2))
  )

celltypes.source.ho.pl  =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = SOURCE)) +
  scattermore::geom_scattermore(pointsize = 17, color="black")+
  scattermore::geom_scattermore(pointsize = 15, color="white") +
  ggrastr::geom_jitter_rast(shape = ".", raster.dpi = 300, scale = .5) +
  cust.theme +
  guides(colour = guide_legend(
    title = NULL, ncol = 1, override.aes = list(shape = 16, size = 3)
  )) + scale_color_manual(values = c("#DDAA33", "#BB5566", "#004488"))

car.reduc.ho.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = MUTATION)) +
  scattermore::geom_scattermore(pointsize = 17, color="black")+
  scattermore::geom_scattermore(pointsize = 15, color="white") +
  ggrastr::geom_jitter_rast(
    data = pd[pd$MUTATION == "CAR-", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[pd$MUTATION == "CAR+", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[pd$MUTATION == "CAR+ S51N-", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[pd$MUTATION == "CAR+ S51N+", ],
    shape = ".", raster.dpi = 300, scale = pt.scale
  ) +
  cust.theme +
  guides(colour = guide_legend(
    title = NULL,  ncol = 1, override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(values = c(
    "CAR+ S51N-" = "#6699CC", "CAR+ S51N+" = "#CC6677", "CAR-" = "#BBBBBB",
    "CAR+" = "#997700"
  )) # #CC3311

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DimReduc Plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ggsave2(
  filename="figures/paper/figure_3_supps_2.png",
  plot_grid(
    plot_grid(
      plot_grid(NULL, celltypes.reduc.pl, nrow = 2, rel_heights = c(.05, 1)),
      NULL,
      plot_grid(NULL, time.reduc.pl, nrow = 2, rel_heights = c(.05, 1)),
      NULL,
      plot_grid(NULL, top.clones.reduc.pl, nrow = 2, rel_heights = c(.05, 1)),
      NULL,
      plot_grid(NULL, car.reduc.pl, nrow = 2, rel_heights = c(.05, 1)),
      NULL,
      plot_grid(NULL, epi.reduc.pl, nrow = 2, rel_heights = c(.05, 1)),
      ncol = 9, rel_widths = c(1.1, .1, .85, .1, 1, .1, 1, .1, 1),
      labels = c("A", "", "B", "", "C", "", "D", "", "E"),
      label_fontface = "bold", label_size = 11, vjust = 1.1
    ),
    NULL,
    plot_grid(
      plot_grid(NULL, celltypes.reduc.braun.pl, nrow = 2, rel_heights = c(.05, 1)),
      NULL,
      plot_grid(NULL, time.reduc.braun.pl, nrow = 2, rel_heights = c(.05, 1)),
      NULL,
      plot_grid(NULL, top.clones.reduc.braun.pl, nrow = 2, rel_heights = c(.05, 1)),
      NULL,
      plot_grid(NULL, car.reduc.braun.pl, nrow = 2, rel_heights = c(.05, 1)),
      NULL,
      plot_grid(NULL, epi.reduc.braun.pl, nrow = 2, rel_heights = c(.05, 1)),
      ncol = 9, rel_widths = c(1.1, .1, .85, .1, 1, .1, 1, .1, 1),
      labels = c("F", "", "G", "", "H", "", "I", "", "J"),
      label_fontface = "bold", label_size = 11, vjust = 1.1
    ),
    NULL,
    plot_grid(
      plot_grid(NULL, celltypes.reduc.hoso.pl, nrow = 2, rel_heights = c(.05, 1)),
      NULL,
      plot_grid(NULL, celltypes.source.hoso.pl, nrow = 2, rel_heights = c(.05, 1)),
      NULL,
      plot_grid(NULL, car.reduc.hoso.pl, nrow = 2, rel_heights = c(.05, 1)),
      NULL,
      NULL,
      ncol = 7, rel_widths = c(1.1, .1, .85, .1, 1, .1, 2.1),
      labels = c("K", "", "L", "", "M"),
      label_fontface = "bold", label_size = 11, vjust = 1.1
    ),
    NULL,
    plot_grid(
      plot_grid(NULL, celltypes.reduc.ho.pl, nrow = 2, rel_heights = c(.05, 1)),
      NULL,
      plot_grid(NULL, celltypes.source.ho.pl, nrow = 2, rel_heights = c(.05, 1)),
      NULL,
      plot_grid(NULL, car.reduc.ho.pl, nrow = 2, rel_heights = c(.05, 1)),
      NULL,
      NULL,
      ncol = 7, rel_widths = c(1.1, .1, .85, .1, 1, .1, 2.1),
      labels = c("O", "", "P", "", "Q"),
      label_fontface = "bold", label_size = 11, vjust = 1.1
    ),
    nrow = 7, rel_heights = c(1, .25, 1, .25, 1, .25, 1)
  )
  ,
  width = 180, height = 80, dpi = 300, bg = "white", units = "mm", scale = 1.6,
  device = png, type = "cairo"
)
