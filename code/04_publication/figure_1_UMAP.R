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
.bioc_packages = c(
  "dittoSeq", "GenomicAlignments", "Rsamtools", "GenomicRanges", "DESeq2"
)

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
  obj@meta.data = droplevels(obj@meta.data)
  obj
}

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
  obj@meta.data$CONDITION = pd$CONDITION[
    match(obj@meta.data$orig.ident, pd$SAMPLE_NAME)
  ]

  obj@meta.data <- obj@meta.data %>% dplyr::mutate_if(is.character,as.factor)
  obj@meta.data = droplevels(obj@meta.data)
  obj
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA | Post Anno
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

seurat.pd = openxlsx::read.xlsx("data/pub_S51A/single_cell_seurat_metadata.xlsx")

se.merz = readRDS(
  paste0(manifest$merz$work, "seurat/03_seurat_anno_t.Rds")
)
se.braun = readRDS(paste0(manifest$braun$work, "seurat/03_seurat_anno_t.Rds"))
se.ho = readRDS(paste0(manifest$ho$work, "seurat/03_seurat_anno_t.Rds"))
se.hoso = readRDS(paste0(manifest$hosoya$work, "seurat/03_seurat_anno_t.Rds"))
se.hoso$orig.ident[se.hoso$orig.ident == "Duodenum"] = "GSM9314767"
se.hoso$orig.ident[se.hoso$orig.ident == "Stomach"] = "GSM9314768"
se.hoso$orig.ident[se.hoso$orig.ident == "Ileum"] = "GSM9314769"
se.alem = readRDS(paste0(manifest$aleman$work, "seurat/03_seurat_anno_t.Rds"))
se.oeke = readRDS(paste0(manifest$oekelen$work, "seurat/03_seurat_anno_t.Rds"))

se.merz = celltype_coarse(se.merz) # Rade
se.braun = celltype_coarse(se.braun) # Braun
se.ho = celltype_coarse(se.ho) # Ho
se.hoso = celltype_coarse(se.hoso) # Hosoya
se.alem = celltype_coarse(se.alem) # Aleman
se.oeke = celltype_coarse(se.oeke) # Aleman

se.merz = add_clin_metadata(obj = se.merz, pd = seurat.pd)
se.merz$STUDY = "Merz"
se.merz$SOURCE = "PB"
se.merz@meta.data$GROUP = ifelse(
  is.na(se.merz@meta.data$GROUP), "Control", "L-IEC"
)
se.braun = add_clin_metadata(obj = se.braun, pd = seurat.pd)
se.ho = add_clin_metadata(obj = se.ho, pd = seurat.pd)
se.hoso = add_clin_metadata(obj = se.hoso, pd = seurat.pd)
se.alem = add_clin_metadata(obj = se.alem, pd = seurat.pd)
se.oeke = add_clin_metadata(obj = se.oeke, pd = seurat.pd)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Merge | Post Anno
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
s.l = list(
  a = se.merz, b = se.braun, c = se.ho, d = se.hoso, e = se.alem, f = se.oeke
)
se.t = merge(s.l[[1]], y = s.l[2:length(s.l)])
se.t@meta.data$orig.ident = factor(se.t@meta.data$orig.ident)
se.t[["RNA"]] <- JoinLayers(se.t[["RNA"]])

se.t = integration(
  obj = DietSeurat(se.t, layers = "counts"),
  no.ftrs = 2000,
  threads = 30,
  .nbr.dims = 15,
  run.integration = T,
  harmony.group.vars = c("orig.ident", "SOURCE"),
  do.cluster = T
)

saveRDS(
  se.t@meta.data %>% dplyr::select(T_LIN_WORK),
  file = "data/t.Rds"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot | Post Anno
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cust.theme = theme(
  legend.position = "right",
  legend.margin = margin(l=-3),
  legend.key.spacing.y= unit(3, "pt"),
  legend.text = element_text(margin = margin(l = 3, unit = "pt")),
  legend.justification = c(0,.5),
  legend.title = element_text(margin = margin(b = 1)),
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  plot.title = element_blank()
)

pd = get_metadata(se.t)

set.seed(1234)
pd = pd[sample(nrow(pd), nrow(pd), replace = F), ]

pd$GROUP = gsub("Neurotoxicity", "irAE", pd$GROUP)
pd$GROUP = ifelse(pd$STUDY == "Aleman", "L-IEC without colitis", pd$GROUP)
pd$GROUP = ifelse(
  pd$STUDY != "Aleman" & pd$GROUP == "L-IEC", "L-IEC with colitis", pd$GROUP
)
table(pd$STUDY, pd$GROUP)

base.size = 20
theme_set(mytheme(base_size = base.size))

reduc.group.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = GROUP)) +
  scattermore::geom_scattermore(pointsize = 10, color="black")+
  scattermore::geom_scattermore(pointsize = 8.5, color="white") +
  ggrastr::geom_jitter_rast(shape = ".", raster.dpi = 300, scale = .5) +
  mytheme(base_size = base.size) +
  cust.theme +
  guides(colour = guide_legend(
    title = NULL,  override.aes = list(shape = 16, size = 6)
  )) +
  scale_color_manual(
    values = c("#DDAA33", "#BB5566", "#004488", "#000000")
  )

reduc.source.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = SOURCE)) +
  scattermore::geom_scattermore(pointsize = 10, color="black")+
  scattermore::geom_scattermore(pointsize = 8.5, color="white") +
  ggrastr::geom_jitter_rast(
    data = pd[pd$SOURCE == "PB", ], shape = ".", raster.dpi = 300, scale = .5
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[pd$SOURCE != "PB", ], shape = ".", raster.dpi = 300, scale = .5
  ) +

  mytheme(base_size = base.size) +
  cust.theme +
  guides(colour = guide_legend(
    title = NULL,  override.aes = list(shape = 16, size = 6)
  )) +
  scale_color_manual(
    values = c("#EECC66", "#997700", "#6699CC", "#004488", "#EE99AA", "#994455", "#000000")
  )

reduc.car.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = CAR_BY_EXPRS)) +
  scattermore::geom_scattermore(pointsize = 10, color="black")+
  scattermore::geom_scattermore(pointsize = 8.5, color="white") +
  ggrastr::geom_jitter_rast(shape = ".", raster.dpi = 300, scale = .5) +
  mytheme(base_size = base.size) +
  cust.theme +
  guides(colour = guide_legend(
    title = NULL,  override.aes = list(shape = 16, size = 6)
  )) +
  scale_color_manual(
    values = c("TRUE" = "#994455", "FALSE" = "#BBBBBB"),
    labels = c("TRUE" = "CAR+", "FALSE" = "CAR-")
  )

ggsave2(
  filename="figures/paper/figure_1_umap_2.png",
  plot_grid(
    NULL,
    plot_grid(
      reduc.group.pl, NULL, reduc.car.pl,
      align = "v",
      nrow = 3, rel_heights = c(1, .15, 1)
    ),
    ncol= 2, rel_widths = c(2.5, 1)
  ),
  width = 180, height = 52.5, dpi = 300, bg = "white", units = "mm", scale = 4,
  device = png, type = "cairo"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA | Pre Anno
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

seurat.pd = openxlsx::read.xlsx("data/pub_S51A/single_cell_seurat_metadata.xlsx")

se.merz.pre.anno = readRDS(
  paste0(manifest$merz$work, "seurat/02_seurat_pre.Rds")
)

se.braun.pre.anno = readRDS(paste0(manifest$braun$work, "seurat/02_seurat_pre.Rds"))
se.ho.pre.anno = readRDS(paste0(manifest$ho$work, "seurat/02_seurat_pre.Rds"))
se.hoso.pre.anno = readRDS(paste0(manifest$hosoya$work, "seurat/02_seurat_pre.Rds"))
se.hoso.pre.anno$orig.ident = as.character(se.hoso.pre.anno$orig.ident)
se.hoso.pre.anno$orig.ident[se.hoso.pre.anno$orig.ident == "Duodenum"] = "GSM9314767"
se.hoso.pre.anno$orig.ident[se.hoso.pre.anno$orig.ident == "Stomach"] = "GSM9314768"
se.hoso.pre.anno$orig.ident[se.hoso.pre.anno$orig.ident == "Ileum"] = "GSM9314769"
se.alem.pre.anno = readRDS(
  paste0(manifest$aleman$work, "seurat/02_seurat_pre.Rds")
)
se.oeke.pre.anno = readRDS(
  paste0(manifest$oekelen$work, "seurat/02_seurat_pre.Rds")
)

se.merz.pre.anno = add_clin_metadata(obj = se.merz.pre.anno, pd = seurat.pd)
se.merz.pre.anno$STUDY = "Merz"
se.merz.pre.anno$SOURCE = "PB"
se.merz.pre.anno@meta.data$GROUP = ifelse(
  is.na(se.merz.pre.anno@meta.data$GROUP), "Control", "L-IEC"
)
se.braun.pre.anno = add_clin_metadata(obj = se.braun.pre.anno, pd = seurat.pd)
se.ho.pre.anno = add_clin_metadata(obj = se.ho.pre.anno, pd = seurat.pd)
se.hoso.pre.anno = add_clin_metadata(obj = se.hoso.pre.anno, pd = seurat.pd)
se.alem.pre.anno = add_clin_metadata(obj = se.alem.pre.anno, pd = seurat.pd)
se.oeke.pre.anno = add_clin_metadata(obj = se.oeke.pre.anno, pd = seurat.pd)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Merge | Pre Anno
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
s.l = list(
  a = se.merz.pre.anno,
  b = se.braun.pre.anno,
  c = se.ho.pre.anno,
  d = se.hoso.pre.anno,
  e = se.alem.pre.anno,
  f = se.oeke.pre.anno
)
se.meta = merge(s.l[[1]], y = s.l[2:length(s.l)])
se.meta@meta.data$orig.ident = factor(se.meta@meta.data$orig.ident)
se.meta[["RNA"]] <- JoinLayers(se.meta[["RNA"]])


tt = se.meta@meta.data
tt = tt[!duplicated(se.meta$orig.ident), ]
tt$PATIENT_ID[tt$orig.ident == "MXMERZ002A_20"] = "P08"
tt$PATIENT_ID = ifelse(
  is.na(tt$PATIENT_ID),
  gsub("_.+", "", tt$orig.ident),
  tt$PATIENT_ID
)
nrow(tt)
length(table(tt$PATIENT_ID))




se.meta = se.meta[, se.meta$DOUBLETS_CONSENSUS == "FALSE"]
t = readRDS("data/t.Rds")
t$T_LIN_WORK = "T"
se.meta$T_LIN = t$T_LIN_WORK[match(colnames(se.meta), rownames(t))]
se.meta$T_LIN = ifelse(
  se.meta$T_LIN == "T", "T cells", NA
)

se.meta = integration(
  obj = DietSeurat(se.meta, layers = "counts"),
  no.ftrs = 2000,
  threads = 30,
  .nbr.dims = 30,
  run.integration = T,
  harmony.group.vars = c("orig.ident", "STUDY"),
  do.cluster = F
)

cust.theme = theme(
  legend.position = "none",
  legend.margin = margin(l=-3),
  legend.key.spacing.y= unit(0, "pt"),
  legend.text = element_text(margin = margin(l = 3, unit = "pt"), size = 18),
  legend.justification = c(0,.5),
  legend.title = element_text(margin = margin(b = 1)),
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  plot.title = element_blank()
)
pd = get_metadata(se.meta)

reduc.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = T_LIN)) +
  scattermore::geom_scattermore(pointsize = 10, color="black")+
  scattermore::geom_scattermore(pointsize = 8.5, color="white") +
  ggrastr::geom_jitter_rast(shape = ".", raster.dpi = 300, scale = .5) +
  mytheme(base_size = base.size) +
  cust.theme +
  guides(colour = guide_legend(
    title = NULL,  override.aes = list(shape = 16, size = 8)
  )) +
  scale_color_manual(values = c("#444444"), na.value = "#BBBBBB", limits = c("T cells")) +
  scale_x_reverse()

ggsave2(
  filename="figures/paper/figure_1_umap_1.png",
  plot_grid(
    NULL,
    reduc.pl,
    ncol= 2, rel_widths = c(2.2, 1)
  ),
  width = 180, height = 55, dpi = 300, bg = "white", units = "mm", scale = 4,
  device = png, type = "cairo"
)

# DimPlot(se.meta, group.by = "CD3_BY_EXPRS")
# DimPlot(se.meta, group.by = "CD4CD8_BY_EXPRS")
# DimPlot(se.meta, group.by = "STUDY")
# FeaturePlot_scCustom(se.meta, features = "KLRF1")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA | Post Anno
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.merz = readRDS(
  paste0(manifest$merz$work, "seurat/03_seurat_anno_t.Rds")
)
se.braun = readRDS(paste0(manifest$braun$work, "seurat/03_seurat_anno_t.Rds"))
se.ho = readRDS(paste0(manifest$ho$work, "seurat/03_seurat_anno_t.Rds"))
se.hoso = readRDS(paste0(manifest$hosoya$work, "seurat/03_seurat_anno_t.Rds"))
se.hoso$orig.ident[se.hoso$orig.ident == "Duodenum"] = "GSM9314767"
se.hoso$orig.ident[se.hoso$orig.ident == "Stomach"] = "GSM9314768"
se.hoso$orig.ident[se.hoso$orig.ident == "Ileum"] = "GSM9314769"
se.alem = readRDS(paste0(manifest$aleman$work, "seurat/03_seurat_anno_t.Rds"))
se.oeke = readRDS(paste0(manifest$oekelen$work, "seurat/03_seurat_anno_t.Rds"))

se.merz = celltype_coarse(se.merz) # Rade
se.braun = celltype_coarse(se.braun) # Braun
se.ho = celltype_coarse(se.ho) # Ho
se.hoso = celltype_coarse(se.hoso) # Hosoya
se.alem = celltype_coarse(se.alem) # Aleman
se.oeke = celltype_coarse(se.oeke) # Aleman

se.merz = add_clin_metadata(obj = se.merz, pd = seurat.pd)
se.merz$STUDY = "Merz"
se.merz$SOURCE = "PB"
se.merz@meta.data$GROUP = ifelse(
  is.na(se.merz@meta.data$GROUP), "Control", "L-IEC"
)
se.braun = add_clin_metadata(obj = se.braun, pd = seurat.pd)
se.ho = add_clin_metadata(obj = se.ho, pd = seurat.pd)
se.hoso = add_clin_metadata(obj = se.hoso, pd = seurat.pd)
se.alem = add_clin_metadata(obj = se.alem, pd = seurat.pd)
se.oeke = add_clin_metadata(obj = se.oeke, pd = seurat.pd)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Merge | Post Anno
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
s.l = list(
  a = se.merz, b = se.braun, c = se.ho, d = se.hoso, e = se.alem, f = se.oeke
)
se.t = merge(s.l[[1]], y = s.l[2:length(s.l)])
se.t@meta.data$orig.ident = factor(se.t@meta.data$orig.ident)
se.t[["RNA"]] <- JoinLayers(se.t[["RNA"]])

se.t = integration(
  obj = DietSeurat(se.t, layers = "counts"),
  no.ftrs = 2000,
  threads = 30,
  .nbr.dims = 20,
  run.integration = T,
  harmony.group.vars = c("orig.ident", "SOURCE"),
  do.cluster = T
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot | Post Anno
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cust.theme = theme(
  legend.position = "right",
  legend.margin = margin(l=-3),
  legend.key.spacing.y= unit(0, "pt"),
  legend.text = element_text(margin = margin(l = 3, unit = "pt")),
  legend.justification = c(0,.5),
  legend.title = element_text(margin = margin(b = 1)),
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  plot.title = element_blank()
)

pd = get_metadata(se.meta)

set.seed(1234)
pd = pd[sample(nrow(pd), nrow(pd), replace = F), ]

unique(pd$GROUP)
pd$GROUP = gsub("Neurotoxicity", "irAE", pd$GROUP)
pd$GROUP = ifelse(pd$STUDY == "Aleman", "L-IEC without colitis", pd$GROUP)
pd$GROUP = ifelse(pd$STUDY != "Aleman" & pd$GROUP == "L-IEC", "L-IEC with colitis", pd$GROUP)

table(pd$STUDY, pd$GROUP)

base.size = 14
theme_set(mytheme(base_size = base.size))


reduc.study.pl =
ggplot(data = pd, aes(x = umap_1, y = umap_2, col = STUDY)) +
  scattermore::geom_scattermore(pointsize = 10, color="black")+
  scattermore::geom_scattermore(pointsize = 8.5, color="white") +
  ggrastr::geom_jitter_rast(shape = ".", raster.dpi = 300, scale = .5) +
  mytheme(base_size = base.size) +
  cust.theme +
  guides(colour = guide_legend(
    title = NULL,  override.aes = list(shape = 16, size = 6)
  )) +
  scale_color_manual(values = colors.pal.10, na.value = "#FFFFFF")

reduc.group.pl =
ggplot(data = pd, aes(x = umap_1, y = umap_2, col = GROUP)) +
  scattermore::geom_scattermore(pointsize = 10, color="black")+
  scattermore::geom_scattermore(pointsize = 8.5, color="white") +
  ggrastr::geom_jitter_rast(shape = ".", raster.dpi = 300, scale = .5) +
  mytheme(base_size = base.size) +
  cust.theme +
  guides(colour = guide_legend(
    title = NULL,  override.aes = list(shape = 16, size = 6)
  )) +
  scale_color_manual(values = colors.pal.10, na.value = "#FFFFFF")

# DimPlot(se.meta, group.by = "celltype_coarse", cols = t.coarse.col)

ggsave2(
  filename="figures/paper/figure_1_umap.png",
  plot_grid(
    reduc.study.pl, NULL, reduc.group.pl, NULL,
    ncol = 4, rel_widths = c(1, .1, 1.1, 2)
  ),
  width = 180, height = 30, dpi = 300, bg = "white", units = "mm", scale = 4,
  device = png, type = "cairo"
)

# pd = se.merz@meta.data
# pd = pd[!duplicated(pd$orig.ident), ]
# pd = metadata_samples_merz[metadata_samples_merz$orig.ident %in% pd$orig.ident, ]
# pd = droplevels(pd)
# pd.late = subset(pd, TIMEPOINT == "Late")
# min(pd.late$TIME_CAR_DAY_30)
# max(pd.late$TIME_CAR_DAY_30)
# mean(pd.late$TIME_CAR_DAY_30)
# pd.verylate = subset(pd, TIMEPOINT == "Very Late")
# min(pd.verylate$TIME_CAR_DAY_100)
# max(pd.verylate$TIME_CAR_DAY_100)
# mean(pd.verylate$TIME_CAR_DAY_100)







