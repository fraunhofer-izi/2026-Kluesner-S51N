# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Library and functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "HGNChelper", "forcats", "cowplot",
  "remotes", "patchwork", "openxlsx", "scCustomize", "ggpubr", "tidyverse",
  "ggh4x", "ggrepel", "anndata", "scico", "matrixStats", "ggalluvial"
)
.bioc_packages = c(
  "dittoSeq", "GenomicAlignments", "Rsamtools", "GenomicRanges", "Nebulosa"
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

source("code//helper/styles.R")
source("code//helper/functions.R")
source("code//helper/ora.R")
source("code//helper/functions_plots.R")
base.size = 8
base_size = 8
theme_set(theme_custom())
# theme_set(mytheme(base_size = base.size))

patient.col = c(
  "P1" = "#4A7BB7",
  "P2" = "#44BB99",
  "P3" = "#EEDD88",
  "P4" = "#EE8866",
  "P5" = "#FFAABB",
  "P6" = "#99DDFF"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

t.mrkr = openxlsx::read.xlsx("data/signatures/mitch_t_marker.xlsx")
t.mrkr$Gene_symbol[t.mrkr$Gene_symbol == "PDCD1"] = "PD1"
t.mrkr$Gene_symbol[t.mrkr$Gene_symbol == "HAVCR2"] = "TIM3"
t.mrkr$Gene_symbol[t.mrkr$Gene_symbol == "MKI67"] = "Ki-67"
t.mrkr$Gene_symbol[t.mrkr$Gene_symbol == "SELL"] = "CD62L"
t.mrkr$Group[grepl("Activation", t.mrkr$Group)] = "Effector,\nCytotoxicity"
t.mrkr$Group[grepl("Exhaustion", t.mrkr$Group)] = "Activation,\nChronic stimulation"
t.mrkr$Group[grepl("Stem", t.mrkr$Group)] = "Stem-like,\nMemory"

seurat.pd = openxlsx::read.xlsx("data/pub_S51A/single_cell_seurat_metadata.xlsx")

se.merz = readRDS(paste0(manifest$merz$work_p65, "seurat/03_seurat_anno_t_p65.Rds"))
se.braun = readRDS(paste0(manifest$braun$work, "seurat/03_seurat_anno_t.Rds"))
se.ho = readRDS(paste0(manifest$ho$work, "seurat/03_seurat_anno_t.Rds"))
se.hoso = readRDS(paste0(manifest$hosoya$work, "seurat/03_seurat_anno_t.Rds"))
se.hoso$orig.ident[se.hoso$orig.ident == "Duodenum"] = "GSM9314767"
se.hoso$orig.ident[se.hoso$orig.ident == "Stomach"] = "GSM9314768"
se.hoso$orig.ident[se.hoso$orig.ident == "Ileum"] = "GSM9314769"
# se.alem = readRDS(paste0(manifest$aleman$work, "seurat/03_seurat_anno_t.Rds"))

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
se.ho = celltype_coarse(se.ho) # Ho
se.hoso = celltype_coarse(se.hoso) # Hosoya
# se.alem = celltype_coarse(se.alem) # Aleman

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
    "Aph", "LP", "2w", "1mo", "2mo", "4mo", "7mo", "8mo",
    "12mo", "13mo", "14-128d"
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
se.ho = add_clin_metadata(obj = se.ho, pd = seurat.pd)
se.hoso = add_clin_metadata(obj = se.hoso, pd = seurat.pd)
# se.alem = add_clin_metadata(obj = se.alem, pd = seurat.pd)

smpls.lbls = setNames(
  paste0(
    seurat.pd$PATIENT_ID, "|", seurat.pd$SOURCE, "\n",
    seurat.pd$TIME_AFTER_CAR
  ),
  seurat.pd$SAMPLE_NAME
)
smpls.lbls = gsub(" NA", "", smpls.lbls)

smpls.lbls.2 = setNames(
  paste0(seurat.pd$SOURCE, "\n", seurat.pd$TIME_AFTER_CAR),
  seurat.pd$SAMPLE_NAME
)
smpls.lbls.2 = gsub("Duodenum", "Duo.", smpls.lbls.2)
smpls.lbls.2 = gsub("Stomach", "Stom.", smpls.lbls.2)

se.ho = se.ho[, se.ho$GROUP == "L-IEC"]
se.ho@meta.data = droplevels(se.ho@meta.data)

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
  se.obj
}

se.merz = pre_process(se.merz)
se.braun = pre_process(se.braun)
se.ho = pre_process(se.ho, parse_vdj = F)
se.hoso = pre_process(se.hoso, parse_vdj = F)
# se.alem = pre_process(se.alem)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Get Barcodes with Ser51Asn muatation
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
variants.path = paste0(
  manifest$work, "project_s51n/DWMCT/Analysis/detectVariants_out/",
  "detectVariants_out_2026_07_21/"
)

read_variants = function(path = NULL){
  variants_files = list.files(
    path, full.names = T, recursive = T, pattern = "cb_mutation.tsv"
  )
  res = parallel::mclapply(variants_files,function(x){
    print(x)
    df = readr::read_tsv(x, id = "path") %>% as.data.frame()
    if(nrow(df) == 0){return()}
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
variants = variants[variants$mut == "G222A", ]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
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
  # obj = obj[, obj$GROUP != "CAR+"]
  # obj = obj[
  #   , !(obj$CAR_BY_EXPRS == "FALSE" & obj$GROUP == "CAR+ S51N+")
  # ]
  obj
}

se.merz = add_mut(obj = se.merz)
se.braun = add_mut(obj = se.braun)
se.ho = add_mut(obj = se.ho)
se.hoso = add_mut(obj = se.hoso)
# se.alem = add_mut(obj = se.alem)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Merge
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
s.l = list(a = se.merz, b = se.braun, d = se.ho, e = se.hoso)
se.meta = merge(s.l[[1]], y = s.l[2:length(s.l)])
se.meta@meta.data$orig.ident = factor(se.meta@meta.data$orig.ident)
se.meta[["RNA"]] <- JoinLayers(se.meta[["RNA"]])

se.meta = add_mut(obj = se.meta)
se.meta = se.meta[ , se.meta$MUTATION %in% c("CAR+ S51N-", "CAR+ S51N+")]

se.meta = integration(
  obj = DietSeurat(se.meta, layers = "counts"),
  no.ftrs = 1000,
  threads = 30,
  .nbr.dims = 15,
  run.integration = T,
  harmony.group.vars = c("orig.ident"),
  do.cluster = F,
  min.dist = .5
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DimReduc
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cust.theme = theme(
  legend.position = "right",
  legend.margin = margin(l=-3),
  # legend.key.spacing.y= unit(-5, "pt"),
  legend.text = element_text(margin = margin(l = 0, unit = "pt"), size = rel(6.5/8)),
  legend.justification = c(0,.5),
  legend.title = element_text(margin = margin(b = 1)),
  legend.key.spacing.y = unit(-5, "pt"),
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  plot.title = element_blank()
)

pd = get_metadata(se.meta)

reduc.mut.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = MUTATION)) +
  scattermore::geom_scattermore(pointsize = 12, color="black")+
  scattermore::geom_scattermore(pointsize = 10.5, color="white") +
  ggrastr::geom_jitter_rast(
    data = pd[pd$MUTATION == "CAR+ S51N-", ],
    shape = ".", raster.dpi = 300, scale = .5
  ) +
  ggrastr::geom_jitter_rast(
    data = pd[pd$MUTATION == "CAR+ S51N+", ],
    shape = ".", raster.dpi = 300, scale = .5
  ) +
  mytheme(base_size = base.size) +
  cust.theme +
  guides(colour = guide_legend(
    title = NULL,  override.aes = list(shape = 16, size = 3)
  )) +
  scale_color_manual(
    values = c("CAR+ S51N-" = "#CCCCCC", "CAR+ S51N+" = "#4A7BB7"),
    na.value = "#FFFFFF"
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Merz | Patient 1
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Alluvial ---------------------------------------------------------------------
con.df = clonalCompare(
  se.merz,
  top.clones = 25,
  group.by = "orig.ident",
  cloneCall = "CLONE_PSEUDO_ID",
  samples = c("P65_1_PB_LP", "P65_1_PB_L", "P65_1_PB_VL", "P065-PB"),
  order.by = c("P65_1_PB_LP", "P65_1_PB_L", "P65_1_PB_VL", "P065-PB"),
  graph = "alluvial",
  relabel.clones = T,
  exportTable = T
)
con.df$clones = gsub("one:", "", con.df$clones)
t = naturalsort(unique(con.df$clones))
con.df$clones = factor(con.df$clones, levels = t)

con.df$Sample = factor(
  con.df$Sample, levels = c("P65_1_PB_LP", "P65_1_PB_L", "P65_1_PB_VL", "P065-PB")
)

cl.c = c(
  "Cl 1" = "#004488", "Cl 2" = "#0077BB", "Cl 3" = "#88CCEE", "Other clones" = "#BBBBBB"
)
cl.lbls = c(
  "Cl 1" = "Cl 1.1", "Cl 2" = "Cl 1.2", "Cl 3" = "Cl 1.3", "Other clones" = "#BBBBBB"
)

top.clones.pl =
  ggplot(con.df, aes(
    x = Sample, fill = clones, group = clones,  stratum = clones,
    alluvium = clones, y = Proportion,  label = clones)
  ) +
  scale_fill_manual(
    values = cl.c, na.value = "#BBBBBB", labels = cl.lbls
  ) +
  theme(
    axis.title.x = element_blank(),
    # axis.text.x = element_text(size = 7),
    panel.border = element_rect(colour = NA),
    axis.line = element_line(colour="black"),
    # legend.key.size = unit(3,"mm"),
    legend.key.spacing.y = unit(3, "pt"),
    legend.margin = margin(l = -5),
    plot.title = element_text(hjust = 0.5)
  ) +
  geom_stratum(linewidth = .1) +
  geom_flow(stat = "alluvium") +
  labs(fill = "Top\nclones") +
  ylab("% of cells") +
  scale_x_discrete(labels = smpls.lbls.2)


# Barplot | Ser51Asn -----------------------------------------------------------
se.merz$CLONE_PSEUDO_ID_TOPS = ifelse(
  se.merz$CLONE_PSEUDO_ID %in% c("Cl 1", "Cl 2", "Cl 3"),
  se.merz$CLONE_PSEUDO_ID,
  "Other clones"
)

pd = se.merz@meta.data

pd = subset(pd, MUTATION != "CAR+")
pd = subset(pd, MUTATION != "CAR-")
pd = droplevels(subset(pd, TIME_AFTER_CAR == "12mo" | TIME_AFTER_CAR == "4mo"))

t = prop.table(table(pd$CLONE_PSEUDO_ID_TOPS, pd$MUTATION), margin = 1) * 100
df = data.frame(t)
df$Var1 = gsub("Other clones", "Other\nclones", df$Var1)

brpl.mut =
  ggplot(df, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat="identity") +
  geom_text(
    aes(label = ifelse(Freq > 10, paste0(round(Freq, 1), "%"), '')),
    position = position_stack(vjust = 0.5), size = 2
  ) +
  ggtitle("P1") +
  ylab("% of cells") + xlab(NULL) + labs(fill = NULL) +
  scale_fill_manual(
    values = c("CAR+ S51N-" = "#BBBBBB", "CAR+ S51N+" = "#4A7BB7")
  ) +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  scale_x_discrete(labels = cl.lbls) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.key.size = unit(2.5,"mm"),
    legend.key.spacing.y = unit(3, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(margin = margin(l = 2)),
    legend.margin = margin(l = 0, r = 42),
    legend.position = "bottom"
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Braun | Patient 2
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Alluvial ---------------------------------------------------------------------
con.df = clonalCompare(
  se.braun,
  top.clones = 25,
  group.by = "orig.ident",
  cloneCall = "CLONE_PSEUDO_ID",
  samples = c("AphNB", "P2248", "P2249", "P2264"),
  order.by = c("AphNB", "P2248", "P2249", "P2264"),
  graph = "alluvial",
  relabel.clones = T,
  exportTable = T
)
con.df$clones = gsub("one:", "", con.df$clones)
t = naturalsort(unique(con.df$clones))
con.df$clones = factor(con.df$clones, levels = t)

con.df$clones = forcats::fct_relevel(con.df$clones, "Cl 3", after = 1)

con.df$Sample = factor(
  con.df$Sample, levels = c("AphNB", "P2248", "P2249", "P2264")
)

cl.c = c(
  "Cl 1" = "#555555", "Cl 2" = "#004488",
  "Cl 3" = "#000000", "Other clones" = "#BBBBBB"
)

cl.lbls = c(
  "Cl 1" = "Cl 1.1", "Cl 2" = "Cl 2", "Cl 3" = "Cl 1.2", "Other clones" = "#BBBBBB"
)

top.clones.braun.pl =
  ggplot(con.df, aes(
    x = Sample, fill = clones, group = clones,  stratum = clones,
    alluvium = clones, y = Proportion,  label = clones)
  ) +
  scale_fill_manual(values = cl.c, na.value = "#BBBBBB", labels = cl.lbls) +
  theme(
    axis.title.x = element_blank(),
    # axis.text.x = element_text(size = 7),
    panel.border = element_rect(colour = NA),
    axis.line = element_line(colour="black"),
    # legend.key.size = unit(3,"mm"),
    legend.key.spacing.y = unit(3, "pt"),
    legend.margin = margin(l = -5),
    plot.title = element_text(hjust = 0.5)
  ) +
  geom_stratum(linewidth = .1) +
  geom_flow(stat = "alluvium") +
  labs(fill = "Top\nclones") +
  ylab("% of cells") +
  scale_x_discrete(labels = smpls.lbls.2)

###

# Barplot | Ser51Asn -----------------------------------------------------------
se.braun$CLONE_PSEUDO_ID_TOPS = ifelse(
  se.braun$CLONE_PSEUDO_ID %in% c("Cl 1", "Cl 2", "Cl 3"),
  se.braun$CLONE_PSEUDO_ID,
  "Other clones"
)

pd = se.braun@meta.data

pd = subset(pd, MUTATION != "CAR+")
pd = subset(pd, MUTATION != "CAR-")

t = prop.table(table(pd$CLONE_PSEUDO_ID_TOPS, pd$MUTATION), margin = 1) * 100
df = data.frame(t)
df$Var1 = gsub("Other clones", "Other\nclones", df$Var1)
df$Var1 = factor(df$Var1, levels = c("Cl 1", "Cl 3", "Cl 2", "Other\nclones"))

brpl.braun.mut =
  ggplot(df, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat="identity") +
  geom_text(
    aes(label = ifelse(Freq > 10, paste0(round(Freq, 1), "%"), '')),
    position = position_stack(vjust = 0.5), size = 2
  ) +
  scale_fill_brewer(palette="Paired") +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  scale_x_discrete(labels = cl.lbls) +
  ggtitle("P2") +
  ylab("% of cells") + xlab(NULL) + labs(fill = NULL) +
  scale_fill_manual(
    values = c("CAR+ S51N-" = "#DDDDDD", "CAR+ S51N+" = "#4A7BB7")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    # legend.key.size = unit(3,"mm"),
    legend.key.spacing.y = unit(3, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(margin = margin(l = 2)),
    legend.margin = margin(l = 0, r = 42),
    legend.position = "bottom"
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Hosoya | Patient 3
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Barplot | Ser51Asn -----------------------------------------------------------
pd = se.hoso@meta.data
pd = subset(pd, MUTATION != "CAR+")
pd = subset(pd, MUTATION != "CAR-")
t = prop.table(table(pd$orig.ident, pd$MUTATION), margin = 1) * 100
df = data.frame(t)
df$Var1 = factor(df$Var1, levels = c("GSM9314768", "GSM9314767", "GSM9314769"))

brpl.hoso.mut =
ggplot(df, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat="identity") +
  geom_text(
    aes(label = ifelse(Freq > 10, paste0(round(Freq, 1), "%"), '')),
    position = position_stack(vjust = 0.5), size = 2
  ) +
  ggtitle("P3") +
  ylab("% of cells") + xlab(NULL) + labs(fill = NULL) +
  scale_fill_manual(
    values = c("CAR+ S51N-" = "#DDDDDD", "CAR+ S51N+" = "#4A7BB7")
  ) +
  scale_x_discrete(labels = smpls.lbls.2) +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  theme(
    plot.title = element_text(hjust = 0.5),
    # legend.key.size = unit(3,"mm"),
    legend.key.spacing.y = unit(3, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(margin = margin(l = 2)),
    legend.position = "bottom"
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Ho | Patient 4
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Barplot | Ser51Asn -----------------------------------------------------------
pd = se.ho@meta.data
pd = subset(pd, MUTATION != "CAR+")
pd = subset(pd, MUTATION != "CAR-")
t = prop.table(table(pd$orig.ident, pd$MUTATION), margin = 1) * 100
df = data.frame(t)

brpl.ho.mut =
  ggplot(df, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat="identity") +
  geom_text(
    aes(label = ifelse(Freq > 10, paste0(round(Freq, 1), "%"), '')),
    position = position_stack(vjust = 0.5), size = 2
  ) +
  ggtitle("P4") +
  ylab("% of cells") + xlab(NULL) + labs(fill = NULL) +
  scale_fill_manual(
    values = c("CAR+ S51N-" = "#DDDDDD", "CAR+ S51N+" = "#4A7BB7")
  ) +
  scale_x_discrete(labels = smpls.lbls.2) +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  theme(
    plot.title = element_text(hjust = 0.5),
    # legend.key.size = unit(3,"mm"),
    legend.key.spacing.y = unit(3, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(margin = margin(l = 2)),
    legend.position = "none"
  )

brpl.ho.mut.zoom =
  ggplot(df, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat="identity") +
  geom_text(
    aes(label = ifelse(Freq > 0, paste0(round(Freq, 1), "%"), '')),
    position = position_stack(vjust = 0.5), size = 2
  ) +
  ggtitle("P4 - Zoomed") +
  ylab(NULL) + xlab(NULL) + labs(fill = NULL) +
  scale_fill_manual(
    values = c("CAR+ S51N-" = "#DDDDDD", "CAR+ S51N+" = "#4A7BB7")
  ) +
  scale_x_discrete(labels = smpls.lbls.2) +
  scale_y_continuous(breaks = c(0, 1, 2)) +
  theme(
    plot.title = element_text(hjust = 0.5),
    # legend.key.size = unit(3,"mm"),
    legend.key.spacing.y = unit(3, "pt"),
    legend.key.spacing.x = unit(5, "pt"),
    legend.text = element_text(margin = margin(l = 2)),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(0, 2.5))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.t.pb = merge(se.merz, y = list(se.braun, se.hoso))
se.t.pb[["RNA"]] <- JoinLayers(se.t.pb[["RNA"]])
se.t.pb$SAMPLE = se.t.pb$orig.ident
se.t.pb$celltype = se.t.pb$celltype_short_3

rownames(se.t.pb)[rownames(se.t.pb) == "PDCD1"] = "PD1"
rownames(se.t.pb)[rownames(se.t.pb) == "HAVCR2"] = "TIM3"
rownames(se.t.pb)[rownames(se.t.pb) == "MKI67"] = "Ki-67"
rownames(se.t.pb)[rownames(se.t.pb) == "SELL"] = "CD62L"

pseudo = Seurat::AggregateExpression(
  se.t.pb, assays = "RNA", return.seurat = T,
  group.by = c("orig.ident", "SAMPLE", "celltype", "MUTATION")
)
pseudo = pseudo[, pseudo$celltype == "CD8 T-Cell"]

pseudo = pseudo[
  , !pseudo$SAMPLE %in% c("GSM9314768", "P65-1-PB-LP", "P65-1-PB-L", "AphNB")
]

pseudo$MUTATION = gsub(
  "\\+ ", "plus", gsub("\\-", "neg", pseudo@meta.data$MUTATION)
)
pseudo$MUTATION = gsub(
  "\\+", "plus", gsub("\\-", "neg", pseudo@meta.data$MUTATION)
)
pseudo@meta.data = pseudo@meta.data %>% dplyr::mutate_if(is.character,as.factor)
pseudo$celltype = gsub(" T-Cell", "", pseudo$celltype)

ctrs.l = list(
  "CAR+ S51N+ vs CAR+ S51N-" = c("CARplusS51Nplus", "CARplusS51Nneg"),
  "CAR+ S51N+ vs CAR-" = c("CARplusS51Nplus", "CARneg"),
  "CAR+ S51N- vs CAR-" = c("CARplusS51Nneg", "CARneg"),
  "CAR+ S51N- vs CAR+" = c("CARplusS51Nneg", "CARplus")
)
dg.res = parallel::mclapply(names(ctrs.l), function(ctrst){
  res = dgea_gex(
    obj = pseudo,
    group = "MUTATION",
    ctrs.grp1 = ctrs.l[[ctrst]][1], ctrs.grp2 = ctrs.l[[ctrst]][2],
    .target = "celltype",
    deseq_paired = T, dsgn = "~SAMPLE + MUTATION",
    logfc.threshold = log2(1.25), min.pct = .8,
    .min.cells = 6
  )
  res$res.dgea$CNTRST = ctrst
  res$res.dgea.sign$CNTRST = ctrst
  res
}, mc.cores = length(ctrs.l))
names(dg.res) = unlist(lapply(dg.res, function(x){x$res.dgea$CNTRST[1]}))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DE | Volcano
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
d = dg.res$`CAR+ S51N+ vs CAR+ S51N-`$res.dgea
if(d[1, ]$significant == "FALSE"){
  dg.res$`CAR+ S51N+ vs CAR+ S51N-`$res.dgea = d[2:nrow(d), ]
}
vlc.pl =
  dgea_volcano(
    dgea.res = dg.res$`CAR+ S51N+ vs CAR+ S51N-`$res.dgea,
    dgea.res.sig = dg.res$`CAR+ S51N+ vs CAR+ S51N-`$res.dgea.sign,
    nbr.tops = 10,
    x.axis.sym = T,
    logFC.column = "avg_log2FC",
    label.size = 1.5,
    x.axis.ext = 1.25,
    box.padding = .3,
    label.padding = .1,
    nudge_x = 4.75,
    pt.size = .25,
    geom.hline = log2(1.25),
    facet.scales = "free_x",
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.box.spacing = unit(1, 'pt'),
    # axis.text = element_text(size = 6),
    axis.title.y = element_text(vjust = + 4)
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DE | Selected DE Genes
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
df = dg.res$`CAR+ S51N+ vs CAR+ S51N-`$res.dgea.sign
m = t.mrkr[t.mrkr$Gene_symbol %in% df$feature, ]
df = df[df$feature %in% m$Gene_symbol, ]
df$Type = t.mrkr$Group[match(df$feature, t.mrkr$Gene_symbol)]

df = df[
  df$feature %in% c("PRF1", "GZMB", "PD1", "TIM3", "Ki-67", "PCNA", "CCR7", "CD62L"),
]

bxpl.lbls = c(
  "CARplusS51Nplus" = "CAR+ S51N+",
  "CARplusS51Nneg" = "CAR+ S51N-"
)

de.expl =
  bxpl_dgea_paired(
    obj = pseudo[, pseudo$MUTATION %in% c("CARplusS51Nplus", "CARplusS51Nneg")],
    ftr = df$feature,
    # ftr = c("TIGIT", "PDCD1", "LAG3", "HAVCR2"),
    facet.ncol = 8,
    facet.scale = "free_y",
    add.facet.groups = t.mrkr,
    panel.space = .5,
    .base.size = 12
  ) +
  scale_fill_manual(values = c("#DDDDDD", "#4A7BB7"), labels = bxpl.lbls) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 4),
    legend.key.size = unit(base_size * 1.2, "pt"),
    #    legend.box.spacing = unit(10, "pt")
  ) +
  ggtitle("")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ggsave2(
  filename="figures/paper/figure_3.png",
  plot_grid(
    plot_grid(
      plot_grid(NULL, reduc.mut.pl, NULL, nrow = 3, rel_heights = c(.05, 1, .1)),
      NULL,
      plot_grid(
        ggdraw() + draw_label(
          "Clonality", hjust = 0.5, size = 8, fontface = "plain"
        ),
        plot_grid(
          top.clones.pl + ggtitle("P1"),
          NULL,
          top.clones.braun.pl + ggtitle("P2"),
          ncol = 3, rel_widths = c(1, .05, 1)
        ),
        NULL,
        nrow = 3, rel_heights = c(.1, 1, .11)
      ),
      labels = c("A", "", "B"),
      label_fontface = "plain", label_size = 12,
      ncol = 3, rel_widths = c(1, .1, 2)
    ),
    NULL,
    plot_grid(
      ggdraw() + draw_label(
        "Mutational composition", hjust = 0.5, size = 8, fontface = "plain"
      ),
      plot_grid(
        brpl.mut + theme(legend.position = "none"),
        NULL,
        brpl.braun.mut + theme(legend.position = "none"),
        NULL,
        brpl.hoso.mut + theme(legend.position = "none"),
        NULL,
        plot_grid(brpl.ho.mut, brpl.ho.mut.zoom, rel_widths = c(1.25, 1)),
        ncol = 7, rel_widths = c(1, .05, 1.2, .05, 1, .05, 1.3),
        # align = "vh",
        labels = c("C", "", "", "", "", "", ""),
        label_fontface = "plain", label_size = 12, vjust = 0
      ),
      ggdraw(
        get_legend(brpl.mut + theme(legend.box.margin=margin(0,0,0,60)))
      ),
      nrow = 3, rel_heights = c(.1, 1, .1)
    ),
    NULL,
    plot_grid(
      plot_grid(
        NULL,
        plot_grid(NULL, vlc.pl, ncol = 2, rel_widths = c(.025, 1)),
        nrow = 2, rel_heights = c(.005, 1)
      ),
      NULL,
      plot_grid(NULL, de.expl, NULL, nrow = 3, rel_heights = c(-.075, 1, .065)),
      ncol = 3, rel_widths = c(1, .1, 2.43),
      labels = c("D", "", "E"),
      label_fontface = "plain", label_size = 12, vjust = .5
    ),
    nrow = 5, rel_heights = c(1, .05, 1, .1, 1.2)
  ),
  width = 180, height = 150, dpi = 300, bg = "white", units = "mm", scale = 1,
  device = png, type = "cairo"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Supp | DGEA for all 4 contrasts
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
theme_set(mytheme(base_size = 12))

dg.meta = rbind(
  dg.res$`CAR+ S51N+ vs CAR+ S51N-`$res.dgea %>%
    dplyr::mutate(cluster = "CAR+ S51N+ vs CAR+ S51N-"),
  dg.res$`CAR+ S51N+ vs CAR-`$res.dgea %>%
    dplyr::mutate(cluster = "CAR+ S51N+ vs CAR-"),
  dg.res$`CAR+ S51N- vs CAR-`$res.dgea %>%
    dplyr::mutate(cluster = "CAR+ S51N- vs CAR-"),
  dg.res$`CAR+ S51N- vs CAR+`$res.dgea %>%
    dplyr::mutate(cluster = "CAR+ S51N- vs CAR+")
)

dg.sign.meta = rbind(
  dg.res$`CAR+ S51N+ vs CAR+ S51N-`$res.dgea.sign %>%
    dplyr::mutate(cluster = "CAR+ S51N+ vs CAR+ S51N-"),
  dg.res$`CAR+ S51N+ vs CAR-`$res.dgea.sign %>%
    dplyr::mutate(cluster = "CAR+ S51N+ vs CAR-"),
  dg.res$`CAR+ S51N- vs CAR-`$res.dgea.sign %>%
    dplyr::mutate(cluster = "CAR+ S51N- vs CAR-"),
  dg.res$`CAR+ S51N- vs CAR+`$res.dgea.sign %>%
    dplyr::mutate(cluster = "CAR+ S51N- vs CAR+")
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DE | Volcano
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
vlc.pl =
  dgea_volcano(
    dgea.res = rbind(
      dg.res$`CAR+ S51N+ vs CAR-`$res.dgea %>%
        dplyr::mutate(cluster = "CAR+ S51N+ vs CAR-"),
      dg.res$`CAR+ S51N- vs CAR-`$res.dgea %>%
        dplyr::mutate(cluster = "CAR+ S51N- vs CAR-"),
      dg.res$`CAR+ S51N- vs CAR+`$res.dgea %>%
        dplyr::mutate(cluster = "CAR+ S51N- vs CAR+")
    ),
    dgea.res.sig = rbind(
      dg.res$`CAR+ S51N+ vs CAR-`$res.dgea.sign %>%
        dplyr::mutate(cluster = "CAR+ S51N+ vs CAR-"),
      dg.res$`CAR+ S51N- vs CAR-`$res.dgea.sign %>%
        dplyr::mutate(cluster = "CAR+ S51N- vs CAR-"),
      dg.res$`CAR+ S51N- vs CAR+`$res.dgea.sign %>%
        dplyr::mutate(cluster = "CAR+ S51N- vs CAR+")
    ),
    nbr.tops = 10,
    x.axis.sym = T,
    logFC.column = "avg_log2FC",
    label.size = 4,
    x.axis.ext = -8,
    box.padding = 0.35,
    label.padding = .1,
    nudge_x = 5,
    geom.hline = log2(1.25),
    facet.scales = "free_x", cut.y.thresh = 1e-40
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.box.spacing = unit(1, 'pt')
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# T cell related DE genes
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
df = dg.meta

m = t.mrkr[t.mrkr$Gene_symbol %in% df[df$significant == "TRUE", ]$feature, ]
df = df[df$feature %in% m$Gene_symbol, ]
df$Type = t.mrkr$Group[match(df$feature, t.mrkr$Gene_symbol)]
df = df[df$Type != "T cell identity", ]
r = colSums(table(df$feature, df$Type) > 0)
df = df[!df$Type %in% names(r[r == 1]), ]

df$Type[grepl("Resident", df$Type)] = "Resident,\nhoming"

max.v = max(abs(df$avg_log2FC))
t.mrkr.tile.pl =
  ggplot(df, aes(feature, CNTRST,fill = avg_log2FC)) +
  geom_tile(color = "white") +
  geom_point(
    data = df[df$significant == "TRUE", ], aes(size = p_val_adj),
    fill = "white", pch=21, show.legend = T
  ) +
  scale_size(range = c(3, 1)) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1, size = 10),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.ticks.length = unit(0.075, 'cm'),
    legend.title = element_text(margin = margin(r = 3, l = 5, unit = "pt")),
    legend.margin = margin(t = 0),
    legend.text = element_text(margin = margin(l = -1, t = 2, b = 2)),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.box.spacing = unit(1, 'pt')
  ) +
  scico::scale_fill_scico(
    palette = "vik", midpoint = 0, direction = 1, limits = c(-max.v, max.v),
    na.value="white"
  ) +
  guides(
    fill = guide_colorbar(
      title = "Log2 Fold Change", title.vjust = 1.1,
      barwidth = unit(6.5, 'lines'), barheight = unit(.4, 'lines'),
      ticks.linewidth = 1/.pt, ticks = T, frame.colour="black",
      frame.linewidth = 0.5/.pt
    ),
    size = guide_legend(title = "FDR <0.05")
  ) +
  facet_grid( ~ Type, scale = "free", space = "free")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ORA | GO
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
df.ora = dg.sign.meta

ftrs.l = df.ora
ftrs.l = split(ftrs.l, ftrs.l$CNTRST)
ftrs.l = lapply(ftrs.l, function(x){
  ftrs = x$feature
  names(ftrs) = x$avg_log2FC
  ftrs
})
lengths(ftrs.l)

ora.t = parallel::mclapply(ftrs.l, function(x){
  run_nmf_ora(genes = x, universe = rownames(pseudo), category = "CD8")
}, mc.cores = 1)

ora.pl =
  ora_bubble(
    gsea.res = ora.t,
    ftrs.list = ftrs.l,
    nbr.tops = 10,
    min.genes = 5,
    dot.range = c(1, 5.5),
    term.length = 100,
    sort.by.padj = F,
    facet.spit = F,
    export.data = F,
    font.size = 12,
    font.size.y = 12,
  ) +
  theme(
    legend.position = "bottom",
    legend.box.spacing = unit(-1, "pt"),
    legend.box.margin = margin(0, 250, 0, 0)
  ) +
  guides(
    fill = guide_colorbar(
      title =  "Pathway\ndirection",
      barheight = unit(.35, 'lines'),
      barwidth = unit(5, 'lines'),
      order = 1, ticks.linewidth = .75/.pt,  frame.linewidth = 0.5/.pt
    ),
    size = guide_legend(title = "-Log10(FDR)", order = 2)
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ggsave2(
  filename="figures/paper/figure_3_supps_3.png",
  plot_grid(
    vlc.pl,
    NULL,
    plot_grid(
      plot_grid(
        NULL, t.mrkr.tile.pl, NULL, nrow = 3, rel_heights = c(-.125, 1.1, 1)
      ),
      NULL,
      # NULL,
      plot_grid(NULL, ora.pl, nrow = 2, rel_heights = c(-.1, 1)),
      ncol = 3, rel_widths = c(2.7, .05, 1),
      labels = c("D", "", "E"),
      label_fontface = "bold", label_size = 16, vjust = -.9
    ),
    nrow = 3, rel_heights = c(1, .1, 1)
  )
  ,
  width = 180, height = 135, dpi = 300, bg = "white", units = "mm", scale = 2,
  device = png, type = "cairo"
)

theme_set(theme_custom())

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA | Supp Table
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
xlsx.filename = "figures/paper/DGEA_Supplementary_Tables.xlsx"
wb <- createWorkbook()
first.row.style = createStyle(fgFill = "lightgrey", textDecoration = "Bold")

sheet = "Table 1"
addWorksheet(wb, sheet)
writeData(
  wb, sheet, "CD8 T cells | DE genes comparing CARplusS51Nplus with CARplusS51Nneg",
  startCol = 1, startRow = 1
)
addStyle(
  wb, sheet, style = first.row.style, rows = 1, cols = 1:7
)
writeData(
  wb, sheet,
  dg.res$`CAR+ S51N+ vs CAR+ S51N-`$res.dgea.sign %>%
    dplyr::mutate(across(c('avg_log2FC'), round, 3)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(p_val_adj, .by_group = T) %>%
    dplyr::select(
      gene = feature, log2FC = avg_log2FC, p_val, p_val_adj,
      cell_type = cluster
    ) %>%
    data.frame(),
  startRow = 3, startCol = 1
)
saveWorkbook(wb, xlsx.filename, overwrite = T)
