# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310992
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries and some Functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "Seurat", "yaml", "dplyr", "stringr", "naturalsort", "data.table", "ggplot2",
  "scales", "cowplot", "scGate"
)
.bioc_packages = c("clustifyr", "scds", "UCell")

# Install CRAN packages (if not already installed)
.inst = .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}

# Install bioconductor packages (if not already installed)
.inst <- .bioc_packages %in% installed.packages()
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

if (any(!"ProjecTILs" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github("carmonalab/ProjecTILs")
}
library(ProjecTILs)

source("code/helper/functions.R")
source("code/helper/styles.R")
source("code/helper/functions_plots.R")
theme_set(mytheme(base_size = 12))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load objects and phenodata
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

output.file = paste0(manifest$hosoya$work, "seurat/02_seurat_pre.Rds")
output.file.2 = paste0(manifest$hosoya$work, "seurat/03_seurat_anno_t.Rds")

se.meta = readRDS(paste0(manifest$hosoya$seurat, "01_seurat_ori.Rds"))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Harmonize CAR construct gene name
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mat.counts = GetAssayData(object = se.meta, assay = "RNA", layer = "counts")
car.exprs = mat.counts["ciltacel-mined", , drop = F]
rownames(car.exprs) = "CAR-BCMA"

new.counts = rbind(
  car.exprs,
  mat.counts[!rownames(se.meta) %in% c("ciltacel-mined"), ]
)

se.meta <- CreateSeuratObject(counts = new.counts, meta.data = se.meta@meta.data)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Add column in metadata wether CD4/CD8, CAR are present
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta = cd4cd8_car_present(obj = se.meta)
print(table(se.meta$CAR_BY_EXPRS))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Normalize
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta = NormalizeData(
  se.meta, assay = 'RNA', normalization.method = "LogNormalize"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Add %MT, %Ribosomal and complexity values
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta[["Perc_of_mito_genes"]] = Seurat::PercentageFeatureSet(
  se.meta, pattern = "^MT-"
)
se.meta[["Perc_of_ribosomal_genes"]] = Seurat::PercentageFeatureSet(
  se.meta, pattern = "^RPL|^RPS"
)
se.meta@meta.data$log10GenesPerUMI = log10(se.meta$nFeature_RNA) / log10(se.meta$nCount_RNA)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cell filtering
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
nFeature_low_cutoff = 250
nFeature_high_cutoff = 8000
nCount_low_cutoff = 1000
nCount_high_cutoff = 100000
mt_cutoff = 65
complx_cutoff = .8

p1 = qc_vln_plot_cell(
  se.meta, low_cutoff = nFeature_low_cutoff, high_cutoff = nFeature_high_cutoff,
)
p2 = qc_vln_plot_cell(
  se.meta, .features = "nCount_RNA", plot_title = "UMIs per Cell",
  y_axis_label = "UMIs", low_cutoff = nCount_low_cutoff, high_cutoff = nCount_high_cutoff,
)
p3 = qc_vln_plot_cell(
  se.meta, .features = "Perc_of_mito_genes", plot_title = "Mito Gene % per Cell",
  y_axis_label = "% Mito Gene Counts", high_cutoff = mt_cutoff,
)
p4 = qc_vln_plot_cell(
  se.meta,  .features = "log10GenesPerUMI", plot_title = "Cell Complexity",
  y_axis_label = "log10(Genes) / log10(UMIs)", high_cutoff = complx_cutoff
)

ggsave2(
  filename="figures/hosoya_et_al/qc_pre_processing/stats_tech_per_cell.pdf",
  plot = cowplot::plot_grid(
    p1, p2, p3, p4, scale = .9, nrow = 1,
    labels = "AUTO", label_fontface = "bold", label_size = 14
  ),
  width = 160,
  height = 42,
  dpi = 100,
  bg = "white",
  units = "mm",
  scale = 2
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cell filtering
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
label_cells_rm = function(obj) {
  obj@meta.data = obj@meta.data %>% mutate(
    KEEP_CELL = case_when(
      (nFeature_RNA < nFeature_low_cutoff) | (nFeature_RNA > nFeature_high_cutoff) |
        (nCount_RNA < nCount_low_cutoff) | (nCount_RNA > nCount_high_cutoff) |
        (Perc_of_mito_genes > mt_cutoff)  | (log10GenesPerUMI < complx_cutoff) ~ FALSE,
      TRUE ~ TRUE
    )
  )
  print(table(obj$orig.ident, obj$KEEP_CELL))
  obj
}

cell.track = count_cells_per_sample(c(se.meta))
se.meta = label_cells_rm(se.meta)
se.meta = subset(se.meta, subset = KEEP_CELL == TRUE)
se.meta@meta.data = droplevels(se.meta@meta.data)
se.meta$KEEP_CELL = NULL
cell.track = count_cells_per_sample(c(se.meta), cell.track, "n1")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("scds")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
obj.l = Split_Object(se.meta, split.by = "orig.ident", threads = 4)
scds_res = parallel::mclapply(obj.l, function(se){
  scds_doublets(se)
}, mc.cores = 1)

scds_res = do.call("rbind", scds_res)
rownames(scds_res) = scds_res$barcode
scds_res$barcode = NULL
se.meta = AddMetaData(se.meta, scds_res)
se.meta$DOUBLETS_CONSENSUS = se.meta$scDblFinder_class == "doublet" &
  se.meta$hybrid_call == TRUE

cell.track = count_cells_per_sample(
  c(subset(se.meta, DOUBLETS_CONSENSUS == FALSE)), cell.track, "n2"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CellCycleScoring
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
source("data/signatures/cellCycleMarkers.R")
se.meta = CellCycleScoring(
  se.meta, s.features = s.genes, g2m.features = g2m.genes,
  assay = 'RNA', search = TRUE
)

obj.l = Split_Object(se.meta, split.by = "orig.ident", threads = 4)
cc.res = parallel::mclapply(obj.l, function(se){
  estimate_cc(se)
}, mc.cores = 4)
cc.res = do.call("rbind", cc.res)
rownames(cc.res) = cc.res$barcode
se.meta = AddMetaData(se.meta, cc.res)
se.meta$barcode = NULL

se.meta@meta.data = se.meta@meta.data %>%
  dplyr::mutate_if(is.character, as.factor)

# Annotate cells that were missed by the cluster approach
cutoff = 0.20
se.meta$CellCycle_Phase = dplyr::case_when(
  se.meta$S.Score > cutoff & se.meta$G2M.Score <= cutoff ~ "S",
  se.meta$G2M.Score > cutoff & se.meta$S.Score <= cutoff ~ "G2M",
  (se.meta$S.Score > cutoff & se.meta$G2M.Score > cutoff) &
    (se.meta$S.Score > se.meta$G2M.Score) ~ "S",
  (se.meta$S.Score > cutoff & se.meta$G2M.Score > cutoff) &
    (se.meta$G2M.Score > se.meta$S.Score) ~ "G2M",
  TRUE ~ se.meta$CellCycle_Phase
)
se.meta$CellCycle_Phase = factor(
  se.meta$CellCycle_Phase, levels = c("G1M", "S", "G2M")
)
se.meta@meta.data$CellCycle = factor(
  ifelse(se.meta$CellCycle_Phase == "G1M", "FALSE", "TRUE")
)

data(cell.cycle.obj) # ProjecTILs package
se.meta = AddModuleScore_UCell(
  se.meta,
  features = list("CellCycle_SCORE" = cell.cycle.obj$human$cycling),
  assay = "RNA", slot = "counts",
  ncores = 25, force.gc = T
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Save")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print(cell.track)

DefaultAssay(se.meta) = "RNA"
Idents(se.meta) = "orig.ident"

saveRDS(se.meta, output.file)
# se.meta = readRDS(output.file)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Split
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta.l = Split_Object(se.meta, split.by = "orig.ident", threads = 10)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("scGate_models")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
scGate_models_DB = get_scGateDB("data/signatures/scGateDB", force_update = F)

models.list.1 <- scGate_models_DB$human$PBMC[
  c("Platelet", "Erythrocyte", "NK", "PlasmaCell", "Bcell", "gdT")
]
models.list.2 <- scGate_models_DB$human$generic[
  c("Immune", "Myeloid", "MoMacDC", "Megakaryocyte", "Tcell", "Macrophage", "Mast")
]

sc_gating = function(obj, obj.l, model) {
  suppressWarnings({
    suppressMessages({
      # x = obj.l[[6]]
      obj.l = parallel::mclapply(obj.l, function(x) {
        x = scGate::scGate(
          x, model = model, assay = "RNA", slot = "data",
          output.col.name = "SCGATE", ncores = 5
        )
        df = x@meta.data[, grepl("^SCGATE", colnames(x@meta.data)), drop = F]
        colnames(df) = toupper(colnames(df))
        df$barcode = rownames(df)
        df
      }, mc.cores = length(obj.l))
      df = do.call("rbind", obj.l)
      rownames(df) = df$barcode
      df$barcode = NULL
      obj = AddMetaData(obj, df)
    })
  })
  obj
}
se.meta = sc_gating(obj = se.meta, obj.l = se.meta.l, model = models.list.1)
se.meta = sc_gating(obj = se.meta, obj.l = se.meta.l, model = models.list.2)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd = se.meta@meta.data; print(nrow(pd))

pd = pd[pd$DOUBLETS_CONSENSUS == FALSE, ]; print(nrow(pd))
pd = pd[pd$SCGATE_IMMUNE == "Pure", ]; print(nrow(pd))
pd = pd[pd$SCGATE_TCELL == "Pure", ]; print(nrow(pd))
pd = pd[pd$SCGATE_MACROPHAGE == "Impure", ]; print(nrow(pd))
pd = pd[pd$SCGATE_MAST == "Impure", ]; print(nrow(pd))
pd = pd[pd$SCGATE_PLATELET == "Impure", ]; print(nrow(pd))
pd = pd[pd$SCGATE_ERYTHROCYTE == "Impure", ]; print(nrow(pd))
pd = pd[pd$SCGATE_MYELOID == "Impure", ]; print(nrow(pd))
pd = pd[pd$SCGATE_MOMACDC == "Impure", ]; print(nrow(pd))
pd = pd[pd$SCGATE_MEGAKARYOCYTE == "Impure", ]; print(nrow(pd))
pd = pd[pd$SCGATE_BCELL == "Impure", ]; print(nrow(pd))
pd = pd[pd$SCGATE_PLASMACELL == "Impure", ]; print(nrow(pd))
pd = pd[pd$SCGATE_GDT == "Impure", ]; print(nrow(pd))
pd = pd[pd$SCGATE_NK == "Impure", ]; print(nrow(pd))
pd = droplevels(pd[pd$CD4CD8_BY_EXPRS != "CD4+CD8+", ]); print(nrow(pd))

se.t = se.meta[, rownames(se.meta@meta.data) %in% rownames(pd)]
se.t@meta.data = droplevels(se.t@meta.data)
table(se.t$CAR_BY_EXPRS, se.t$CD4CD8_BY_EXPRS)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("CD4/CD8 Imputation")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
exprs.smoothed.l = parallel::mclapply(
  as.character(unique(se.t$orig.ident)),function(x){

    se = subset(se.t, orig.ident == x)
    se = se %>%
      FindVariableFeatures(verbose = F) %>%
      ScaleData(verbose = F) %>%
      RunPCA(verbose = F)

    exprs.smoothed = UCell::SmoothKNN(
      obj = se,
      signature.names = c("CD4", "CD8A"),
      assay="RNA", reduction="pca",
      k = 10, suffix = "_CD4CD8_smooth"
    )
    DefaultAssay(exprs.smoothed) = "RNA_CD4CD8_smooth"
    exprs.smoothed = FetchData(
      exprs.smoothed, vars = c("CD4", "CD8A"), layer = "data"
    )
    colnames(exprs.smoothed) = c("CD4", "CD8")
    exprs.smoothed$GRP = se$CD4CD8_BY_EXPRS[
      match(rownames(exprs.smoothed), rownames(se@meta.data))
    ]
    exprs.smoothed$SAMPLE = as.character(se$orig.ident[1])
    exprs.smoothed

  }, mc.cores = 10)

exprs.smoothed = do.call("rbind", exprs.smoothed.l)
exprs.smoothed$T_LIN = se.t$CD4CD8_BY_EXPRS[
  match(rownames(exprs.smoothed), colnames(se.t))
]

exprs.smoothed = exprs.smoothed %>% dplyr::mutate(
  T_LIN = dplyr::case_when(
    T_LIN == "CD4-CD8+" ~ "CD8",
    T_LIN == "CD4+CD8-" ~ "CD4",
    T_LIN == "CD4-CD8-" ~ "T"
  )
)

thres.x = .2
thres.y = .3

ggplot(exprs.smoothed, aes(CD4, CD8)) +
  geom_point(size = .1, alpha = .7) +
  guides(alpha = 'none') +
  scico::scale_color_scico(palette = "acton", direction = -1) +
  mytheme() +
  theme(
    aspect.ratio = 1,
    panel.spacing = unit(1, "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
  ) +
  geom_hline(yintercept = thres.y, lwd = .2, linetype = "dashed") +
  geom_vline(xintercept = thres.x, lwd = .2, linetype = "dashed") +
  facet_wrap( ~ T_LIN)

# Annotate T cell as CD4/8 cell based on knn smoothing
df = exprs.smoothed %>% mutate(
  T_LIN_WORK = case_when(
    (CD4 > thres.x & CD8 < thres.y) & T_LIN == "T" ~ "CD4",
    (CD4 < thres.x & CD8 > thres.y)  & T_LIN == "T" ~ "CD8",
    TRUE ~ T_LIN
  )
)
table(df$T_LIN_WORK, df$T_LIN)

add.meta = df %>% dplyr::select(
  T_LIN_WORK, CD4_SMOOTH = CD4, CD8_SMOOTH = CD8
)
se.t = AddMetaData(se.t, add.meta)
# se.t = subset(se.t, T_LIN_WORK != "T")

se.t@meta.data$T_LIN_WORK = ifelse(
  se.t@meta.data$T_LIN_WORK == "T",
  "dnT", se.t@meta.data$T_LIN_WORK
)

table(se.t@meta.data$T_LIN_WORK, useNA = "always")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("ProjecTILs (T-cells)")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print(paste0("Nbr of T-Cells: ", ncol(se.t)))
ref.maps <- get.reference.maps(
  collection = "human", directory = paste0(manifest$work, "references/atlases/"),
  update = FALSE
)
ref.cd8 <- ref.maps$human$CD8
ref.cd4 <- ref.maps$human$CD4

ref.cd8 = ref.cd8[, !grepl("EX$", ref.cd8$functional.cluster)]

# Classify CD8 T subtypes
se.cd8 <- ProjecTILs.classifier(
  subset(se.t, T_LIN_WORK == "CD8"), ref.cd8, ncores = 20,
  split.by = "orig.ident", filter.cells = FALSE, min.confidence = 0
)

# Classify CD4 T subtypes
se.cd4 <- ProjecTILs.classifier(
  subset(se.t, T_LIN_WORK == "CD4"), ref.cd4, ncores = 20,
  split.by = "orig.ident", filter.cells = FALSE, min.confidence = 0
)

se.tmp = merge(se.cd4, se.cd8)
se.tmp = se.tmp[, !is.na(se.tmp$functional.cluster)]
se.t = merge(se.tmp, se.t[, se.t$T_LIN_WORK == "dnT" | se.t$T_LIN_WORK == "gdT"])
se.t[["RNA"]] <- JoinLayers(se.t[["RNA"]])
# se.t = se.tmp

colnames(se.t@meta.data)[
  colnames(se.t@meta.data) == "functional.cluster"
] = "SPICA_TCELL"
colnames(se.t@meta.data)[
  colnames(se.t@meta.data) == "functional.cluster.conf"
] = "SPICA_TCELL_CONF"

se.t@meta.data = se.t@meta.data[, !colnames(se.t@meta.data) %in% c(
  "CTgene", "CTnt", "CTaa", "CTstrict", "clonalProportion",
  "clonalFrequency", "cloneSize", "clonePseudoID",
  colnames(se.t@meta.data)[grepl("SCGATE", colnames(se.t@meta.data))],
  colnames(se.t@meta.data)[grepl("CD4_SMOOTH|CD8_SMOOTH", colnames(se.t@meta.data))]
)]

se.t@meta.data = se.t@meta.data %>% mutate(
  celltype = case_when(
    grepl("^CD4", SPICA_TCELL) & CellCycle == T ~ "CD4.Cycling",
    grepl("^CD8", SPICA_TCELL) & CellCycle == T ~ "CD8.Cycling",
    SPICA_TCELL == "CD4.CTL_EOMES" ~ "CD4.CTL EOMES+",
    SPICA_TCELL == "CD4.CTL_GNLY" ~ "CD4.CTL GNLY+",
    T_LIN_WORK == "dnT" ~ "dnT",
    TRUE ~ as.character(SPICA_TCELL)
  )
)
se.t@meta.data$celltype = gsub("CD4.", "CD4 ", se.t@meta.data$celltype)
se.t@meta.data$celltype = gsub("CD8.", "CD8 ", se.t@meta.data$celltype)

se.t@meta.data = se.t@meta.data %>%
  mutate(
    celltype_short_3 = case_when(
      grepl("^CD4", celltype) ~ "CD4 T-Cell",
      grepl("^CD8", celltype) ~ "CD8 T-Cell",
      grepl("gdT", celltype) ~ "gd T-Cell",
      grepl("dpT", celltype) ~ "dp T-Cell",
      TRUE ~ celltype
    )
  )

se.t@meta.data = droplevels(se.t@meta.data)
se.t@meta.data$celltype = factor(se.t@meta.data$celltype)
se.t@meta.data$celltype_short_3 = factor(se.t@meta.data$celltype_short_3)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Save")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DefaultAssay(se.t) = "RNA"
Idents(se.t) = "orig.ident"
saveRDS(se.t, output.file.2)

