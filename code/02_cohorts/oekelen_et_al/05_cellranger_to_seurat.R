# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182527

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("Seurat", "yaml", "dplyr", "doParallel", "parallel", "data.table")
.bioc_packages = c("SingleCellExperiment", "scDblFinder", "CiteFuse")

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

source("code/helper/functions.R")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load Rawcounts and create a merged Seurat object
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")
work.path = paste0(manifest$oekelen$work, "cellranger/")
seurat.path = paste0(manifest$oekelen$work, "seurat/")
dir.create(seurat.path, recursive = T)

obj.path = paste0(paste0(seurat.path, "01_seurat_ori.Rds"))

cellranger.dirs = list.dirs(path = work.path, full.names = T, recursive = T)

fltrd.dirs = cellranger.dirs[
  grepl("sample_filtered_feature_bc_matrix", cellranger.dirs)
]
tmp = gsub(work.path, "", fltrd.dirs)
tmp = gsub("/", "", gsub("out.*", "", tmp))
tmp = gsub("multi_", "", tmp)
names(fltrd.dirs) = tmp
print(length(names(fltrd.dirs)))

fltrd.counts = Read10X(data.dir = fltrd.dirs, gene.column = 2)
seurat = CreateSeuratObject(counts = fltrd.counts, project = "oekelen")
seurat$barcode = gsub(".+_", "", rownames(seurat@meta.data))
DefaultAssay(seurat) = "RNA"
Idents(seurat) = "orig.ident"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Demultiplexing
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mat = fltrd.counts = Read10X(
  data.dir = paste0(manifest$oekelen$work, "filtered_feature_bc_matrix"),
  gene.column = 2
)

se_hto = CreateSeuratObject(counts = mat[[2]][grepl("HTO", rownames(mat[[2]])), ])

sce_citeseq = preprocessing(
  list(
    RNA = mat[[1]],
    ADT = mat[[2]][!grepl("HTO", rownames(mat[[2]])), ],
    HTO = mat[[2]][grepl("HTO", rownames(mat[[2]])), ]
  )
)
sce_citeseq = normaliseExprs(
  sce = sce_citeseq, altExp_name = "HTO", transform = "log"
)
sce_citeseq = scater::runTSNE(
  sce_citeseq,altexp = "HTO", name = "TSNE_HTO", pca = TRUE
)
visualiseDim(sce_citeseq, dimNames = "TSNE_HTO")

sce_citeseq <- crossSampleDoublets(sce_citeseq)
visualiseDim(
  sce_citeseq, dimNames = "TSNE_HTO", colour_by = "doubletClassify_between_label"
)
sce_citeseq <- withinSampleDoublets(sce_citeseq, minPts = 10)

visualiseDim(
  sce_citeseq, dimNames = "TSNE_HTO", colour_by = "doubletClassify_within_label"
)

# Filter out doublets
sce_citeseq = sce_citeseq[
  , sce_citeseq$doubletClassify_between_class != "doublet/multiplet"
]
sce_citeseq = sce_citeseq[
  , sce_citeseq$doubletClassify_between_class != "negative"
]

####

sce_citeseq.pd = sce_citeseq@colData
se_hto@meta.data$SAMPLE_CLUSTER =  sce_citeseq.pd$doubletClassify_between_label[match(
 rownames( se_hto@meta.data), rownames(sce_citeseq.pd)
)]
se_hto = se_hto[, !is.na(se_hto$SAMPLE_CLUSTER)]

df = AverageExpression(
  se_hto, assays = "RNA", layer = "counts", group.by = "SAMPLE_CLUSTER"
)[[1]]
df = reshape2::melt(as.matrix(df))
df = df %>%
  dplyr::group_by(Var1) %>%
  slice_max(value, n = 1) %>% data.frame()
colnames(df) = c("HTO", "Cluster", "AVE")

df$Cluster = gsub("g", "", df$Cluster)
pd = read.csv(paste0(manifest$oekelen$work, "metadata/meta_sample.csv"))
pd$Hashing.BC = gsub("_.+", "", pd$Hashing.BC)
pd$GROUP = gsub(".+_", "", pd$X)
pd$CLUSTER = df$Cluster[match(pd$Hashing.BC, df$HTO)]

seurat@meta.data$SAMPLE_CLUSTER = se_hto@meta.data$SAMPLE_CLUSTER[match(
  seurat@meta.data$barcode, rownames(se_hto@meta.data)
)]
seurat@meta.data$orig.ident = pd$X[match(
  seurat@meta.data$SAMPLE_CLUSTER, pd$CLUSTER
)]
seurat@meta.data$SOURCE = pd$Tissue[match(
  seurat@meta.data$SAMPLE_CLUSTER, pd$CLUSTER
)]
seurat@meta.data$GROUP = pd$GROUP[match(
  seurat@meta.data$SAMPLE_CLUSTER, pd$CLUSTER
)]

table(is.na(seurat@meta.data$SAMPLE_CLUSTER))
seurat = seurat[, !is.na(seurat@meta.data$SAMPLE_CLUSTER)]

seurat@meta.data <- seurat@meta.data %>% mutate_if(is.character,as.factor)

# seurat = seurat %>%
#   NormalizeData() %>%
#   FindVariableFeatures() %>%
#   ScaleData() %>%
#   RunPCA() %>%
#   RunUMAP(dims = 1:20)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# scDblFinder
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
obj.l = Split_Object(seurat, split.by = "orig.ident", threads = 10)
scDb.res = parallel::mclapply(obj.l, function(se){
  set.seed(1234)
  sce = scDblFinder(GetAssayData(se, slot="counts"))
  df = data.frame(sce@colData) %>% dplyr::select(scDblFinder.score, scDblFinder.class)
  colnames(df) = c("scDblFinder_score", "scDblFinder_class")
  df$barcode = rownames(df)
  df
}, mc.cores = 1)
scDb.res = do.call("rbind", scDb.res)
rownames(scDb.res) = scDb.res$barcode
scDb.res$barcode = NULL
seurat = AddMetaData(seurat, scDb.res)
DefaultAssay(seurat) = "RNA"
Idents(seurat) = "orig.ident"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Save
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Idents(seurat) = "orig.ident"

saveRDS(seurat, file = obj.path)
