# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310992

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("Seurat", "yaml", "dplyr", "doParallel", "parallel", "data.table")
.bioc_packages = c("SingleCellExperiment", "scDblFinder")

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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load Rawcounts and create a merged Seurat object
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")
work.path = paste0(manifest$hosoya$work, "cellranger/")
seurat.path = manifest$hosoya$seurat
dir.create(seurat.path, recursive = T)

obj.path = paste0(paste0(seurat.path, "01_seurat_ori.Rds"))

cellranger.dirs = list.dirs(path = work.path, full.names = T, recursive = T)

fltrd.dirs = cellranger.dirs[grepl("sample_filtered_feature_bc_matrix", cellranger.dirs)]
tmp = gsub(work.path, "", fltrd.dirs)
tmp = gsub("/", "", gsub("out.*", "", tmp))
tmp = gsub("multi_", "", tmp)
names(fltrd.dirs) = tmp
print(length(names(fltrd.dirs)))

# i = names(fltrd.dirs)[1]
bpparam = BiocParallel::MulticoreParam(workers = 1)
seurat.l = BiocParallel::bplapply(names(fltrd.dirs), function(i) {

  id = i
  fltrd.counts = Read10X(data.dir = fltrd.dirs[names(fltrd.dirs) == id], gene.column = 2)
  seu.obj = CreateSeuratObject(counts = fltrd.counts, project = id)
  seu.obj@meta.data$orig.ident = id

  set.seed(1234)
  sce = scDblFinder(GetAssayData(seu.obj, slot="counts"))
  df = data.frame(sce@colData) %>%
    dplyr::select(scDblFinder.score, scDblFinder.class)
  colnames(df) = c("scDblFinder_score", "scDblFinder_class")
  seu.obj = AddMetaData(seu.obj, df)

  seu.obj

}, BPPARAM = bpparam)

seurat = merge(
  seurat.l[[1]], y = seurat.l[2:length(seurat.l)],
  add.cell.ids = names(seurat.l), project = "hosoya"
)
seurat[["RNA"]] <- JoinLayers(seurat[["RNA"]])
DefaultAssay(seurat) = "RNA"
Idents(seurat) = "orig.ident"

seurat@meta.data$GSM = seurat@meta.data$orig.ident
seurat@meta.data = seurat@meta.data %>% dplyr::mutate(
  orig.ident = dplyr::case_when(
    orig.ident == "GSM9314767" ~ "Duodenum",
    orig.ident == "GSM9314768" ~ "Stomach",
    orig.ident == "GSM9314769" ~ "Ileum"
  )
)
seurat@meta.data$orig.ident = factor(seurat@meta.data$orig.ident)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Save
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
saveRDS(seurat, file = obj.path)
