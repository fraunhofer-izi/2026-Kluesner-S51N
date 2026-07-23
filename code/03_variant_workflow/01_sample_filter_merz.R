# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries and some Functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "Seurat", "yaml", "dplyr", "stringr", "naturalsort", "cowplot", "data.table",
  "ggplot2", "ggthemes", "patchwork", "devtools", "scCustomize"
)
.bioc_packages = c()

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

source("code/helper/styles.R")
theme_set(mytheme(base_size = 12))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load objects and phenodata
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")
# se.t = readRDS(paste0(
#   manifest$merz$objects, "/meta/integration/06_seurat_harmony_t_all.Rds"
# ))
se.t = readRDS(
  paste0(manifest$merz$work, "seurat/03_seurat_anno_t.Rds")
)
se.t@meta.data = droplevels(se.t@meta.data)

clin.pd = readRDS("/homes/olymp/michael.rade/projects/2025-LT-Persister/data/metadata_samples/MasterTabelleCART_04112025.Rds")

pdata = readRDS("data/metadata_samples_merz.Rds")

# se.t@meta.data = se.t@meta.data  %>% mutate(across(where(is.factor), as.character))

unique(se.t$orig.ident) %in% pdata$orig.ident
se.t$PRODUCT = as.character(
  pdata$PRODUCT[match(se.t@meta.data$orig.ident, pdata$orig.ident)]
)
se.t$TIMEPOINT = as.character(
  pdata$TIMEPOINT[match(se.t@meta.data$orig.ident, pdata$orig.ident)]
)
se.t$SAMPLE_ID = as.character(
  pdata$SAMPLE_ID[match(se.t@meta.data$orig.ident, pdata$orig.ident)]
)
se.t$PATIENT_ID = as.character(
  pdata$PATIENT_ID[match(se.t@meta.data$orig.ident, pdata$orig.ident)]
)
se.t@meta.data$PRODUCT[se.t@meta.data$orig.ident == "P065-PB"] = "cilta"
se.t@meta.data$TIMEPOINT[se.t@meta.data$orig.ident == "P065-PB"] = "Very Very Late"
se.t@meta.data$SAMPLE_ID[se.t@meta.data$orig.ident == "P065-PB"] = "Patient065_1"
se.t@meta.data$PATIENT_ID[se.t@meta.data$orig.ident == "P065-PB"] = "Patient065"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FC Data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pdata.facs = clin.pd$pdata.immune.subsets
pdata.facs = pdata.facs[pdata.facs$DAY == "Day 30" | pdata.facs$DAY == "Day 100", ]
pdata.facs$DAY = ifelse(pdata.facs$DAY == "Day 30", "Late", "Very Late")
pdata.facs$ID = paste0(pdata.facs$SAMPLE_ID, "_", pdata.facs$DAY)
pdata.facs = pdata.facs %>%
  dplyr::select(SAMPLE_ID, PATIENT_ID, DAY, ID, CD3_CAR_PERC)
pdata.facs = pdata.facs[!is.na(pdata.facs$PATIENT_ID), ]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Merge FACS and RNA Data | Overview plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd = se.t@meta.data
df = data.frame(
  prop.table(table(pd$orig.ident, pd$CAR_BY_EXPRS), margin = 1) * 100
)
df = df[df$Var2 == "TRUE", ][, c(1,3)]
colnames(df) = c("orig.ident", "CD3_CAR_PERC_SC")

df = df %>% dplyr::left_join(
  pd %>%
    dplyr::group_by(orig.ident) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(orig.ident, TIMEPOINT, PRODUCT, STUDY, SAMPLE_ID, PATIENT_ID),
  by = "orig.ident"
) %>%
  # dplyr::filter(TIMEPOINT != "LP") %>%
  dplyr::arrange(desc(CD3_CAR_PERC_SC))

df$CD3_CAR_PERC_FC = pdata.facs$CD3_CAR_PERC[
  match(paste0(df$SAMPLE_ID, "_", df$TIMEPOINT), pdata.facs$ID)
]

df.pl = rbind(
  df %>% dplyr::select(
    PERC = CD3_CAR_PERC_SC, TIMEPOINT, PRODUCT, orig.ident
  ) %>%
    dplyr::mutate(GROUP = "RNA"),
  df %>% dplyr::select(
    PERC = CD3_CAR_PERC_FC, TIMEPOINT, PRODUCT, orig.ident
  ) %>%
    dplyr::mutate(GROUP = "FC")
)
df.pl = df.pl[df.pl$TIMEPOINT != "LP", ]
df.pl$orig.ident = factor(df.pl$orig.ident, levels = unique(df.pl$orig.ident))

ggplot(df.pl, aes(x=orig.ident, y=PERC, fill = GROUP)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(
    axis.text.x = element_text(angle=45, vjust=1, hjust=1)
  ) +
  geom_hline(yintercept = 1, linewidth = .3, linetype = "dashed") +
  # facet_wrap(~ STUDY) +
  ylab("Proportion") +
  xlab("Sample") +
  ggtitle("Estimated percentages of CAR+ cells for cilta-cel treated patients") +
  facet_wrap(~ PRODUCT, nrow = 2, scales = "free")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Filter for samples with CAR+ cells estimated by FACS or RNA-Seq
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
car.pos.thres = 0

tops = subset(df, CD3_CAR_PERC_SC > car.pos.thres | CD3_CAR_PERC_FC > car.pos.thres) %>%
  droplevels()
tops.lp = droplevels(df[df$TIMEPOINT == "LP",])
tops.lp = tops.lp[tops.lp$SAMPLE_ID %in% names(table(tops$SAMPLE_ID)), ]
tops = rbind(tops.lp, tops)
tops = tops[order(tops$TIMEPOINT, tops$CD3_CAR_PERC_FC), ]
tops$orig.ident.old = pd$orig.ident.old[match(
  tops$orig.ident, pd$orig.ident
)]
tops = dplyr::rename(tops, orig.ident.new = orig.ident)
tops = tops %>% dplyr::relocate(CD3_CAR_PERC_SC , .after = CD3_CAR_PERC_FC) %>%
  droplevels()

tops.cilta = droplevels(subset(tops, PRODUCT == "cilta"))
tops.ide = droplevels(subset(tops, PRODUCT == "ide"))

tops.cilta = tops.cilta[tops.cilta$TIMEPOINT != "LP", ]
tops.cilta = droplevels(tops.cilta)
table(tops.cilta$PATIENT_ID, tops.cilta$TIMEPOINT)
table(tops.cilta$STUDY)
nrow(tops.cilta)

saveRDS(tops.cilta, "data/pub_S51A/cellranger_samples_cilta_merz.Rds")
# tops.cilta = readRDS("data/pub_S51A/cellranger_samples_cilta.Rds")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Make Seurat object for demux | demultiplexOnSeurat.py
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta.t.sub = se.t[c("CAR-BCMA", "CD4", "CD8A", "CD8B"), ]
se.meta.t.sub = se.meta.t.sub[, grepl("CD4|CD8", se.meta.t.sub$celltype_short_3)]

se.meta.t.sub = se.meta.t.sub[
  , se.meta.t.sub$orig.ident %in% tops.cilta$orig.ident.new
]

smpls.run.1to4 = openxlsx::read.xlsx(
  "/homes/olymp/michael.rade/projects/2023-BCMA-CAR/publication/supplementary_info/table_qc_demux.xlsx",
  sheet = 1, rowNames = F, startRow = 1, detectDates = T
)

se.meta.t.sub$Multiplexed = smpls.run.1to4$Multiplexed[match(
  se.meta.t.sub$orig.ident, smpls.run.1to4$Sample.ID
)]
se.meta.t.sub$FASTQ_FILE_NAME = smpls.run.1to4$FASTQ.file.name[match(
  se.meta.t.sub$orig.ident, smpls.run.1to4$Sample.ID
)] %>% as.character()
se.meta.t.sub$Multiplexed = ifelse(
  is.na(se.meta.t.sub$Multiplexed),
  "no", se.meta.t.sub$Multiplexed
)
se.meta.t.sub$FASTQ_FILE_NAME = ifelse(
  is.na(as.character(se.meta.t.sub$FASTQ_FILE_NAME)),
  as.character(se.meta.t.sub$orig.ident.old),
  as.character(se.meta.t.sub$FASTQ_FILE_NAME)
)

se.meta.t.sub@meta.data = droplevels(se.meta.t.sub@meta.data)
se.meta.t.sub$barcodes = gsub(".*_", "", colnames(se.meta.t.sub))
se.meta.t.sub$GROUP = "mined"
se.meta.t.sub$SAMPLE_NEW = se.meta.t.sub$orig.ident
se.meta.t.sub$ID = paste0(
  se.meta.t.sub$PATIENT_ID, "_", se.meta.t.sub$TIMEPOINT
)
se.meta.t.sub@meta.data = se.meta.t.sub@meta.data %>%
  dplyr::select(
    orig.ident, barcodes, GROUP, SAMPLE_NEW, PATIENT_ID, TIMEPOINT,
    celltype_short_3, CAR_BY_EXPRS, ID, Multiplexed, FASTQ_FILE_NAME
  )

se.meta.t.sub = RenameCells(
  se.meta.t.sub, new.names = paste0("mined_", colnames(se.meta.t.sub))
)
se.meta.t.sub@meta.data$FASTQ_FILE_NAME[
  se.meta.t.sub@meta.data$orig.ident == "P065-PB"
] = "P065-PB"

saveRDS(se.meta.t.sub, "data/pub_S51A/seurat_t_cilta_merz.Rds")

