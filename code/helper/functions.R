# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Normalize, Harmony, Clustering
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
integration = function(
    obj,
    obj.l = NULL,
    no.ftrs = 500,
    .assay = "RNA",
    max.cl = 1,
    min.cells.per.sample = 25,
    threads = 5,
    .nbr.dims = 15,
    do.cluster = T,
    do.dimreduc = T,
    hvg.union = T,
    custom.features = NULL,
    run.integration = T,
    harmony.group.vars = NULL,
    obj.split.by = "orig.ident",
    .perp = 50,
    dmreduc.dims = NULL,
    n.neighbors = 30,
    min.dist = 0.3) {

  library(SignatuR)
  data(SignatuR)
  library(parallel)
  library(BiocParallel)
  library(harmony)

  if (any(!"readgmt" %in% installed.packages())) {
    Sys.unsetenv("GITHUB_PAT")
    devtools::install_github("jhrcook/readgmt")
  }
  library(readgmt)

  set.seed(1234)

  start.time <- Sys.time()

  if(is.null(dmreduc.dims)) {
    dmreduc.dims = .nbr.dims
  }

  DefaultAssay(obj) = "RNA"
  obj@meta.data = droplevels(obj@meta.data)
  obj = DietSeurat(obj, counts = TRUE, data = TRUE)
  obj = NormalizeData(obj)


  # Gene categories to exclude from variable genes
  bl <- c(
    SignatuR::GetSignature(SignatuR$Hs$Compartments$Mito)[[1]],
    SignatuR::GetSignature(SignatuR$Hs$Compartments$TCR)[[1]]
  )
  bl = c(bl, c("RPS4Y1", "EIF1AY", "DDX3Y", "KDM5D", "XIST"))
  bl <- unique(bl)

  if (hvg.union == T) {

    if(is.null(obj.l)) {
      print("Split object")
      obj.l = Split_Object(obj, split.by = obj.split.by, threads = threads)
    }

    select.bool = unlist(
      lapply(obj.l, function(x){ncol(x) >= min.cells.per.sample})
    )
    print(table(select.bool))
    obj.l = obj.l[select.bool]
    length(obj.l)

    print("HVG")
    obj.l = parallel::mclapply(obj.l, function(x) {
      x = x[!rownames(x) %in% bl, ]
      x = FindVariableFeatures(
        x, selection.method = "vst",  assay = .assay, verbose = FALSE
      )
      x
    }, mc.cores = threads)

    features = SelectIntegrationFeatures(
      object.list = obj.l, nfeatures = no.ftrs
    )

    VariableFeatures(obj) = features

    rm(obj.l); gc()

  } else if (!is.null(custom.features)) {
    VariableFeatures(obj) = custom.features
  } else {
    obj = FindVariableFeatures(
      obj, selection.method = "vst", nfeatures = no.ftrs, assay = .assay
    )
  }

  obj = ScaleData(obj, assay = .assay)
  obj = RunPCA(obj, assay = .assay)
  # plot(ElbowPlot(obj, ndims = 50))

  if(run.integration == T){
    obj = RunHarmony(
      obj, group.by.vars = harmony.group.vars, # theta = c(2,3),
      reduction.use ='pca', dims.use = 1:.nbr.dims,
      max_iter = 15, ncores = threads, verbose = T
    )
    comp.wrk = 'harmony'
  } else {
    comp.wrk = 'pca'
  }

  if (do.cluster == T) {
    print("Find clusters")
    obj = FindNeighbors(
      obj, reduction = comp.wrk, dims = 1:.nbr.dims, verbose = F
    )
    reso = seq(0,max.cl,.1)
    names(reso) = reso
    suppressWarnings({
      suppressMessages({
        findclusters.res = parallel::mclapply(reso, function(x) {
          FindClusters(obj, resolution = x, verbose = F)@meta.data[, "seurat_clusters", drop = F]
        }, mc.cores = length(reso))
      })
    })
    res.names = names(findclusters.res)
    findclusters.res = do.call("cbind", findclusters.res)
    colnames(findclusters.res) = paste0("RNA_snn_res.", res.names)
    stopifnot(identical(rownames(obj@meta.data), rownames(findclusters.res)))
    obj = AddMetaData(obj, findclusters.res)
  }
  if (do.dimreduc == T) {
    set.seed(1234)
    print("UMAP")
    obj = RunUMAP(
      obj, reduction = comp.wrk, dims = 1:dmreduc.dims, seed.use = 1234,
      min.dist = min.dist, n.neighbors = n.neighbors, verbose = F
    )
  }

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)

  return(obj)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# % of present cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
present_cells_per_ftr = function(
    obj,
    features,
    group1,
    group2,
    .assay = "RNA",
    perc_expr_thres = 0
){

  count_mat = GetAssayData(obj, assay = .assay, layer = "counts")
  meta<- obj@meta.data %>%
    tibble::rownames_to_column(var = "cell")

  # get percentage of positive cells matrix
  count_df<- as.matrix(count_mat) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="gene") %>%
    tidyr::pivot_longer(!gene, names_to = "cell", values_to = "count") %>%
    left_join(meta) %>%
    group_by(gene, .data[[group1]], .data[[group2]]) %>%
    summarise(percentage = mean(count > perc_expr_thres)) %>%
    tidyr::pivot_wider(
      names_from = c(.data[[group1]], .data[[group2]]),
      values_from= percentage,
      names_sep="|"
    )

  percent_mat<- count_df[, -1] %>% as.matrix()
  rownames(percent_mat)<- count_df$gene

  percent_mat
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA: GEX
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dgea_gex = function(
    obj,
    .target = "celltype_short_3",
    group = "GROUP",
    ctrs.grp1 = NULL,
    ctrs.grp2 = NULL,
    split.by.tp = F,
    min.pct = .25,
    logfc.threshold = log2(1.25),
    .min.cells = 10,
    subsample = F,
    subsample.n = 250,
    test.method = "MAST",
    latent.vars = NULL,
    threads = 20,
    min.de.genes = NULL,
    deseq_paired = F,
    dsgn = NULL,
    shrink_and_rank = F
){

  se.w = obj
  se.w@meta.data = droplevels(se.w@meta.data)

  se.w@meta.data[[.target]] = gsub("_", "\\.", se.w@meta.data[[.target]])
  if(split.by.tp == T) {
    se.w$TMP  = paste0(se.w@meta.data[[.target]], "_", se.w@meta.data$TIMEPOINT)
    tbl = as.data.frame.matrix(table(se.w$TMP, se.w@meta.data[[group]]))
    tbl = tbl[, colnames(tbl) %in% c(ctrs.grp1, ctrs.grp2)]
    tbl = tbl[matrixStats::rowSums2(tbl > .min.cells) == 2, ]
    print(tbl)

    # pd = se.w@meta.data
    # pd = pd[pd$TMP %in% rownames(tbl), ]
    # pd = pd[pd[[group]] %in% c(ctrs.grp1, ctrs.grp2), ]
    # for (i in naturalsort(unique(pd$TMP))) {
    #   d = droplevels(pd[pd$TMP %in% i, ])
    #   p = lapply(split(d, d[[group]]), function(x){
    #     x = droplevels(x)
    #     unname(table(x$orig.ident))
    #   })
    #   print("--------------------")
    #   print(i)
    #   print(p)
    # }

  } else {
    se.w$TMP = se.w@meta.data[[.target]]
    print(table(se.w@meta.data[[.target]], se.w@meta.data[[group]]))
  }

  res.dgea = run_de(
    obj = se.w, target = "TMP", min.cells = .min.cells,
    lfc.thresh = logfc.threshold, min.pct.thres = min.pct,
    contrast.group = group, contrast = c(ctrs.grp1, ctrs.grp2),
    subsample = subsample, subsample.n = subsample.n,
    test.method = test.method, latent.vars = latent.vars,
    threads = threads, deseq_paired = deseq_paired, dsgn = dsgn,
    shrink_and_rank = shrink_and_rank
  )

  res.dgea$celltype = gsub("_.+", "", res.dgea$cluster)
  res.dgea$timepoint = gsub(".+_", "", res.dgea$cluster)
  res.dgea.sign = subset(res.dgea, significant == T)
  res.dgea.sign = res.dgea.sign[
    abs(res.dgea.sign$avg_log2FC) > logfc.threshold,
  ]
  res.dgea.sign = res.dgea.sign[
    res.dgea.sign$pct.1 > min.pct | res.dgea.sign$pct.2 > min.pct,
  ]

  if(!is.null(min.de.genes)){
    ct.keep = names(table(res.dgea.sign$cluster)[table(res.dgea.sign$cluster) > min.de.genes])
    res.dgea.sign = res.dgea.sign[res.dgea.sign$cluster %in% ct.keep, ]
  }

  if(nrow(res.dgea.sign) == 0){
    return(print("No DE genes found"))
  }

  res.dgea$ID = paste0(res.dgea$cluster, "_", res.dgea$feature)
  res.dgea.sign$ID = paste0(res.dgea.sign$cluster, "_", res.dgea.sign$feature)
  res.dgea$significant = ifelse(
    res.dgea$ID %in% res.dgea.sign$ID, TRUE, FALSE
  )
  gc()
  print(table(res.dgea.sign$cluster, res.dgea.sign$group))
  list(res.dgea = res.dgea, res.dgea.sign = res.dgea.sign)
}

run_de = function(
    obj = NULL,
    target = NULL,
    contrast.group = NULL,
    contrast = NULL,
    min.cells = 10,
    slot = "data",
    assay = "RNA",
    padj.thresh = 0.05,
    lfc.thresh = .1,
    min.pct.thres = 0,
    test.method = "wilcox",
    latent.vars = NULL,
    subsample = F,
    subsample.n = 200,
    threads = 20,
    deseq_paired = F,
    dsgn = NULL,
    shrink_and_rank = F
){

  DefaultAssay(obj) = assay

  obj = DietSeurat(
    obj,
    counts = TRUE,
    data = T,
    scale.data = FALSE,
    features = NULL,
    assays = assay,
    dimreducs = F,
    graphs = F
  )

  obj@meta.data = droplevels(obj@meta.data)
  obj@meta.data$RM = is.na(obj@meta.data[[target]])
  obj = subset(obj, subset = RM == F)

  obj = obj[, obj@meta.data[[contrast.group]] %in% contrast]
  obj@meta.data = droplevels(obj@meta.data)

  print("Split object")
  obj.l = Split_Object(obj, split.by = target, threads = threads)

  tbl = table(obj@meta.data[[contrast.group]], obj@meta.data[[target]])

  # Filter out celltypes with less than x cells in one group
  obj.l = obj.l[colnames(tbl)[colSums(tbl >= min.cells) == 2]]

  options(warn = 1)
  bpparam = BiocParallel::MulticoreParam(workers = threads)

  if(deseq_paired == F){print(paste0("Run DGEA: ", test.method))}

  dgea.res = BiocParallel::bplapply(names(obj.l), function(x) {

    print(x)

    o = obj.l[[x]]
    if(subsample == T){
      Idents(o) = "orig.ident"
      o = subset(o, downsample = subsample.n)
    }

    if(deseq_paired) {

      library(DESeq2)

      print("DESeq2 | Paired mode")

      p = o@meta.data %>% mutate_if(is.character,as.factor)
      dds = DESeqDataSetFromMatrix(
        countData = GetAssayData(o, layer = "counts"),
        colData = p,
        design = ~ 0
      )
      design(dds) = formula(dsgn)
      dds[[contrast.group]] <- relevel(dds[[contrast.group]], ref = contrast[2])

      # keep.exprs = edgeR::filterByExpr(y = counts(dds), group = dds$GROUP)
      # dds = DESeq(dds[keep.exprs,])
      dds = DESeq(dds)

      # resultsNames(dds)
      ctrst = paste0(contrast.group, "_", contrast[1], "_vs_", contrast[2])
      print(paste0("contrast: ", ctrst))
      markers = results(dds, name=ctrst) %>% data.frame()

      if(shrink_and_rank == T) {
        library(topconfects)
        dconfects <- topconfects::deseq2_confects(dds, name=ctrst)
        markers$confect = dconfects$table$confect[match(
          rownames(markers), dconfects$table$name
        )]
        markers$rank = dconfects$table$rank[match(
          rownames(markers), dconfects$table$name
        )]

        lfc.shrink = lfcShrink(dds, coef=ctrst, type="apeglm")
        markers$avg_log2FC_Shrink = lfc.shrink$log2FoldChange[match(
          rownames(markers), rownames(lfc.shrink)
        )]
      }

      frctn = present_cells_per_ftr(
        obj = o, features = markers$feature,
        group1 = contrast.group, group2 = "TMP"
      )
      colnames(frctn) = gsub("\\|.+", "", colnames(frctn))
      markers$pct.1 = frctn[, contrast[1]][
        match(rownames(markers), rownames(frctn))
      ]
      markers$pct.2 = frctn[, contrast[2]][
        match(rownames(markers), rownames(frctn))
      ]

      # print(markers %>% head())

      if(shrink_and_rank == T) {
        markers = markers %>%
          dplyr::select(
            p_val = pvalue,
            avg_log2FC = log2FoldChange,
            avg_log2FC_Shrink = avg_log2FC_Shrink,
            pct.1, pct.2,
            p_val_adj = padj,
            stat = stat,
            baseMean,
            confect,
            rank
          )
      } else {
        markers = markers %>%
          dplyr::select(
            p_val = pvalue,
            avg_log2FC = log2FoldChange,
            pct.1, pct.2,
            p_val_adj = padj,
            stat = stat,
            baseMean
          )
      }

    } else {
      markers = Seurat::FindMarkers(
        object = o,
        ident.1 = contrast[1], ident.2 = contrast[2],
        group.by = contrast.group, assay = assay,
        logfc.threshold = lfc.thresh, min.pct = min.pct.thres,
        test.use = test.method, latent.vars = latent.vars,
        min.cells.group = 2
      )
    }

    markers$feature = rownames(markers)
    markers$cluster = x
    markers$group.1 = contrast[1]
    rownames(markers) = paste0(markers$feature, "_", markers$cluster)
    markers

  }, BPPARAM = bpparam)

  dgea.res = do.call("rbind", dgea.res)

  dgea.res = dgea.res[order(dgea.res$p_val_adj, decreasing = F), ]
  dgea.res$significant = (dgea.res$p_val_adj < padj.thresh) &
    (abs(dgea.res$avg_log2FC) > lfc.thresh) &
    (dgea.res$pct.1 > min.pct.thres | dgea.res$pct.2 > min.pct.thres)

  dgea.res
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# get metadata from Seurat
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
get_metadata <- function(obj, ..., embedding = names(obj@reductions), nbr.dim = 2) {

  res = as_tibble(obj@meta.data, rownames = "cell")

  if (!is.null(embedding)) {
    if (any(!embedding %in% names(obj@reductions))) {
      stop(paste0(embedding, " not found in seurat object\n"), call. = FALSE)
    }
    embed_dat = purrr::map(names(obj@reductions), ~obj@reductions[[.x]]@cell.embeddings[, 1:nbr.dim]) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell")

    res = dplyr::left_join(res, embed_dat, by = "cell")
  }

  if (length(list(...)) > 0) {
    cols_to_get <- setdiff(..., colnames(obj@meta.data))
    if (length(cols_to_get) > 0) {
      res = Seurat::FetchData(obj, vars = cols_to_get) %>%
        tibble::rownames_to_column("cell") %>%
        dplyr::left_join(res, ., by = "cell")
    }
  }
  res
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Split Seurat object (BiocParallel)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Split_Object = function(object, split.by = "orig.ident", threads = 5) {

  library(parallel)
  library(BiocParallel)

  groupings <- FetchData(object = object, vars = split.by)[, 1]
  groupings <- unique(x = as.character(x = groupings))
  names(groupings) = groupings

  if (is.null(threads)) {
    bpparam = BiocParallel::MulticoreParam(workers = length(groupings))
  } else {
    bpparam = BiocParallel::MulticoreParam(workers = threads)
  }

  obj.list = BiocParallel::bplapply(groupings, function(grp) {
    cells <- which(x = object[[split.by, drop = TRUE]] == grp)
    cells <- colnames(x = object)[cells]
    se = subset(x = object, cells = cells)
    se@meta.data = droplevels(se@meta.data)
    se
  }, BPPARAM = bpparam)

  return(obj.list)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Assign VDJ
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
assign_vdj = function(
    obj,
    vdj = "vdj_t",
    batch = NULL,
    present.bool = TRUE,
    export.table = FALSE,
    .filterMulti = NULL,
    .remove.na = NULL
){

  stopifnot(vdj %in% c("vdj_t", "vdj_b"))

  library(scRepertoire)

  l = list()
  for (i in batch) {
    print(i)
    cellranger.dirs = list.dirs(
      path = i, full.names = T, recursive = F
    )
    cellranger.samples = basename(cellranger.dirs)
    fltrd.vdj = paste0(
      cellranger.dirs, "/outs/per_sample_outs/", cellranger.samples, "/", vdj,
      "/filtered_contig_annotations.csv"
    )
    names(fltrd.vdj) = gsub("multi_", "", cellranger.samples)
    print(length(fltrd.vdj))
    l[[i]] = fltrd.vdj
  }
  names(l) = NULL
  fltrd.vdj = unlist(l, use.names = T)

  ###

  contig_list <- lapply(fltrd.vdj, function(x) {
    tryCatch(read.csv(x), error=function(e) NULL)
  })
  length(contig_list)
  print(paste0(
    "Empty file: ", names(lengths(contig_list)[lengths(contig_list) == 0])
  ))
  contig_list = contig_list[lengths(contig_list) != 0]

  contig_list = lapply(contig_list, function(x){
    x$sample = NULL
    x
  })

  if(vdj == "vdj_t"){
    combined <- scRepertoire::combineTCR(
      contig_list,
      filterMulti = .filterMulti,
      remove.na = .remove.na,
      filter.nonproductive = T,
      samples = paste0(names(contig_list))
    )
  } else {
    combined <- scRepertoire::combineBCR(
      contig_list,
      filterMulti = .filterMulti,
      remove.na = .remove.na,
      filter.nonproductive = T,
      samples = paste0(names(contig_list))
    )
  }

  combined = data.table::rbindlist(combined) %>% data.frame()

  if(present.bool == TRUE){
    if(vdj == "vdj_t"){
      obj$VDJ_T_AVAIL = combined$CTnt[
        match(rownames(obj@meta.data), combined$barcode)
      ]
      obj$VDJ_T_AVAIL = ifelse(is.na(obj$VDJ_T_AVAIL), FALSE, TRUE)
      print(table(obj$VDJ_T_AVAIL))
    } else {
      obj$VDJ_B_AVAIL = combined$CTnt[
        match(rownames(obj@meta.data), combined$barcode)
      ]
      obj$VDJ_B_AVAIL = ifelse(is.na(obj$VDJ_B_AVAIL), FALSE, TRUE)
      print(table(obj$VDJ_B_AVAIL))
    }
    obj
  } else {

    combined$sample = obj$orig.ident[
      match(combined$barcode, rownames(obj@meta.data))
    ]
    combined = combined[!is.na(combined$sample), ]
    combined = split(combined, combined$sample)

    combined.keep = lapply(combined, function(x){nrow(x)}) > 0
    combined = combined[combined.keep]

    max.clonotypes = max(
      unlist(lapply(combined, function(x){max(unname(table(x$CTstrict)))}))
    )
    obj <- combineExpression(
      combined, obj,
      cloneCall = "strict",
      group.by = "sample",
      proportion = F,
      cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=max.clonotypes)
    )
    obj$cloneSize = droplevels(obj$cloneSize)
    # obj
    if(export.table == T){
      obj@meta.data[, (length(colnames(obj@meta.data))-6):length(colnames(obj@meta.data))]
    } else{
      obj
    }
  }
}

combineVDJ = function(
    vdj = "vdj_t",
    batch = NULL,
    .filterMulti = NULL,
    .remove.na = NULL
){

  stopifnot(vdj %in% c("vdj_t", "vdj_b"))

  library(scRepertoire)

  l = list()
  for (i in batch) {
    print(i)
    cellranger.dirs = list.dirs(
      path = i, full.names = T, recursive = F
    )
    cellranger.samples = basename(cellranger.dirs)
    fltrd.vdj = paste0(
      cellranger.dirs, "/outs/per_sample_outs/", cellranger.samples, "/", vdj,
      "/filtered_contig_annotations.csv"
    )
    names(fltrd.vdj) = gsub("multi_", "", cellranger.samples)
    print(length(fltrd.vdj))
    l[[i]] = fltrd.vdj
  }
  names(l) = NULL
  fltrd.vdj = unlist(l, use.names = T)

  ###

  contig_list <- lapply(fltrd.vdj, function(x) {
    tryCatch(read.csv(x), error=function(e) NULL)
  })
  length(contig_list)
  print(paste0(
    "Empty file: ", names(lengths(contig_list)[lengths(contig_list) == 0])
  ))
  contig_list = contig_list[lengths(contig_list) != 0]

  contig_list = lapply(contig_list, function(x){
    x$sample = NULL
    x
  })

  if(vdj == "vdj_t"){
    combined <- scRepertoire::combineTCR(
      contig_list,
      filterMulti = .filterMulti,
      remove.na = .remove.na,
      samples = paste0(names(contig_list))
    )
  } else {
    combined <- scRepertoire::combineBCR(
      contig_list,
      filterMulti = .filterMulti,
      remove.na = .remove.na,
      samples = paste0(names(contig_list))
    )
  }

  combined = data.table::rbindlist(combined) %>% data.frame()
  combined

}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# clean Clonotype: Clonotypes must not have CD4 and CD8 cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
clean_clonotypes = function(
    obj = NULL,
    celltype = "celltype_short_3",
    clone_id = "CTstrict"
){

  obj$celltype = obj[[celltype]]
  obj$clone_id = obj[[clone_id]]

  pd = obj@meta.data %>%
    dplyr::select(celltype, clone_id)
  pd$barcode = rownames(pd)
  cl.lin.diff =
    pd %>%
    dplyr::filter(!is.na(clone_id)) %>%
    dplyr::filter(grepl("CD4|CD8", celltype)) %>%
    dplyr::group_by(clone_id) %>%
    dplyr::mutate(clone_size = n()) %>%
    dplyr::mutate(cd4 = grepl("CD4", celltype)) %>%
    dplyr::mutate(cd8 = grepl("CD8", celltype)) %>%
    dplyr::mutate(nbr_cd4 = sum(cd4)) %>%
    dplyr::mutate(nbr_cd8 = sum(cd8)) %>%
    dplyr::filter(nbr_cd4 > 0 & nbr_cd8 > 0)

  rm.cl.1 = cl.lin.diff %>%
    dplyr::filter(clone_size <= 5) %>%
    dplyr::pull(barcode) %>%
    unique()

  cl.lin.diff = cl.lin.diff[!cl.lin.diff$barcode %in% rm.cl.1, ]

  if(nrow(cl.lin.diff) == 0){
    return(c(rm.cl.1))
  }

  cl.lin.diff$ratio = NA
  for (i in unique(cl.lin.diff$clone_id)) {
    cl.sub = cl.lin.diff[cl.lin.diff$clone_id == i, ]
    cl.sub = cl.sub[!duplicated(cl.sub$clone_id), ]
    .max = names(which.max(cl.sub[, c("nbr_cd4", "nbr_cd8")]))
    .min = names(which.min(cl.sub[, c("nbr_cd4", "nbr_cd8")]))
    cl.lin.diff[cl.lin.diff$clone_id == i, ]$ratio = cl.sub[[.max]] / cl.sub[[.min]]
  }

  rm.cl.2 = cl.lin.diff %>%
    dplyr::filter(ratio <= 3) %>%
    data.frame() %>%
    dplyr::pull(barcode) %>%
    unique()

  cl.lin.diff = cl.lin.diff[!cl.lin.diff$barcode %in% rm.cl.2, ]

  if(nrow(cl.lin.diff) == 0){
    return(c(rm.cl.1, rm.cl.2))
  }

  rm.cl.3 = cl.lin.diff %>%
    dplyr::group_by(clone_id, celltype) %>%
    dplyr::mutate(lin = n()) %>%
    dplyr::group_by(clone_id) %>%
    dplyr::filter(lin == min(lin)) %>%
    dplyr::pull(barcode) %>%
    unique()

  return(c(rm.cl.1, rm.cl.2, rm.cl.3))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
count_cells_per_sample = function(
    obj = NULL,
    count.base = NULL,
    col.name = NULL
){

  l = lapply(obj, function(x){
    x@meta.data %>% dplyr::select(orig.ident)
  })
  df = do.call("rbind", l)
  df = df %>%  dplyr::count(orig.ident)
  if(is.null(count.base)) {
    return(df)
  } else {
    count.base[[col.name]] = df$n[match(count.base$orig.ident, df$orig.ident)]
    return(count.base)
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Doublet detection | scds
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
scds_doublets = function(se){
  suppressWarnings({
    suppressMessages({
      sce = as.SingleCellExperiment(se)
      set.seed(1234)
      sce = scds::cxds(sce)
      sce = scds::bcds(sce)
      sce = scds::cxds_bcds_hybrid(sce, estNdbl = T)
      CD  = data.frame(sce@colData)
      CD = CD %>% dplyr::select(
        cxds_score, bcds_score, hybrid_score, cxds_call, bcds_call, hybrid_call
      )
      CD$barcode = rownames(CD)
      CD
    })
  })
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cd4cd8_car_present = function(
    obj
){
  cd8cd4 = FetchData(obj, c("CD8A", "CD8B", "CD4"), layer = "counts")
  ct.cd8cd4 = cd8cd4 %>% mutate(
    CD4CD8_BY_EXPRS = case_when(
      CD4 > 0 & CD8A == 0 & CD8B == 0 ~ "CD4+CD8-",
      CD4 == 0 & (CD8A > 0 | CD8B > 0) ~ "CD4-CD8+",
      CD4 == 0 & CD8A == 0 & CD8B == 0 ~ "CD4-CD8-",
      CD4 > 0 & (CD8A > 0 | CD8B > 0) ~ "CD4+CD8+",
      TRUE ~ "unresolved"
    )) %>%
    dplyr::select(CD4CD8_BY_EXPRS)
  rownames(ct.cd8cd4) = rownames(cd8cd4)
  obj = AddMetaData(obj, ct.cd8cd4)

  cd3 = FetchData(obj, c("CD3D", "CD3E", "CD3G"), slot = "counts")
  cd3 = cd3 %>% mutate(
    CD3_BY_EXPRS = case_when(
      CD3D > 0 | CD3E > 0 | CD3G > 0 ~ "CD3",
      TRUE ~ "unresolved"
    )) %>%
    dplyr::select(CD3_BY_EXPRS)
  obj = AddMetaData(obj, cd3)

  if("CAR-BCMA" %in% rownames(obj)){
    car.ftr = FetchData(obj, c("CAR-BCMA"), layer = "counts")
    obj$CAR_BY_EXPRS = as.factor(car.ftr[[1]] > 0)
  } else {
    obj$CAR_BY_EXPRS = FALSE
    obj$CAR_BY_EXPRS = as.factor(obj$CAR_BY_EXPRS)
  }

  obj
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cell Cycle | Gene set enrichment
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
estimate_cc = function(se){

  library(ProjecTILs)

  suppressWarnings({
    suppressMessages({

      data(cell.cycle.obj) # ProjecTILs package
      DefaultAssay(se) = "RNA"

      se = se %>%
        FindVariableFeatures(verbose = F) %>%
        ScaleData(verbose = F) %>%
        RunPCA(npcs = 20, verbose = F) %>%
        FindNeighbors(reduction = "pca", dims = 1:20, verbose = F)

      tmp =  tryCatch(
        FindClusters(se, resolution = 2, verbose = F),
        error=function(e) "error"
      )
      if(class(tmp) != "Seurat") {
        se = FindClusters(se, resolution = 1, verbose = F)
      } else {
        se = tmp
      }

      quiet <- function(x) {
        sink(tempfile())
        on.exit(sink())
        invisible(force(x))
      }

      cc.phase = quiet(
        clustifyr::run_gsea(
          GetAssayData(object = se, assay = "RNA", slot = "data"),
          query_genes = cell.cycle.obj$human$cycling,
          cluster_ids =  se@meta.data[["seurat_clusters"]], n_perm = 1000
        )
      )

      cc.phase$pval_adj = p.adjust(cc.phase$pval, method = "BH")
      cc.cl = rownames(cc.phase[cc.phase$pval < 0.05, ])

      cc.cells = (se@meta.data[["seurat_clusters"]] %in% cc.cl)
      se@meta.data$CellCycle = factor(cc.cells)

      pd.cc = se@meta.data[cc.cells, ]
      pd.cc = pd.cc %>% dplyr::mutate(
        CellCycle_Phase = dplyr::case_when(
          G2M.Score > S.Score ~ "G2M",
          TRUE ~ "S"
        )
      )
      se$CellCycle_Phase = pd.cc$CellCycle_Phase[
        match(rownames(se@meta.data), rownames(pd.cc))
      ]
      se$CellCycle_Phase[is.na(se$CellCycle_Phase)] = "G1M"
      se$CellCycle_Phase = factor(
        se$CellCycle_Phase, levels = c("G1M", "S", "G2M")
      )

      res = se@meta.data %>% dplyr::select(CellCycle, CellCycle_Phase)
      res$barcode = rownames(res)
      res
    })
  })
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DecoupleR
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
run_decoupler = function(obj = NULL, cores = 10){

  net.cyto = read.csv2(
    "data/signatures/cytosig_signature.tsv", sep = "\t",
    check.names = F
  )
  net.cyto = reshape2::melt(as.matrix(net.cyto)) %>% as_tibble()
  colnames(net.cyto) = c("target", "source", "weight")
  net.cyto = net.cyto %>% dplyr::select(source, target, weight)
  net.cyto$weight = as.numeric(net.cyto$weight)
  net.cyto$source = as.character(net.cyto$source)
  net.cyto$target = as.character(net.cyto$target)

  if(!file.exists("data/signatures/tf_db.Rds")){
    net.tf <- decoupleR::get_collectri(organism = 'human')
    net.tf = dplyr::rename(net.tf, weight = mor)
    saveRDS(net.tf, file = "data/signatures/tf_db.Rds")
  } else {
    net.tf = readRDS("data/signatures/tf_db.Rds")
  }

  if(!file.exists("data/signatures/progeny_db.Rds")){
    net.progeny = decoupleR::get_progeny(organism = 'human', top = 500)
    saveRDS(net.progeny, file = "data/signatures/progeny_db.Rds")
  } else {
    net.progeny = readRDS("data/signatures/progeny_db.Rds")
  }

  run = function(df, net.obj = NULL) {
    res.l = parallel::mclapply(unique(df$celltype), function(x) {

      deg.ct = df[df$celltype == x, ]
      rownames(deg.ct) = deg.ct$feature

      res <- decoupleR::run_ulm(
        mat = deg.ct[, 't', drop = FALSE],
        net = net.obj,
        .source = 'source',
        .target = 'target',
        .mor = "weight",
        minsize = 5
      )

      res$FDR = p.adjust(res$p_value, method = "BH")
      res$celltype = x
      res

    }, mc.cores = 11)

    res = do.call("rbind", res.l) %>% data.frame()
    res[order(res$FDR, decreasing = F), ]
  }


  deg = obj$res.dgea
  deg$t = deg$stat
  deg = deg[!is.na(deg$t), ]

  dcplr.cyto = run(df = deg, net.obj = net.cyto)
  dcplr.tf = run(df = deg, net.obj = net.tf)
  dcplr.progeny = run(df = deg, net.obj = net.progeny)

  deg.sign = obj$res.dgea.sign

  dcplr.cyto$DE = deg.sign$feature[match(
    paste0(dcplr.cyto$source, "_", dcplr.cyto$celltype), rownames(deg.sign)
  )]
  dcplr.cyto$DE = ifelse(!is.na(dcplr.cyto$DE), "Y", NA)

  dcplr.tf$DE = deg.sign$feature[match(
    paste0(dcplr.tf$source, "_", dcplr.tf$celltype), rownames(deg.sign)
  )]
  dcplr.tf$DE = ifelse(!is.na(dcplr.tf$DE), "Y", NA)

  l = list(
    dcplr.cyto = dcplr.cyto,
    dcplr.tf = dcplr.tf,
    dcplr.progeny = dcplr.progeny
  )
  l
}

