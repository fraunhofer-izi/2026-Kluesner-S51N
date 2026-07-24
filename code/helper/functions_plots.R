# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Make grandient for ggplot (background) | for DGEA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
make_gradient <- function(deg = 45, n = 100, cols = blues9) {

  .cran_packages = c("grid", "ggplot2","RColorBrewer")
  .inst = .cran_packages %in% installed.packages()
  if (any(!.inst)) {
    install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
  }
  for (pack in .cran_packages) {
    suppressMessages(library(
      pack,
      quietly = TRUE,
      verbose = FALSE,
      character.only = TRUE
    ))
  }

  cols <- colorRampPalette(cols)(n + 1)
  rad <- deg / (180 / pi)
  mat <- matrix(
    data = rep(seq(0, 1, length.out = n) * cos(rad), n),
    byrow = TRUE,
    ncol = n
  ) +
    matrix(
      data = rep(seq(0, 1, length.out = n) * sin(rad), n),
      byrow = FALSE,
      ncol = n
    )
  mat <- mat - min(mat)
  mat <- mat / max(mat)
  mat <- 1 + mat * n
  mat <- matrix(data = cols[round(mat)], ncol = n)
  grid::rasterGrob(
    image = mat,
    width = unit(1, "npc"),
    height = unit(1, "npc"),
    interpolate = TRUE
  )
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA Volcano
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dgea_volcano = function(
    dgea.res = NULL,
    dgea.res.sign = NULL,
    nbr.tops = 7,
    cl.label = NULL, # e.g.  cl.label = setNames(c("CD4", "CD8"), c("CD4 T-Cell", "CD8 T-Cell"))
    sort.by.p = T,
    pval.column = "p_val",
    logFC.column = "avg_log2FC",
    facet.scales = "free_y",
    facet.rows = 1,
    nudge_x = 2,
    x.axis.sym = F,
    x.axis.ext = 0,
    geom.hline = log2(1.5),
    box.padding = 0.3,
    label.padding = .12,
    label.size = 2,
    leg.title = "FDR <0.05",
    cut.y.thresh = NULL,
    cut.x.thresh = NULL,
    pt.size = .1
){

  dgea.res$ID = paste0(dgea.res$cluster, "_", dgea.res$feature)
  dgea.res.sign$ID = paste0(dgea.res.sign$cluster, "_", dgea.res.sign$feature)
  dgea.res = dgea.res[dgea.res$cluster %in% unique(dgea.res.sign$cluster), ]
  dgea.res = dgea.res[!is.na(dgea.res$p_val), ]

  if(sort.by.p == T) {
    tops_up = dgea.res.sign %>%
      dplyr::filter(.data[[logFC.column]] > 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::arrange(.data[[pval.column]], -abs(.data[[logFC.column]])) %>%
      dplyr::slice_head(n=nbr.tops)
    tops_up$GENE_CL = paste0(tops_up$feature, "_", tops_up$cluster)

    tops_down = dgea.res.sign %>%
      dplyr::filter(.data[[logFC.column]] < 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::arrange(.data[[pval.column]], -abs(.data[[logFC.column]])) %>%
      dplyr::slice_head(n=nbr.tops)
    tops_down$GENE_CL = paste0(tops_down$feature, "_", tops_down$cluster)
  } else {
    tops_up = dgea.res.sign %>%
      dplyr::filter(.data[[logFC.column]] > 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(.data[[logFC.column]], n = nbr.tops)
    tops_up$GENE_CL = paste0(tops_up$feature, "_", tops_up$cluster)

    tops_down = dgea.res.sign %>%
      dplyr::filter(.data[[logFC.column]] < 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::top_n(n = nbr.tops, wt = -.data[[logFC.column]])
    tops_down$GENE_CL = paste0(tops_down$feature, "_", tops_down$cluster)
  }

  dgea.res$GENE_CL = paste0(dgea.res$feature, "_", dgea.res$cluster)
  dgea.res = dgea.res %>% dplyr::mutate(
    label_up = ifelse(GENE_CL %in% tops_up$GENE_CL, feature, "")
  )
  dgea.res = dgea.res %>% dplyr::mutate(
    label_down = ifelse(GENE_CL %in% tops_down$GENE_CL, feature, "")
  )
  dgea.res$SIGNIFICANT = dgea.res$ID %in% dgea.res.sign$ID
  dgea.res$SIGNIFICANT = factor(dgea.res$SIGNIFICANT, levels = c(TRUE, FALSE))

  dgea.res = dgea.res %>%
    mutate(p_val = ifelse(p_val == 0, 1e-300, dgea.res[[pval.column]]))
  if(!is.null(cut.y.thresh)){
    dgea.res$p_val[dgea.res$p_val < cut.y.thresh] = cut.y.thresh
  }
  if(!is.null(cut.x.thresh)){
    dgea.res$avg_log2FC[dgea.res$avg_log2FC < -(cut.x.thresh)] = -(cut.x.thresh)
    dgea.res$avg_log2FC[dgea.res$avg_log2FC > cut.x.thresh] = cut.x.thresh
  }

  axis.max = max(abs(dgea.res[[logFC.column]]), na.rm = T)

  if(is.null(cl.label)) {
    tbl = table(
      dgea.res.sign$cluster,
      (dgea.res.sign$significant & dgea.res.sign$avg_log2FC > 0)
    )
    if (length(colnames(tbl)) == 1) {
      if (colnames(tbl) == TRUE) {
        cl.label = setNames(
          paste0(rownames(tbl), " | ", "up: ", tbl[, 1], " | down: 0"),
          rownames(tbl)
        )
      } else {
        cl.label = setNames(
          paste0(rownames(tbl), " | ", "up: 0", " | down: ", tbl[, 1]),
          rownames(tbl)
        )
      }
    } else {
      cl.label = setNames(
        paste0(rownames(tbl), " | ", "up: ", tbl[, 2], " | down: ", tbl[, 1]),
        rownames(tbl)
      )
    }
  }

  g <- make_gradient(
    deg = 180, n = 500,
    cols = scico::scico(
      9, palette = 'vik', begin = .3, end = .7, direction = -1,
    )
  )
  set.seed(42)

  pl =
    ggplot(dgea.res, aes(x = .data[[logFC.column]], y = -log10(p_val))) +
    # annotation_custom(
    #   grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    # ) +
    geom_point(
      data = subset(dgea.res, SIGNIFICANT == F), aes(color = SIGNIFICANT),
      size = pt.size
    ) +
    geom_point(
      data = subset(dgea.res, SIGNIFICANT == T), aes(color = SIGNIFICANT),
      size = pt.size
    ) +
    theme(
      panel.spacing = unit(1.5, "lines"),
      panel.grid.minor = element_blank(),
      legend.position="bottom",
      strip.text = element_text(size = rel(1), face = "plain")
    ) +
    facet_wrap(
      ~ cluster, scales = facet.scales, labeller = labeller(cluster = cl.label),
      nrow = facet.rows
    ) +
    geom_label_repel(
      data = subset(dgea.res, SIGNIFICANT == T & label_up != ""),
      label = subset(dgea.res, SIGNIFICANT == T & label_up != "")$label_up,
      segment.colour = "black",
      size = label.size,
      direction = "y",
      hjust = .5,
      nudge_x = nudge_x,
      # nudge_y = -2,
      segment.size = .2,
      box.padding = box.padding,
      label.padding = label.padding,
      min.segment.length = 0,
      max.overlaps = 50
    ) +
    geom_label_repel(
      data = subset(dgea.res, SIGNIFICANT == T & label_down != ""),
      label = subset(dgea.res, SIGNIFICANT == T & label_down != "")$label_down,
      segment.colour = "black",
      size = label.size,
      direction = "y",
      hjust = .5,
      # xlim = c(NA, 0),
      nudge_x = -nudge_x,
      segment.size = .2,
      box.padding = box.padding,
      label.padding = label.padding,
      min.segment.length = 0,
      max.overlaps = 50
    ) +
    geom_vline(xintercept = -geom.hline, linetype = "dashed", linewidth = .2) +
    geom_vline(xintercept = geom.hline, linetype = "dashed", linewidth = .2) +
    scale_color_manual(values = c("TRUE" = "#555555", "FALSE" = "#BBBBBB")) +
    labs(y = "-Log10(p-value)", x = "Log2 fold change", colour = leg.title)  +
    guides(colour = guide_legend(override.aes = list(size=3)))
  if(x.axis.sym == T) {
    pl = pl + xlim(-axis.max - x.axis.ext, axis.max + x.axis.ext)
  }
  pl
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA | Paired Boxplots
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bxpl_dgea_paired = function(
    obj = NULL,
    dgea.res.df = NULL,
    ftr = c("NFKB1", "TBX21", "EOMES"),
    facet.scale = "free_y",
    facet.ncol = 4,
    vst.norm = F,
    add.facet.groups = NULL,
    panel.space = .5,
    .base.size = 8
){

  df = FetchData(obj, vars = c("MUTATION", "SAMPLE", "celltype", ftr))
  df = reshape2::melt(df)

  if(!is.null(add.facet.groups)){
    df$FACET_GROUP = add.facet.groups$Group[match(
      df$variable, add.facet.groups$Gene_symbol
    )]
  }

  if(vst.norm == T){
    mat = GetAssayData(obj, layer = "counts") %>% as.matrix()
    dds = DESeq2::DESeqDataSetFromMatrix(
      countData = mat,
      colData = FetchData(obj, vars = c("MUTATION", "SAMPLE")),
      design = ~ 0
    )
    vsd <-  DESeq2::varianceStabilizingTransformation(dds, blind=T)
    mat = assay(vsd)[ftr, ]
    mat.melt = reshape2::melt(mat)
    mat.melt$SAMPLE = colData(dds)$SAMPLE[
      match(mat.melt$Var2, rownames(colData(dds)))
    ]
    mat.melt$MUTATION = colData(dds)$MUTATION[
      match(mat.melt$Var2, rownames(colData(dds)))
    ]
    df$value = mat.melt$value[match(
      paste0(df$variable, df$MUTATION, df$SAMPLE),
      paste0(mat.melt$Var1, mat.melt$MUTATION, mat.melt$SAMPLE)
    )]
  }

  if(!is.null(dgea.res.df)){
    df$FDR = dgea.res.df$p_val_adj[match(df$variable, dgea.res.df$feature)]
    df = rstatix::add_significance(df, p.col = "FDR")

    lbls = df[!duplicated(df$variable), ]
    lbls = setNames(
      paste0(lbls$variable, " | ", lbls$FDR.signif),
      as.character(lbls$variable)
    )
  } else {
    lbls = setNames(
      as.character(unique(df$variable)),
      as.character(unique(df$variable))
    )
  }

  p = ggplot(df, aes(MUTATION, value, fill = MUTATION)) +
    geom_boxplot(outliers = F) +
    geom_point(size = .1) +
    geom_line(aes(group = SAMPLE), linewidth = .3) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(panel.space, "lines"),
      legend.position = "right",
      ggh4x.facet.nestline = element_line(colour = "black", linewidth = .2)
    ) +
    ylab("Norm. Expression") +
    labs(fill = NULL)

  if(is.null(add.facet.groups)){
    p = p + facet_wrap(
      ~ variable, scale = facet.scale, ncol = facet.ncol,
      labeller = labeller(variable = lbls)
    )
  } else {
    p = p + facet_nested_wrap(
      ~ FACET_GROUP + variable, scale = facet.scale,
      labeller = labeller(variable = lbls), ncol = facet.ncol
    )
  }
  p

}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DecoupleR results | Heamtmap
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
decoupler_tileplot = function(
    obj = NULL,
    alpha = 0.05,
    title = NULL,
    max.cutoff = .99,
    shape.de = T,
    dot.range = c(3, 1),
    ftr.fltr = 0,
    split.by.group = F,
    nbr.tops = NULL,
    barwidth = 7,
    barheight = .35,
    ftrs.clustering = F
){

  ftrs = obj %>%
    dplyr::filter(FDR < alpha) %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::filter(n > ftr.fltr) %>%
    dplyr::pull(var = source)

  obj = obj[obj$source %in% ftrs, ]

  if(!is.null(nbr.tops)){
    ftrs.up = obj %>%
      dplyr::filter(FDR < alpha) %>%
      dplyr::filter(score > 0) %>%
      dplyr::group_by(celltype) %>%
      dplyr::slice_min(FDR , n = nbr.tops) %>%
      dplyr::pull(source)
    ftrs.down = obj %>%
      dplyr::filter(FDR < alpha) %>%
      dplyr::filter(score < 0) %>%
      dplyr::group_by(celltype) %>%
      dplyr::slice_min(FDR , n = nbr.tops) %>%
      dplyr::pull(source)
    obj = obj[obj$source %in% c(ftrs.up, ftrs.down), ]
    # obj$source = factor(obj$source, levels = c(ftrs.down, rev(ftrs.up)))
  }


  if(ftrs.clustering == T){
    if(split.by.group == T){

      df = obj %>%
        dplyr::mutate(celltype_grp = paste0(celltype, GRP)) %>%
        tidyr::pivot_wider(
          id_cols = 'celltype_grp', names_from = 'source', values_from = 'score'
        ) %>%
        tibble::column_to_rownames('celltype_grp')

      df[is.na(df)] = Inf
      hclust_ftr <- hclust(dist(t(df)), method = "complete")

    } else {
      df = obj %>%
        tidyr::pivot_wider(
          id_cols = 'celltype', names_from = 'source', values_from = 'score'
        ) %>%
        tibble::column_to_rownames('celltype')

      df[is.na(df)] = 0
      hclust_ftr <- hclust(dist(t(df)), method = "complete")
    }
    obj$source = factor(
      obj$source, levels = hclust_ftr$labels[hclust_ftr$order]
    )
  }

  p = quantile(obj$score[obj$score > 0], max.cutoff)
  n = -quantile(abs(obj$score[obj$score < 0]), max.cutoff)
  if(p > abs(n)){
    obj$score[obj$score > p] = p
    obj$score[obj$score < -p] = -p
  } else if (abs(n) > p) {
    obj$score[obj$score < n] = n
    obj$score[obj$score > abs(n)] = abs(n)
  }
  max.v = max(abs(obj$score))

  pl =
    ggplot(obj, aes(source,  y = celltype)) +
    geom_tile(
      aes(fill = score), color = "white",  linewidth = .3
    )
  if(shape.de == T){
    pl = pl + geom_point(
      data = obj[obj$FDR < alpha, ],
      aes(size = FDR),
      fill = "white", pch=21, show.legend = T
    ) +
      scale_size(range = dot.range) +
      scale_shape_manual(
        values = c(8), name = "test", labels = c("FDR < 0.05", "")
      )
  }
  pl = pl + scico::scale_fill_scico(
    palette = "vik", midpoint = 0, begin = .05, end = .95,
    limits = c(-max.v, max.v)
  ) +
    xlab(NULL) + ylab(NULL) +
    theme(
      panel.spacing = unit(.75, "lines"),
      axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = rel(1)),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom",
      legend.ticks.length = unit(0.075, 'cm'),
      # legend.title = element_text(margin = margin(r = 3, l = 5, unit = "pt")),
      # legend.margin = margin(t = 0),
      # legend.text = element_text(margin = margin(l = -1, t = 2, b = 2)),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    guides(
      fill = guide_colorbar(
        title = "Activity score", title.vjust = 1.1,
        barwidth = unit(barwidth, 'lines'), barheight = unit(barheight, 'lines'),
        ticks.linewidth = 1/.pt, ticks = T, frame.colour="black",
        frame.linewidth = 0.5/.pt, order = 1
      )
    ) +
    ggtitle(title)
  pl
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# QC violing plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
qc_vln_plot_cell = function(
    obj = se.meta,
    .features = "nFeature_RNA",
    .group.by = "orig.ident",
    plot_title = "Genes Per Cell",
    x_axis_label = NULL,
    y_axis_label = "Features",
    low_cutoff = NULL,
    high_cutoff = NULL
){
  library(ggplot2)
  library(ggthemes)

  if(!is.list(obj)) {
    df = obj@meta.data %>% dplyr::select(.data[[.group.by]], .data[[.features]])
    df
  }
  if(class(obj) == "list") {
    l = lapply(obj, function(x){
      df = x@meta.data %>% dplyr::select(.data[[.group.by]], .data[[.features]])
    })
    df = do.call("rbind", l)
    df
  }
  if(class(obj) == "data.frame") {
    df = obj %>% dplyr::select(.data[[.group.by]], .data[[.features]])
    df
  }

  ggplot(data = df, mapping = aes(x = .data[[.group.by]], y = .data[[.features]])) +
    geom_violin(
      size = .1,
      width = 1,
      scale = "area",
      na.rm = TRUE
    ) +
    stat_summary(
      fun.min = function(z) { quantile(z,0.25) },
      fun.max = function(z) { quantile(z,0.75) },
      fun = median, colour = "#0077BB", size = .2) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle=45, vjust=1, hjust=1),
      axis.ticks.x = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = rel(1.2))
    ) +
    geom_hline(
      yintercept = c(low_cutoff, high_cutoff), linetype = "dashed",
      color = "#BB5566", size = .3
    ) +
    xlab(x_axis_label) +
    ylab(y_axis_label) +
    ggtitle(plot_title)
}

