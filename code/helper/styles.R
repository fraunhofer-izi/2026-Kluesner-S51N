# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Colors
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.inst = c("ggthemes", "scales") %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}

colors.pal.10 = ggthemes::tableau_color_pal("Tableau 10")(10)
colors.pal.20 = ggthemes::tableau_color_pal("Tableau 20")(20)
colors_stata =  ggthemes::stata_pal("s2color")(15)

ct.col = c(
  "Plasma(blast)" = "#7C2529",
  "Plasmablast" = "#7C2529",
  "Plasma cell" = "#7C2529",
  "B-Cell" = "#E18A8D",
  "NK" = "#BC8400",
  "NK CD56bright" = "#e8d725",
  "CD56 bright NK" = "#BC8400",
  "CD4 T-Cell" = "#e8d725",
  "CD8 T-Cell" = "#AFA10D",
  "gd T-Cell" = "#20581C",
  "dp T-Cell" = "#CCDDAA",
  "T-Cell (cycling)" = "#CC3311",
  "Mono CD14" = "#93aeba",
  "CD14 Mono" = "#93aeba",
  "Mono CD16" = "#4B859F",
  "CD16 Mono" = "#4B859F",
  "cDC" = "#9C9BDB",
  "pDC" = "#194573",
  "other DC" = "#D1BBD7",
  "Macrophage" = "black",
  "Erythrocyte" = "#B281A6",
  "Platelet" = "#AA4499",
  "Progenitor" = "#555555",
  "HSPC" = "#555555",
  "Other" = "black",
  "Cycling" = "grey",
  "Not Estimable" = "black"
)


cd8.col = c(
  "CD8.NaiveLike" = "#0077BB",
  "CD8.CM" = "#33BBEE",
  "CD8.EM" = "#0c6e63",
  "CD8.TEMRA" = "#72b28a",
  "CD8.EMRA" = "#72b28a",
  "CD8.EMRA.1" = "#72b28a",
  "CD8.EMRA.2" = "#CCDDAA",
  "CD8 EMRA 1" = "#72b28a",
  "CD8 EMRA 2" = "#CCDDAA",
  "CD8 TEMRA" = "#CCDDAA",
  "CD8.MAIT" = "#997700",
  "CD8.TPEX" = "#EE3377",
  "CD8.TEX" = "#EE6677",
  "CD8.Cycling" = "#4f2535",
  "CD8 NaiveLike" = "#0077BB",
  "CD8 CM" = "#33BBEE",
  "CD8 EM" = "#0c6e63",
  "CD8 EMRA" = "#72b28a",
  "CD8 EMRA KLRC2+" = "#CCDDAA",
  "CD8 MAIT" = "#997700",
  "CD8 TPEX" = "#EE3377",
  "CD8 TEX" = "#EE6677",
  "CD8 Cycling" = "#4f2535"
)

cd4.col = c(
  "CD4.NaiveLike" = "#b7d2e0",
  "CD4 Memory" = "#216385",
  "CD4.CTL_EOMES" = "#da6f6f",
  "CD4.CTL_GNLY" = "#e5bfaf",
  "CD4 CTL_Exh" = "#EE7733",
  "CD4 CTL-Exh" = "#EE7733",
  "CD4.Cycling" = "#994455",
  "CD4.Tfh" = "#aca6e0",
  "CD4.Th17" = "#f5d39f",
  "CD4.Treg" = "#fdbfd4",
  "CD4 NaiveLike" = "#b7d2e0",
  "CD4 CTL EOMES+" = "#da6f6f",
  "CD4 CTL GNLY+" = "#e5bfaf",
  "CD4 CTL Exh" = "#EE7733",
  "CD4 Cycling" = "#994455",
  "CD4 Tfh" = "#aca6e0",
  "CD4 Th17" = "#f5d39f",
  "CD4 Treg" = "#fdbfd4"
)

t.coarse.col = c(
  "CD4 NaiveLike" = "#b7d2e0",
  "CD4 Exhausted" = "#da6f6f",
  "CD4 CTL" = "#e5bfaf",
  "CD4 Cycling" = "#994455",
  "CD4 Helper" = "#aca6e0",
  "CD4 Treg" = "#fdbfd4",
  "CD4 NaiveLike" = "#b7d2e0",
  "CD8 Effector" = "#72b28a",
  "CD8 EM" = "#72b28a",
  "CD8 TEMRA" = "#CCDDAA",
  "CD8 NaiveLike" = "#216385",
  "CD8 NaiveLike/CM" = "#216385",
  "CD8 MAIT" = "#997700",
  "CD8 Exhausted" = "#EE6677",
  "CD8 Cycling" = "#4f2535",
  "dnT" = "black"
)

til.col = c(cd8.col, cd4.col, "dnT" = "#666633")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ggplot theme
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mytheme = function(base_size = 8, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      axis.ticks = element_line(colour = "black"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(size = rel(1), colour = "black"),
      strip.text.y = element_text(size = rel(1), colour = "black"),
      strip.text = element_text(size = rel(1), colour = "black"),
      axis.text = element_text(size = rel(1), colour = "black"),
      axis.title = element_text(size = rel(1), colour = "black"),
      legend.title = element_text(colour = "black", size = rel(1)),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = .3),
      legend.key.size = unit(1, "lines"),
      legend.text = element_text(size = rel(1), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(hjust = 0, face = "plain", colour = "black", size = rel(1)),
      plot.subtitle = element_text(colour = "black", size = rel(.85))
    )
}

mytheme_grid = function(base_size = 8, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      axis.ticks = element_line(colour = "black"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(size = rel(1), colour = "black"),
      strip.text.y = element_text(size = rel(1), colour = "black"),
      strip.text = element_text(size = rel(1), colour = "black"),
      axis.text = element_text(size = rel(1), colour = "black"),
      axis.title = element_text(size = rel(1), colour = "black"),
      legend.title = element_text(colour = "black", size = rel(1)),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = .3),
      legend.key.size = unit(1, "lines"),
      legend.text = element_text(size = rel(1), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(hjust = 0, face = "plain", colour = "black", size = rel(1)),
      plot.subtitle = element_text(colour = "black", size = rel(.85))
    )
}

rm.axis = theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  plot.title = element_blank(),
  axis.line = element_blank()
)

pt_to_mm <- function(x) x / ggplot2::.pt
theme_custom = function(base.size = 8){theme_bw(base_size = base.size) +
    theme(
      axis.text.x = element_text(size = rel(1), color = "black"),
      axis.text.y = element_text(size = rel(1), color = "black"),
      axis.title.x = element_text(size = rel(1), color = "black"),
      axis.title.y = element_text(size = rel(1), color = "black"),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = pt_to_mm(0.75)
      ),
      strip.background = element_rect(
        color = "black",
        fill = "grey90",
        linewidth = pt_to_mm(0.75)
      ),
      legend.text  = element_text(size = rel(6.5/8)),
      legend.title = element_text(size = rel(6.5/8)),
      legend.key.size = unit(base_size * 0.7, "pt"),
      legend.spacing.y = unit(base_size * 0.3, "pt"),
      legend.spacing.x = unit(base_size * 0.5, "pt"),
      strip.text = element_text(
        size = rel(0.875),
        margin = margin(t = 1, r = 2, b = 1, l = 2)
      ),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = rel(1)),
      plot.subtitle = element_text(size = rel(1)),
      plot.caption = element_text(size = rel(1)),
      text = element_text(color = "black"),
      axis.line = element_line(linewidth = pt_to_mm(0.5)),
      axis.ticks = element_line(linewidth = pt_to_mm(0.5)) #,
      # panel.grid.major = element_line(linewidth = pt_to_mm(0.25)),
      # plot.margin = margin(
      #   t = 15,
      #   r = 15,
      #   b = 15,
      #   l = 15,
      #   unit = "mm"
      # )
    )}
