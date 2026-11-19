#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

qq_data <- map(commandArgs(trailingOnly = TRUE), \(f) {
  x <- readRDS(f)
  if (!"analysis" %in% names(x$int)) {
    x$int$analysis <- "GxC"
  }
  if (!"analysis" %in% names(x$main)) {
    x$main$analysis <- "GxC"
  }
  x
})

plot_qq <- function(data) {
  data |>
    mutate(
      perm = recode(perm, Global = "Genotype"),
      perm = factor(perm, levels = c("Genotype", "Interaction", "Within-sample")),
      analysis = factor(analysis, levels = c("GxC", "Grouped"))
    ) |>
    ggplot(aes(log_x_bin_mid, log_y_pvalue,
               ymin = log_lower_ci, ymax = log_upper_ci,
               colour = perm, fill = perm)) +
    geom_point(alpha = 0.8) +
    geom_abline(linetype = "dashed") +
    geom_ribbon(linetype = 2, alpha = 0.1, colour = NA) +
    facet_grid(
      rows = vars(analysis),
      cols = vars(perm),
      labeller = labeller(perm = \(x) paste(x, "permutation"))
    ) +
    labs(
      x = "Expected -log10(p-value)",
      y = "Observed -log10(p-value)"
    ) +
    theme_jp() +
    theme(legend.position = "none")
}

int_p <- plot_qq(map_dfr(qq_data, "int"))

ggsave(
  "perm-sc-all-int-plot.pdf",
  int_p,
  width = 12,
  height = 8
)
saveRDS(int_p, "perm-sc-all-int-plot.rds")

ggsave(
  "perm-sc-all-main-plot.pdf",
  plot_qq(map_dfr(qq_data, "main")),
  width = 12,
  height = 8
)
