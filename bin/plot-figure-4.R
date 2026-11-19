#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

cell_frac <- readRDS(args[1])
prop_n_cells <- readRDS(args[2])

figure <- wrap_plots(
  cell_frac,
  prop_n_cells,
  ncol = 1,
  heights = c(10, 7)
) +
  plot_annotation(tag_levels = "a")

ggsave(
  "figure-4.pdf",
  figure,
  width = 14,
  height = 17,
  limitsize = FALSE
)
