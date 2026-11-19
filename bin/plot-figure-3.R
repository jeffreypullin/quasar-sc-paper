#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

power <- readRDS(args[1])
concordance <- readRDS(args[2])

figure <- wrap_plots(
  concordance,
  wrap_elements(full = power),
  ncol = 1,
  heights = c(4.5, 8)
) +
  plot_annotation(tag_levels = "a")

ggsave(
  "figure-3.pdf",
  figure,
  width = 14,
  height = 13.5,
  limitsize = FALSE
)
