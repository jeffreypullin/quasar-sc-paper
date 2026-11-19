#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

scatter <- readRDS(args[1])
power <- readRDS(args[2])
time <- readRDS(args[3])

figure <- wrap_plots(
  scatter,
  wrap_elements(full = power),
  time,
  ncol = 1
) +
  plot_annotation(tag_levels = "a")

ggsave(
  "sc-pglmm-vs-lmm-figure.pdf",
  figure,
  width = 14,
  height = 22,
  limitsize = FALSE
)
