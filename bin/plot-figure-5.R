#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggupset)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

interaction <- readRDS(args[1])
time <- readRDS(args[2])
upset <- readRDS(args[3])

figure <- wrap_plots(
  interaction,
  time,
  wrap_elements(full = upset),
  design = "AA\nBC",
  heights = c(1.4, 1)
) +
  plot_annotation(tag_levels = "a")

ggsave(
  "figure-5.pdf",
  figure,
  width = 16,
  height = 14,
  limitsize = FALSE
)
