#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

perm <- readRDS(args[1])
concordance <- readRDS(args[2])
power <- readRDS(args[3])
time <- readRDS(args[4])

figure <- wrap_plots(
  perm,
  concordance,
  wrap_elements(full = power),
  time,
  ncol = 2
) +
  plot_annotation(tag_levels = "a")

ggsave("figure-2.pdf", figure, width = 18, height = 11)
