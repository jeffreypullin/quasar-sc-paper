#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

cells_per_indiv <- read_tsv(args[1], show_col_types = FALSE)

p <- cells_per_indiv |>
  ggplot(aes(n_cells)) +
  geom_histogram(bins = 30, fill = "#44AA99", color = "white", linewidth = 0.2) +
  facet_wrap(~ cell_label, scales = "free") +
  labs(
    x = "Number of cells per individual",
    y = "Number of individuals"
  ) +
  theme_jp()

ggsave(
  "cells-per-indiv-plot.pdf",
  plot = p,
  width = 16,
  height = 12
)
