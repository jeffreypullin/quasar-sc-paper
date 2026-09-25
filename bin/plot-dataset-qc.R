#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(purrr)
  library(readr)
  library(stringr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

label_dataset <- function(x) {
  case_when(
    x == "onek1k" ~ "OneK1K",
    x == "tremor" ~ "Tremor",
    str_starts(x, "randolph") ~ "Randolph",
    TRUE ~ x
  )
}

count_files <- args[str_detect(basename(args), "cell-counts")]
cells_files <- args[str_detect(basename(args), "cells-per-indiv")]

cell_counts <- map_dfr(count_files, read_tsv, show_col_types = FALSE) |>
  mutate(dataset = label_dataset(dataset))

cells_per_indiv <- map_dfr(cells_files, read_tsv, show_col_types = FALSE) |>
  mutate(dataset = label_dataset(dataset)) |>
  summarise(n_cells = sum(n_cells), .by = c(dataset, individual))

dataset_levels <- c("OneK1K", "Tremor", "Randolph") |>
  intersect(unique(c(cell_counts$dataset, cells_per_indiv$dataset)))
dataset_levels <- c(
  dataset_levels,
  sort(setdiff(unique(c(cell_counts$dataset, cells_per_indiv$dataset)), dataset_levels))
)
dataset_cols <- setNames(
  rep(c("#009988", "#0077BB", "#EE7733", "#CC3311", "#33BBEE", "#AA4499"), length.out = length(dataset_levels)),
  dataset_levels
)

cell_counts <- cell_counts |>
  mutate(dataset = factor(dataset, levels = dataset_levels))
cells_per_indiv <- cells_per_indiv |>
  mutate(dataset = factor(dataset, levels = dataset_levels))

p_counts <- cell_counts |>
  ggplot(aes(log10(n_counts), colour = dataset, fill = dataset)) +
  geom_density(alpha = 0.2, linewidth = 0.8) +
  scale_colour_manual(values = dataset_cols) +
  scale_fill_manual(values = dataset_cols) +
  labs(
    x = expression(log[10]("UMI counts per cell")),
    y = "Density"
  ) +
  theme_jp()

p_cells <- cells_per_indiv |>
  ggplot(aes(n_cells, weight = n_cells, fill = dataset)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    bins = 30,
    colour = "white",
    linewidth = 0.2,
    position = "identity",
    alpha = 0.7
  ) +
  scale_fill_manual(values = dataset_cols) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(
    x = "Number of cells per individual",
    y = "Proportion of cells"
  ) +
  theme_jp()

ggsave(
  "dataset-qc.pdf",
  plot = p_counts + p_cells +
    plot_layout(guides = "collect") &
    theme(legend.position = "top"),
  width = 12,
  height = 5
)
