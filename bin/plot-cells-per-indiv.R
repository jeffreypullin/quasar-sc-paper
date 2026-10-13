#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

main_types <- unname(cell_type_lookup[c(
  "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_ET", "CD8_NC", "CD8_S100B",
  "NK", "NK_R", "Plasma", "B_Mem", "B_IN", "Mono_C", "Mono_NC", "DC"
)])

cells_per_indiv <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(cell_label %in% main_types)

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

# CV = sd / mean, so types with very different typical counts can be compared.
cv_data <- cells_per_indiv |>
  summarise(
    cv = sd(n_cells) / mean(n_cells),
    .by = cell_label
  ) |>
  transmute(cell_type = cell_label, cv) |>
  arrange(cell_type)

write_tsv(cv_data, "cells-per-indiv-cv.tsv")

cv_data <- cv_data |>
  mutate(cell_label = reorder(cell_type, cv))

p_cv <- cv_data |>
  ggplot(aes(cell_label, cv)) +
  geom_col(fill = "#44AA99") +
  coord_flip() +
  labs(
    x = "Cell type",
    y = "CV of cells per donor"
  ) +
  theme_jp_vgrid()

ggsave(
  "cells-per-indiv-cv-plot.pdf",
  plot = p_cv,
  width = 8,
  height = 8
)
