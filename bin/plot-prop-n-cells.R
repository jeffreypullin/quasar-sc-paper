#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(tidyr)
  library(data.table)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)
perm_sc_data_files <- read_tsv(args[1], show_col_types = FALSE)
cells_per_indiv <- read_tsv(args[2], show_col_types = FALSE)

main_types <- c("CD4_NC", "CD4_ET", "CD4_SOX4", "B_IN", "Mono_NC", "Plasma")

p_threshold <- 1e-4

tail_summary <- function(variant_file) {
  pvalue <- fread(variant_file, select = "pvalue")$pvalue
  pvalue <- pvalue[!is.na(pvalue)]
  tibble(n_below = sum(pvalue < p_threshold), n = length(pvalue))
}

tail_data <- perm_sc_data_files |>
  filter(
    model == "p_glmm_sc",
    int_cov == "none",
    as.character(k) == "none",
    indiv_frac == 1,
    count_frac == 1,
    as.numeric(min_n_cells) < 0,
    cell_type %in% main_types
  ) |>
  mutate(summary = map(variant_file, tail_summary)) |>
  unnest(summary) |>
  summarise(
    n_below = sum(n_below),
    n_tests = sum(n),
    obs_exp = n_below / (n_tests * p_threshold),
    .by = c(cell_type, cell_frac)
  )

# Mirrors subsample_cells_within_individuals: floor(frac * n), at least 1 cell.
cells_data <- cells_per_indiv |>
  filter(n_cells > 0) |>
  inner_join(
    tibble(cell_type = main_types, cell_label = unname(cell_type_lookup[main_types])),
    by = "cell_label"
  ) |>
  inner_join(distinct(tail_data, cell_type, cell_frac), by = "cell_type", relationship = "many-to-many") |>
  mutate(n_cells_ds = pmax(1, floor(cell_frac * n_cells))) |>
  summarise(
    mean_cells_per_donor = mean(n_cells_ds),
    median_cells_per_donor = median(n_cells_ds),
    .by = c(cell_type, cell_frac)
  )

plot_data <- inner_join(tail_data, cells_data, by = c("cell_type", "cell_frac")) |>
  arrange(cell_type, cell_frac)

write_tsv(plot_data, "prop-n-cells.tsv")

ct_levels <- unname(cell_type_lookup[main_types[main_types %in% plot_data$cell_type]])

p <- plot_data |>
  mutate(cell_type = factor(unname(cell_type_lookup[cell_type]), levels = ct_levels)) |>
  ggplot(aes(mean_cells_per_donor, obs_exp, colour = cell_type)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_path() +
  geom_point(size = 2) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 100, 1000)) +
  scale_colour_manual(
    values = setNames(cell_type_cols, cell_type_lookup[names(cell_type_cols)]),
    drop = TRUE
  ) +
  labs(
    x = "Mean cells per donor",
    y = "Observed / expected permutation p < 1e-4",
    colour = "Cell type"
  ) +
  theme_jp_vgrid()

ggsave("prop-n-cells-plot.pdf", p, width = 10, height = 7)
