#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)

cell_types <- c("Plasma", "B_IN", "CD4_NC")

pb_manifest <- pb_data_files |>
  filter(
    model == "lm",
    int_cov == "none",
    indiv_frac == 1,
    count_frac == 1,
    cell_type %in% cell_types
  )

sc_models <- c("p_glmm_sc", "lmm_sc")
method_levels <- unname(method_lookup[c("lm", sc_models)])

sc_manifest <- sc_data_files |>
  filter(
    model %in% sc_models,
    cov_spec == "bulk_pca",
    k == "none",
    indiv_frac == 1,
    count_frac == 1,
    cell_type %in% cell_types
  )

group_cols <- c("cell_type", "method", "min_n_cells", "cell_frac")

read_egenes <- function(manifest) {
  manifest |>
    rowwise() |>
    mutate(gene_tbl = list(read_tsv(region_file, show_col_types = FALSE))) |>
    unnest(gene_tbl) |>
    ungroup() |>
    mutate(
      region_bh = p.adjust(pvalue, method = "BH"),
      .by = all_of(group_cols)
    ) |>
    filter(region_bh < 0.05) |>
    select(all_of(group_cols), feature_id) |>
    distinct()
}

egenes <- bind_rows(
  pb_manifest |> mutate(method = "lm"),
  sc_manifest |> mutate(method = model)
) |>
  read_egenes()

min_n_cells_levels <- sort(unique(c(pb_manifest$min_n_cells, sc_manifest$min_n_cells)))
filter_labels <- if_else(
  min_n_cells_levels < 0,
  "All donors",
  sprintf("Donors with >= %d cells", min_n_cells_levels)
)

grid <- bind_rows(
  pb_manifest |> distinct(cell_type, min_n_cells, cell_frac) |> mutate(method = "lm"),
  sc_manifest |> distinct(cell_type, min_n_cells, cell_frac, model) |> rename(method = model)
) |>
  mutate(
    method_label = factor(method_lookup[method], levels = method_levels),
    filter_label = factor(
      filter_labels[match(min_n_cells, min_n_cells_levels)],
      levels = filter_labels
    ),
    cell_type_label = factor(
      cell_type_lookup[cell_type],
      levels = unname(cell_type_lookup[cell_types])
    )
  )

count_data <- grid |>
  left_join(
    egenes |> summarise(n_egene = n(), .by = all_of(group_cols)),
    by = group_cols
  ) |>
  mutate(n_egene = replace_na(n_egene, 0L))

method_cols <- method_col_lookup[method_levels]

cell_frac_breaks <- sort(unique(grid$cell_frac))
cell_frac_labels <- percent(cell_frac_breaks, accuracy = NULL, drop0trailing = TRUE)

plot_cell_frac <- function(data, y, y_lab) {
  data |>
    ggplot(aes(
      cell_frac, {{ y }},
      colour = method_label,
      linetype = filter_label,
      shape = filter_label,
      group = interaction(method_label, filter_label)
    )) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_x_log10(breaks = cell_frac_breaks, labels = cell_frac_labels) +
    scale_colour_manual(values = method_cols) +
    labs(
      x = "Cell fraction",
      y = y_lab,
      colour = "Method",
      linetype = "Donor filter",
      shape = "Donor filter"
    ) +
    theme_jp() +
    theme(legend.position = "right")
}

p_counts <- plot_cell_frac(count_data, n_egene, "Number of eGenes") +
  facet_wrap(~cell_type_label, scales = "free_y")

ggsave(
  "power-cell-frac-counts-plot.pdf",
  p_counts,
  width = 12,
  height = 5
)

truth <- egenes |>
  filter(method == "lm", cell_frac == 1, min_n_cells < 0) |>
  select(cell_type, feature_id) |>
  distinct()

truth_n <- truth |>
  summarise(n_truth = n(), .by = cell_type)

recall_data <- grid |>
  left_join(
    egenes |>
      inner_join(truth, by = c("cell_type", "feature_id")) |>
      summarise(n_recovered = n(), .by = all_of(group_cols)),
    by = group_cols
  ) |>
  left_join(truth_n, by = "cell_type") |>
  mutate(
    n_recovered = replace_na(n_recovered, 0L),
    recall = if_else(n_truth > 0, n_recovered / n_truth, NA_real_)
  )

p_recall <- plot_cell_frac(recall_data, recall, "Proportion of full LM eGenes recovered") +
  scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
  facet_wrap(~cell_type_label)

ggsave(
  "power-cell-frac-recall-plot.pdf",
  p_recall,
  width = 12,
  height = 5
)
