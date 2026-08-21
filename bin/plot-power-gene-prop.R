#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(stringr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)

cell_types <- c("Plasma", "B_IN", "CD4_NC")

props <- c(
  "sc_mean" = "Single-cell mean",
  "sc_var" = "Single-cell variance",
  "sc_cv" = "Single-cell CV",
  "sc_non_zero_frac" = "Single-cell non-zero fraction",
  "pb_mean" = "Pseudobulk mean",
  "pb_var" = "Pseudobulk variance",
  "pb_cv" = "Pseudobulk CV",
  "pb_non_zero_frac" = "Pseudobulk non-zero fraction"
)

log_props <- c("sc_mean", "sc_var", "sc_cv", "pb_mean", "pb_var", "pb_cv")

pb_manifest <- pb_data_files |>
  filter(
    model == "lm",
    int_cov == "none",
    indiv_frac == 1,
    cell_frac == 1,
    cell_type %in% cell_types
  )

sc_manifest <- sc_data_files |>
  filter(
    model == "p_glmm_sc",
    cov_spec == "bulk_pca",
    k == "none",
    indiv_frac == 1,
    cell_frac == 1,
    cell_type %in% cell_types
  )

read_egenes <- function(manifest, method) {
  manifest |>
    mutate(method = method) |>
    rowwise() |>
    mutate(gene_tbl = list(read_tsv(region_file, show_col_types = FALSE))) |>
    unnest(gene_tbl) |>
    ungroup() |>
    mutate(
      region_bh = p.adjust(pvalue, method = "BH"),
      .by = c(cell_type, method)
    ) |>
    filter(region_bh < 0.05) |>
    select(cell_type, method, feature_id) |>
    distinct()
}

egenes <- bind_rows(
  read_egenes(pb_manifest, "lm"),
  read_egenes(sc_manifest, "p_glmm_sc")
)

overlap_levels <- c(
  "Unique to Pseudobulk LM",
  "Common to both",
  "Unique to Single-cell P-GLMM"
)

overlap_cols <- c(
  "Unique to Pseudobulk LM" = unname(method_col_lookup[["Pseudobulk LM"]]),
  "Common to both" = "#4477AA",
  "Unique to Single-cell P-GLMM" = unname(method_col_lookup[["Single-cell P-GLMM"]])
)

egenes_overlap <- egenes |>
  mutate(present = TRUE) |>
  pivot_wider(
    names_from = method,
    values_from = present,
    values_fill = list(present = FALSE)
  ) |>
  mutate(
    overlap = case_when(
      lm & p_glmm_sc ~ "Common to both",
      lm & !p_glmm_sc ~ "Unique to Pseudobulk LM",
      !lm & p_glmm_sc ~ "Unique to Single-cell P-GLMM"
    ),
    overlap = factor(overlap, levels = overlap_levels)
  ) |>
  select(cell_type, feature_id, overlap)

gene_props <- bind_rows(
  pb_manifest |> distinct(cell_type, gene_prop_file),
  sc_manifest |> distinct(cell_type, gene_prop_file)
) |>
  distinct(cell_type, gene_prop_file) |>
  rowwise() |>
  mutate(props_tbl = list(read_tsv(gene_prop_file, show_col_types = FALSE))) |>
  unnest(props_tbl) |>
  ungroup() |>
  distinct(cell_type, feature_id, .keep_all = TRUE)

plot_prop_power <- function(prop_col, prop_label) {
  plot_data <- egenes_overlap |>
    inner_join(gene_props, by = c("cell_type", "feature_id")) |>
    filter(is.finite(.data[[prop_col]])) |>
    mutate(
      prop_value = .data[[prop_col]],
      cell_type_label = factor(
        cell_type_lookup[cell_type],
        levels = unname(cell_type_lookup[cell_types])
      )
    )

  p <- ggplot(plot_data, aes(overlap, prop_value, fill = overlap)) +
    geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3, linewidth = 0.4) +
    scale_fill_manual(values = overlap_cols, guide = "none") +
    facet_wrap(~cell_type_label, scales = "free_y") +
    labs(
      x = NULL,
      y = prop_label
    ) +
    theme_jp() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1)
    )

  if (prop_col %in% log_props) {
    p <- p + scale_y_log10()
  }

  p
}

walk2(names(props), unname(props), function(prop_col, prop_label) {
  out_file <- paste0(
    "power-gene-prop-",
    str_replace_all(prop_col, "_", "-"),
    "-plot.pdf"
  )
  ggsave(
    out_file,
    plot_prop_power(prop_col, prop_label),
    width = 12,
    height = 5
  )
})
