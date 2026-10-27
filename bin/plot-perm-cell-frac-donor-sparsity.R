#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(tidyr)
  library(ggh4x)
  library(data.table)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)
sc_data_files <- read_tsv(args[1], show_col_types = FALSE)

cell_types <- c("B_IN", "CD4_NC")
sparsity_levels <- c("0-5%", "5-10%", "10-25%", "25-50%", "50-75%", "75-100%")

compute_qq_data <- function(pvalue) {
  n <- length(pvalue)
  m <- (1:n) / (n + 1)
  c <- abs(qnorm(0.05 / 2))
  v <- (1:n) * (n - (1:n) + 1) / (n + 1)^2 / (n + 2)
  s <- sqrt(v)

  tibble(
    x_bin = cut(-log10(m), breaks = seq(0, 6, by = 0.1), include.lowest = TRUE),
    y_pvalue = sort(pvalue),
    lower_ci = m - c * s,
    upper_ci = m + c * s
  ) |>
    summarise(
      log_y_pvalue = -log10(mean(y_pvalue)),
      log_lower_ci = -log10(mean(lower_ci)),
      log_upper_ci = -log10(mean(upper_ci)),
      .by = x_bin
    ) |>
    mutate(log_x_bin_mid = seq(0.05, 5.95, by = 0.1)[as.numeric(x_bin)])
}

read_binned_pvalues <- function(files, prop_file) {
  props <- read_tsv(prop_file, show_col_types = FALSE) |>
    select(feature_id, pb_non_zero_frac)

  map(files, \(x) fread(x, select = c("feature_id", "pvalue"))) |>
    bind_rows() |>
    filter(!is.na(pvalue)) |>
    inner_join(props, by = "feature_id") |>
    mutate(
      donor_sparsity = cut(
        pb_non_zero_frac,
        breaks = c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 1),
        labels = sparsity_levels,
        right = FALSE,
        include.lowest = TRUE
      )
    ) |>
    select(feature_id, pb_non_zero_frac, donor_sparsity, pvalue)
}

pvalues <- sc_data_files |>
  filter(
    model == "p_glmm_sc",
    cov_spec == "bulk_pca",
    as.character(k) == "none",
    as.numeric(indiv_frac) == 1,
    as.numeric(count_frac) == 1,
    as.numeric(min_n_cells) < 0,
    cell_type %in% cell_types
  ) |>
  summarise(
    file_list = list(variant_file),
    prop_file = first(prop_file),
    .by = c(cell_type, cell_frac)
  ) |>
  mutate(pvalues = map2(file_list, prop_file, read_binned_pvalues)) |>
  select(cell_type, cell_frac, pvalues) |>
  unnest(pvalues)

frac_levels <- sort(unique(pvalues$cell_frac))

to_facets <- function(df) {
  df |>
    mutate(
      cell_type = factor(cell_type_lookup[cell_type], levels = cell_type_lookup[c("CD4_NC", "B_IN")]),
      cell_frac = factor(cell_frac, levels = frac_levels)
    )
}

qq_by <- function(df, by) {
  df |>
    nest(.by = all_of(by)) |>
    mutate(qq_data = map(data, \(d) compute_qq_data(d$pvalue))) |>
    select(-data) |>
    unnest(qq_data) |>
    to_facets()
}

plot_data <- pvalues |>
  qq_by(c("cell_type", "cell_frac", "donor_sparsity")) |>
  mutate(donor_sparsity = factor(donor_sparsity, levels = sparsity_levels))

p <- plot_data |>
  ggplot(aes(
    log_x_bin_mid, log_y_pvalue,
    ymin = log_lower_ci, ymax = log_upper_ci,
    colour = donor_sparsity, fill = donor_sparsity
  )) +
  geom_point(alpha = 0.8, size = 0.8) +
  geom_abline(linetype = "dashed") +
  geom_ribbon(linetype = 2, alpha = 0.05, colour = NA) +
  scale_colour_viridis_d(end = 0.9, drop = FALSE) +
  scale_fill_viridis_d(end = 0.9, drop = FALSE) +
  facet_grid2(
    vars(cell_type), vars(cell_frac),
    scales = "free_y",
    independent = "y"
  ) +
  labs(
    title = sprintf("%s by donor-level non-zero fraction", method_lookup[["p_glmm_sc"]]),
    x = "Expected -log10(p-value)",
    y = "Observed -log10(p-value)",
    colour = "Donors with\nnon-zero counts",
    fill = "Donors with\nnon-zero counts"
  ) +
  theme_jp()

ggsave(
  "perm-cell-frac-p_glmm_sc-donor-sparsity-plot.pdf",
  p,
  width = max(12, 2.2 * length(frac_levels)),
  height = 8,
  limitsize = FALSE
)

filt_pvalues <- pvalues |>
  filter(pb_non_zero_frac > 0.1)

gene_counts <- pvalues |>
  summarise(
    n_kept = n_distinct(feature_id[pb_non_zero_frac > 0.1]),
    n_total = n_distinct(feature_id),
    .by = c(cell_type, cell_frac)
  ) |>
  to_facets()

p_filt <- filt_pvalues |>
  qq_by(c("cell_type", "cell_frac")) |>
  ggplot(aes(
    log_x_bin_mid, log_y_pvalue,
    ymin = log_lower_ci, ymax = log_upper_ci,
    colour = cell_type
  )) +
  geom_point(alpha = 0.8, size = 0.8) +
  geom_abline(linetype = "dashed") +
  geom_ribbon(linetype = 2, alpha = 0.1) +
  geom_text(
    data = gene_counts,
    aes(x = 0, y = Inf, label = sprintf("%d / %d genes", n_kept, n_total)),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.5,
    colour = "grey40",
    size = 3
  ) +
  scale_colour_manual(
    values = setNames(cell_type_cols, cell_type_lookup[names(cell_type_cols)])
  ) +
  facet_grid2(
    vars(cell_type), vars(cell_frac),
    scales = "free_y",
    independent = "y"
  ) +
  labs(
    title = sprintf("%s, genes with non-zero counts in > 10%% of donors", method_lookup[["p_glmm_sc"]]),
    x = "Expected -log10(p-value)",
    y = "Observed -log10(p-value)"
  ) +
  theme_jp() +
  theme(legend.position = "none")

ggsave(
  "perm-cell-frac-p_glmm_sc-donor-nz10-plot.pdf",
  p_filt,
  width = max(12, 2.2 * length(frac_levels)),
  height = 8,
  limitsize = FALSE
)
