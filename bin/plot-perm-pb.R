#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(qvalue)
  library(purrr)
  library(stringr)
  library(tidyr)
  library(ggh4x)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

pb_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(int_cov == "none")

compute_qq_data_pb <- function(variant_files, prop_file) {

  variant_data <- bind_rows(!!!map(
    variant_files,
    function(x) {
      read_tsv(x, show_col_types = FALSE) |>
        select(feature_id, pvalue, maf)
    })
  )

  prop_data <- read_tsv(prop_file, show_col_types = FALSE)

  pvalue_data <- left_join(
    variant_data,
    prop_data,
    by = "feature_id"
  ) |>
    filter(!is.na(pvalue))

  pvalue <- pvalue_data$pvalue
  n <- length(pvalue)
  m <- (1:n) / (n + 1)
  c <- abs(qnorm(0.05 / 2))
  v <- (1:n) * (n - (1:n) + 1) / (n + 1)^2 / (n + 2)
  s <- sqrt(v)
  lower_ci <- m - c * s
  upper_ci <- m + c * s

  log_x_pvalue <- -log10(m)
  y_pvalue <- sort(pvalue)
  x_bin <- cut(
    log_x_pvalue,
    breaks = seq(0, 6, by = 0.1),
    include.lowest = TRUE
  )

  tibble(x_bin, y_pvalue, lower_ci, upper_ci) |>
    summarise(
      log_y_pvalue = -log10(mean(y_pvalue)),
      log_lower_ci = -log10(mean(lower_ci)),
      log_upper_ci = -log10(mean(upper_ci)),
      .by = c(x_bin)
    ) |>
    mutate(log_x_bin_mid = seq(0.05, 5.95, by = 0.1)[as.numeric(x_bin)])
}

pb_qq_data <- pb_data_files |>
  filter(cell_frac == 1, indiv_frac == 1, count_frac == 1) |>
  summarise(
    file_list = list(variant_file),
    prop_file = first(prop_file),
    .by = c(dataset, cell_type, model)
  ) |>
  rowwise() |>
  mutate(qq_data = list(compute_qq_data_pb(file_list, prop_file))) |>
  ungroup() |>
  select(-file_list) |>
  unnest(cols = qq_data) |>
  filter(cell_type %in% names(cell_type_lookup))

ct_levels <- unname(cell_type_lookup[
  names(cell_type_lookup) %in% unique(pb_qq_data$cell_type)
])

for (model_name in unique(pb_qq_data$model)) {
  plot_data <- pb_qq_data |>
    filter(model == model_name) |>
    mutate(cell_type = factor(unname(cell_type_lookup[cell_type]), levels = ct_levels))
  n_ct <- n_distinct(plot_data$cell_type)
  n_col <- min(n_ct, 4)

  p <- plot_data |>
    ggplot(aes(log_x_bin_mid, log_y_pvalue,
               ymin = log_lower_ci, ymax = log_upper_ci,
               colour = cell_type)) +
    geom_point(alpha = 0.8) +
    geom_abline(linetype = "dashed") +
    geom_ribbon(linetype = 2, alpha = 0.1) +
    scale_colour_manual(
      values = setNames(cell_type_cols, cell_type_lookup[names(cell_type_cols)]),
      drop = TRUE
    ) +
    facet_wrap(vars(cell_type), ncol = n_col) +
    labs(
      title = coalesce(unname(method_lookup[model_name]), model_name),
      x = "Expected -log10(p-value)",
      y = "Observed -log10(p-value)"
    ) +
    theme_jp() +
    theme(legend.position = "none")

  ggsave(
    paste0("perm-pb-", model_name, "-plot.pdf"),
    p,
    width = 4 * n_col,
    height = 4 * ceiling(n_ct / n_col) + 1
  )
}

#pb_overall_int_plot_data <- pb_data_files |>
#  filter(cell_frac == 1) |>
#  filter(indiv_frac == 1) |>
#  summarise(file_list = list(variant_file), .by = cell_type) |>
#  rowwise() |>
#  mutate(qq_data = list(compute_qq_data_pb_int(file_list, "age"))) |>
#  ungroup() |>
#  select(-file_list) |>
#  unnest(cols = qq_data)

#pb_overall_int_p <- pb_overall_int_plot_data |>
#  ggplot(aes(log_x_bin_mid, log_y_pvalue,
#             ymin = log_lower_ci, ymax = log_upper_ci,
#             colour = cell_type)) +
#  geom_point(alpha = 0.8) +
#  geom_abline(linetype = "dashed") +
#  geom_ribbon(linetype = 2, alpha = 0.1) +
#  facet_wrap(~cell_type) +
#  labs(
#    x = "Expected -log10(p-value)",
#    y = "Observed -log10(p-value)",
#    colour = "Cell type"
#  )

#pb_cell_frac_plot_data <- pb_data_files |>
#  filter(indiv_frac == 1) |>
#  mutate(cell_frac = factor(cell_frac)) |>
#  summarise(
#    file_list = list(variant_file),
#    prop_file = first(prop_file),
#    .by = c(cell_frac, cell_type)
#  ) |>
#  rowwise() |>
#  mutate(qq_data = list(compute_qq_data_pb(file_list, prop_file))) |>
#  ungroup() |>
#  select(-file_list) |>
#  unnest(cols = qq_data)

#pb_cell_frac_p <- pb_cell_frac_plot_data |>
#  mutate(cell_type = factor(cell_type,
#                            levels = c("CD4_NC", "B_IN", "Plasma"))) |>
#  ggplot(aes(log_x_bin_mid, log_y_pvalue,
#             ymin = log_lower_ci, ymax = log_upper_ci)) +
#  geom_point(alpha = 0.8) +
#  geom_abline(linetype = "dashed") +
#  geom_ribbon(linetype = 2, alpha = 0.1) +
#  facet_grid2(vars(cell_type), vars(cell_frac),
#              scales = "free_y", independent = "y") +
#  labs(
#    x = "Expected -log10(p-value)",
#    y = "Observed -log10(p-value)"
#  )