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

sc_data_files <- read_tsv(args[1], show_col_types = FALSE)

compute_qq_data <- function(files) {

  pvalue <- unlist(map(
    files,
    function(x) read_tsv(x, show_col_types = FALSE)$pvalue)
  )
  pvalue <- pvalue[!is.na(pvalue)]
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
      .by = x_bin
    ) |>
    mutate(log_x_bin_mid = seq(0.05, 5.95, by = 0.1)[as.numeric(x_bin)])
}

compute_qq_data_by_prop <- function(variant_files, prop_file, prop) {

  variant_data <- bind_rows(!!!map(
    variant_files,
    function(x) {
      read_tsv(x, show_col_types = FALSE) |>
        select(feature_id, pvalue)
    })
  )

  prop_data <- read_tsv(prop_file, show_col_types = FALSE) |>
    select(feature_id, all_of(prop))

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

  prop <- pvalue_data[[prop]]
  tibble(x_bin, y_pvalue, prop, lower_ci, upper_ci) |>
    summarise(
      log_y_pvalue = -log10(mean(y_pvalue)),
      log_lower_ci = -log10(mean(lower_ci)),
      log_upper_ci = -log10(mean(upper_ci)),
      .by = c(x_bin, prop)
    ) |>
    mutate(log_x_bin_mid = seq(0.05, 5.95, by = 0.1)[as.numeric(x_bin)])
}

compute_qq_data_by_maf <- function(variant_files) {

  variant_data <- bind_rows(!!!map(
    variant_files,
    function(x) {
      read_tsv(x, show_col_types = FALSE) |>
        select(maf, pvalue)
    })
  ) |>
    filter(!is.na(pvalue))

  pvalue <- variant_data$pvalue
  maf <- variant_data$maf
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

  maf_bin <- cut(
    maf,
    breaks = c(0, 0.05, 0.10, 0.50),
    labels = c("<5%", "5-10%", "10-50%"),
    right = FALSE,
    include.lowest = TRUE
  )

  tibble(x_bin, y_pvalue, maf_bin, lower_ci, upper_ci) |>
    summarise(
      log_y_pvalue = -log10(mean(y_pvalue)),
      log_lower_ci = -log10(mean(lower_ci)),
      log_upper_ci = -log10(mean(upper_ci)),
      .by = c(x_bin, maf_bin)
    ) |>
    mutate(log_x_bin_mid = seq(0.05, 5.95, by = 0.1)[as.numeric(x_bin)])
}

overall_plot_data <- sc_data_files |>
  filter(cell_frac == 1) |>
  filter(indiv_frac == 1) |>
  filter(count_frac == 1) |>
  summarise(file_list = list(variant_file), .by = c(cell_type, model)) |>
  rowwise() |>
  mutate(qq_data = list(compute_qq_data(file_list))) |>
  ungroup() |>
  select(-file_list) |>
  unnest(cols = qq_data)

overall_p <- overall_plot_data |>
  filter(cell_type %in% names(cell_type_lookup)) |>
  filter(model %in% names(method_lookup)) |>
  mutate(
    cell_type = factor(
      cell_type_lookup[cell_type],
      levels = cell_type_lookup[c("CD4_NC", "B_IN", "Plasma")]
    ),
    model = factor(
      method_lookup[model],
      levels = unique(unname(method_lookup[c("p_glmm_sc", "lmm_sc")]))
    )
  ) |>
  ggplot(aes(log_x_bin_mid, log_y_pvalue,
             ymin = log_lower_ci, ymax = log_upper_ci,
             colour = cell_type)) +
  geom_point(alpha = 0.8) +
  geom_abline(linetype = "dashed") +
  geom_ribbon(linetype = 2, alpha = 0.1) +
  scale_colour_manual(
    values = setNames(cell_type_cols, cell_type_lookup[names(cell_type_cols)])
  ) +
  facet_grid(vars(model), vars(cell_type)) +
  labs(
    x = "Expected -log10(p-value)",
    y = "Observed -log10(p-value)",
    colour = "Cell type"
  ) +
  theme_jp() +
  theme(legend.position = "none")

ggsave(
  "perm-overall-plot.pdf",
  overall_p,
  width = 16,
  height = 10
)

#prop_plot_data <- sc_data_files |>
#  filter(cell_frac == 1) |>
#  filter(indiv_frac == 1) |>
#  summarise(
#    file_list = list(variant_file),
#    prop_file = first(prop_file),
#    .by = cell_type
#  ) |>
#  rowwise() |>
#  mutate(qq_data = list(compute_qq_data_by_prop(file_list, prop_file, "sc_non_zero_frac"))) |>
#  ungroup() |>
#  select(-file_list) |>
#  unnest(cols = qq_data) |>
#  mutate(
#    prop_fct = cut(
#      prop,
#      breaks = c(0.005, 0.01, 0.05, 0.1, 0.5, Inf),
#      labels = c("0.5-1%", "1-5%", "5-10%", "10-50%", "50%-100%"),
#      right = FALSE,
#      include.lowest = TRUE
#    )
#  )

#prop_p <- prop_plot_data |>
#  ggplot(aes(log_x_bin_mid, log_y_pvalue,
#             ymin = log_lower_ci, ymax = log_upper_ci,
#             colour = cell_type)) +
#  geom_point(alpha = 0.8) +
#  geom_abline(linetype = "dashed") +
#  geom_ribbon(linetype = 2, alpha = 0.1) +
#  facet_grid(vars(cell_type), vars(prop_fct)) + 
#  labs(
#    x = "Expected -log10(p-value)",
#    y = "Observed -log10(p-value)",
#    colour = "Cell type"
#  )

#ggsave(
#  "perm-prop-plot.pdf",
#  prop_p,
#  width = 12,
#  height = 10
#)

#maf_plot_data <- sc_data_files |>
#  filter(cell_frac == 1) |>
#  filter(indiv_frac == 1) |>
#  summarise(
#    file_list = list(variant_file),
#    prop_file = first(prop_file),
#    .by = cell_type
#  ) |>
#  rowwise() |>
#  mutate(qq_data = list(compute_qq_data_by_maf(file_list))) |>
#  ungroup() |>
#  select(-file_list) |>
#  unnest(cols = qq_data)

#maf_p <- maf_plot_data |>
#  ggplot(aes(log_x_bin_mid, log_y_pvalue,
#             ymin = log_lower_ci, ymax = log_upper_ci,
#             colour = cell_type)) +
#  geom_point(alpha = 0.8) +
#  geom_abline(linetype = "dashed") +
#  geom_ribbon(linetype = 2, alpha = 0.1) +
#  facet_grid(vars(cell_type), vars(maf_bin)) + 
#  labs(
#    x = "Expected -log10(p-value)",
#    y = "Observed -log10(p-value)",
#    colour = "Cell type"
#  )

#ggsave(
#  "perm-maf-plot.pdf",
#  maf_p,
#  width = 12,
#  height = 10
#)

#cell_frac_plot_data <- sc_data_files |>
#  filter(indiv_frac == 1) |>
#  mutate(cell_frac = factor(cell_frac)) |>
#  summarise(file_list = list(variant_file), .by = c(cell_frac, cell_type)) |>
#  rowwise() |>
#  mutate(qq_data = list(compute_qq_data(file_list))) |>
#  ungroup() |>
#  select(-file_list) |>
#  unnest(cols = qq_data)

#cell_frac_p <- cell_frac_plot_data |>
#  mutate(cell_type = factor(cell_type, levels = c("CD4_NC", "B_IN", "Plasma"))) |>
#  ggplot(aes(log_x_bin_mid, log_y_pvalue,
#             ymin = log_lower_ci, ymax = log_upper_ci)) +
#  geom_point(alpha = 0.8) +
#  geom_abline(linetype = "dashed") +
#  geom_ribbon(linetype = 2, alpha = 0.1) +
#  facet_grid2(vars(cell_type), vars(cell_frac), scales = "free_y", independent = "y") + 
#  labs(
#    x = "Expected -log10(p-value)",
#    y = "Observed -log10(p-value)"
#  )

#ggsave(
#  "perm-cell-frac-plot.pdf",
#  cell_frac_p,
#  width = 14,
#  height = 10
#)

#indiv_frac_plot_data <- sc_data_files |>
#  filter(cell_type == "CD4_NC") |>
#  filter(cell_frac == 1) |>
#  mutate(indiv_frac = factor(indiv_frac)) |>
#  summarise(file_list = list(variant_file), .by = indiv_frac) |>
#  rowwise() |>
#  mutate(qq_data = list(compute_qq_data(file_list))) |>
#  ungroup() |>
#  select(-file_list) |>
#  unnest(cols = qq_data)

#indiv_frac_p <- indiv_frac_plot_data |>
#  ggplot(aes(log_x_bin_mid, log_y_pvalue,
#             ymin = log_lower_ci, ymax = log_upper_ci)) +
#  geom_point(alpha = 0.8) +
#  geom_abline(linetype = "dashed") +
#  geom_ribbon(linetype = 2, alpha = 0.1) +
#  facet_wrap(~indiv_frac) +
#  labs(
#    x = "Expected -log10(p-value)",
#    y = "Observed -log10(p-value)"
#  )

#ggsave(
#  "perm-indiv-frac-plot.pdf",
#  indiv_frac_p,
#  width = 12,
#  height = 10
#)