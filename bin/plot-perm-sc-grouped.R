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

args <- commandArgs(trailingOnly = TRUE)

sc_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(k != "none")

compute_qq_data <- function(variant_files, prop_file, pvalue_col) {

  variant_data <- bind_rows(!!!map(
    variant_files,
    function(x) {
      read_tsv(x, show_col_types = FALSE) |>
        select(feature_id, maf, all_of(pvalue_col))
    })
  )

  prop_data <- read_tsv(prop_file, show_col_types = FALSE)

  pvalue_data <- left_join(
    variant_data,
    prop_data,
    by = "feature_id"
  ) |>
    filter(maf > 0.1) |>
    filter(pb_non_zero_frac > 0.1)

  pvalue <- pvalue_data[[pvalue_col]]
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
      .by = c(x_bin)
    ) |>
    mutate(log_x_bin_mid = seq(0.05, 5.95, by = 0.1)[as.numeric(x_bin)])
}

compute_qq_data_het <- function(variant_files, prop_file) {
  compute_qq_data(variant_files, prop_file, "group_linear_pvalue")
}

compute_qq_data_group1 <- function(variant_files, prop_file, int_cov) {
  col <- paste0(int_cov, "_q1_pvalue")
  compute_qq_data(variant_files, prop_file, col)
}

compute_qq_data_pvalue <- function(variant_files, prop_file) {
  compute_qq_data(variant_files, prop_file, "pvalue")
}

build_qq_plot_data <- function(data_files, qq_fn) {
  data_files |>
    filter(cell_frac == 1) |>
    filter(indiv_frac == 1) |>
    summarise(
      file_list = list(variant_file),
      prop_file = first(prop_file),
      int_cov = first(int_cov),
      .by = c(cell_type, int_cov, k)
    ) |>
    rowwise() |>
    mutate(qq_data = list(qq_fn(file_list, prop_file, int_cov))) |>
    ungroup() |>
    select(-file_list) |>
    unnest(cols = qq_data)
}

het_plot_data <- build_qq_plot_data(sc_data_files, function(files, prop, cov) compute_qq_data_het(files, prop)) |>
  mutate(k_label = paste0("K=", k))

group1_plot_data <- build_qq_plot_data(sc_data_files, compute_qq_data_group1) |>
  mutate(k_label = paste0("K=", k))

pvalue_plot_data <- build_qq_plot_data(sc_data_files, function(files, prop, cov) compute_qq_data_pvalue(files, prop)) |>
  mutate(k_label = paste0("K=", k))

het_p <- het_plot_data |>
  ggplot(aes(log_x_bin_mid, log_y_pvalue,
             ymin = log_lower_ci, ymax = log_upper_ci,
             colour = k_label)) +
  geom_point(alpha = 0.8) +
  geom_abline(linetype = "dashed") +
  geom_ribbon(linetype = 2, alpha = 0.1) +
  facet_grid(vars(cell_type, int_cov), vars(k_label)) +
  labs(
    x = "Expected -log10(het_pvalue)",
    y = "Observed -log10(het_pvalue)",
    colour = "Grouping"
  )

group1_p <- group1_plot_data |>
  ggplot(aes(log_x_bin_mid, log_y_pvalue,
             ymin = log_lower_ci, ymax = log_upper_ci,
             colour = k_label)) +
  geom_point(alpha = 0.8) +
  geom_abline(linetype = "dashed") +
  geom_ribbon(linetype = 2, alpha = 0.1) +
  facet_grid(vars(cell_type, int_cov), vars(k_label)) +
  labs(
    x = "Expected -log10(q1_pvalue)",
    y = "Observed -log10(q1_pvalue)",
    colour = "Grouping"
  )

pvalue_p <- pvalue_plot_data |>
  ggplot(aes(log_x_bin_mid, log_y_pvalue,
             ymin = log_lower_ci, ymax = log_upper_ci,
             colour = k_label)) +
  geom_point(alpha = 0.8) +
  geom_abline(linetype = "dashed") +
  geom_ribbon(linetype = 2, alpha = 0.1) +
  facet_grid(vars(cell_type, int_cov), vars(k_label)) +
  labs(
    x = "Expected -log10(pvalue)",
    y = "Observed -log10(pvalue)",
    colour = "Grouping"
  )

ggsave(
  "perm-sc-grouped-het-plot.pdf",
  het_p,
  width = 12,
  height = 10
)

ggsave(
  "perm-sc-grouped-q1-plot.pdf",
  group1_p,
  width = 12,
  height = 10
)

ggsave(
  "perm-sc-grouped-pvalue-plot.pdf",
  pvalue_p,
  width = 12,
  height = 10
)
