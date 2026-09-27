#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(tidyr)
  library(ggh4x)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

sc_data_files <- read_tsv(args[1], show_col_types = FALSE)
pb_data_files <- read_tsv(args[2], show_col_types = FALSE)

models <- c("lm", "p_glmm_sc", "lmm_sc")
nz_breaks <- c(0.01, 0.05, 0.10, 0.50, 1)
nz_labels <- c("1-5%", "5-10%", "10-50%", "50-100%")

compute_qq_data <- function(pvalue) {
  pvalue <- pvalue[!is.na(pvalue)]
  n <- length(pvalue)
  if (n == 0) {
    return(tibble(
      log_x_bin_mid = numeric(),
      log_y_pvalue = numeric(),
      log_lower_ci = numeric(),
      log_upper_ci = numeric()
    ))
  }

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

compute_qq_by_nz <- function(variant_files, prop_file) {
  variant_data <- map(
    variant_files,
    function(x) {
      read_tsv(x, show_col_types = FALSE) |>
        select(feature_id, pvalue)
    }
  ) |>
    list_rbind()

  prop_data <- read_tsv(prop_file, show_col_types = FALSE) |>
    select(feature_id, sc_non_zero_frac)

  pvalue_data <- variant_data |>
    left_join(prop_data, by = "feature_id") |>
    filter(!is.na(pvalue), !is.na(sc_non_zero_frac)) |>
    mutate(
      nz_bin = as.character(cut(
        sc_non_zero_frac,
        breaks = nz_breaks,
        labels = nz_labels,
        include.lowest = TRUE,
        right = FALSE
      ))
    ) |>
    filter(!is.na(nz_bin))

  pvalue_data |>
    group_split(nz_bin, .keep = TRUE) |>
    map(function(df) {
      compute_qq_data(df$pvalue) |>
        mutate(nz_bin = df$nz_bin[[1]])
    }) |>
    list_rbind()
}

plot_perm_nz_frac <- function(plot_data, model_label, out_file) {
  ct_levels <- unname(cell_type_lookup[
    names(cell_type_lookup) %in% unique(plot_data$cell_type)
  ])

  p <- plot_data |>
    mutate(
      cell_type = factor(
        cell_type_lookup[cell_type],
        levels = ct_levels
      ),
      nz_bin = factor(nz_bin, levels = nz_labels)
    ) |>
    ggplot(aes(
      log_x_bin_mid, log_y_pvalue,
      ymin = log_lower_ci, ymax = log_upper_ci,
      colour = cell_type
    )) +
    geom_point(alpha = 0.8, size = 0.8) +
    geom_abline(linetype = "dashed") +
    geom_ribbon(linetype = 2, alpha = 0.1) +
    scale_colour_manual(
      values = setNames(cell_type_cols, cell_type_lookup[names(cell_type_cols)])
    ) +
    facet_grid2(
      vars(cell_type), vars(nz_bin),
      scales = "free_y",
      independent = "y"
    ) +
    labs(
      title = model_label,
      x = "Expected -log10(p-value)",
      y = "Observed -log10(p-value)"
    ) +
    theme_jp() +
    theme(legend.position = "none")

  n_ct <- n_distinct(plot_data$cell_type)
  n_bin <- n_distinct(plot_data$nz_bin)
  ggsave(
    out_file,
    p,
    width = max(12, 2.2 * n_bin),
    height = max(10, 2.5 * n_ct),
    limitsize = FALSE
  )
}

build_sc_plot_data <- function(model_name) {
  sc_data_files |>
    filter(
      dataset == "tremor",
      model == model_name,
      cov_spec == "bulk_pca",
      k == "none",
      cell_frac == 1,
      indiv_frac == 1,
      count_frac == 1
    ) |>
    summarise(
      file_list = list(variant_file),
      prop_file = first(prop_file),
      .by = cell_type
    ) |>
    rowwise() |>
    mutate(qq_data = list(compute_qq_by_nz(file_list, prop_file))) |>
    ungroup() |>
    select(-file_list, -prop_file) |>
    unnest(cols = qq_data)
}

build_pb_plot_data <- function(model_name) {
  pb_data_files |>
    filter(
      dataset == "tremor",
      model == model_name,
      int_cov == "none",
      cell_frac == 1,
      indiv_frac == 1,
      count_frac == 1
    ) |>
    summarise(
      file_list = list(variant_file),
      prop_file = first(prop_file),
      .by = cell_type
    ) |>
    rowwise() |>
    mutate(qq_data = list(compute_qq_by_nz(file_list, prop_file))) |>
    ungroup() |>
    select(-file_list, -prop_file) |>
    unnest(cols = qq_data)
}

for (model_name in models) {
  plot_data <- if (model_name %in% c("p_glmm_sc", "lmm_sc")) {
    build_sc_plot_data(model_name)
  } else {
    build_pb_plot_data(model_name)
  }

  plot_perm_nz_frac(
    plot_data,
    method_lookup[[model_name]],
    sprintf("tremor-perm-nz-frac-%s-plot.pdf", model_name)
  )
}
