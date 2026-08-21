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

cell_types <- c("Plasma", "B_IN", "CD4_NC")
models <- c("lm", "nb_glm", "p_glmm_sc", "lmm_sc")

compute_qq_data <- function(files) {
  pvalue <- unlist(map(
    files,
    function(x) {
      read_tsv(x, show_col_types = FALSE)$pvalue
    }
  ))
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

plot_perm_cell_frac <- function(plot_data, model_label, out_file) {
  p <- plot_data |>
    mutate(
      cell_type = factor(
        cell_type_lookup[cell_type],
        levels = cell_type_lookup[c("CD4_NC", "B_IN", "Plasma")]
      ),
      cell_frac = factor(cell_frac, levels = sort(unique(cell_frac)))
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
      vars(cell_type), vars(cell_frac),
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

  n_frac <- n_distinct(plot_data$cell_frac)
  ggsave(
    out_file,
    p,
    width = max(12, 2.2 * n_frac),
    height = 10,
    limitsize = FALSE
  )
}

build_sc_plot_data <- function(model_name) {
  sc_data_files |>
    filter(
      model == model_name,
      cov_spec == "bulk_pca",
      k == "none",
      indiv_frac == 1,
      cell_type %in% cell_types
    ) |>
    summarise(file_list = list(variant_file), .by = c(cell_type, cell_frac)) |>
    rowwise() |>
    mutate(qq_data = list(compute_qq_data(file_list))) |>
    ungroup() |>
    select(-file_list) |>
    unnest(cols = qq_data)
}

build_pb_plot_data <- function(model_name) {
  pb_data_files |>
    filter(
      model == model_name,
      int_cov == "none",
      indiv_frac == 1,
      cell_type %in% cell_types
    ) |>
    summarise(file_list = list(variant_file), .by = c(cell_type, cell_frac)) |>
    rowwise() |>
    mutate(qq_data = list(compute_qq_data(file_list))) |>
    ungroup() |>
    select(-file_list) |>
    unnest(cols = qq_data)
}

for (model_name in models) {
  plot_data <- if (model_name %in% c("p_glmm_sc", "lmm_sc")) {
    build_sc_plot_data(model_name)
  } else {
    build_pb_plot_data(model_name)
  }

  plot_perm_cell_frac(
    plot_data,
    method_lookup[[model_name]],
    sprintf("perm-cell-frac-%s-plot.pdf", model_name)
  )
}
