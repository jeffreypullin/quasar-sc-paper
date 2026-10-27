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
skipped_combos <- read_tsv(args[2], show_col_types = FALSE)

cell_types <- c("B_IN", "CD4_NC")

compute_qq_data <- function(files) {
  pvalue <- unlist(map(files, function(x) fread(x, select = "pvalue")$pvalue))
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

plot_perm_cell_frac <- function(plot_data, skipped, model_label, out_file) {
  if (nrow(plot_data) == 0 && nrow(skipped) == 0) {
    p <- ggplot() +
      labs(title = model_label) +
      theme_jp()
    ggsave(out_file, p, width = 8, height = 4)
    return(invisible())
  }

  frac_levels <- sort(unique(c(plot_data$cell_frac, skipped$cell_frac)))
  n_frac <- length(frac_levels)
  n_type <- n_distinct(c(plot_data$cell_type, skipped$cell_type))

  to_facets <- function(df) {
    df |>
      mutate(
        cell_type = factor(
          cell_type_lookup[cell_type],
          levels = cell_type_lookup[c("CD4_NC", "B_IN")]
        ),
        cell_frac = factor(cell_frac, levels = frac_levels)
      )
  }

  p <- plot_data |>
    to_facets() |>
    ggplot(aes(
      log_x_bin_mid, log_y_pvalue,
      ymin = log_lower_ci, ymax = log_upper_ci,
      colour = cell_type
    )) +
    geom_point(alpha = 0.8, size = 0.8) +
    geom_abline(linetype = "dashed") +
    geom_ribbon(linetype = 2, alpha = 0.1) +
    geom_blank(
      data = to_facets(skipped) |> cross_join(tibble(x = c(0, 3), y = c(0, 3))),
      aes(x, y),
      inherit.aes = FALSE
    ) +
    geom_text(
      data = to_facets(skipped),
      aes(x = 0, y = 3, label = sprintf("not run\n%d donors", n_donors)),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      colour = "grey40",
      size = 3.5
    ) +
    scale_colour_manual(
      values = setNames(cell_type_cols, cell_type_lookup[names(cell_type_cols)])
    ) +
    facet_grid2(
      vars(cell_type), vars(cell_frac),
      scales = "free_y",
      independent = "y",
      drop = TRUE
    ) +
    labs(
      title = model_label,
      x = "Expected -log10(p-value)",
      y = "Observed -log10(p-value)"
    ) +
    theme_jp() +
    theme(legend.position = "none")

  ggsave(
    out_file,
    p,
    width = max(8, 2.2 * n_frac),
    height = max(4, 3.5 * n_type),
    limitsize = FALSE
  )
}

build_sc_plot_data <- function(min_n_cells_val) {
  sc_data_files |>
    filter(
      model == "p_glmm_sc",
      cov_spec == "bulk_pca",
      as.character(k) == "none",
      as.numeric(indiv_frac) == 1,
      as.numeric(count_frac) == 1,
      cell_type %in% cell_types,
      as.numeric(min_n_cells) == min_n_cells_val
    ) |>
    summarise(file_list = list(variant_file), .by = c(cell_type, cell_frac)) |>
    mutate(qq_data = map(file_list, compute_qq_data)) |>
    select(cell_type, cell_frac, qq_data) |>
    unnest(cols = qq_data)
}

for (min_cells in c(5L, 10L)) {
  skipped <- skipped_combos |>
    filter(
      dataset == "onek1k",
      cell_type %in% cell_types,
      indiv_frac == 1,
      count_frac == 1,
      min_n_cells == min_cells
    ) |>
    select(cell_type, cell_frac, n_donors)
  plot_perm_cell_frac(
    build_sc_plot_data(min_cells),
    skipped,
    sprintf(
      "%s (>= %d cells per donor)",
      method_lookup[["p_glmm_sc"]],
      min_cells
    ),
    sprintf("perm-cell-frac-p_glmm_sc-min%d-plot.pdf", min_cells)
  )
}
