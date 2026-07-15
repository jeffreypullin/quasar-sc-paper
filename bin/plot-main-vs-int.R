#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(stringr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

set.seed(300)

args <- commandArgs(trailingOnly = TRUE)

to_snake <- function(x) {
  out <- tolower(x)
  out <- gsub("[^a-z0-9]", "_", out)
  out <- gsub("_+", "_", out)
  out <- gsub("_$", "", out)
  out
}

neg_log10_p <- function(p) {
  p <- pmin(pmax(p, 1e-300), 1 - 1e-300)
  -log10(p)
}

sc_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(int_cov != "none") |>
  filter(k == "none") |>
  filter(cell_frac == 1, indiv_frac == 1)

int_covs <- sc_data_files |>
  distinct(int_cov) |>
  pull(int_cov)

for (int_cov in int_covs) {
  int_p_col <- paste0("snp_x_", to_snake(int_cov), "_pvalue")
  snake_int_cov <- to_snake(int_cov)

  analysis_files <- sc_data_files |>
    filter(int_cov == .env$int_cov)

  cell_type <- unique(analysis_files$cell_type)
  if (length(cell_type) != 1) {
    stop("Expected one cell type per interaction analysis, got: ", paste(cell_type, collapse = ", "))
  }
  cell_type <- cell_type[[1]]

  variant_data <- analysis_files |>
    pull(variant_file) |>
    unique() |>
    map(function(path) {
      read_tsv(path, show_col_types = FALSE) |>
        select(
          feature_id,
          snp_id,
          chrom,
          snp_pvalue,
          all_of(int_p_col)
        )
    }) |>
    bind_rows() |>
    filter(
      is.finite(snp_pvalue),
      is.finite(.data[[int_p_col]]),
      snp_pvalue > 0,
      .data[[int_p_col]] > 0
    ) |>
    mutate(
      neg_log10_main = neg_log10_p(snp_pvalue),
      neg_log10_int = neg_log10_p(.data[[int_p_col]])
    )

  n_sample <- min(50000, nrow(variant_data))
  plot_data <- variant_data |>
    slice_sample(n = n_sample)

  p <- plot_data |>
    ggplot(aes(x = neg_log10_main, y = neg_log10_int, color = chrom)) +
    geom_point(alpha = 0.2, size = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3) +
    labs(
      title = paste0("Main vs interaction p-values (", int_cov, ")"),
      subtitle = paste0(cell_type, "; n = ", n_sample, " SNP tests"),
      x = expression(-log[10](main~pvalue)),
      y = expression(-log[10](interaction~pvalue)),
      color = "Chromosome"
    ) +
    theme_jp()

  ggsave(
    paste0("main-vs-int-", snake_int_cov, "-plot.pdf"),
    p,
    width = 7,
    height = 6
  )
}
