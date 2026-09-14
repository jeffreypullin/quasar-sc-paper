#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(tidyr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

sc_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(int_cov != "none")
if ("k" %in% names(sc_data_files)) {
  sc_data_files <- filter(sc_data_files, k == "none")
}
if ("cell_frac" %in% names(sc_data_files)) {
  sc_data_files <- filter(sc_data_files, cell_frac == 1)
}
if ("indiv_frac" %in% names(sc_data_files)) {
  sc_data_files <- filter(sc_data_files, indiv_frac == 1)
}
if (!"data_type" %in% names(sc_data_files)) {
  sc_data_files$data_type <- "counts"
}
sc_data_files <- sc_data_files |>
  mutate(
    data_type = ifelse(is.na(data_type) | data_type == "", "counts", data_type)
  ) |>
  filter(data_type %in% c("counts", "sct_counts"))

int_pvalue_col <- function(int_cov) {
  paste0("snp_x_", to_snake(int_cov), "_perm_pvalue")
}

read_variant_pvalues <- function(variant_files, int_cov, type) {
  pvalue_col <- if (type == "main") {
    "snp_pvalue"
  } else {
    int_pvalue_col(int_cov)
  }

  bind_rows(!!!map(
    variant_files,
    function(x) {
      read_tsv(x, show_col_types = FALSE) |>
        select(feature_id, any_of(c("snp_id", "maf", pvalue_col)))
    }
  )) |>
    filter(!is.na(.data[[pvalue_col]]))
}

compute_qq_data <- function(pvalue) {
  pvalue <- pvalue[!is.na(pvalue) & is.finite(pvalue) & pvalue > 0]
  n <- length(pvalue)
  if (n == 0) {
    return(tibble())
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

build_qq_plot_data <- function(data_files, type) {
  data_files |>
    summarise(
      file_list = list(variant_file),
      .by = c(cell_type, int_cov, model, data_type)
    ) |>
    rowwise() |>
    mutate(
      qq_data = list({
        pvals <- read_variant_pvalues(file_list, int_cov, type)
        pcol <- if (type == "main") "snp_pvalue" else int_pvalue_col(int_cov)
        compute_qq_data(pvals[[pcol]])
      })
    ) |>
    ungroup() |>
    select(-file_list) |>
    unnest(cols = qq_data)
}

qq_plot <- function(plot_data, ylab_suffix) {
  plot_data |>
    mutate(
      data_type = factor(
        data_type,
        levels = c("counts", "sct_counts")
      )
    ) |>
    ggplot(aes(
      log_x_bin_mid, log_y_pvalue,
      ymin = log_lower_ci, ymax = log_upper_ci,
      colour = data_type,
      fill = data_type
    )) +
    geom_ribbon(linetype = 2, alpha = 0.1, colour = NA) +
    geom_abline(linetype = "dashed") +
    geom_point(alpha = 0.8) +
    facet_wrap(~cell_type + int_cov + model) +
    scale_colour_manual(
      values = c(counts = "#2166AC", sct_counts = "#B2182B"),
      labels = c(counts = "counts", sct_counts = "sct_counts")
    ) +
    scale_fill_manual(
      values = c(counts = "#2166AC", sct_counts = "#B2182B"),
      guide = "none"
    ) +
    labs(
      x = "Expected -log10(p-value)",
      y = paste0("Observed -log10(p-value)", ylab_suffix),
      colour = "Phenotype"
    ) +
    theme_jp()
}

int_plot_data <- build_qq_plot_data(sc_data_files, "int")
int_p <- qq_plot(int_plot_data, "")

ggsave(
  "perm-sc-int-logcount-bins-plot.pdf",
  int_p,
  width = 12,
  height = 10
)

main_plot_data <- build_qq_plot_data(sc_data_files, "main")
main_p <- qq_plot(main_plot_data, " (main effect)")
ggsave(
  "perm-sc-int-logcount-bins-main-plot.pdf",
  main_p,
  width = 12,
  height = 10
)
