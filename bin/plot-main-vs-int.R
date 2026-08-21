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

neg_log10_p <- function(p) {
  p <- pmin(pmax(p, 1e-300), 1 - 1e-300)
  -log10(p)
}

method_levels <- unique(unname(method_lookup))

sc_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(int_cov != "none") |>
  filter(k == "none") |>
  filter(cell_frac == 1, indiv_frac == 1)

if ("data_type" %in% names(sc_data_files)) {
  sc_data_files <- filter(sc_data_files, data_type == "counts")
}
if ("sc_type" %in% names(sc_data_files)) {
  sc_data_files <- filter(sc_data_files, sc_type == "counts")
}

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
    distinct(model, variant_file) |>
    pmap(function(model, variant_file) {
      model_id <- model
      read_tsv(variant_file, show_col_types = FALSE) |>
        select(snp_pvalue, all_of(int_p_col)) |>
        rename(int_pvalue = all_of(int_p_col)) |>
        mutate(model = model_id)
    }) |>
    bind_rows() |>
    filter(
      is.finite(snp_pvalue),
      is.finite(int_pvalue),
      snp_pvalue > 0,
      int_pvalue > 0
    ) |>
    mutate(
      method = factor(method_lookup[model], levels = method_levels),
      neg_log10_main = neg_log10_p(snp_pvalue),
      neg_log10_int = neg_log10_p(int_pvalue)
    )

  # The Spearman correlation is the quantity of interest: how much main-effect
  # signal leaks into the conditional interaction test. Computed on every test,
  # not just the subsample that gets plotted.
  model_summary <- variant_data |>
    group_by(method) |>
    summarise(
      n_tests = n(),
      spearman = cor(snp_pvalue, int_pvalue, method = "spearman"),
      .groups = "drop"
    ) |>
    mutate(
      label = paste0(
        method, "\n", format(n_tests, big.mark = ","), " tests; ",
        "Spearman rho = ", sprintf("%.3f", spearman)
      )
    )

  print(model_summary |> select(method, n_tests, spearman))

  n_per_model <- min(50000, min(model_summary$n_tests))
  plot_data <- variant_data |>
    group_by(method) |>
    slice_sample(n = n_per_model) |>
    ungroup() |>
    left_join(select(model_summary, method, label), by = "method") |>
    mutate(label = factor(label, levels = model_summary$label))

  p <- plot_data |>
    ggplot(aes(x = neg_log10_main, y = neg_log10_int, color = method)) +
    geom_point(alpha = 0.2, size = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6,
                color = "black") +
    facet_wrap(vars(label)) +
    scale_color_manual(values = method_col_lookup, guide = "none") +
    labs(
      title = paste0("Main vs interaction p-values (", int_cov, ")"),
      subtitle = paste0(cell_type, "; n = ", n_per_model, " SNP tests per model"),
      x = expression(-log[10](main~pvalue)),
      y = expression(-log[10](interaction~pvalue))
    ) +
    theme_jp()

  ggsave(
    paste0("main-vs-int-", snake_int_cov, "-plot.pdf"),
    p,
    width = 5.5 * nrow(model_summary),
    height = 6
  )
}
