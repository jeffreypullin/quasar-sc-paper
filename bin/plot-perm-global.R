#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(tidyr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

pb_data_files <- read_tsv(args[[1]], show_col_types = FALSE) |>
  filter(int_cov != "none", cell_frac == 1, indiv_frac == 1)

if ("count_frac" %in% names(pb_data_files)) {
  pb_data_files <- filter(pb_data_files, count_frac == 1)
}

empty_qq <- tibble(
  log_x_bin_mid = double(),
  log_y_pvalue = double(),
  log_lower_ci = double(),
  log_upper_ci = double()
)

qq_bins <- function(pvalue) {
  pvalue <- pvalue[!is.na(pvalue)]
  n <- length(pvalue)
  if (n == 0) {
    return(empty_qq)
  }
  m <- (1:n) / (n + 1)
  c <- abs(qnorm(0.05 / 2))
  v <- (1:n) * (n - (1:n) + 1) / (n + 1)^2 / (n + 2)
  s <- sqrt(v)
  log_x_pvalue <- -log10(m)
  x_bin <- cut(
    log_x_pvalue,
    breaks = seq(0, 6, by = 0.1),
    include.lowest = TRUE
  )
  tibble(
    x_bin,
    y_pvalue = sort(pvalue),
    lower_ci = m - c * s,
    upper_ci = m + c * s
  ) |>
    summarise(
      log_y_pvalue = -log10(mean(y_pvalue)),
      log_lower_ci = -log10(mean(lower_ci)),
      log_upper_ci = -log10(mean(upper_ci)),
      .by = c(x_bin)
    ) |>
    mutate(log_x_bin_mid = seq(0.05, 5.95, by = 0.1)[as.numeric(x_bin)])
}

compute_qq_data_pb_int <- function(variant_files, prop_file, int_cov, type) {
  variant_data <- bind_rows(!!!map(
    variant_files,
    function(x) {
      read_tsv(x, show_col_types = FALSE) |>
        select(feature_id, maf, starts_with("snp"))
    }
  ))

  pvalue_col <- if (type == "main") {
    "snp_pvalue"
  } else {
    paste0("snp_x_", to_snake(int_cov), "_pvalue")
  }

  prop_data <- read_tsv(prop_file, show_col_types = FALSE)
  pvalue_data <- left_join(variant_data, prop_data, by = "feature_id") |>
    filter(maf > 0.1, pb_non_zero_frac > 0.1)

  if (!pvalue_col %in% names(pvalue_data)) {
    return(empty_qq)
  }
  qq_bins(pvalue_data[[pvalue_col]])
}

build_qq_plot_data <- function(data_files, type) {
  if (nrow(data_files) == 0) {
    return(empty_qq |> mutate(
      cell_type = character(),
      int_cov = character(),
      model = character()
    ))
  }
  data_files |>
    summarise(
      file_list = list(variant_file),
      prop_file = first(prop_file),
      .by = c(cell_type, int_cov, model)
    ) |>
    rowwise() |>
    mutate(qq_data = list(compute_qq_data_pb_int(file_list, prop_file, int_cov, type))) |>
    ungroup() |>
    select(-file_list) |>
    unnest(cols = qq_data)
}

qq_plot <- function(plot_data, title) {
  if (nrow(plot_data) == 0) {
    return(
      ggplot() +
        labs(title = title) +
        theme_jp()
    )
  }
  plot_data |>
    ggplot(aes(
      log_x_bin_mid, log_y_pvalue,
      ymin = log_lower_ci, ymax = log_upper_ci,
      colour = cell_type
    )) +
    geom_point(alpha = 0.8) +
    geom_abline(linetype = "dashed") +
    geom_ribbon(linetype = 2, alpha = 0.1) +
    facet_wrap(~cell_type + int_cov + model) +
    labs(
      x = "Expected -log10(p-value)",
      y = "Observed -log10(p-value)",
      colour = "Cell type",
      title = title
    ) +
    theme_jp()
}

ggsave(
  "perm-global-pb-overall-int-plot.pdf",
  qq_plot(build_qq_plot_data(pb_data_files, "int"), "Genotype-permuted interaction p-values"),
  width = 12,
  height = 10
)

ggsave(
  "perm-global-pb-overall-int-main-plot.pdf",
  qq_plot(build_qq_plot_data(pb_data_files, "main"), "Genotype-permuted main-effect p-values"),
  width = 12,
  height = 10
)
