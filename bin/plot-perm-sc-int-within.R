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
  filter(int_cov == "pseudotime")

compute_qq_data_int <- function(variant_files, prop_file, int_cov, type) {

  variant_data <- bind_rows(!!!map(
    variant_files,
    function(x) {
      read_tsv(x, show_col_types = FALSE) |>
        select(feature_id, maf, starts_with("snp"))
    })
  )

  if (type == "main") {
    pvalue_col <- "snp_pvalue"
  } else if (type == "int") {
    pvalue_col <- paste0("snp_x_", int_cov, "_perm_pvalue")
  }

  prop_data <- read_tsv(prop_file, show_col_types = FALSE)

  pvalue_data <- left_join(
    variant_data,
    prop_data,
    by = "feature_id"
  )

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

build_qq_plot_data <- function(data_files, type) {
  data_files |>
    filter(cell_frac == 1) |>
    filter(indiv_frac == 1) |>
    summarise(
      file_list = list(variant_file),
      prop_file = first(prop_file),
      .by = c(cell_type, int_cov)
    ) |>
    rowwise() |>
    mutate(qq_data = list(compute_qq_data_int(file_list, prop_file, int_cov, type))) |>
    ungroup() |>
    select(-file_list) |>
    unnest(cols = qq_data)
}

sc_int_plot_data <- build_qq_plot_data(sc_data_files, "int") |>
  mutate(source = "sc (within-sample perm)")

int_plot_data <- bind_rows(sc_int_plot_data)

int_p <- int_plot_data |>
  ggplot(aes(log_x_bin_mid, log_y_pvalue,
             ymin = log_lower_ci, ymax = log_upper_ci,
             colour = source)) +
  geom_point(alpha = 0.8) +
  geom_abline(linetype = "dashed") +
  geom_ribbon(linetype = 2, alpha = 0.1) +
  facet_wrap(~source + cell_type + int_cov) +
  labs(
    x = "Expected -log10(p-value)",
    y = "Observed -log10(p-value)",
    colour = "Source",
    title = "Within-sample pseudotime permutation null (interaction)"
  )

ggsave(
  "perm-sc-int-within-plot.pdf",
  int_p,
  width = 12,
  height = 10
)

sc_main_plot_data <- build_qq_plot_data(sc_data_files, "main") |>
  mutate(source = "sc (within-sample perm)")
main_plot_data <- bind_rows(sc_main_plot_data)

main_p <- main_plot_data |>
  ggplot(aes(log_x_bin_mid, log_y_pvalue,
             ymin = log_lower_ci, ymax = log_upper_ci,
             colour = source)) +
  geom_point(alpha = 0.8) +
  geom_abline(linetype = "dashed") +
  geom_ribbon(linetype = 2, alpha = 0.1) +
  facet_wrap(~source + cell_type + int_cov) +
  labs(
    x = "Expected -log10(p-value)",
    y = "Observed -log10(p-value)",
    colour = "Source",
    title = "Within-sample pseudotime permutation null (main effect)"
  )

ggsave(
  "perm-sc-int-within-main-plot.pdf",
  main_p,
  width = 12,
  height = 10
)
