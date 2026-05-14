#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(qvalue)
  library(purrr)
  library(stringr)
  library(tidyr)
  library(ggh4x)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

sc_gwas_data_files <- read_tsv(args[1], show_col_types = FALSE)
pb_gwas_data_files <- read_tsv(args[2], show_col_types = FALSE)

gwas_data <- bind_rows(
  sc_gwas_data_files |>
    mutate(method = "single-cell"),
  pb_gwas_data_files |>
    mutate(method = paste0("pseudobulk-", model))
) |>
  rowwise() |>
  mutate(
    n_total = as.numeric(read_lines(n_variants, n_max = 1)),
    pvals   = list(read_tsv(sig_variant_file, show_col_types = FALSE)[[10]])
  ) |>
  ungroup() |>
  expand_grid(threshold = c(5e-8, 5e-11)) |>
  rowwise() |>
  mutate(
    expected = n_total * threshold,
    observed = sum(pvals < threshold, na.rm = TRUE)
  ) |>
  ungroup()

all_counts <- gwas_data |>
  group_by(method, model, cell_type, threshold) |>
  summarise(
    expected = sum(expected),
    observed = sum(observed),
    .groups = "drop"
  ) |>
  pivot_longer(
    cols = c(expected, observed),
    names_to = "type",
    values_to = "count"
  ) |>
  mutate(
    type = factor(
      type,
      levels = c("expected", "observed"),
      labels = c("Expected", "Observed")
    ),
    threshold = formatC(threshold, format = "e", digits = 0),
  )

p <- ggplot(all_counts, aes(x = threshold, y = count, fill = type)) +
  geom_col(position = "dodge") +
  facet_wrap(~method) +
  labs(
    x = "Significance threshold",
    y = "Number of significant variants",
    fill = NULL
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("perm-gwas-overall-plot.pdf", p, width = 8, height = 5)
