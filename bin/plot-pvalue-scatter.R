#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(tidyr)
})

set.seed(300)

args <- commandArgs(trailingOnly = TRUE)

pb_files_data <- read_tsv(args[1], show_col_types = FALSE)
sc_files_data <- read_tsv(args[2], show_col_types = FALSE)

pb_variant_data <- pb_files_data |>
  filter(model == "nb_glm") |>
  filter(chr == "chr21") |>
  filter(cell_frac == 1, indiv_frac == 1) |>
  filter(int_cov == "none") |>
  rowwise() |>
  mutate(var_data = list(read_tsv(variant_file, show_col_types = FALSE))) |>
  ungroup() |>
  unnest(cols = c(var_data))

sc_variant_data <- sc_files_data |>
  filter(chr == "chr21") |>
  filter(cell_frac == 1, indiv_frac == 1) |>
  rowwise() |>
  mutate(var_data = list(read_tsv(variant_file, show_col_types = FALSE))) |>
  ungroup() |>
  unnest(cols = c(var_data))

plot_data <- inner_join(
  sc_variant_data,
  pb_variant_data,
  by = c("snp_id", "cell_type", "feature_id"),
  suffix = c("_sc", "_pb")
) |>
  group_by(cell_type, feature_id) |>
  slice_sample(n = 100) |>
  ungroup() |>
  mutate(
    neg_log10_sc = -log10(pvalue_sc),
    neg_log10_pb = -log10(pvalue_pb)
  )

p <- plot_data |>
  ggplot(aes(x = neg_log10_pb, y = neg_log10_sc)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3) +
  facet_wrap(~cell_type, scales = "free") +
  labs(
    x = "Pseudobulk -log10 pvalue",
    y = "Single-cell -log10 pvalue",
  ) +
  theme_bw()

ggsave("pvalue-scatter.pdf", p, width = 7, height = 6)
