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

sc_data_files |>
  rowwise() |>
  mutate(var_data = list(read_tsv(variant_file, show_col_types = FALSE))) |>
  unnest(var_data) |>
  select(cell_type, feature_id, snp_id, snp_pvalue, snp_x_pseudotime_beta, snp_x_pseudotime_se, snp_x_pseudotime_pvalue) |>
  arrange(snp_x_pseudotime_pvalue) |>
  distinct(feature_id, .keep_all = TRUE) |>
  filter(snp_x_pseudotime_pvalue < 5e-6) |>
  select(feature_id, snp_pvalue, snp_x_pseudotime_pvalue) |>
  print(n = 30)

2 + "fdjlks"
