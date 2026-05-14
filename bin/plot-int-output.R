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

pb_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(int_cov != "none")

pb_data_files |>
  filter(int_cov == "age") |>
  rowwise() |>
  mutate(gene_tbl = list(read_tsv(region_file, show_col_types = FALSE))) |>
  unnest(gene_tbl) |>
  ungroup() |>
  mutate(
    bh = p.adjust(int_acat_pvalue, method = "BH"),
    .by = c(cell_type)
  ) |>
  filter(bh < 0.25) |>
  select(cell_type, feature_id, bh) |>
  arrange(bh)

2 + "fdjlks"

pb_data_files |>
  filter(int_cov == "age") |>
  rowwise() |>
  mutate(var_tbl = list(read_tsv(variant_file, show_col_types = FALSE))) |>
  unnest(var_tbl) |>
  filter(snp_x_age_pvalue < 5e-6) |>
  select(cell_type, feature_id, snp_id, snp_pvalue, snp_x_age_pvalue) |>
  distinct(feature_id, .keep_all = TRUE) |>
  arrange(snp_x_age_pvalue) |>
  print(n = 50)

2 + "fdskl"

pb_sig_sex_genes <- pb_data_files |>
  filter(int_cov == "sex") |>
  rowwise() |>
  mutate(gene_tbl = list(read_tsv(region_file, show_col_types = FALSE))) |>
  unnest(gene_tbl) |>
  ungroup() |>
  mutate(
    bh = p.adjust(int_acat_pvalue, method = "BH"),
    .by = c(cell_type)
  ) |>
  filter(bh < 0.01) |>
  select(cell_type, feature_id, bh) |>
  arrange(bh)

pb_data_files |>
  filter(int_cov == "age") |>
  rowwise() |>
  mutate(gene_tbl = list(read_tsv(region_file, show_col_types = FALSE))) |>
  unnest(gene_tbl) |>
  ungroup() |>
  mutate(
    bh = p.adjust(main_acat_pvalue, method = "BH"),
    .by = c(cell_type)
  ) |>
  filter(bh < 0.01) |>
  select(cell_type, feature_id, bh) |>
  arrange(bh) |>
  print(n = 20)


2 + "fdjsl"
