#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

int_cov <- "starcat_cytotoxic"
int_beta_col <- paste0("snp_x_", int_cov, "_beta")
int_se_col <- paste0("snp_x_", int_cov, "_se")
int_p_col <- paste0("snp_x_", int_cov, "_pvalue")

sc_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(k == "none", int_cov == "starcat_Cytotoxic")

all_variants <- sc_data_files |>
  rowwise() |>
  mutate(var_data = list(read_tsv(variant_file, show_col_types = FALSE))) |>
  unnest(var_data) |>
  ungroup()

output_cols <- c(
  "hit_type",
  "cell_type",
  "feature_id",
  "snp_id",
  "snp_beta",
  "snp_se",
  "snp_pvalue",
  int_beta_col,
  int_se_col,
  int_p_col
)

int_hits <- all_variants |>
  group_by(cell_type, feature_id) |>
  arrange(.data[[int_p_col]], .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  filter(.data[[int_p_col]] < 5e-8) |>
  mutate(hit_type = "interaction") |>
  select(all_of(output_cols)) |>
  arrange(.data[[int_p_col]])

main_only_hits <- all_variants |>
  group_by(cell_type, feature_id) |>
  arrange(snp_pvalue, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  filter(snp_pvalue < 5e-8, .data[[int_p_col]] > 0.5) |>
  anti_join(
    int_hits |> select(cell_type, feature_id),
    by = c("cell_type", "feature_id")
  ) |>
  arrange(snp_pvalue) |>
  slice_head(n = 5) |>
  mutate(hit_type = "main_only") |>
  select(all_of(output_cols))

sig_hits <- bind_rows(int_hits, main_only_hits)

write_tsv(sig_hits, "sc-int-eqtl-hits.tsv")
