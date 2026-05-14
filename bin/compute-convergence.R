#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
chrom <- args[1]
cell_type <- args[2]
int_cov <- args[3]
variant_data <- fread(args[4])
gene_props_data <- read_tsv(args[5], show_col_types = FALSE)

pvalue_col <- if (int_cov != "none") "snp_pvalue" else "pvalue"

if ("pvalue" %in% colnames(variant_data)) {
  pvalue_col <- "pvalue"
}

na_data <- variant_data |>
  group_by(feature_id) |>
  filter(any(is.na(.data[[pvalue_col]]))) |>
  ungroup() |>
  left_join(gene_props_data, by = "feature_id")

infinte_data <- variant_data |>
  group_by(feature_id) |>
  filter(any(is.infinite(.data[[pvalue_col]]))) |>
  ungroup() |>
  left_join(gene_props_data, by = "feature_id")

zero_data <- variant_data |>
  group_by(feature_id) |>
  filter(any(.data[[pvalue_col]] == 0)) |>
  ungroup() |>
  left_join(gene_props_data, by = "feature_id")

convergence_data <- bind_rows(
  na_data,
  infinte_data,
  zero_data
)

out_file <- paste0(chrom, "-", cell_type, "-problem-variants.tsv")
write_tsv(convergence_data, out_file)
