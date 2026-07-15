#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]
chrom <- args[2]
cell_type <- args[3]
int_cov <- args[4]
variant_data <- fread(args[5])
gene_props_data <- read_tsv(args[6], show_col_types = FALSE)

pvalue_col <- if (int_cov != "none") "snp_pvalue" else "pvalue"

if ("pvalue" %in% colnames(variant_data)) {
  pvalue_col <- "pvalue"
}

n_unique_feature_ids <- length(unique(variant_data$feature_id))

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

convergence_data <- bind_rows(
  na_data,
  infinte_data,
) |>
  mutate(n_unique_feature_ids = n_unique_feature_ids)

out_file <- paste0(dataset, "-", chrom, "-", cell_type, "-problem-variants.tsv")
write_tsv(convergence_data, out_file)
