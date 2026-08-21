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
variant_data <- read_tsv(args[4], show_col_types = FALSE)
gene_props_data <- read_tsv(args[5], show_col_types = FALSE)

pvalue_col <- if (int_cov != "none" && int_cov != "null") "snp_pvalue" else "pvalue"

if ("pvalue" %in% colnames(variant_data)) {
  pvalue_col <- "pvalue"
}

n_sig_variant_data <- variant_data |>
  left_join(
    gene_props_data |>
      select(feature_id, pb_non_zero_frac),
    by = "feature_id"
  ) |>
  summarise(n_sig_variant = sum(.data[[pvalue_col]] < 5e-6, na.rm = TRUE))

out_file <- paste0(chrom, "-", cell_type, "-n-sig-variants.tsv")
write_tsv(n_sig_variant_data, out_file)
